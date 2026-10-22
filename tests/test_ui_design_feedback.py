"""Display qualifications and navigation preserve the underlying results."""
from copy import deepcopy
from types import SimpleNamespace
import re

from gabes_ui.diagnostics import diagnostic_sections
from gabes_ui.readout import caption_html, partition_readout, present_metrics, scheme_label, strip_html
from gabes_ui.shell import apply_scheme_link


def test_noise_scope_is_attached_and_numeric_units_are_not_duplicated():
    metrics = [dict(label="Squeezing Degree [dB]", value="-8.19 dB", tier="hero",
                    delta="conditional passive completion; physical squeezing unavailable"),
               dict(label="Seed gain G_s", value="394.88", tier="hero"),
               dict(label="Conjugate gain G_c", value="397.11")]
    original = deepcopy(metrics)
    presented = present_metrics(metrics, scheme_name="fwm", raw={"claim_gate": {"quantitative_gain_supported": False}})
    h, s, status, _ = partition_readout(presented)
    html = strip_html(h, s, status)
    assert "Squeezing Degree [dB]" not in html
    assert "Squeezing Degree</div>" in html and ">dB</span>" in html
    assert re.sub(r"<[^>]+>", "", html).count("physical squeezing unavailable") == 1
    assert "g-cell-status" in html and "g-chip--warn" in html
    assert "Seed gain G<sub>s</sub>" in html
    assert [m["value"] for m in presented[1:]] == ["395", "397"]
    assert metrics == original


def test_supported_and_measured_reference_gain_values_keep_their_precision():
    metrics = [dict(label="Seed gain G_s", value="394.88"),
               dict(label="Probe gain, pump-off reference", value="15.51")]
    result = present_metrics(metrics, scheme_name="fwm", raw={"claim_gate": {"quantitative_gain_supported": True}})
    assert result == metrics
    assert present_metrics(metrics, scheme_name="sas") == metrics


def test_symbols_are_restored_after_escaping_markup():
    assert caption_html("eta_ext = .869, G_s <b>") == "η<sub>ext</sub> = .869, G<sub>s</sub> &lt;b&gt;"
    assert caption_html("η_probe") == "η<sub>probe</sub>"


def test_short_scheme_labels_do_not_change_identity_and_links_are_validated():
    title = "Four-wave mixing (Squeezing / Biphoton)"
    schemes = [SimpleNamespace(name="fwm", title=title)]
    assert scheme_label(title) == "Four-wave mixing"
    state = {}
    apply_scheme_link(state, {"scheme": "unknown"}, schemes)
    assert not state
    apply_scheme_link(state, {"scheme": "fwm"}, schemes)
    assert state["_scheme_choice"] == title
    state["_scheme_choice"] = "another manual choice"
    apply_scheme_link(state, {"scheme": "fwm"}, schemes)
    assert state["_scheme_choice"] == "another manual choice"


def test_grouped_diagnostics_keep_every_row_and_full_reason_text():
    rows = ["| Model scope | diagnostic |", "| Atomic solver | Floquet |",
            "| 2-D reference velocity geometry | (v_z, v_x) |",
            "| Passive completion minimum eigenvalue | -1.234e-12 |",
            "| Unfamiliar future diagnostic | 123.456 |"]
    reasons = "missing atomic diffusion; production uses 1-D velocities; no independent validation"
    table = "| Check | Value |\n|---|---|\n" + "\n".join(rows) + f"\n| Claim-gate reasons | {reasons} |"
    sections, separated = diagnostic_sections(table)
    assert sorted(row for group in sections.values() for row in group) == sorted(rows)
    assert "; ".join(separated) == reasons
