"""Readout strip, chips and captions (gabes_ui.readout)."""
from copy import deepcopy
from xml.etree import ElementTree

import pytest

from gabes_ui.readout import (
    caption_html,
    clean_choice_label,
    metrics_table_markdown,
    magneto_status_summary,
    partition_readout,
    status_tone,
    strip_html,
)


MAGNETO_STATUSES = [
    dict(label="Width status", value="sample-grid consistent", kind="status"),
    dict(label="Absolute contrast status", value="external validation required", kind="status",
         help="Absolute contrast requires an independent measurement."),
    dict(label="Model scope", value="closed hyperfine pair; open-channel leakage and repump omitted",
         kind="status", help="Other ground-F leakage and repumping are omitted."),
]


def test_magneto_statuses_become_one_closed_disclosure_without_losing_details():
    before = deepcopy(MAGNETO_STATUSES)
    html = strip_html([], [], MAGNETO_STATUSES, compact_statuses=True)
    root = ElementTree.fromstring(html)
    details, = root.findall(".//details")
    assert "open" not in details.attrib
    assert "Model estimate" in "".join(details.find("summary").itertext())
    assert len(details.findall(".//dt")) == len(MAGNETO_STATUSES)
    text = "".join(details.itertext())
    for metric in MAGNETO_STATUSES:
        assert metric["label"] in text and metric["value"] in text
        if metric.get("help"):
            assert metric["help"] in text
    assert "not experimentally validated" in text
    assert "g-chip" not in html
    assert MAGNETO_STATUSES == before


@pytest.mark.parametrize("status, label", [
    ({"label": "Status", "value": "invalid transition"}, "Check transition"),
    ({"label": "Rotation readout status", "value": "undefined for circular input"}, "Check polarization"),
    ({"label": "Rotation readout status", "value": "unsupported for elliptical input"}, "Check polarization"),
    ({"label": "Width status", "value": "physical zero outside scan"}, "Check field range"),
    ({"label": "Slope status", "value": "nonconverged on sampled grid"}, "Check resolution"),
    ({"label": "Peak shot-noise SNR", "value": "unavailable"}, "Check detection"),
    ({"label": "Atomic solver", "value": "steady-state validation failed"}, "Check atomic solve"),
    ({"label": "Optimum status", "value": "no valid sensitivity samples"}, "Check sensitivity"),
    ({"label": "Optimum status", "value": "sampled model minimum at scan edge; expand optical scan"}, "Expand optical scan"),
    ({"label": "Signal linearity", "value": "signal peak exceeds tested local field interval"}, "Reduce field signal"),
    ({"label": "Optical excitation", "value": "optical elimination requires full-OBE validation"}, "Check model range"),
    ({"label": "Propagation scope", "value": "OD > 0.1: homogeneous intensity"}, "Check propagation"),
])
def test_compact_summary_surfaces_actionable_conditions(status, label):
    assert magneto_status_summary([status]) == (label, "warn")
    root = ElementTree.fromstring(strip_html([], [], [status], compact_statuses=True))
    assert label in "".join(root.find(".//summary").itertext())


def test_compact_status_details_escape_user_content_and_do_not_change_numeric_cells():
    metric = dict(label="Peak shot-noise SNR", value="1172 ×")
    status = dict(label="<scope>", value="<unsafe>", help="<script>alert('x')</script>")
    compact = ElementTree.fromstring(strip_html([metric], [], [status], compact_statuses=True))
    normal = ElementTree.fromstring(strip_html([metric], [], [status]))
    assert ElementTree.tostring(compact[0]) == ElementTree.tostring(normal[0])
    assert compact.find(".//script") is None
    assert "<unsafe>" in "".join(compact.itertext())

SAS_LIKE = [
    dict(label="SAS resolution", value="resolution-limited", tier="hero", kind="status",
         help="Fewer than six sample intervals span the FWHM."),
    dict(label="Sub-Doppler FWHM", value="18.87 MHz"),
    dict(label="Half-height edges", value="-1304.16 to -1285.29 MHz"),
    dict(label="Samples / FWHM", value="2.3"),
    dict(label="Scan-edge distance", value="3505.6 MHz"),
    dict(label="Peak OD", value="0.93"),
    dict(label="Gaussian Doppler FWHM", value="518.7 MHz"),
    dict(label="Peak phase shift", value="688.3 mrad"),
]


def test_status_never_takes_a_value_cell():
    heroes, secondary, statuses, overflow = partition_readout(SAS_LIKE, hero_count=1)
    assert [m["label"] for m in statuses] == ["SAS resolution"]
    assert [m["label"] for m in heroes] == ["Sub-Doppler FWHM"]
    assert len(secondary) == 4
    assert [m["label"] for m in overflow] == ["Gaussian Doppler FWHM", "Peak phase shift"]


def test_every_metric_lands_somewhere_exactly_once():
    parts = partition_readout(SAS_LIKE, hero_count=1)
    placed = [m["label"] for part in parts for m in part]
    assert sorted(placed) == sorted(m["label"] for m in SAS_LIKE)


def test_hero_tier_is_honoured_among_values():
    metrics = [dict(label="a", value="1"), dict(label="b", value="2", tier="hero")]
    heroes, secondary, _, _ = partition_readout(metrics, hero_count=1)
    assert heroes[0]["label"] == "b" and secondary[0]["label"] == "a"


def test_strip_markup_splits_number_and_unit_and_escapes():
    html = strip_html([dict(label="<FWHM>", value="18.87 MHz", help="a 'b'")], [], [])
    assert "&lt;FWHM&gt;" in html
    assert "18.87<span class='g-unit'>MHz</span>" in html
    assert "g-num" in html and "&#x27;b&#x27;" in html


def test_text_values_are_not_monospaced_and_status_becomes_a_chip():
    html = strip_html([dict(label="Zero-field feature", value="crossover")], [],
                      [SAS_LIKE[0]])
    assert "g-cell-value'>crossover" in html
    assert "g-chip g-chip--warn" in html and ">resolution-limited<" in html


def test_status_tone_reads_caveats():
    assert status_tone("resolution-limited") == "warn"
    assert status_tone("external validation required") == "warn"
    assert status_tone("resolved") == "neutral"


def test_delta_becomes_a_subline():
    html = strip_html([dict(label="G", value="15.50", delta="semi-empirical estimate")], [], [])
    assert "<div class='g-cell-sub'>semi-empirical estimate</div>" in html


def test_metrics_table_keeps_notes_and_escapes_pipes():
    md = metrics_table_markdown([dict(label="a|b", value="1", delta="d", help="h")])
    assert "| a\\|b | 1 | d h |" in md


def test_clean_choice_label_drops_leading_emoji_only():
    assert clean_choice_label("🕳️ EIT dip") == "EIT dip"
    assert clean_choice_label("✖️ Buffer LCA") == "Buffer LCA"
    assert clean_choice_label("(CPT) scan") == "(CPT) scan"
    assert clean_choice_label("⁸⁵Rb") == "⁸⁵Rb"


def test_caption_restores_physics_symbols():
    assert caption_html("85Rb Rydberg-EIT: Omega_c = 3.00 MHz, probe = 6.0 uW") == (
        "85Rb Rydberg-EIT: Ω<sub>c</sub> = 3.00 MHz, probe = 6.0 µW")
    assert caption_html("Delta = 0.9 GHz,  T = 121 C,  eta = 0.869") == (
        "Δ = 0.9 GHz, T = 121 °C, η = 0.869")
    assert caption_html("F=2 -> F'=1,  QWP=0.0 deg, I=0.80 mW/cm^2") == (
        "F=2 → F&#x27;=1, QWP=0.0°, I=0.80 mW/cm²")
    assert caption_html("beta decay <b>") == "beta decay &lt;b&gt;"
    assert caption_html(None) == ""
