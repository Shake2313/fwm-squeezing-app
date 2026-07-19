"""
Reference checks for the Rydberg-EIT electrometry scheme.

The public UI shows the static spectrum only. Experimental sensitivity values
from arXiv:2606.04354 are kept as internal constants and tested here.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import schemes  # noqa: E402
from gabes.schemes.rydberg import RydbergEITScheme  # noqa: E402


def _metric_value(view, label):
    for metric in view["metrics"]:
        if metric["label"] == label:
            return float(metric["value"].split()[0])
    raise AssertionError(f"missing metric {label!r}: {view['metrics']}")


def test_reference_defaults_match_rydberg_eit_paper():
    sc = schemes.get("rydberg_eit")
    sets = sc.recommended_defaults(sc.defaults())
    ref = sets["AT electrometry"]
    assert ref["probe_power_uw"] == 6.0
    assert ref["coupling_power_mw"] == 30.0
    assert ref["beam_diameter_mm"] == 0.15
    assert ref["cell_mm"] == 50.0
    assert ref["coupling_rabi_mhz"] == 3.0
    assert ref["lo_rabi_mhz"] == 3.7
    assert ref["mw_frequency_ghz"] == 37.0
    assert ref["if_khz"] == 40.0


def test_reference_eit_linewidth_near_experiment():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["EIT"]
    view = sc.observables(sc.compute(params), params)
    linewidth = _metric_value(view, "EIT linewidth")
    plt.close(view["figure"])
    assert 1.3 <= linewidth <= 1.9


def test_microwave_at_splitting_tracks_lo_rabi():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    view = sc.observables(sc.compute(params), params)
    split = _metric_value(view, "RF AT splitting")
    plt.close(view["figure"])
    # The transmission-peak separation tracks Ω_LO but sits a few % inside it
    # (the dressed peaks are pulled toward line centre by the absorptive
    # background), so it is a constant fraction just below 1, not exactly 1.
    assert 0.88 <= split / params["lo_rabi_mhz"] <= 1.0


def test_sensitivity_reference_constants_are_internal_only():
    sc = schemes.get("rydberg_eit")
    assert isinstance(sc, RydbergEITScheme)
    assert sc.REFERENCE_SENSITIVITY_NV_CM_SQRT_HZ == 12.5
    assert sc.REFERENCE_PSN_LIMIT_NV_CM_SQRT_HZ == 11.2
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    view = sc.observables(sc.compute(params), params)
    labels = " ".join(metric["label"].lower() for metric in view["metrics"])
    tables = " ".join(table["markdown"].lower() for table in view["tables"])
    plt.close(view["figure"])
    assert "sensitivity" not in labels
    assert "psn" not in labels
    assert "12.5" not in tables and "11.2" not in tables


def test_if_proxy_metric_and_temperature_dephasing_are_opt_in():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["EIT"]
    view = sc.observables(sc.compute(params), params)
    labels = {metric["label"] for metric in view["metrics"]}
    plt.close(view["figure"])
    assert "IF discriminator" in labels
    assert "IF optimum detuning" in labels

    warm = dict(params, temp_c=60.0)
    assert sc.compute(warm)["temperature_dephasing_mhz"] == 0.0
    broadened = sc.compute(dict(warm, temp_dephasing_mhz_per_c=0.01))
    assert abs(broadened["temperature_dephasing_mhz"] - 0.40) < 1e-12


def test_at_heroes_prioritize_split_and_if_discriminator():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    view = sc.headless_observables(sc.compute(params), params)
    heroes = [m for m in view["metrics"] if m.get("tier") == "hero"]
    assert [m["label"] for m in heroes] == [
        "RF AT splitting", "IF discriminator"]


def test_detuned_at_promotes_informative_center_shift():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    params["mw_detuning_mhz"] = 4.0
    view = sc.headless_observables(sc.compute(params), params)
    heroes = [m for m in view["metrics"] if m.get("tier") == "hero"]
    assert [m["label"] for m in heroes] == [
        "RF AT splitting", "AT center shift"]


def test_weak_lo_reports_unresolved_at_status():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    params["lo_rabi_mhz"] = 0.1
    view = sc.headless_observables(sc.compute(params), params)
    heroes = [m for m in view["metrics"] if m.get("tier") == "hero"]
    assert [m["label"] for m in heroes] == ["AT status", "IF discriminator"]
    assert heroes[0].get("kind") == "status"


def test_unresolved_eit_width_uses_status_hero():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["EIT"]
    raw = sc.compute(params)
    raw["chi_bar"] = np.zeros_like(raw["chi_bar"])
    view = sc.headless_observables(raw, params)
    heroes = [m for m in view["metrics"] if m.get("tier") == "hero"]
    assert [m["label"] for m in heroes] == [
        "Transmission at resonance", "EIT status"]
    assert heroes[1].get("kind") == "status"
    assert all("nan" not in str(m["value"]).lower() for m in heroes)


def test_coupling_power_and_waist_drive_rabi():
    """481 nm coupling power/waist set Ω_c via √(P/d²), anchored at reference."""
    sc = schemes.get("rydberg_eit")
    base = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    # Reference operating point reproduces the fitted anchor exactly.
    assert abs(sc.compute(base)["coupling_rabi_mhz"] - 3.0) < 1e-9
    # Doubling power scales Ω_c by √2 (intensity ∝ power at fixed waist).
    hi_p = dict(base, coupling_power_mw=2 * base["coupling_power_mw"])
    assert abs(sc.compute(hi_p)["coupling_rabi_mhz"] - 3.0 * 2 ** 0.5) < 1e-6
    # Doubling the beam diameter halves Ω_c (intensity ∝ 1/d²).
    wide = dict(base, beam_diameter_mm=2 * base["beam_diameter_mm"])
    assert abs(sc.compute(wide)["coupling_rabi_mhz"] - 1.5) < 1e-6


def test_at_center_shift_tracks_microwave_detuning():
    """The dressed-transparency centre is ~0 on resonance and shifts by ≈ −Δ_mw/2
    (the dressed-doublet midpoint) when the microwave is detuned."""
    sc = schemes.get("rydberg_eit")
    base = sc.recommended_defaults(sc.defaults())["AT electrometry"]

    on_res = sc.observables(sc.compute(base), base)
    center0 = _metric_value(on_res, "AT center shift")
    plt.close(on_res["figure"])
    assert abs(center0) < 0.1

    pos = dict(base, mw_detuning_mhz=4.0)
    vp = sc.observables(sc.compute(pos), pos)
    cp = _metric_value(vp, "AT center shift")
    plt.close(vp["figure"])
    neg = dict(base, mw_detuning_mhz=-4.0)
    vn = sc.observables(sc.compute(neg), neg)
    cn = _metric_value(vn, "AT center shift")
    plt.close(vn["figure"])
    # Dressed centre moves to −Δ_mw/2: positive detuning -> negative shift.
    assert cp < -0.1 and cn > 0.1
    assert abs(cp + cn) < 1e-6      # antisymmetric in detuning


def test_doppler_on_broadens_eit_linewidth():
    """Residual two-photon Doppler (per-level k) washes out the narrow EIT
    feature, so Doppler-on is broader than the suppressed static model."""
    sc = schemes.get("rydberg_eit")
    eit = sc.recommended_defaults(sc.defaults())["EIT"]
    off = dict(eit, doppler="off")
    on = dict(eit, doppler="on")
    w_off = _metric_value(sc.observables(sc.compute(off), off), "EIT linewidth")
    w_on = _metric_value(sc.observables(sc.compute(on), on), "EIT linewidth")
    plt.close("all")
    assert w_on > w_off


def test_per_level_doppler_ratio_is_backward_compatible():
    """A doppler_ratios entry of 1.0 reproduces the plain doppler_levels S_v, so
    existing schemes are unchanged; the Rydberg ladder carries a residual ratio."""
    import numpy as np
    from gabes import atoms
    plain = atoms.AtomModel(
        name="t", n_levels=2, labels=("g", "e"), ground=(0,), excited=(1,),
        decay=((1, 0, 1.0),), dephasing=(), doppler_levels=(1,))
    explicit = atoms.AtomModel(
        name="t", n_levels=2, labels=("g", "e"), ground=(0,), excited=(1,),
        decay=((1, 0, 1.0),), dephasing=(), doppler_levels=(1,),
        doppler_ratios=((1, 1.0),))
    assert np.allclose(plain.S_v, explicit.S_v)

    ryd = RydbergEITScheme()._atom(1.0e6, 1.0e6)
    assert not np.allclose(ryd.S_v, 0.0)   # residual two-photon Doppler is carried


# --- Ju et al. (arXiv:2606.04354) Fig. 2 matching ---

def test_default_view_is_eit_fig2a():
    """The scheme opens on the EIT regime so the landing figure is Fig. 2(a)."""
    sc = schemes.get("rydberg_eit")
    view_spec = next(s for s in sc.param_schema() if s.name == "view")
    assert view_spec.default == "EIT"


def test_reference_eit_compensation_pair():
    """Fig. 2(a): with B-field compensation ≈1.6 MHz, without ≈1.9 MHz, and the
    compensated transparency peak is the taller of the two."""
    sc = schemes.get("rydberg_eit")
    eit = sc.recommended_defaults(sc.defaults())["EIT"]
    raw = sc.compute(eit)
    x, T_c, _ = sc._transmission(raw["chi_bar"], raw, eit)
    _, T_u, _ = sc._transmission(raw["chi_bar_uncomp"], raw, eit)
    w_c, _ = sc._eit_features(x, T_c)
    w_u, _ = sc._eit_features(x, T_u)
    assert 1.45 <= w_c <= 1.75            # compensated ~1.6 MHz
    assert 1.8 <= w_u <= 2.1              # uncompensated ~1.9 MHz
    assert w_u > w_c                       # compensation narrows the line
    ic = int(np.argmin(np.abs(x)))
    assert T_c[ic] > T_u[ic]               # and raises the transparency peak


def test_probe_power_broadens_eit():
    """Fig. 2(b): raising the probe power power-broadens the EIT linewidth."""
    sc = schemes.get("rydberg_eit")
    eit = sc.recommended_defaults(sc.defaults())["EIT"]
    widths = []
    for p_uw in (1.0, 6.0, 10.0):
        p = dict(eit, probe_power_uw=p_uw)
        raw = sc.compute(p)
        x, T, _ = sc._transmission(raw["chi_bar"], raw, p)
        widths.append(sc._eit_features(x, T)[0])
    assert widths[0] < widths[1] < widths[2]


def test_fig2b_extra_view_runs():
    """The Fig. 2(b) probe-power panel computes a picklable sweep and renders a
    two-panel figure with monotonically rising peak amplitude."""
    import matplotlib
    matplotlib.use("Agg")
    sc = schemes.get("rydberg_eit")
    eit = sc.recommended_defaults(sc.defaults())["EIT"]
    ev = sc.extra_views()[0]
    s = ev.compute(eit)
    assert len(s["powers"]) > 4
    assert all(b >= a - 1e-9 for a, b in zip(s["comp"]["amp"], s["comp"]["amp"][1:]))
    fig = ev.render(s)
    assert len(fig.axes) == 2
    plt.close(fig)


if __name__ == "__main__":
    test_reference_defaults_match_rydberg_eit_paper()
    test_reference_eit_linewidth_near_experiment()
    test_microwave_at_splitting_tracks_lo_rabi()
    test_sensitivity_reference_constants_are_internal_only()
    test_if_proxy_metric_and_temperature_dephasing_are_opt_in()
    test_coupling_power_and_waist_drive_rabi()
    test_at_center_shift_tracks_microwave_detuning()
    test_doppler_on_broadens_eit_linewidth()
    test_per_level_doppler_ratio_is_backward_compatible()
    test_default_view_is_eit_fig2a()
    test_reference_eit_compensation_pair()
    test_probe_power_broadens_eit()
    test_fig2b_extra_view_runs()
    print("Rydberg-EIT reference checks OK.")
