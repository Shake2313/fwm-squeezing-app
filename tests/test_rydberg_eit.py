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

from gabes import constants, schemes  # noqa: E402
from gabes.schemes.rydberg import MHZ, RydbergEITScheme  # noqa: E402


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


def test_absolute_sensitivity_is_computed_without_reference_injection():
    sc = schemes.get("rydberg_eit")
    assert isinstance(sc, RydbergEITScheme)
    assert sc.REFERENCE_SENSITIVITY_NV_CM_SQRT_HZ == 12.5
    assert sc.REFERENCE_PSN_LIMIT_NV_CM_SQRT_HZ == 11.2
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    view = sc.observables(sc.compute(params), params)
    by_label = {metric["label"]: metric for metric in view["metrics"]}
    tables = " ".join(table["markdown"].lower() for table in view["tables"])
    plt.close(view["figure"])
    total = float(by_label["Total sensitivity"]["value"].split()[0])
    psn = float(by_label["PSN-limited sensitivity"]["value"].split()[0])
    assert total >= psn > 0.0
    # The public result comes from the declared dipole/detector chain.  It is not
    # forced to either experimental reference constant.
    assert abs(total - sc.REFERENCE_SENSITIVITY_NV_CM_SQRT_HZ) > 0.1
    assert abs(psn - sc.REFERENCE_PSN_LIMIT_NV_CM_SQRT_HZ) > 0.1
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


def test_at_heroes_prioritize_total_and_psn_sensitivity():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    view = sc.headless_observables(sc.compute(params), params)
    heroes = [m for m in view["metrics"] if m.get("tier") == "hero"]
    assert [m["label"] for m in heroes] == [
        "Total sensitivity", "PSN-limited sensitivity"]


def test_detuned_at_promotes_informative_center_shift():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    params["mw_detuning_mhz"] = 4.0
    view = sc.headless_observables(sc.compute(params), params)
    heroes = [m for m in view["metrics"] if m.get("tier") == "hero"]
    assert [m["label"] for m in heroes] == [
        "Total sensitivity", "AT center shift"]


def test_weak_lo_reports_unresolved_at_status():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    params["lo_rabi_mhz"] = 0.1
    view = sc.headless_observables(sc.compute(params), params)
    heroes = [m for m in view["metrics"] if m.get("tier") == "hero"]
    assert [m["label"] for m in heroes] == ["Total sensitivity", "AT status"]
    assert heroes[1].get("kind") == "status"


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


def test_at_field_conversion_and_peak_rms_conventions_are_explicit():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    view = sc.headless_observables(sc.compute(params), params)
    split_mhz = _metric_value(view, "RF AT splitting")
    field_mvm = _metric_value(view, "RF field from fitted Ω_RF")
    split_estimate_mvm = _metric_value(view, "AT split field estimate")
    coupling = sc._rf_coupling(params)
    expected_mvm = coupling.field_from_cyclic_rabi_hz(
        params["lo_rabi_mhz"] * 1e6) * 1e3
    assert np.isclose(field_mvm, expected_mvm, rtol=5e-3)
    assert np.isclose(
        split_estimate_mvm,
        coupling.field_from_at_splitting_hz(split_mhz * 1e6) * 1e3,
        rtol=5e-3,
    )

    peak = dict(params, rf_field_convention="Peak")
    peak_view = sc.headless_observables(sc.compute(peak), peak)
    peak_field_mvm = _metric_value(peak_view, "RF field from fitted Ω_RF")
    assert np.isclose(peak_field_mvm / field_mvm, np.sqrt(2.0), rtol=5e-3)

    detuned = dict(params, mw_detuning_mhz=4.0)
    detuned_view = sc.headless_observables(sc.compute(detuned), detuned)
    detuned_field_mvm = _metric_value(detuned_view, "RF field from fitted Ω_RF")
    assert np.isclose(detuned_field_mvm, field_mvm, rtol=5e-3)


def test_superhet_low_if_matches_static_lo_derivative():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    params["if_khz"] = 0.001
    raw = sc.compute(params)
    x, transmission, _ = sc._transmission(raw["chi_bar"], raw, params)
    result = sc._superheterodyne_readout(raw, params, x, transmission)
    detuning = result["optimum_detuning_mhz"]
    response_index = int(np.argmin(np.abs(
        result["response_detuning_mhz"] - detuning)))
    dynamic = result["transmission_phasor_per_angular_rabi"][response_index]

    epsilon_mhz = 1e-4
    static_values = []
    for sign in (-1.0, 1.0):
        shifted = dict(
            params,
            lo_rabi_mhz=params["lo_rabi_mhz"] + sign * epsilon_mhz,
        )
        shifted_raw = sc.compute(shifted)
        shifted_x, shifted_t, _ = sc._transmission(
            shifted_raw["chi_bar"], shifted_raw, shifted)
        static_values.append(float(np.interp(detuning, shifted_x, shifted_t)))
    static_derivative = (
        (static_values[1] - static_values[0])
        / (2.0 * epsilon_mhz * MHZ)
    )
    assert np.isclose(dynamic.real, static_derivative, rtol=3e-3)
    assert abs(dynamic.imag) < abs(dynamic.real) * 1e-3


def test_electronics_noise_raises_total_but_not_psn_sensitivity():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    base = sc._readout(sc.compute(params), params)["superhet"]["sensitivity"]
    noisy = dict(params, detector_electronic_noise_pa_sqrt_hz=5.0)
    with_electronics = sc._readout(
        sc.compute(noisy), noisy)["superhet"]["sensitivity"]
    assert np.isclose(
        with_electronics.psn_field_asd_v_m_per_sqrt_hz,
        base.psn_field_asd_v_m_per_sqrt_hz,
    )
    assert (with_electronics.total_field_asd_v_m_per_sqrt_hz
            > base.total_field_asd_v_m_per_sqrt_hz)


def test_superhet_optimizes_noise_equivalent_field_not_only_response():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    raw = sc.compute(params)
    result = sc._readout(raw, params)["superhet"]
    finite = np.asarray(result["total_field_asd_v_m_per_sqrt_hz"], dtype=float)
    local_best = int(np.argmin(finite))
    assert np.isclose(
        result["sensitivity"].total_field_asd_v_m_per_sqrt_hz,
        finite[local_best],
    )
    assert np.isclose(
        result["optimum_detuning_mhz"],
        result["response_detuning_mhz"][local_best],
    )


def test_separated_temperature_uses_effective_motion_and_cold_spot_pressure():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    separated = dict(
        params,
        temperature_model="Separated",
        heater_setpoint_c=50.0,
        effective_temp_c=44.0,
        cold_spot_temp_c=42.0,
    )
    raw = sc.compute(separated)
    cold_density = sc._temperature_state(separated).saturated_cold_spot_density_m3()
    expected_density = cold_density * (42.0 + 273.15) / (44.0 + 273.15)
    assert raw["heater_setpoint_c"] == 50.0
    assert raw["effective_temp_c"] == 44.0
    assert raw["cold_spot_temp_c"] == 42.0
    assert np.isclose(raw["N"], expected_density)
    # Heater metadata alone must not change the atomic solve.
    other_setpoint = sc.compute(dict(separated, heater_setpoint_c=70.0))
    assert np.array_equal(other_setpoint["chi_bar"], raw["chi_bar"])


def test_density_dephasing_and_effective_atom_number_are_exposed():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    warm = dict(
        params, temp_c=50.0,
        density_dephasing_mhz_per_1e16_m3=0.02,
        atom_participation_fraction=0.5,
        beam_overlap_efficiency=0.8,
    )
    raw = sc.compute(warm)
    readout = sc._readout(raw, warm)
    assert raw["density_dephasing_mhz"] > 0.0
    full = sc._effective_atom_number(
        raw, dict(warm, atom_participation_fraction=1.0,
                  beam_overlap_efficiency=1.0))
    assert np.isclose(readout["effective_atoms"].atoms, full.atoms * 0.4)


def test_scheme_sam_calibration_reports_uncertainty_and_far_field_warning():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    params.update({
        "sam_enabled": "On",
        "sam_source_power_dbm": -40.0,
        "sam_antenna_gain_dbi": 10.0,
        "sam_distance_m": 0.10,
        "sam_source_power_std_db": 0.2,
        "sam_distance_std_m": 0.002,
        "sam_antenna_max_dimension_m": 0.10,
    })
    readout = sc._readout(sc.compute(params), params)
    calibration = readout["sam"]
    assert calibration.field_v_m > 0.0
    assert calibration.standard_uncertainty_v_m > 0.0
    assert calibration.far_field_ratio < 1.0
    assert calibration.warning
    labels = {metric["label"] for metric in readout["metrics"]}
    assert "SAM RF field" in labels
    assert "Fitted-RF/SAM field ratio" in labels


def test_cell_heating_extra_view_reports_absolute_sensitivity():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    event = sc.extra_views()[1]
    sweep = event.compute(params)
    assert len(sweep["temperatures_c"]) == 10
    assert len(sweep["total_sensitivity_nv_cm_sqrt_hz"]) == 10
    assert np.all(np.isfinite(sweep["total_sensitivity_nv_cm_sqrt_hz"]))
    assert sweep["best_total_sensitivity_nv_cm_sqrt_hz"] == min(
        sweep["total_sensitivity_nv_cm_sqrt_hz"])
    figure = event.render(sweep)
    assert len(figure.axes) >= 4
    plt.close(figure)


def test_rydberg_view_declares_experimental_csv_overlay():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["EIT"]
    view = sc.headless_observables(sc.compute(params), params)
    assert view["comparison"]["axis_index"] == 0
    assert view["comparison"]["x_unit"] == "MHz"


if __name__ == "__main__":
    test_reference_defaults_match_rydberg_eit_paper()
    test_reference_eit_linewidth_near_experiment()
    test_microwave_at_splitting_tracks_lo_rabi()
    test_absolute_sensitivity_is_computed_without_reference_injection()
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
