"""
Generic SFWM / biphoton checks.

    python tests/test_fwm_generic.py   # or: pytest tests/test_fwm_generic.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import hyperfine, observables, zeeman  # noqa: E402
from gabes.schemes import fwm  # noqa: E402


def _params(**updates):
    scheme = fwm.FWMScheme()
    params = scheme.defaults()
    params.update({
        "mode": fwm.MODE_BIPHOTON,
        "biphoton_velocity_step": 10.0,
        "tau_max_ns": 12.0,
    })
    params.update(updates)
    return params


def _recommended_params(**updates):
    scheme = fwm.FWMScheme()
    params = _params(**updates)
    defaults = scheme.recommended_defaults(params)[fwm.MODE_BIPHOTON]
    defaults.update({
        "biphoton_velocity_step": params["biphoton_velocity_step"],
        "tau_max_ns": params["tau_max_ns"],
    })
    defaults.update(updates)
    return defaults


def _stats(raw, params, target=True):
    return observables.biphoton_stats(
        raw["tau_axis_ns"], raw["psi_tau"], raw["pair_rate_cps"],
        signal_eff=params["signal_eff_pct"] / 100.0,
        idler_eff=params["idler_eff_pct"] / 100.0,
        dark_signal_cps=params["dark_signal_cps"],
        dark_idler_cps=params["dark_idler_cps"],
        coincidence_window_ns=params["coincidence_window_ns"],
        timing_jitter_ns=params["timing_jitter_ns"],
        filter_bandwidth_mhz=params["filter_bandwidth_mhz"],
        target_g2_peak=raw["topology"].target_g2_peak if target else None,
    )


def test_topology_energy_and_roles():
    for topo in (fwm.TOPOLOGY_RB87_TELECOM, fwm.TOPOLOGY_CS_BTW,
                 fwm.TOPOLOGY_DIAMOND):
        spec = fwm.topology_from_params(_params(topology=topo))
        roles = {field.role for field in spec.fields}
        assert {"pump", "coupling", "signal", "idler"} == roles
        assert abs(fwm.energy_mismatch_hz(spec.fields)) < 2e8


def test_phase_matching_reference_angle_is_maximum():
    params = _recommended_params(topology=fwm.TOPOLOGY_RB87_TELECOM)
    spec = fwm.topology_from_params(params)
    L = spec.default_cell_mm * 1e-3
    angles = np.linspace(0.0, 2.5, 101)
    weights = np.array([
        fwm.phase_matching_weight(
            np.array([fwm.phase_mismatch_vector(
                spec.fields,
                signal_angle_deg=params["signal_angle_deg"],
                idler_angle_deg=a,
                reference_delta_k=spec.reference_delta_k,
            )["delta_k_vector"]]), L)[0]
        for a in angles
    ])
    best = angles[int(np.argmax(weights))]
    assert abs(best - params["idler_angle_deg"]) <= 0.05
    exact = fwm.phase_matching_weight(
        np.array([fwm.phase_mismatch_vector(
            spec.fields,
            signal_angle_deg=params["signal_angle_deg"],
            idler_angle_deg=params["idler_angle_deg"],
            reference_delta_k=spec.reference_delta_k,
        )["delta_k_vector"]]), L)[0]
    assert np.isclose(exact, 1.0)


def test_rb87_default_vector_phase_match_has_positive_rate():
    scheme = fwm.FWMScheme()
    params = _recommended_params(topology=fwm.TOPOLOGY_RB87_TELECOM)
    raw = scheme.compute(params)
    assert raw["phase_match_weight"] > 0.99
    assert raw["pair_rate_cps"] > 0
    assert params["signal_angle_deg"] == 1.5
    assert abs(params["idler_angle_deg"] - 0.77) < 0.02


def test_rb87_equal_angles_are_transversely_suppressed():
    scheme = fwm.FWMScheme()
    matched = _recommended_params(topology=fwm.TOPOLOGY_RB87_TELECOM)
    good = scheme.compute(matched)
    bad = dict(matched, signal_angle_deg=1.5, idler_angle_deg=1.5,
               signal_side=fwm.SIDE_PLUS, idler_side=fwm.SIDE_PLUS)
    raw = scheme.compute(bad)
    assert raw["phase_match_weight"] < 1e-4
    assert raw["pair_rate_cps"] < good["pair_rate_cps"] * 1e-4


def test_side_flip_suppresses_matched_geometry():
    scheme = fwm.FWMScheme()
    params = _recommended_params(topology=fwm.TOPOLOGY_RB87_TELECOM)
    params["idler_side"] = fwm.SIDE_MINUS
    raw = scheme.compute(params)
    assert raw["phase_match_weight"] < 1e-4
    assert raw["pair_rate_cps"] < 1e-2


def test_detector_background_and_window_reduce_car():
    tau = np.linspace(0.0, 10.0, 201)
    wave = np.exp(-tau / 2.0)
    base = observables.biphoton_stats(
        tau, wave, 5_000.0, signal_eff=0.2, idler_eff=0.2,
        dark_signal_cps=100.0, dark_idler_cps=100.0,
        coincidence_window_ns=1.0)
    noisy = observables.biphoton_stats(
        tau, wave, 5_000.0, signal_eff=0.2, idler_eff=0.2,
        dark_signal_cps=20_000.0, dark_idler_cps=20_000.0,
        coincidence_window_ns=50.0)
    assert noisy["CAR"] < base["CAR"]
    assert noisy["g2_peak"] < base["g2_peak"]


def test_timing_jitter_broadens_waveform():
    tau = np.linspace(0.0, 10.0, 401)
    wave = np.exp(-tau / 0.8)
    sharp = observables.biphoton_stats(tau, wave, 1_000.0, timing_jitter_ns=0.0)
    broad = observables.biphoton_stats(tau, wave, 1_000.0, timing_jitter_ns=0.8)
    assert broad["fwhm_ns"] > sharp["fwhm_ns"]
    assert broad["g2_SI_tau"].shape == tau.shape


def test_long_jitter_kernel_preserves_axis_length():
    tau = np.linspace(0.0, 1.0, 481)
    wave = np.exp(-tau / 0.1)
    stats = observables.biphoton_stats(tau, wave, 1_000.0, timing_jitter_ns=0.55)
    assert stats["g2_SI_tau"].shape == tau.shape
    assert stats["tau_axis_ns"].shape == tau.shape


def test_reference_g2_uses_explicit_added_accidentals():
    tau = np.linspace(0.0, 10.0, 401)
    wave = np.exp(-tau / 0.8)
    stats = observables.biphoton_stats(
        tau, wave, 1_000.0, signal_eff=0.1, idler_eff=0.1,
        coincidence_window_ns=1.0, target_g2_peak=44.0)
    assert np.isclose(stats["g2_peak"], 44.0)
    assert stats["raw_g2_peak"] > stats["g2_peak"]
    assert stats["added_accidental_cps"] > 0


def test_rb87_telecom_preset_smoke():
    scheme = fwm.FWMScheme()
    # Calibrated mode reproduces the injected reference waveform (FWHM ~0.56 ns).
    params = _recommended_params(topology=fwm.TOPOLOGY_RB87_TELECOM,
                                 biphoton_model=fwm.BIPHOTON_CALIBRATED)
    raw = scheme.compute(params)
    stats = _stats(raw, params)
    assert np.isfinite(stats["g2_peak"]) and stats["g2_peak"] > 2
    assert stats["pair_rate_cps"] > 0
    assert stats["fwhm_ns"] < 1.0
    # Predictive mode: waveform is solved (absolute width approximate for the
    # extreme 780/1529 nm ratio), so only finiteness/positivity is asserted here.
    pred = _recommended_params(topology=fwm.TOPOLOGY_RB87_TELECOM,
                               biphoton_model=fwm.BIPHOTON_PREDICTIVE)
    praw = scheme.compute(pred)
    pstats = _stats(praw, pred, target=False)
    assert np.all(np.isfinite(praw["psi_tau"]))
    assert 0.0 < pstats["fwhm_ns"] < 50.0
    assert praw["regime"] in ("group-delay", "damped-Rabi")


def test_cs_btw_channels_have_different_widths():
    scheme = fwm.FWMScheme()
    p917 = _recommended_params(topology=fwm.TOPOLOGY_CS_BTW,
                                cs_channel=fwm.CS_CHANNEL_917)
    p795 = _recommended_params(topology=fwm.TOPOLOGY_CS_BTW,
                                cs_channel=fwm.CS_CHANNEL_795)
    s917 = _stats(scheme.compute(p917), p917)
    s795 = _stats(scheme.compute(p795), p795)
    assert s917["pair_rate_cps"] > 0
    assert s795["pair_rate_cps"] > 0
    assert abs(s917["fwhm_ns"] - s795["fwhm_ns"]) > 0.05


def test_cs_btw_predictive_width_ordering():
    """Predictive: the 852-917 nm channel BTW is narrower than 852-795 nm — the
    wavelength-dependent collective two-photon-coherence ordering of Kim et al.
    (the absolute ns-widths are approximate; only the ordering is asserted)."""
    scheme = fwm.FWMScheme()
    p917 = _recommended_params(topology=fwm.TOPOLOGY_CS_BTW,
                               cs_channel=fwm.CS_CHANNEL_917,
                               biphoton_model=fwm.BIPHOTON_PREDICTIVE)
    p795 = _recommended_params(topology=fwm.TOPOLOGY_CS_BTW,
                               cs_channel=fwm.CS_CHANNEL_795,
                               biphoton_model=fwm.BIPHOTON_PREDICTIVE)
    s917 = _stats(scheme.compute(p917), p917, target=False)
    s795 = _stats(scheme.compute(p795), p795, target=False)
    assert s917["fwhm_ns"] < s795["fwhm_ns"]


def test_predictive_coupling_rabi_broadens_two_photon_resonance():
    """The Ω_c² Autler-Townes term lives in the two-photon denominator (vs the old
    weak-coupling drive in the numerator), so raising the coupling drive changes
    the source bandwidth — a coupling-power dependence the calibrated model lacks."""
    scheme = fwm.FWMScheme()
    base = _recommended_params(topology=fwm.TOPOLOGY_CS_BTW,
                               cs_channel=fwm.CS_CHANNEL_795,
                               biphoton_model=fwm.BIPHOTON_PREDICTIVE)
    weak = scheme.compute(dict(base, coupling_mw=0.25))
    strong = scheme.compute(dict(base, coupling_mw=8.0))
    assert np.isfinite(weak["source_bandwidth_mhz"])
    assert np.isfinite(strong["source_bandwidth_mhz"])
    assert weak["source_bandwidth_mhz"] != strong["source_bandwidth_mhz"]


def test_biphoton_ui_render_modes():
    scheme = fwm.FWMScheme()
    for topo in (fwm.TOPOLOGY_RB87_TELECOM, fwm.TOPOLOGY_CS_BTW,
                 fwm.TOPOLOGY_DIAMOND):
        params = _recommended_params(topology=topo)
        raw = scheme.compute(params)
        view = scheme.observables(raw, params)
        assert view.get("figure") is not None
        assert view.get("metrics")


def test_squeezing_hides_twin_beam_coincidence_figure():
    scheme = fwm.FWMScheme()
    params = scheme.defaults()
    params["resolution"] = "Fast  (~3 s)"
    raw = scheme.compute(params)
    view = scheme.observables(raw, params)
    assert view.get("figure") is not None
    assert not view.get("figures", [])
    assert all("Twin-beam coincidence" not in table["title"]
               for table in view.get("tables", []))


def test_squeezing_observables_tolerate_rapid_tpd_changes():
    import matplotlib.pyplot as plt

    scheme = fwm.FWMScheme()
    params = scheme.defaults()
    params["resolution"] = "Fast  (~3 s)"
    raw = scheme.compute(params)
    for tpd in (-480.0, -8.0, 0.0, 245.0, 500.0):
        params["tpd"] = tpd
        view = scheme.observables(raw, params)
        fig = view["figure"]
        fig.canvas.draw()
        labels = [ax.get_ylabel() for ax in fig.axes] + [fig.axes[-1].get_xlabel()]
        assert all("$" not in label for label in labels)
        plt.close(fig)


def test_seeded_phase_detail_modes_are_gated_by_resolution():
    center = fwm.branch_center_GHz(0.9, -1)
    common = dict(
        T=394.15, P_pump=0.6, P_probe=8e-6, line_strength=1.0,
        coarse_points=11, fine_points=0, scan_min=center - 0.02,
        scan_max=center + 0.02, velocity_step=20.0, velocity_cutoff=1.0,
        branch=-1,
    )
    legacy = fwm.compute_spectrum(0.9, **common)
    balanced = fwm.compute_spectrum(
        0.9, phase_detail=fwm.PHASE_BALANCED,
        pump_probe_angle_deg=fwm.SEEDED_PHASE_ANGLE_DEG, **common)
    fine = fwm.compute_spectrum(
        0.9, phase_detail=fwm.PHASE_FINE,
        pump_probe_angle_deg=fwm.SEEDED_PHASE_ANGLE_DEG, **common)

    assert legacy["delta_k_z"] is None
    assert balanced["delta_k_z"] is not None
    assert balanced["phase_segments"] == 1
    assert fine["delta_k_z"] is not None
    assert fine["phase_segments"] > 1


def test_fidelity_alias_and_ultra_tiny_grid():
    assert fwm.normalize_fidelity("Fine  (~20 s)") == fwm.FIDELITY_HIGH

    center = fwm.branch_center_GHz(0.9, -1)
    raw = fwm.compute_spectrum(
        0.9, T=394.15, P_pump=0.6, P_probe=8e-6, line_strength=1.0,
        coarse_points=7, fine_points=0, scan_min=center - 0.01,
        scan_max=center + 0.01, velocity_step=30.0, velocity_cutoff=0.4,
        phase_detail=fwm.PHASE_ULTRA, model_fidelity=fwm.FIDELITY_ULTRA,
        branch=-1)
    assert raw["delta_k_z"] is not None
    assert raw["phase_segments"] == fwm.ULTRA_PROPAGATION_SEGMENTS
    assert raw["ultra_phase_iterations"] == fwm.ULTRA_PHASE_ITERATIONS
    assert raw["ultra_dynamic_depletion"] is True
    assert raw["ultra_in_cell_loss_noise"] is True
    assert np.all(np.isfinite(raw["G_s"]))
    assert np.all(np.isfinite(raw["G_c"]))
    assert np.all(np.isfinite(raw["S_dB"]))
    assert np.nanmax(raw["G_s"]) <= raw["pump_depletion_cap"] * (1.0 + 1e-9)


def test_loss_noise_never_improves_squeezing():
    Gs = np.array([2.0, 10.0])
    Gc = np.array([1.0, 9.0])
    ideal = observables.intensity_difference_squeezing_dB(Gs, Gc, 0.9)
    lossy = observables.segmented_loss_noise_squeezing_dB(
        Gs, Gc, 0.9, in_cell_loss_frac=0.1)
    assert np.all(lossy >= ideal)


def test_rb85_fwm_zeeman_cg_sum_rules_match_lumped_strengths():
    atom = zeeman.rb85_d1_double_lambda_zeeman()
    assert atom.n_levels == 24
    for key, cf2 in hyperfine.CF2.items():
        assert np.isclose(atom.lumped_strengths[key], 3.0 * cf2, rtol=1e-12)
    assert np.isclose(atom.lumped_strength_correction, 1.0, rtol=1e-12)


def test_biphoton_fine_phase_adds_absolute_and_2d_map():
    scheme = fwm.FWMScheme()
    params = _recommended_params(
        topology=fwm.TOPOLOGY_RB87_TELECOM,
        phase_detail="Fine",
        biphoton_velocity_step=20.0,
    )
    raw = scheme.compute(params)
    assert np.isfinite(raw["delta_k_absolute"])
    assert np.isfinite(raw["phase_match_weight_absolute"])
    assert raw["phase_matching_2d"] is not None
    assert raw["phase_matching_2d"].ndim == 2

    idx = np.unravel_index(np.argmax(raw["phase_matching_2d"]),
                           raw["phase_matching_2d"].shape)
    best_signal = raw["signal_angle_axis_deg"][idx[0]]
    best_idler = raw["idler_angle_axis_2d_deg"][idx[1]]
    assert abs(best_signal - params["signal_angle_deg"]) <= 0.06
    assert abs(best_idler - params["idler_angle_deg"]) <= 0.06


def test_fwm_default_buttons_are_squeezing_and_contextual_biphoton():
    scheme = fwm.FWMScheme()
    defaults = scheme.recommended_defaults(scheme.defaults())
    assert set(defaults) == {fwm.MODE_SEEDED, fwm.MODE_BIPHOTON}
    assert defaults[fwm.MODE_SEEDED]["mode"] == fwm.MODE_SEEDED
    schema = {spec.name: spec for spec in scheme.param_schema()}
    assert schema["mode"].applies_defaults
    assert schema["topology"].applies_defaults
    assert schema["cs_channel"].applies_defaults

    cs_defaults = scheme.recommended_defaults(_params(
        topology=fwm.TOPOLOGY_CS_BTW,
        cs_channel=fwm.CS_CHANNEL_795,
    ))[fwm.MODE_BIPHOTON]
    assert cs_defaults["mode"] == fwm.MODE_BIPHOTON
    assert cs_defaults["topology"] == fwm.TOPOLOGY_CS_BTW
    assert cs_defaults["cs_channel"] == fwm.CS_CHANNEL_795
    assert cs_defaults["biphoton_temp_c"] == 75.0
    cs917_defaults = scheme.recommended_defaults(_params(
        topology=fwm.TOPOLOGY_CS_BTW,
        cs_channel=fwm.CS_CHANNEL_917,
    ))[fwm.MODE_BIPHOTON]
    assert abs(cs917_defaults["idler_angle_deg"] - 1.39) < 0.03
    assert abs(cs_defaults["idler_angle_deg"] - 1.61) < 0.03
    assert cs_defaults["signal_side"] == fwm.SIDE_PLUS
    assert cs_defaults["idler_side"] == fwm.SIDE_PLUS

    diamond_defaults = scheme.recommended_defaults(_params(
        topology=fwm.TOPOLOGY_DIAMOND,
    ))[fwm.MODE_BIPHOTON]
    assert diamond_defaults["topology"] == fwm.TOPOLOGY_DIAMOND
    assert diamond_defaults["diamond_idler_nm"] == 761.702
    assert abs(diamond_defaults["idler_angle_deg"] - 1.92) < 0.03
    assert diamond_defaults["idler_side"] == fwm.SIDE_MINUS


def test_cs_btw_short_window_render_no_shape_error():
    scheme = fwm.FWMScheme()
    params = _recommended_params(
        topology=fwm.TOPOLOGY_CS_BTW,
        cs_channel=fwm.CS_CHANNEL_917,
        tau_max_ns=1.0,
        timing_jitter_ns=0.55,
    )
    raw = scheme.compute(params)
    view = scheme.observables(raw, params)
    assert view.get("figure") is not None
    assert any(table["title"].startswith("Reference ")
               for table in view.get("tables", []))


if __name__ == "__main__":
    test_topology_energy_and_roles()
    test_phase_matching_reference_angle_is_maximum()
    test_rb87_default_vector_phase_match_has_positive_rate()
    test_rb87_equal_angles_are_transversely_suppressed()
    test_side_flip_suppresses_matched_geometry()
    test_detector_background_and_window_reduce_car()
    test_timing_jitter_broadens_waveform()
    test_long_jitter_kernel_preserves_axis_length()
    test_reference_g2_uses_explicit_added_accidentals()
    test_rb87_telecom_preset_smoke()
    test_cs_btw_channels_have_different_widths()
    test_biphoton_ui_render_modes()
    test_squeezing_hides_twin_beam_coincidence_figure()
    test_seeded_phase_detail_modes_are_gated_by_resolution()
    test_fidelity_alias_and_ultra_tiny_grid()
    test_loss_noise_never_improves_squeezing()
    test_rb85_fwm_zeeman_cg_sum_rules_match_lumped_strengths()
    test_biphoton_fine_phase_adds_absolute_and_2d_map()
    test_fwm_default_buttons_are_squeezing_and_contextual_biphoton()
    test_cs_btw_short_window_render_no_shape_error()
    print("Generic SFWM / biphoton checks OK.")
