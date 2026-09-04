"""
Generic SFWM / biphoton checks.

    python tests/test_fwm_generic.py   # or: pytest tests/test_fwm_generic.py
"""
import sys
from pathlib import Path
from unittest.mock import patch

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

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
    bad = dict(matched, signal_angle_deg=1.5,
               idler_angle_offset_deg=1.5 - matched["idler_angle_deg"],
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
    assert np.isclose(sharp["source_fwhm_ns"], sharp["detected_fwhm_ns"])
    assert np.isclose(sharp["fwhm_ns"], sharp["detected_fwhm_ns"])
    assert np.isclose(broad["source_fwhm_ns"], sharp["source_fwhm_ns"])
    assert broad["detected_fwhm_ns"] > broad["source_fwhm_ns"]
    assert broad["fwhm_ns"] == broad["detected_fwhm_ns"]
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
    # Predictive mode: the waveform is solved, but the absolute per-source width
    # remains approximate, so only finiteness/positivity is asserted here.
    pred = _recommended_params(topology=fwm.TOPOLOGY_RB87_TELECOM,
                               biphoton_model=fwm.BIPHOTON_PREDICTIVE)
    praw = scheme.compute(pred)
    pstats = _stats(praw, pred, target=False)
    assert np.all(np.isfinite(praw["psi_tau"]))
    assert 0.0 < pstats["fwhm_ns"] < 50.0
    assert pstats["fwhm_ns"] == pstats["detected_fwhm_ns"]
    assert pstats["source_fwhm_ns"] < pstats["detected_fwhm_ns"]
    assert praw["regime"] in ("group-delay", "damped-Rabi")
    pview = scheme.headless_observables(praw, pred)
    rate_metric = next(m for m in pview["metrics"] if m["label"] == "Pair rate")
    assert "Reference-anchored" in rate_metric["help"]
    metric_labels = {metric["label"] for metric in pview["metrics"]}
    assert {"Intrinsic biphoton width", "Detected biphoton width"} <= metric_labels
    unconverged_view = scheme.headless_observables(
        dict(praw, velocity_converged=False), pred)
    unconverged_metrics = {
        metric["label"]: metric["value"]
        for metric in unconverged_view["metrics"]
    }
    assert unconverged_metrics["Intrinsic biphoton width"] == "unconverged"
    assert unconverged_metrics["Detected biphoton width"] == "unconverged"
    unconverged_detection = next(
        table["markdown"] for table in unconverged_view["tables"]
        if table["title"] == "Detection estimates")
    assert "| Intrinsic biphoton width | unconverged |" in unconverged_detection
    assert "| Detected biphoton width | unconverged |" in unconverged_detection


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
    assert s917["source_fwhm_ns"] < s795["source_fwhm_ns"]
    assert s917["detected_fwhm_ns"] < s795["detected_fwhm_ns"]


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


def test_fwm_headless_observables_skip_figure_generation():
    import matplotlib.pyplot as plt

    scheme = fwm.FWMScheme()
    assert scheme.supports_headless_observables is True

    seeded_params = scheme.defaults()
    seeded_params["resolution"] = fwm.FIDELITY_FAST
    seeded_raw = scheme.compute(seeded_params)
    biphoton_params = _recommended_params(
        topology=fwm.TOPOLOGY_RB87_TELECOM,
        biphoton_model=fwm.BIPHOTON_CALIBRATED,
        phase_detail="Fine",
    )
    biphoton_raw = scheme.compute(biphoton_params)
    assert biphoton_raw["phase_matching_2d"] is None

    original_subplots = plt.subplots

    def _fail_subplots(*_args, **_kwargs):
        raise AssertionError("headless observables must not build figures")

    try:
        plt.subplots = _fail_subplots
        for raw, params in (
            (seeded_raw, seeded_params),
            (biphoton_raw, biphoton_params),
        ):
            view = scheme.headless_observables(raw, params)
            direct = scheme.observables(raw, params, include_figures=False)
            assert view.get("figure") is None
            assert not view.get("figures", [])
            assert direct.get("figure") is None
            assert not direct.get("figures", [])
            assert view.get("metrics")
            assert view.get("tables")
            assert [m["label"] for m in view["metrics"]] == [
                m["label"] for m in direct["metrics"]
            ]
            assert [t["title"] for t in view["tables"]] == [
                t["title"] for t in direct["tables"]
            ]
    finally:
        plt.subplots = original_subplots


def test_biphoton_basic_ui_keeps_only_lab_controls():
    scheme = fwm.FWMScheme()
    specs = scheme.param_schema()
    biphoton_only = {"mode": fwm.MODE_BIPHOTON}
    basic = {
        sp.name for sp in specs
        if sp.visible_if == biphoton_only and not sp.advanced and not sp.hidden
    }
    assert basic == {
        "topology",
        "biphoton_temp_c",
        "pump_biphoton_uw",
        "coupling_mw",
        "pump_detuning_mhz",
        "two_photon_detuning_mhz",
        "signal_angle_deg",
        "idler_angle_offset_deg",
        "coincidence_window_ns",
        "filter_bandwidth_mhz",
    }
    hidden = {sp.name for sp in specs if sp.hidden}
    assert {"coupling_detuning_mhz", "idler_angle_deg",
            "signal_side", "idler_side"} <= hidden

    params = _recommended_params(
        topology=fwm.TOPOLOGY_RB87_TELECOM,
        pump_detuning_mhz=120.0,
        two_photon_detuning_mhz=25.0,
        signal_angle_deg=1.8,
    )
    runtime = scheme._biphoton_runtime_params(params)
    expected_idler = fwm.transverse_matched_angle_deg(1529.37, 780.24, 1.8)
    assert runtime["coupling_detuning_mhz"] == -95.0
    assert np.isclose(runtime["idler_angle_deg"], expected_idler)


def test_fwm_ui_uses_short_labels_without_changing_stored_values():
    scheme = fwm.FWMScheme()
    specs = {sp.name: sp for sp in scheme.param_schema()}

    assert scheme.title == "Four-wave mixing (Squeezing / Biphoton)"
    assert "mean-field Squeezing indicator" in scheme.caption
    assert specs["mode"].choices == (fwm.MODE_SEEDED, fwm.MODE_BIPHOTON)
    assert specs["mode"].choice_labels[fwm.MODE_SEEDED] == "Squeezing"
    assert specs["topology"].choice_labels == fwm.TOPOLOGY_LABELS
    assert specs["biphoton_model"].choice_labels == fwm.BIPHOTON_MODEL_LABELS
    assert specs["resolution"].choice_labels == fwm.FIDELITY_LABELS
    assert specs["cs_channel"].choice_labels == fwm.CS_CHANNEL_LABELS

    assert specs["probe_uw"].label == "Seed power"
    assert specs["opd"].label == "One-photon detuning Δ"
    assert specs["tpd"].label == "Two-photon detuning δ"
    assert specs["detection_eff_pct"].label == "Detection efficiency η"
    assert specs["phase_detail"].choice_labels == {
        "Balanced": "1D", "Fine": "1D + 2D"
    }

    hidden = {
        "eom_residual_carrier_uw", "eom_other_sidebands_uw",
        "loss_pct", "qe_pct", "line_strength", "mode_overlap_penalty",
        "polarization_penalty", "zeeman_participation_penalty",
        "floquet_order", "phase_detail",
    }
    assert all(specs[name].hidden for name in hidden)


def test_squeezing_hides_twin_beam_coincidence_figure():
    scheme = fwm.FWMScheme()
    params = scheme.defaults()
    params["resolution"] = fwm.FIDELITY_FAST
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
    params["resolution"] = fwm.FIDELITY_FAST
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
    for result in (balanced, fine):
        np.testing.assert_allclose(result["delta_k_z"],
                                   result["delta_k_z_vacuum"])
        _, expected_probe, expected_conj = fwm.seeded_option_a_wavenumbers(
            0.9, result["probe_axis_GHz"])
        np.testing.assert_allclose(result["k_probe_propagation_per_m"],
                                   expected_probe)
        np.testing.assert_allclose(result["k_conjugate_propagation_per_m"],
                                   expected_conj)


def test_fidelity_alias_and_ultra_tiny_grid():
    assert fwm.normalize_fidelity("Fine  (~20 s)") == fwm.FIDELITY_BALANCED

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
    np.testing.assert_allclose(raw["delta_k_z"], raw["delta_k_z_vacuum"])
    assert raw["propagation_convention"].startswith("Option A")
    assert raw["ultra_dynamic_depletion"] is True
    assert raw["ultra_in_cell_loss_noise"] is False
    assert raw["squeezing_status"].startswith("unavailable: gain-referred diagnostic")
    assert raw["claim_gate"]["level"] == "MEAN_FIELD_DIAGNOSTIC"
    assert not raw["claim_gate"]["quantitative_gain_supported"]
    assert not raw["claim_gate"]["physical_squeezing_prediction"]
    assert np.all(np.isfinite(raw["G_s"]))
    assert np.all(np.isfinite(raw["G_c"]))
    assert np.all(np.isfinite(raw["S_dB"]))
    assert raw["T_field_small_signal"].shape[-2:] == (2, 2)
    assert raw["T_canonical_small_signal"].shape == raw["T_field_small_signal"].shape
    assert raw["Q_photon_flux"].shape == raw["T_field_small_signal"].shape
    assert raw["canonical_mode_status"].startswith("conditional")
    assert np.all(np.isfinite(raw["photon_flux_gap_smallsignal"]))
    assert np.all(np.isfinite(raw["commutator_defect_max_smallsignal"]))
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


def test_biphoton_fine_phase_map_is_lazy_and_figure_only():
    scheme = fwm.FWMScheme()
    params = _recommended_params(
        topology=fwm.TOPOLOGY_RB87_TELECOM,
        phase_detail="Fine",
        biphoton_velocity_step=20.0,
    )
    raw = scheme.compute(params)
    assert np.isfinite(raw["delta_k_absolute"])
    assert np.isfinite(raw["phase_match_weight_absolute"])
    assert raw["phase_matching_2d"] is None
    phase_spec = {sp.name: sp for sp in scheme.param_schema()}["phase_detail"]
    assert phase_spec.recompute is False
    assert phase_spec.hidden is True
    assert phase_spec.choice_labels == {"Balanced": "1D", "Fine": "1D + 2D"}
    assert "phase_detail" not in scheme.recompute_keys()

    balanced_raw = scheme.compute(dict(params, phase_detail="Balanced"))
    for key in ("tau_axis_ns", "psi_tau", "v_grid", "velocity_weights", "source_v"):
        assert np.array_equal(raw[key], balanced_raw[key])
    for key in ("pair_rate_cps", "phase_match_weight", "delta_k_vector"):
        assert raw[key] == balanced_raw[key]

    def unexpected_map(*args, **kwargs):
        raise AssertionError("headless observables must not build the 2-D phase map")

    with patch.object(fwm, "biphoton_phase_matching_map", unexpected_map):
        headless = scheme.headless_observables(raw, params)
    assert headless["figures"] == []
    assert headless["figure_controls"] == ["phase_detail"]

    signal_axis, idler_axis, phase_map = fwm.biphoton_phase_matching_map(
        raw["fields"], raw["cell_length_m"],
        signal_angle_deg=params["signal_angle_deg"],
        idler_angle_deg=params["idler_angle_deg"],
        reference_delta_k=raw["topology"].reference_delta_k)
    assert phase_map.shape == (121, 121)
    inline_map = fwm.phase_matching_weight(
        fwm.phase_mismatch_vector_grid(
            raw["fields"], signal_axis, idler_axis,
            reference_delta_k=raw["topology"].reference_delta_k),
        raw["cell_length_m"])
    assert np.array_equal(phase_map, inline_map)

    idx = np.unravel_index(np.argmax(phase_map), phase_map.shape)
    best_signal = signal_axis[idx[0]]
    best_idler = idler_axis[idx[1]]
    assert abs(best_signal - params["signal_angle_deg"]) <= 0.06
    assert abs(best_idler - params["idler_angle_deg"]) <= 0.06

    view = scheme.observables(raw, params)
    assert "2D phase matching" in {title for title, _ in view["figures"]}

    import matplotlib.pyplot as plt
    plt.close(view["figure"])
    for _, fig in view["figures"]:
        plt.close(fig)


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


def test_fwm_info_keeps_zeeman_diagnostic_scope_explicit():
    scheme = fwm.FWMScheme()
    info = scheme.info()
    assert "Squeezing indicator" in info
    assert "not a physical squeezing spectrum" in info
    assert "Reduced model" in info
    assert "pair-rate scale is anchored" in info
    assert "first principles" not in info
    line_strength = next(
        sp for sp in scheme.param_schema() if sp.name == "line_strength")
    assert line_strength.hidden is True
    assert "has not been refitted" in line_strength.help
    assert "does not anchor measured gain or squeezing" in line_strength.help
    assert "no measurements identify a unique split" in line_strength.help

    factors = {
        sp.name: sp for sp in scheme.param_schema()
        if sp.name in {
            "mode_overlap_penalty", "polarization_penalty",
            "zeeman_participation_penalty"
        }
    }
    assert set(factors) == {
        "mode_overlap_penalty", "polarization_penalty",
        "zeeman_participation_penalty"
    }
    assert all(sp.default == 1.0 for sp in factors.values())
    assert all(sp.vmin == 0.0 and sp.vmax == 1.0 for sp in factors.values())
    assert all(sp.hidden for sp in factors.values())
    assert "count that geometry twice" in factors[
        "mode_overlap_penalty"].help
    assert "not a raw Stokes" in factors["polarization_penalty"].help
    assert "does not turn on" in factors[
        "zeeman_participation_penalty"].help


def test_seeded_coupling_penalties_route_to_both_spectrum_paths():
    scheme = fwm.FWMScheme()
    params = scheme.defaults()
    params.update(
        mode_overlap_penalty=0.9,
        polarization_penalty=0.75,
        zeeman_participation_penalty=0.6,
    )
    marker = {"marker": True}

    with patch.object(fwm, "compute_spectrum", return_value=marker) as solve:
        assert scheme.compute(params) is marker
    kwargs = solve.call_args.kwargs
    assert kwargs["mode_overlap_penalty"] == 0.9
    assert kwargs["polarization_penalty"] == 0.75
    assert kwargs["zeeman_participation_penalty"] == 0.6

    full_view = scheme.extra_views()[0]
    with patch.object(fwm, "full_spectrum", return_value=marker) as full:
        assert full_view.compute(params) is marker
    kwargs = full.call_args.kwargs
    assert kwargs["mode_overlap_penalty"] == 0.9
    assert kwargs["polarization_penalty"] == 0.75
    assert kwargs["zeeman_participation_penalty"] == 0.6


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
    metrics = {metric["label"]: metric for metric in view["metrics"]}
    assert float(metrics["Detected biphoton width"]["value"].split()[0]) > float(
        metrics["Intrinsic biphoton width"]["value"].split()[0])
    detection = next(table["markdown"] for table in view["tables"]
                     if table["title"] == "Detection estimates")
    assert "Intrinsic biphoton width" in detection
    assert "Detected biphoton width" in detection
    assert "Net timing-difference response FWHM" in detection
    assert any(table["title"] == "Literature comparison"
               for table in view.get("tables", []))


def test_beam_geometry_knobs_default_to_the_legacy_constants():
    """Promoting the fixed geometry to knobs must not move the anchored point."""
    scheme = fwm.FWMScheme()
    params = scheme.defaults()
    assert params["pump_waist_um"] == fwm.W_PUMP * 1e6
    assert params["probe_waist_um"] == fwm.W_PROBE * 1e6
    assert params["qe_pct"] == fwm.QE_DETECTOR * 100.0
    assert params["detection_eff_pct"] == pytest.approx(
        fwm.SEEDED_DETECTION_EFFICIENCY_PCT)

    raw = scheme.compute(params)
    center = fwm.branch_center_GHz(params["opd"], -1)
    res = fwm.FWM_FIDELITY[fwm.normalize_fidelity(params["resolution"])]
    legacy = fwm.compute_spectrum(
        params["opd"], T=params["temp_c"] + 273.15,
        P_pump=params["pump_mw"] * 1e-3, P_probe=params["probe_uw"] * 1e-6,
        line_strength=params["line_strength"],
        L=params["cell_mm"] * 1e-3, loss_frac=params["loss_pct"] / 100.0,
        coarse_points=res["coarse_points"], fine_points=0,
        scan_min=center - fwm.WINDOW_GHZ, scan_max=center + fwm.WINDOW_GHZ,
        velocity_step=res["velocity_step"],
        velocity_cutoff=res.get("velocity_cutoff", 3.0),
        phase_detail=res["phase_detail"],
        pump_probe_angle_deg=params["seeded_angle_deg"],
        model_fidelity=fwm.normalize_fidelity(params["resolution"]),
        branch=-1)
    for key in ("probe_axis_GHz", "G_s", "G_c", "S_dB"):
        assert np.array_equal(np.asarray(raw[key]), np.asarray(legacy[key])), \
            f"{key} drifted when the geometry became a knob"
    assert raw["w_pump_m"] == fwm.W_PUMP
    assert raw["w_probe_m"] == fwm.W_PROBE
    assert raw["qe"] == fwm.QE_DETECTOR


def test_beam_geometry_knobs_route_to_both_spectrum_paths():
    scheme = fwm.FWMScheme()
    params = scheme.defaults()
    params.update(pump_waist_um=800.0, probe_waist_um=450.0, qe_pct=90.45)
    marker = {"marker": True}

    with patch.object(fwm, "compute_spectrum", return_value=marker) as solve:
        assert scheme.compute(params) is marker
    kwargs = solve.call_args.kwargs
    assert kwargs["w_pump"] == pytest.approx(800e-6)
    assert kwargs["w_probe"] == pytest.approx(450e-6)
    assert kwargs["qe"] == pytest.approx(0.9045)
    assert kwargs["detection_efficiency"] == pytest.approx(
        0.9045 * (1.0 - params["loss_pct"] / 100.0))

    full_view = scheme.extra_views()[0]
    with patch.object(fwm, "full_spectrum", return_value=marker) as full:
        assert full_view.compute(params) is marker
    kwargs = full.call_args.kwargs
    assert kwargs["w_pump"] == pytest.approx(800e-6)
    assert kwargs["w_probe"] == pytest.approx(450e-6)
    assert kwargs["qe"] == pytest.approx(0.9045)
    assert kwargs["detection_efficiency"] == pytest.approx(
        0.9045 * (1.0 - params["loss_pct"] / 100.0))


def test_detection_efficiency_routes_to_focused_and_full_spectrum():
    scheme = fwm.FWMScheme()
    params = scheme.defaults()
    params["detection_eff_pct"] = 50.0
    marker = {"marker": True}

    with patch.object(fwm, "compute_spectrum", return_value=marker) as solve:
        assert scheme.compute(params) is marker
    assert solve.call_args.kwargs["detection_efficiency"] == pytest.approx(0.5)

    with patch.object(fwm, "full_spectrum", return_value=marker) as full:
        assert scheme.extra_views()[0].compute(params) is marker
    assert full.call_args.kwargs["detection_efficiency"] == pytest.approx(0.5)


def test_full_scan_plot_uses_display_solver_detail():
    axis = np.array([-2.2, -2.1])
    spectrum = {
        "probe_axis_GHz": axis,
        "G_s": np.array([1.0, 1.1]),
        "gain_referred_noise_dB": np.array([0.0, -0.1]),
        "model_fidelity": fwm.FIDELITY_FAST,
        "propagation_convention": "internal propagation token",
    }
    full = {"D_GHz": 0.9, "minus": spectrum, "plus": dict(spectrum)}
    fig = fwm.FWMScheme().extra_views()[0].render(full)
    assert fig.axes[0].get_title() == "Solver detail: Fast"
    assert "internal propagation token" not in fig.axes[0].get_title()

    import matplotlib.pyplot as plt
    plt.close(fig)


def test_qe_knob_sets_detection_eta_and_the_squeezing_floor():
    """eta = QE·(1−loss) is the hard ceiling 10·log10(1−eta) on any readout."""
    scheme = fwm.FWMScheme()
    params = scheme.defaults()
    params.update(qe_pct=90.45, loss_pct=5.5, resolution=fwm.FIDELITY_FAST)
    raw = scheme.compute(params)
    eta = 0.9045 * (1.0 - 0.055)
    assert raw["qe"] == pytest.approx(0.9045)
    assert raw["eta"] == pytest.approx(eta)
    floor_dB = 10.0 * np.log10(1.0 - eta)
    assert np.nanmin(np.asarray(raw["S_dB"])) >= floor_dB - 1e-9


def test_pump_waist_moves_the_gain():
    """A geometry knob that cannot change the solve would be cosmetic."""
    scheme = fwm.FWMScheme()
    base = scheme.defaults()
    wide = dict(base, pump_waist_um=900.0)
    g_base = np.nanmax(np.asarray(scheme.compute(base)["G_s"]))
    g_wide = np.nanmax(np.asarray(scheme.compute(wide)["G_s"]))
    assert not np.isclose(g_base, g_wide, rtol=1e-6)


def test_responsivity_tracks_quantum_efficiency():
    # Display-only conversion R = QE·λ/hc; the legacy caption showed 0.59 A/W.
    assert abs(fwm.responsivity_AW(0.92) - 0.59) < 5e-3
    assert abs(fwm.responsivity_AW(0.9045) - 0.58) < 5e-3


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
    test_fwm_headless_observables_skip_figure_generation()
    test_squeezing_hides_twin_beam_coincidence_figure()
    test_seeded_phase_detail_modes_are_gated_by_resolution()
    test_fidelity_alias_and_ultra_tiny_grid()
    test_loss_noise_never_improves_squeezing()
    test_rb85_fwm_zeeman_cg_sum_rules_match_lumped_strengths()
    test_biphoton_fine_phase_map_is_lazy_and_figure_only()
    test_fwm_default_buttons_are_squeezing_and_contextual_biphoton()
    test_cs_btw_short_window_render_no_shape_error()
    test_beam_geometry_knobs_default_to_the_legacy_constants()
    test_beam_geometry_knobs_route_to_both_spectrum_paths()
    test_qe_knob_sets_detection_eta_and_the_squeezing_floor()
    test_pump_waist_moves_the_gain()
    test_responsivity_tracks_quantum_efficiency()
    print("Generic SFWM / biphoton checks OK.")
