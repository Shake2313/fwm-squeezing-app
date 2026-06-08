"""
Generic SFWM / biphoton checks.

    python tests/test_fwm_generic.py   # or: pytest tests/test_fwm_generic.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import observables  # noqa: E402
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
    params = _params(topology=fwm.TOPOLOGY_RB87_TELECOM)
    spec = fwm.topology_from_params(params)
    L = spec.default_cell_mm * 1e-3
    angles = np.linspace(0.0, 4.0, 81)
    weights = np.array([
        fwm.phase_matching_weight(
            np.array([fwm.phase_mismatch(
                spec.fields, idler_angle_deg=a,
                reference_delta_k=spec.reference_delta_k)]), L)[0]
        for a in angles
    ])
    best = angles[int(np.argmax(weights))]
    assert abs(best - 1.5) <= 0.1
    assert np.isclose(weights.max(), 1.0)


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
    params = _params(topology=fwm.TOPOLOGY_RB87_TELECOM)
    raw = scheme.compute(params)
    stats = _stats(raw, params)
    assert np.isfinite(stats["g2_peak"]) and stats["g2_peak"] > 2
    assert stats["pair_rate_cps"] > 0
    assert stats["fwhm_ns"] < 1.0


def test_cs_btw_channels_have_different_widths():
    scheme = fwm.FWMScheme()
    p917 = _params(topology=fwm.TOPOLOGY_CS_BTW, cs_channel=fwm.CS_CHANNEL_917)
    p795 = _params(topology=fwm.TOPOLOGY_CS_BTW, cs_channel=fwm.CS_CHANNEL_795)
    s917 = _stats(scheme.compute(p917), p917)
    s795 = _stats(scheme.compute(p795), p795)
    assert abs(s917["fwhm_ns"] - s795["fwhm_ns"]) > 0.05


def test_biphoton_ui_render_modes():
    scheme = fwm.FWMScheme()
    for topo in (fwm.TOPOLOGY_RB87_TELECOM, fwm.TOPOLOGY_CS_BTW,
                 fwm.TOPOLOGY_DIAMOND):
        params = _params(topology=topo)
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


def test_fwm_default_buttons_are_squeezing_and_contextual_biphoton():
    scheme = fwm.FWMScheme()
    defaults = scheme.recommended_defaults(scheme.defaults())
    assert set(defaults) == {fwm.MODE_SEEDED, fwm.MODE_BIPHOTON}
    assert defaults[fwm.MODE_SEEDED]["mode"] == fwm.MODE_SEEDED

    cs_defaults = scheme.recommended_defaults(_params(
        topology=fwm.TOPOLOGY_CS_BTW,
        cs_channel=fwm.CS_CHANNEL_795,
    ))[fwm.MODE_BIPHOTON]
    assert cs_defaults["mode"] == fwm.MODE_BIPHOTON
    assert cs_defaults["topology"] == fwm.TOPOLOGY_CS_BTW
    assert cs_defaults["cs_channel"] == fwm.CS_CHANNEL_795
    assert cs_defaults["biphoton_temp_c"] == 75.0

    diamond_defaults = scheme.recommended_defaults(_params(
        topology=fwm.TOPOLOGY_DIAMOND,
    ))[fwm.MODE_BIPHOTON]
    assert diamond_defaults["topology"] == fwm.TOPOLOGY_DIAMOND
    assert diamond_defaults["diamond_idler_nm"] == 761.702


def test_cs_btw_short_window_render_no_shape_error():
    scheme = fwm.FWMScheme()
    params = _params(
        topology=fwm.TOPOLOGY_CS_BTW,
        cs_channel=fwm.CS_CHANNEL_917,
        tau_max_ns=1.0,
        timing_jitter_ns=0.55,
    )
    raw = scheme.compute(params)
    view = scheme.observables(raw, params)
    assert view.get("figure") is not None
    assert any(table["title"] == "Reference validation (medium model)"
               for table in view.get("tables", []))


if __name__ == "__main__":
    test_topology_energy_and_roles()
    test_phase_matching_reference_angle_is_maximum()
    test_detector_background_and_window_reduce_car()
    test_timing_jitter_broadens_waveform()
    test_long_jitter_kernel_preserves_axis_length()
    test_reference_g2_uses_explicit_added_accidentals()
    test_rb87_telecom_preset_smoke()
    test_cs_btw_channels_have_different_widths()
    test_biphoton_ui_render_modes()
    test_squeezing_hides_twin_beam_coincidence_figure()
    test_fwm_default_buttons_are_squeezing_and_contextual_biphoton()
    test_cs_btw_short_window_render_no_shape_error()
    print("Generic SFWM / biphoton checks OK.")
