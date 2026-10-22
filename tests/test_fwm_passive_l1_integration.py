"""Current Jasperse NF routing over the preserved passive-L1 transfer API."""
from copy import deepcopy

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from gabes import observables
from gabes.passive_readout import passive_l1_from_spectrum
from gabes.schemes import fwm

@pytest.fixture
def passive_spectrum():
    center = fwm.branch_center_GHz(0.9, -1)
    amplifier = np.array([[2.0, np.sqrt(3.0)], [np.sqrt(3.0), 2.0]])
    maps = np.array([
        np.diag([np.sqrt(0.9), 1.0]) @ amplifier,
        np.diag([np.sqrt(1.1), 1.0]) @ amplifier,
        amplifier,
    ], dtype=complex)
    # Displayed powers deliberately differ from the transfer-map flux gains.
    # The passive readout must not reconstruct a covariance from these values.
    raw = {
        "D_GHz": 0.9,
        "probe_axis_GHz": center + np.array([-8.0, 0.0, 8.0]) * 1e-3,
        "raman_center_minus_GHz": center,
        "raman_center_plus_GHz": fwm.branch_center_GHz(0.9, 1),
        "T_canonical_small_signal": maps,
        "input_seed_power_W": 8e-6,
        "G_s": np.array([1.5, 8.0, 2.5]),
        "G_c": np.array([0.3, 7.0, 1.4]),
        "gain_referred_noise_dB": np.array([-0.4, -4.7, -6.0]),
        "S_dB": np.array([-0.4, -4.7, -6.0]),
        "physical_squeezing_dB": None,
        "eta": 0.8,
        "phase_detail": fwm.PHASE_LEGACY,
        "model_fidelity": fwm.FIDELITY_FAST,
        "N_atoms": 1e18,
        "sigma_v": 200.0,
        "n_velocity": 3,
        "Op_A_2pi_GHz": 0.5,
        "Os_2pi_MHz": 0.1,
        "claim_gate": {"physical_squeezing_prediction": False},
        "squeezing_status": "unavailable: microscopic atomic diffusion absent",
    }
    raw.update(passive_l1_from_spectrum(raw))
    return raw


def _current_reference(raw, params):
    _, k, _ = fwm.seeded_option_a_wavenumbers(raw["D_GHz"], raw["probe_axis_GHz"])
    b = 2*fwm.constants.HBAR*fwm.constants.C_LIGHT*k/raw["input_seed_power_W"]
    fano = 1.0 + params.get("seed_classical_rin_scaled", 0.0)*1e-14/b
    eta = params["after_cell_eff_pct"]*params["qe_pct"]/1e4
    return passive_l1_from_spectrum(dict(raw, eta=eta), seed_fano=fano)


def test_focused_headless_current_nf_mask_nonpassive_points_without_cache_mutation(passive_spectrum):
    import matplotlib.pyplot as plt
    scheme = fwm.FWMScheme()
    params = dict(scheme.defaults(), noise_model=fwm.NOISE_PASSIVE_L1, tpd=-8.0, seed_classical_rin_scaled=3.0)
    original = deepcopy(passive_spectrum)
    expected = _current_reference(passive_spectrum, params)
    rendered = scheme.observables(passive_spectrum, params)
    try:
        curve = rendered["figure"].axes[1].lines[0].get_ydata()
        line = rendered["figure"].axes[1].lines[0]
        delta = (passive_spectrum["probe_axis_GHz"]-passive_spectrum["raman_center_minus_GHz"])*1e3
        np.testing.assert_allclose(np.interp(delta, line.get_xdata(), curve),
                                   expected["NF_dB"], equal_nan=True)
        assert np.isnan(np.interp(0.0, line.get_xdata(), curve))
        assert rendered["metrics"][0]["label"] == "Squeezing Degree [dB]"
        assert rendered["metrics"][0]["value"] == f"{expected['NF_dB'][0]:.2f} dB"
        current = rendered["jasperse_readout"]
        np.testing.assert_allclose(current["seed_fano_input"],
                                   expected["seed_fano_input"], rtol=1e-14)
        headless = scheme.headless_observables(passive_spectrum, params)
        assert headless["metrics"] == rendered["metrics"]
        np.testing.assert_allclose(headless["jasperse_readout"]["NF_dB"],
                                   expected["NF_dB"], equal_nan=True)
        for tpd in (0.0, -4.0):
            missing = scheme.headless_observables(passive_spectrum, dict(params, tpd=tpd))
            assert missing["metrics"][0]["value"] == "unavailable"
    finally:
        plt.close(rendered["figure"])
    for key in ("G_s", "G_c", "gain_referred_noise_dB", "S_dB",
                "T_canonical_small_signal", "passive_l1_noise_dB"):
        np.testing.assert_array_equal(passive_spectrum[key], original[key])
    assert passive_spectrum["physical_squeezing_dB"] is None
    assert passive_spectrum["claim_gate"] == original["claim_gate"]


@pytest.mark.parametrize("classical_rin", (0.0, 3.0))
def test_full_scan_uses_current_seed_rin_and_shared_external_efficiency(passive_spectrum, classical_rin):
    import matplotlib.pyplot as plt
    params = dict(fwm.FWMScheme().defaults(), noise_model=fwm.NOISE_PASSIVE_L1, seed_classical_rin_scaled=classical_rin,
                  after_cell_eff_pct=70.0, qe_pct=80.0)
    full = {"D_GHz": 0.9, "minus": passive_spectrum,
            "plus": dict(passive_spectrum, eta=0.5)}
    original = deepcopy(full)
    figure = fwm.FWMScheme().extra_views()[0].render(full, params)
    try:
        for index, branch in enumerate(("minus", "plus")):
            expected = _current_reference(full[branch], params)
            np.testing.assert_allclose(figure.axes[1].lines[index].get_ydata(),
                                       expected["NF_dB"], equal_nan=True)
            assert np.isnan(figure.axes[1].lines[index].get_ydata()[1])
            np.testing.assert_array_equal(full[branch]["S_dB"], original[branch]["S_dB"])
    finally:
        plt.close(figure)


def test_current_nf_independent_of_retired_pump_scatter_and_input_fano_ui(passive_spectrum):
    scheme = fwm.FWMScheme()
    params = dict(scheme.defaults(), noise_model=fwm.NOISE_PASSIVE_L1, tpd=-8.0)
    shifted = dict(passive_spectrum,
        gain_referred_noise_dB=passive_spectrum["gain_referred_noise_dB"]+9.0,
        hardened_noise={"pump_scatter_noise": 50.0})
    baseline = scheme.headless_observables(passive_spectrum, params)
    changed = scheme.headless_observables(shifted, dict(params, seed_fano=80.0,
        excess_noise=12.0, excess_noise_mode=fwm.NOISE_GAIN_PROPORTIONAL,
        excess_noise_slope=9.0))
    assert baseline["metrics"][0] == changed["metrics"][0]
    np.testing.assert_allclose(baseline["jasperse_readout"]["NF"],
                               changed["jasperse_readout"]["NF"], equal_nan=True)


def test_all_invalid_curve_and_point_outside_scan_remain_unavailable(passive_spectrum):
    import matplotlib.pyplot as plt
    raw = deepcopy(passive_spectrum)
    raw["T_canonical_small_signal"][:] = raw["T_canonical_small_signal"][1]
    scheme = fwm.FWMScheme()
    params = dict(scheme.defaults(), noise_model=fwm.NOISE_PASSIVE_L1, tpd=0.0)
    view = scheme.observables(raw, params)
    try:
        assert view["metrics"][0]["value"] == "unavailable"
        assert np.all(np.isnan(view["figure"].axes[1].lines[0].get_ydata()))
    finally:
        plt.close(view["figure"])
    outside = scheme.headless_observables(passive_spectrum, dict(params, tpd=20.0))
    assert outside["metrics"][0]["value"] == "unavailable"
    assert outside["passive_l1_operating_point"]["passive_l1_status"] == "L1_OUTSIDE_SAMPLED_SCAN"


def test_selected_nf_uses_canonical_transfer_instead_of_displayed_capped_gain(passive_spectrum):
    import matplotlib.pyplot as plt
    raw = deepcopy(passive_spectrum)
    raw["T_canonical_small_signal"][1] = raw["T_canonical_small_signal"][2]
    params = dict(fwm.FWMScheme().defaults(), noise_model=fwm.NOISE_PASSIVE_L1, tpd=-4.0, seed_classical_rin_scaled=3.0)
    reference = _current_reference(raw, params)
    expected_nf = 0.5*(reference["NF"][0]+reference["NF"][1])
    rendered = fwm.FWMScheme().observables(raw, params)
    try:
        point = rendered["passive_l1_operating_point"]
        assert point["NF"] == pytest.approx(expected_nf)
        assert point["NF_dB"] == pytest.approx(10*np.log10(expected_nf))
        axis = rendered["figure"].axes[1]
        assert np.interp(-4.0, *axis.lines[0].get_data()) == pytest.approx(point["NF_dB"])
        assert axis.collections[0].get_offsets()[0, 1] == pytest.approx(point["NF_dB"])
        changed = dict(raw, G_s=100*raw["G_s"], G_c=100*raw["G_c"])
        second = fwm.FWMScheme().headless_observables(changed, params)
        assert second["metrics"][0] == rendered["metrics"][0]
    finally:
        plt.close(rendered["figure"])


def test_result_export_separates_current_a_b_nf_from_cached_coherent_reference(passive_spectrum):
    from gabes_ui.export import build_export_payload
    scheme = fwm.FWMScheme()
    params = dict(scheme.defaults(), noise_model=fwm.NOISE_PASSIVE_L1, seed_classical_rin_scaled=3.0, tpd=-8.0)
    view = scheme.headless_observables(passive_spectrum, params)
    payload = build_export_payload(scheme=scheme, raw=passive_spectrum, params=params,
        view=view, plot_payloads=[], comparison_payload=None,
        cache_version=scheme.cache_version, readout_cache_version="test")
    expected = _current_reference(passive_spectrum, params)
    np.testing.assert_allclose(payload["view"]["jasperse_readout"]["seed_fano_input"],
                               expected["seed_fano_input"])
    assert payload["raw"]["seed_fano_input"] == [1.0]*3
    assert payload["view"]["passive_l1_operating_point"]["NF"] == pytest.approx(
        view["passive_l1_operating_point"]["NF"])
    assert payload["view"]["passive_l1_operating_point"]["NF_dB"] == pytest.approx(
        view["passive_l1_operating_point"]["NF_dB"])
    assert payload["view"]["sidebar_readouts"]["seed_shot_noise_b"] == (
        view["sidebar_readouts"]["seed_shot_noise_b"])

def test_compute_adds_separate_l1_without_changing_legacy_gain_math():
    center = fwm.branch_center_GHz(0.9, -1)
    raw = fwm.compute_spectrum(
        0.9, coarse_points=3, fine_points=0,
        scan_min=center - 0.005, scan_max=center + 0.005,
        velocity_step=200.0, velocity_cutoff=0.1,
        phase_detail=fwm.PHASE_BALANCED, excess_noise_model=False)
    legacy = observables.gain_referred_noise_dB(raw["G_s"], raw["G_c"], raw["eta"])
    np.testing.assert_array_equal(raw["gain_referred_noise_dB"], legacy)
    np.testing.assert_array_equal(raw["S_dB"], legacy)
    expected = passive_l1_from_spectrum(raw)
    np.testing.assert_allclose(raw["passive_l1_noise_dB"],
                               expected["passive_l1_noise_dB"], equal_nan=True)
    np.testing.assert_allclose(raw["passive_l1_G_c"],
                               np.abs(raw["T_canonical_small_signal"][:, 1, 0])**2)
    assert raw["physical_squeezing_dB"] is None
    assert not raw["claim_gate"]["physical_squeezing_prediction"]
