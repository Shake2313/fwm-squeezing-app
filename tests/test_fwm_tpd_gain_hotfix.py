"""Laboratory gain acceptance and empirical/physical boundary regressions."""
import json
from pathlib import Path

import numpy as np
import pytest

from gabes import fwm_gain_closure
from gabes.schemes import fwm
from analysis.fwm_gain_hotfix.snapshot_convention import snapshot_directory

RECORD = (Path(__file__).resolve().parents[1] /
          "analysis/fwm_gain_hotfix/tpd_affine_20261001")
FIT = json.loads((RECORD / "fit.json").read_text(encoding="utf-8"))
EOM = np.array(FIT["eom_GHz"])
DELTA_MHZ = (fwm.constants.NU_HF / 1e9 - EOM) * 1000.


def _params(tier, **overrides):
    params = dict(fwm.FWMScheme().defaults(), resolution=tier, temp_c=118.,
                  pump_mw=380., opd=1.06, probe_uw=3.7)
    params.update(overrides)
    return params


def _direct(tier, **overrides):
    kwargs = dict(T=391.15, P_pump=0.380, P_probe=3.7e-6, L=0.0125,
                  w_pump=530e-6, w_probe=330e-6, pump_probe_angle_deg=0.32,
                  transit_rate=2*np.pi*100e3, branch=-1, line_strength=0.74,
                  floquet_order=3, phase_detail=fwm.PHASE_ULTRA,
                  detection_efficiency=0.8694,
                  model_fidelity=tier, response_method=fwm.RESPONSE_POLE,
                  velocity_step=1., velocity_cutoff=4., coarse_points=10,
                  fine_points=0, scan_min=1.06-EOM[0], scan_max=1.06-EOM[-1])
    kwargs.update(overrides)
    return fwm.compute_spectrum(1.06, **kwargs)


@pytest.fixture(scope="module", params=(fwm.FIDELITY_FAST, fwm.FIDELITY_BALANCED))
def laboratory(request):
    tier = request.param
    return tier, _direct(tier), fwm.FWMScheme().compute(_params(tier))


def test_all_ten_measured_gains_match_in_direct_and_frontend_paths(laboratory):
    _, direct, displayed = laboratory
    for raw in (direct, displayed):
        values = [fwm.operating_point(raw, delta) for delta in DELTA_MHZ]
        for key, arm in (("G_s", "probe"), ("G_c", "conjugate")):
            prediction = np.array([value[key] for value in values])
            measured = np.array(FIT[arm]["measured_gain"])
            assert np.max(np.abs(prediction/measured - 1.)) < 0.08
            assert np.all(np.diff(prediction) > 0.)
        calibration = raw["gain_closure"]["output_calibration"]
        assert calibration["applied"]
        assert calibration["statistical_confidence_interval"] is None
        assert not calibration["independently_validated"]
        assert not calibration["transfer_matrix_calibrated"]
        assert raw["physical_squeezing_dB"] is None
        assert not raw["claim_gate"]["quantitative_gain_supported"]


def test_detector_efficiency_changes_noise_but_not_fast_atomic_sampling(laboratory):
    tier, _, displayed = laboratory
    other = fwm.FWMScheme().compute(_params(tier, detection_eff_pct=50.))
    for key in ("probe_axis_GHz", "G_s", "G_c", "G_s_smallsignal", "G_c_smallsignal"):
        np.testing.assert_array_equal(displayed[key], other[key])
    assert not np.array_equal(displayed["S_dB"], other["S_dB"])


def test_noise_configuration_does_not_steer_calibrated_fast_gain():
    kwargs = dict(response_method=fwm.RESPONSE_POLE_ADAPTIVE)
    default = _direct(fwm.FIDELITY_FAST, **kwargs)
    for noise in (False, {"pump_scatter_kappa": 4.0}):
        changed = _direct(fwm.FIDELITY_FAST, excess_noise_model=noise, **kwargs)
        for key in ("probe_axis_GHz", "G_s", "G_c"):
            np.testing.assert_array_equal(default[key], changed[key])


def test_output_fit_preserves_canonical_transfer_and_switch_off_baseline(laboratory):
    tier, direct, _ = laboratory
    baseline_name = ("before_parent.json" if snapshot_directory().name == "before_parent"
                     else "before.json")
    before = json.loads((RECORD / baseline_name).read_text(encoding="utf-8"))
    previous = before["cases"]["exact_" + tier.split()[0]]
    for key in ("G_s_smallsignal", "G_c_smallsignal"):
        np.testing.assert_allclose(direct[key], previous[key], rtol=1e-7, atol=1e-10)
    bypass = _direct(tier, gain_closure_enabled=False)
    previous = before["cases"]["off_" + tier.split()[0]]
    for key in ("G_s", "G_c", "S_dB"):
        np.testing.assert_allclose(bypass[key], previous[key], rtol=1e-7, atol=1e-10)
    assert not bypass["gain_closure"]["output_calibration"]["applied"]


@pytest.mark.parametrize("opd,branch,enabled,eligible", (
    (0.9, -1, True, True), (0.8, -1, True, True), (1.22, -1, True, True),
    (1.06, 1, True, True), (1.06, -1, False, True), (1.06, -1, True, False)))
def test_ineligible_output_map_is_an_exact_bypass(opd, branch, enabled, eligible):
    calibration = fwm_gain_closure.laboratory_provenance(
        enabled=enabled, eligible=eligible, D_GHz=opd, branch=branch)
    gs, gc = np.array([1., 9., 17.]), np.array([0., 8., 16.])
    actual_s, actual_c = fwm_gain_closure.apply_power_gain(
        gs, gc, calibration=calibration, P_pump=0.380, P_seed=3.7e-6)
    assert actual_s is gs and actual_c is gc


def test_active_output_map_preserves_no_gain_monotonicity_and_pump_bounds():
    calibration = fwm_gain_closure.laboratory_provenance(
        enabled=True, eligible=True, D_GHz=1.06, branch=-1)
    gs = np.r_[1., np.geomspace(1.01, 1e9, 100)]
    gc = np.r_[0., np.geomspace(0.01, 1e9, 100)]
    corrected_s, corrected_c = fwm_gain_closure.apply_power_gain(
        gs, gc, calibration=calibration, P_pump=0.380, P_seed=3.7e-6)
    assert corrected_s[0] == 1. and corrected_c[0] == 0.
    assert np.all(np.diff(corrected_s) >= 0.)
    assert np.all(np.diff(corrected_c) >= 0.)
    assert np.all((corrected_s-1.)*3.7e-6 <= 0.190)
    assert np.all(corrected_c*3.7e-6 <= 0.190)
    no_pump_s, no_pump_c = fwm_gain_closure.apply_power_gain(
        gs, gc, calibration=calibration, P_pump=0., P_seed=3.7e-6)
    np.testing.assert_array_equal(no_pump_s, np.ones_like(gs))
    np.testing.assert_array_equal(no_pump_c, np.zeros_like(gc))
