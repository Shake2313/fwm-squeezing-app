"""Current Fast atomic sampling is independent of external detection settings."""
import numpy as np
import pytest

from gabes.schemes import fwm


@pytest.mark.parametrize("detuning_ghz", [0.9, 1.4])
def test_current_fast_source_is_identical_across_detector_efficiencies(detuning_ghz):
    center = fwm.branch_center_GHz(detuning_ghz, -1)
    settings = dict(
        T=394.15, P_pump=0.6, P_probe=8e-6, line_strength=0.74,
        coarse_points=401, fine_points=0, velocity_step=1.0,
        velocity_cutoff=4.0, branch=-1, phase_detail=fwm.PHASE_ULTRA,
        scan_min=center-fwm.WINDOW_GHZ, scan_max=center+fwm.WINDOW_GHZ,
        model_fidelity=fwm.FIDELITY_FAST, response_method=fwm.RESPONSE_POLE_ADAPTIVE,
        legacy_gain_calibration_fixture=False, legacy_mean_field_fixture=False)
    outputs = [fwm.compute_spectrum(
        detuning_ghz, detection_efficiency=eta, **settings)
        for eta in (0.0, 0.5, 1.0)]
    for output in outputs[1:]:
        for key in ("probe_axis_GHz", "G_s", "G_c", "G_s_smallsignal",
                    "G_c_smallsignal", "T_field_small_signal", "T_canonical_small_signal"):
            np.testing.assert_array_equal(output[key], outputs[0][key])
        assert output["response_estimator"]["solved_probe_indices"] == (
            outputs[0]["response_estimator"]["solved_probe_indices"])
    assert outputs[0]["response_estimator"]["solved_probe_points"] < 401


def test_current_detector_and_seed_noise_controls_are_outside_atomic_cache():
    scheme = fwm.FWMScheme()
    keys = set(scheme.recompute_keys())
    assert {"after_cell_eff_pct", "qe_pct", "loss_pct", "seed_classical_rin_scaled",
            "measured_gain_calibration"}.isdisjoint(keys)
