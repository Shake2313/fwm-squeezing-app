import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import observables  # noqa: E402


def test_asymmetric_balanced_noise_recovers_symmetric_eta_formula():
    G_s = np.array([1.2, 3.0, 15.0, 80.0])
    G_c = G_s - 1.0
    eta = 0.87

    old = observables.intensity_difference_squeezing_dB(G_s, G_c, eta)
    new = observables.balanced_twin_beam_squeezing_dB(
        G_s, G_c, eta_s=eta, eta_c=eta, reference_weight="raw")

    assert np.allclose(new, old, rtol=1e-12, atol=1e-12)


def test_shot_weight_minimizes_normalized_noise_for_unequal_powers():
    G_s = 12.0
    G_c = 10.8
    eta_s = 0.42
    eta_c = 0.91

    shot = observables.balanced_twin_beam_noise(
        G_s, G_c, eta_s, eta_c, reference_weight="shot")
    raw = observables.balanced_twin_beam_noise(
        G_s, G_c, eta_s, eta_c, reference_weight="raw")
    dc = observables.balanced_twin_beam_noise(
        G_s, G_c, eta_s, eta_c, reference_weight="dc")

    assert shot <= raw
    assert shot <= dc
    assert shot < 1.0


def test_measured_source_noise_anchor_propagates_symmetric_loss():
    source_noise = 10 ** (-8.3 / 10.0)
    eta = 0.72
    got = observables.balanced_twin_beam_noise(
        5.6, 4.6, eta_s=eta, eta_c=eta, reference_weight="raw",
        source_noise=source_noise)
    expected = eta * source_noise + (1.0 - eta)

    assert np.isclose(got, expected)


if __name__ == "__main__":
    test_asymmetric_balanced_noise_recovers_symmetric_eta_formula()
    test_shot_weight_minimizes_normalized_noise_for_unequal_powers()
    test_measured_source_noise_anchor_propagates_symmetric_loss()
    print("balanced readout tests OK")
