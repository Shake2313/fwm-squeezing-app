"""Independent physical-quadrature reference for the conditional L1 readout.

The drift below is written directly for (x_a, p_a, x_b, p_b), with vacuum
variance one. Its Van Loan block convention and coordinates differ from the
handoff's reference. The reference never reconstructs diffusion from M.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from gabes.passive_readout import passive_l1_readout


_OMEGA = np.kron(np.eye(2), [[0.0, 1.0], [-1.0, 0.0]])


def _segment_matrices(segment):
    coupling, phase, delta_a, delta_b, loss_a, loss_b, length = segment
    real = coupling * np.cos(phase)
    imag = coupling * np.sin(phase)
    drift = np.array([
        [-loss_a / 2, -delta_a, real, imag],
        [delta_a, -loss_a / 2, imag, -real],
        [real, imag, -loss_b / 2, -delta_b],
        [imag, -real, delta_b, -loss_b / 2],
    ])
    diffusion = np.diag([loss_a, loss_a, loss_b, loss_b])
    nambu_drift = np.array([
        [-loss_a / 2 + 1j * delta_a, real + 1j * imag],
        [real - 1j * imag, -loss_b / 2 - 1j * delta_b],
    ])
    return drift, diffusion, nambu_drift, length


def _quadrature_reference(segments, fano, phase_fano=1.0, amplitude_phase=0.0):
    covariance = np.diag([fano, phase_fano, 1.0, 1.0])
    covariance[0, 1] = covariance[1, 0] = amplitude_phase
    mean = np.array([1.0, 0.0, 0.0, 0.0])
    nambu_transfer = np.eye(2, dtype=complex)
    for segment in segments:
        drift, diffusion, nambu_drift, length = _segment_matrices(segment)
        transfer = expm(drift * length)
        block = np.block([
            [drift, diffusion], [np.zeros((4, 4)), -drift.T],
        ])
        integral = expm(block * length)[:4, 4:] @ transfer.T
        covariance = transfer @ covariance @ transfer.T + integral
        mean = transfer @ mean
        nambu_transfer = expm(nambu_drift * length) @ nambu_transfer
    probe, conjugate = mean[:2], mean[2:]
    gs, gc = probe @ probe, conjugate @ conjugate
    difference_direction = np.r_[probe, -conjugate]
    return {
        "M": nambu_transfer,
        "V": covariance,
        "mean": mean,
        "Gs": gs,
        "Gc": gc,
        "S": difference_direction @ covariance @ difference_direction / (gs + gc),
        "Fp": probe @ covariance[:2, :2] @ probe / gs,
        "Fc": conjugate @ covariance[2:, 2:] @ conjugate / gc,
    }


def _random_profiles(count=160):
    rng = np.random.default_rng(2810)
    for _ in range(count):
        segments = [
            (rng.uniform(0.03, 0.8), rng.uniform(-np.pi, np.pi),
             *rng.uniform(-1.0, 1.0, 2), rng.uniform(0.0, 1.0), 0.0,
             rng.uniform(0.15, 0.8))
            for _ in range(rng.integers(1, 6))
        ]
        yield segments, rng.uniform(1.0, 10.0)


def test_van_loan_reference_agrees_with_independent_covariance_ode():
    segment = (0.63, 0.71, 0.42, -0.33, 0.8, 0.4, 0.73)
    drift, diffusion, _, length = _segment_matrices(segment)
    initial_covariance = np.diag([3.7, 2.0, 1.0, 1.0])
    initial = np.r_[[1.0, 0.0, 0.0, 0.0], initial_covariance.ravel()]

    def derivative(_, state):
        covariance = state[4:].reshape(4, 4)
        return np.r_[drift @ state[:4],
                     (drift @ covariance + covariance @ drift.T + diffusion).ravel()]

    solution = solve_ivp(derivative, (0.0, length), initial, method="DOP853",
                         rtol=2e-12, atol=2e-13)
    assert solution.success
    reference = _quadrature_reference([segment], 3.7, phase_fano=2.0)
    assert_allclose(solution.y[:4, -1], reference["mean"], atol=3e-12)
    assert_allclose(solution.y[4:, -1].reshape(4, 4), reference["V"], atol=3e-12)


def test_complex_distributed_probe_loss_matches_independent_reference():
    upper_bound_cases = 0
    for segments, fano in _random_profiles():
        reference = _quadrature_reference(segments, fano)
        result = passive_l1_readout(reference["M"], seed_fano=fano, eta=0.713)
        assert result["passive_completion_valid"]
        assert_allclose(
            [result["passive_l1_S_src"], result["passive_l1_F_p"],
             result["passive_l1_F_c"]],
            [reference["S"], reference["Fp"], reference["Fc"]],
            rtol=2e-11, atol=2e-12,
        )
        assert_allclose(result["passive_l1_S_det"], 0.287 + 0.713 * reference["S"],
                        rtol=2e-11, atol=2e-12)
        det_abs2 = np.exp(-sum(segment[4] * segment[-1] for segment in segments))
        assert_allclose(result["passive_l1_det_abs2"], det_abs2, rtol=2e-11)
        assert_allclose(result["passive_l1_Lambda"], 2.0 * det_abs2,
                        rtol=2e-11, atol=2e-12)
        assert np.linalg.eigvalsh(reference["V"] + 1j * _OMEGA)[0] >= -2e-11
        difference, total = result["passive_l1_D"], result["passive_l1_T"]
        assert difference <= 1.0 + 2e-12
        assert abs(reference["M"][0, 1]) ** 2 - abs(reference["M"][1, 1]) ** 2 <= -1 + 2e-12
        assert abs(result["passive_l1_q1"]) ** 2 >= max(0, difference ** 2 - difference) - 2e-12
        if 0.0 <= difference <= 1.0:
            upper_bound_cases += 1
            upper = (2.0 - difference + (fano - 1.0) * difference ** 2) / total
            assert reference["S"] <= upper + 2e-12
    assert upper_bound_cases >= 100


def test_seed_phase_and_amplitude_phase_covariance_cancel():
    for segments, fano in _random_profiles(40):
        baseline = _quadrature_reference(segments, fano)
        for phase_fano, correlation in [(50.0, 0.0), (3.0, 0.7)]:
            changed = _quadrature_reference(segments, fano, phase_fano, correlation)
            assert_allclose([changed["S"], changed["Fp"], changed["Fc"]],
                            [baseline["S"], baseline["Fp"], baseline["Fc"]],
                            rtol=2e-11, atol=2e-12)
            difference = changed["Gs"] - changed["Gc"]
            closure = changed["S"] * (changed["Gs"] + changed["Gc"]) - difference * (
                changed["Fp"] - changed["Fc"] - 1.0)
            assert_allclose(closure, 2.0 * abs(np.linalg.det(changed["M"])) ** 2,
                            rtol=2e-11, atol=2e-12)


def test_conjugate_vacuum_loss_adds_nonnegative_noise_and_lambda():
    rng = np.random.default_rng(3902)
    for segments, fano in _random_profiles(80):
        segments = [(*segment[:5], rng.uniform(0.05, 0.8), segment[-1])
                    for segment in segments]
        reference = _quadrature_reference(segments, fano)
        result = passive_l1_readout(reference["M"], seed_fano=fano)
        assert reference["S"] >= result["passive_l1_S_src"] - 2e-12
        difference, total = reference["Gs"] - reference["Gc"], reference["Gs"] + reference["Gc"]
        closure = reference["S"] * total - difference * (reference["Fp"] - reference["Fc"] - 1)
        assert closure >= 2.0 * abs(np.linalg.det(reference["M"])) ** 2 - 2e-12


def test_zero_mean_difference_can_still_have_positive_sub_sql_variance():
    squeezing = 1.0
    transmission = np.tanh(squeezing) ** 2
    segments = [(squeezing, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0),
                (0.0, 0.0, 0.0, 0.0, -np.log(transmission), 0.0, 1.0)]
    reference = _quadrature_reference(segments, 4.05)
    result = passive_l1_readout(reference["M"], seed_fano=4.05)
    assert abs(result["passive_l1_D"]) < 2e-14
    assert 0.0 < result["passive_l1_S_src"] < 1.0
    assert_allclose(result["passive_l1_S_src"], reference["S"], atol=2e-12)
    assert_allclose(result["passive_l1_S_src"], 1 / np.cosh(squeezing) ** 2, atol=2e-12)


def test_phase_sensitive_squeezed_reservoir_disproves_unconditional_lower_bound():
    transmission, squeezing = 0.4, 1.0
    transfer = np.diag([np.sqrt(transmission), 1.0])
    result = passive_l1_readout(transfer)
    output_x = transmission + (1 - transmission) * np.exp(-2 * squeezing)
    output_p = transmission + (1 - transmission) * np.exp(2 * squeezing)
    quadrature_transfer = np.diag([np.sqrt(transmission), np.sqrt(transmission), 1, 1])
    added_covariance = np.diag([(1 - transmission) * np.exp(-2 * squeezing),
                               (1 - transmission) * np.exp(2 * squeezing), 0, 0])
    cp_matrix = added_covariance + 1j * (
        _OMEGA - quadrature_transfer @ _OMEGA @ quadrature_transfer.T)
    assert np.linalg.eigvalsh(cp_matrix)[0] >= -2e-12
    assert result["passive_completion_valid"]
    assert output_x * output_p >= 1.0  # physical Gaussian output uncertainty
    assert output_x < result["passive_l1_S_src"]


def test_endpoint_completion_does_not_imply_distributed_probe_loss_bound():
    result = passive_l1_readout(np.diag([0.9, 10.0]))
    assert result["passive_completion_valid"]
    assert_allclose(result["passive_l1_det_abs2"], 81.0)
    assert_allclose(result["passive_l1_Lambda"], 162.0)


def test_negative_formal_variance_is_retained_but_never_displayed():
    result = passive_l1_readout(np.array([[0.0, 0.0], [1.0, 0.0]]), eta=0.2)
    assert not result["passive_completion_valid"]
    assert result["passive_l1_S_src"] == -1.0
    assert np.isfinite(result["passive_l1_formal_noise_dB"])
    assert np.isnan(result["passive_l1_noise_dB"])


def test_high_gain_does_not_hide_resolved_negative_completion_eigenvalue():
    squeezing = 8.0
    matrix = 1.0001 * np.array([[np.cosh(squeezing), np.sinh(squeezing)],
                              [np.sinh(squeezing), np.cosh(squeezing)]])
    result = passive_l1_readout(matrix)
    assert result["passive_completion_min_eig"] < -1e-4
    assert not result["passive_completion_valid"]
    assert np.isnan(result["passive_l1_noise_dB"])


@pytest.mark.parametrize("matrix,kwargs", [
    (np.eye(3), {}),
    (np.array([[np.nan, 0.0], [0.0, 1.0]]), {}),
    (np.eye(2), {"seed_fano": -0.1}),
    (np.eye(2), {"seed_fano": 0.0}),
    (np.eye(2), {"seed_fano": np.inf}),
    (np.eye(2), {"eta": -0.1}),
    (np.eye(2), {"eta": 1.1}),
    (np.eye(2), {"eta": np.nan}),
])
def test_invalid_inputs_are_rejected(matrix, kwargs):
    with pytest.raises(ValueError):
        passive_l1_readout(matrix, **kwargs)


def test_parameter_broadcast_preserves_lossless_fano_and_detection_formula():
    squeezing = np.array([0.0, 0.7, 1.2])
    matrices = np.empty((3, 2, 2))
    matrices[:, 0, 0] = matrices[:, 1, 1] = np.cosh(squeezing)
    matrices[:, 0, 1] = matrices[:, 1, 0] = np.sinh(squeezing)
    fano = np.array([[1.0], [4.05]])
    efficiency = np.array([[1.0], [0.7068]])
    result = passive_l1_readout(matrices, seed_fano=fano, eta=efficiency)
    expected_source = fano / np.cosh(2.0 * squeezing)
    assert result["passive_l1_S_src"].shape == (2, 3)
    assert_allclose(result["passive_l1_S_src"], expected_source, atol=2e-12)
    assert_allclose(result["passive_l1_S_det"],
                    1 - efficiency + efficiency * expected_source, atol=2e-12)
    assert result["passive_l1_noise_dB"][1, 0] > 0.0  # above-SQL noise is retained


def test_no_bright_output_has_no_source_normalization_or_display_value():
    result = passive_l1_readout(np.zeros((2, 2)))
    assert result["passive_l1_T"] == 0.0
    assert np.isnan(result["passive_l1_S_src"])
    assert np.isnan(result["passive_l1_noise_dB"])
    assert result["passive_l1_status"] == "L1_NO_BRIGHT_OUTPUT_OR_INVALID_VARIANCE"


def test_masked_l1_exports_as_strict_json_with_null_and_preserved_status():
    import json

    from gabes_ui.export import to_json_primitive

    matrices = np.stack([np.eye(2), np.array([[0.0, 0.0], [1.0, 0.0]])])
    readout = passive_l1_readout(matrices)
    payload = to_json_primitive(readout)
    # An unavailable point needs a JSON null, never the nonstandard NaN token.
    assert payload["passive_l1_noise_dB"] == [0.0, None]
    assert payload["passive_l1_S_src"][1] == -1.0
    assert not payload["passive_completion_valid"][1]
    encoded = json.dumps(payload, allow_nan=False)
    assert json.loads(encoded)["passive_l1_status"] == payload["passive_l1_status"]


def test_strict_result_bundle_keeps_current_fano_and_invalid_point_evidence():
    import json

    from gabes.schemes import fwm
    from gabes_ui.export import build_export_payload

    scheme = fwm.FWMScheme()
    params = dict(scheme.defaults(), seed_classical_rin_scaled=20.0)
    raw = {"T_canonical_small_signal": np.stack([
        np.eye(2), np.array([[0.0, 0.0], [1.0, 0.0]]), np.eye(2)]),
        "G_s": np.ones(3), "eta": 0.8, "probe_axis_GHz": np.array([0.0, 1.0, 2.0])}
    raw.update(passive_l1_readout(raw["T_canonical_small_signal"], eta=raw["eta"]))
    current = fwm._seeded_passive_l1_readout(raw, params)
    point = fwm._passive_l1_selected(current, raw, 0.0)
    payload = build_export_payload(
        scheme=scheme, raw=raw, params=params,
        view={"passive_l1_readout": current, "passive_l1_operating_point": point},
        plot_payloads=[], comparison_payload=None,
        cache_version=scheme.cache_version, readout_cache_version="independent-test")
    decoded = json.loads(json.dumps(payload, allow_nan=False))
    assert decoded["raw"]["seed_fano_input"] == [1.0] * 3
    omega = fwm.seeded_option_a_wavenumbers(params["opd"], raw["probe_axis_GHz"])[1] * fwm.constants.C_LIGHT
    expected_fano = 1.0 + 20e-14 / (2*fwm.constants.HBAR*omega/(params["probe_uw"]*1e-6))
    np.testing.assert_allclose(
        decoded["view"]["passive_l1_readout"]["seed_fano_input"], expected_fano, rtol=1e-14)
    assert decoded["view"]["passive_l1_readout"]["passive_l1_noise_dB"][1] is None
    assert not decoded["view"]["passive_l1_readout"]["passive_completion_valid"][1]
    assert decoded["view"]["passive_l1_readout"]["passive_l1_status"][1] == "L1_COMPLETION_NOT_PASSIVE"
    assert decoded["view"]["passive_l1_operating_point"]["passive_l1_available"]
