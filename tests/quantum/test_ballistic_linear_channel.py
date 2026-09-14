"""Physical checks for the independent passive finite-velocity oscillator bridge."""

from dataclasses import replace
import json

import numpy as np
import pytest
from scipy.integrate import quad_vec, solve_ivp
from scipy.linalg import expm

from analysis.grand_challenge.reference.ballistic_linear_channel import (
    Parameters, Segment, atomic_source_covariance, build_control,
    direct_profile_ode, drift, optical_diagnostics, propagate_segments, segment_channel,
)


ALGEBRAIC = 5e-12
ODE = 2e-10


@pytest.fixture(scope="module")
def control():
    return build_control()


@pytest.mark.parametrize("omega,gamma,kappa", [
    (0.8, 0.35, 0.9+0.35j), (-1.7, 0.8, -0.3+1.1j),
    (0.0, 0.0, 1.2j), (0.9, 0.6, 0.0),
])
def test_gramian_matches_independent_quadrature_and_preserves_commutator(omega, gamma, kappa):
    p = Parameters(omega=omega, gamma=gamma)
    segment = Segment(1.3, kappa)
    channel = segment_channel(segment, p)
    k, q = drift(kappa, p), np.diag([0.0, 2*gamma/p.v])
    np.testing.assert_allclose(k+k.conj().T+q, 0, atol=ALGEBRAIC, rtol=0)

    def integrand(s):
        e = expm(k*s)
        return e @ q @ e.conj().T

    reference, _ = quad_vec(integrand, 0, segment.length, epsabs=1e-13, epsrel=1e-13)
    np.testing.assert_allclose(channel.reservoir, reference, atol=ALGEBRAIC, rtol=0)
    np.testing.assert_allclose(channel.transfer @ channel.transfer.conj().T+channel.reservoir,
                               np.eye(2), atol=ALGEBRAIC, rtol=0)
    assert np.linalg.eigvalsh(channel.reservoir).min() >= -ALGEBRAIC
    assert min(optical_diagnostics(channel)["cp_eigenvalues"]) >= -ALGEBRAIC


def test_lossless_exact_beamsplitter_and_zero_length():
    p = Parameters(omega=0, gamma=0)
    segment = Segment(0.71, 0.4+0.9j)
    theta = abs(segment.kappa)*segment.length
    phase = segment.kappa/abs(segment.kappa)
    expected = np.array([[np.cos(theta), -1j*phase*np.sin(theta)],
                         [-1j*phase.conjugate()*np.sin(theta), np.cos(theta)]])
    channel = segment_channel(segment, p)
    np.testing.assert_allclose(channel.transfer, expected, atol=ALGEBRAIC, rtol=0)
    np.testing.assert_array_equal(channel.reservoir, np.zeros((2, 2)))
    zero = segment_channel(Segment(0, 7j), p)
    np.testing.assert_array_equal(zero.transfer, np.eye(2))
    np.testing.assert_array_equal(zero.reservoir, np.zeros((2, 2)))
    np.testing.assert_array_equal(propagate_segments([]).transfer, np.eye(2))


def test_uncoupled_free_fluxes_and_atomic_attenuation():
    p = Parameters()
    length = 1.7
    channel = segment_channel(Segment(length, 0), p)
    expected_t = np.diag([np.exp(1j*p.omega*length/p.c),
                          np.exp((1j*p.omega-p.gamma)*length/p.v)])
    expected_w = np.diag([0, -np.expm1(-2*p.gamma*length/p.v)])
    np.testing.assert_allclose(channel.transfer, expected_t, atol=ALGEBRAIC, rtol=0)
    np.testing.assert_allclose(channel.reservoir, expected_w, atol=ALGEBRAIC, rtol=0)
    assert optical_diagnostics(channel)["atomic_boundary_commutator"] == 0


def test_ordered_segments_propagate_earlier_noise_and_match_split_ode():
    p = Parameters(omega=-0.7)
    segments = [Segment(0.4, 1.0+0.2j), Segment(0.7, -0.3+0.8j), Segment(0.2, 0.1-0.5j)]
    result = propagate_segments(segments, p)
    local = [segment_channel(s, p) for s in segments]
    e0, e1, e2 = [ch.transfer for ch in local]
    expected_t = e2 @ e1 @ e0
    expected_w = (e2 @ e1 @ local[0].reservoir @ e1.conj().T @ e2.conj().T
                  + e2 @ local[1].reservoir @ e2.conj().T + local[2].reservoir)
    np.testing.assert_allclose(result.transfer, expected_t, atol=ALGEBRAIC, rtol=0)
    np.testing.assert_allclose(result.reservoir, expected_w, atol=ALGEBRAIC, rtol=0)
    assert np.max(abs(result.transfer-e0 @ e1 @ e2)) > 0.05
    assert np.max(abs(result.reservoir-sum(ch.reservoir for ch in local))) > 0.02

    # Restart only at known discontinuities; integrate T,W jointly without exponentials.
    state = np.concatenate([np.eye(2, dtype=complex).ravel(), np.zeros(4, complex)])
    q = np.diag([0, 2*p.gamma/p.v])
    for segment in segments:
        k = drift(segment.kappa, p)

        def rhs(z, values):
            t, w = values[:4].reshape(2, 2), values[4:].reshape(2, 2)
            return np.concatenate([(k @ t).ravel(), (k @ w+w @ k.conj().T+q).ravel()])

        solution = solve_ivp(rhs, (0, segment.length), state, rtol=2e-12, atol=2e-14, method="DOP853")
        assert solution.success
        state = solution.y[:, -1]
    np.testing.assert_allclose(result.transfer, state[:4].reshape(2, 2), atol=ODE, rtol=0)
    np.testing.assert_allclose(result.reservoir, state[4:].reshape(2, 2), atol=ODE, rtol=0)


@pytest.mark.parametrize("na,nr", [(0.0, 0.0), (0.7, 0.0), (0.0, 0.3), (0.7, 0.3)])
def test_optical_thermal_channel_matches_full_covariance_ode(na, nr):
    p = Parameters()
    segment = Segment(1.3, 0.9+0.35j)
    channel = segment_channel(segment, p)
    optical = optical_diagnostics(channel, atomic_occupation=na, reservoir_occupation=nr)
    k = drift(segment.kappa, p)
    q_sym = np.diag([0, (nr+0.5)*2*p.gamma/p.v])

    def rhs(z, values):
        covariance = values.reshape(2, 2)
        return (k @ covariance+covariance @ k.conj().T+q_sym).ravel()

    initial = np.diag([0.5, na+0.5]).astype(complex)
    solution = solve_ivp(rhs, (0, segment.length), initial.ravel(),
                         method="DOP853", rtol=2e-12, atol=2e-14)
    assert solution.success
    expected_variance = solution.y[:, -1].reshape(2, 2)[0, 0].real
    assert optical["optical_vacuum_input_output_variance"] == pytest.approx(expected_variance, abs=ODE)
    assert optical["variance_occupation_identity_error"] < ALGEBRAIC
    assert min(optical["cp_eigenvalues"]) >= -ALGEBRAIC
    assert optical["optical_vacuum_input_output_occupation"] == pytest.approx(
        na*optical["atomic_boundary_commutator"]+nr*optical["reservoir_commutator"], abs=ALGEBRAIC)
    if na == nr == 0:
        assert expected_variance == pytest.approx(0.5, abs=ODE)
        assert optical["optical_vacuum_input_output_occupation"] == 0


def test_constant_segment_matches_independent_direct_ode(control):
    constant = control["constant_segment"]
    assert constant["ode_transfer_error"] < ODE
    assert constant["ode_reservoir_error"] < ODE
    channel = direct_profile_ode(lambda z: 2j, 0)
    np.testing.assert_array_equal(channel.transfer, np.eye(2))


def test_nonconstant_complex_profile_refines_at_second_order(control):
    rows = control["nonconstant_profile"]["refinements"]
    for old, new in zip(rows, rows[1:]):
        assert new["transfer_error"] < old["transfer_error"]/3.7
        assert new["reservoir_error"] < old["reservoir_error"]/3.7
        assert new["observed_order"] > 1.9
    assert max(rows[-1]["transfer_error"], rows[-1]["reservoir_error"]) < 2e-5
    for row in rows:
        assert row["full_commutator_error"] < ALGEBRAIC
        assert row["reservoir_min_eigenvalue"] >= -ALGEBRAIC
        assert min(row["optical"]["cp_eigenvalues"]) >= -ALGEBRAIC
        assert row["optical"]["optical_vacuum_input_output_variance"] == pytest.approx(0.5, abs=ALGEBRAIC)
    assert control["negative_controls"]["reversed_profile_transfer_error"] > 0.05


def test_atomic_source_retains_complex_spatial_cross_covariance():
    p = Parameters()
    z = np.array([0.0, 0.2, 0.7, 1.3])
    kernels = atomic_source_covariance(z, p, atomic_occupation=0.7, reservoir_occupation=0.3)
    lam = (1j*p.omega-p.gamma)/p.v
    for i, zi in enumerate(z):
        for j, zj in enumerate(z):
            def integrand(s):
                return 2*p.gamma/p.v*np.exp(lam*(zi-s)+lam.conjugate()*(zj-s))
            expected, _ = quad_vec(integrand, 0, min(zi, zj), epsabs=1e-13)
            assert kernels["reservoir_commutator"][i, j] == pytest.approx(expected, abs=ALGEBRAIC)
    assert abs(kernels["reservoir_commutator"][1, 3]) > 0.1
    assert abs(kernels["reservoir_commutator"][1, 3].imag) > 0.05
    for name in ("boundary_commutator", "reservoir_commutator", "total_symmetrized"):
        np.testing.assert_allclose(kernels[name], kernels[name].conj().T, atol=ALGEBRAIC, rtol=0)
        assert np.linalg.eigvalsh(kernels[name]).min() >= -ALGEBRAIC
    np.testing.assert_allclose(kernels["total_symmetrized"],
                               1.2*kernels["boundary_commutator"]+0.8*kernels["reservoir_commutator"],
                               atol=ALGEBRAIC, rtol=0)


def test_zero_damping_source_is_shared_boundary_not_independent_noise():
    p = Parameters(gamma=0)
    z = [0, 0.3, 0.9]
    kernels = atomic_source_covariance(z, p)
    np.testing.assert_array_equal(kernels["reservoir_commutator"], np.zeros((3, 3)))
    np.testing.assert_allclose(abs(kernels["boundary_commutator"]), 1, atol=ALGEBRAIC, rtol=0)
    assert np.linalg.matrix_rank(kernels["boundary_commutator"], tol=ALGEBRAIC) == 1


def test_correlated_green_kernel_reproduces_optical_noise_and_diagonalization_fails(control):
    spatial = control["spatial_atomic_source"]
    assert spatial["stationary_kernel_identity_error"] < ALGEBRAIC
    rows = spatial["quadrature_refinements"]
    assert all(row["boundary_error"] < ALGEBRAIC for row in rows)
    for old, new in zip(rows, rows[1:]):
        assert new["reservoir_error"] < old["reservoir_error"]/3.8
    assert rows[-1]["full_commutator_error"] < 5e-6
    assert rows[-1]["diagonalized_bath_commutator_deficit"] > 0.3
    assert rows[-1]["diagonalized_reservoir_commutator"] > 0
    assert "analysis only" in spatial["negative_control_status"]


def test_dropping_atomic_inflow_is_noncanonical_despite_positive_noise(control):
    negative = control["negative_controls"]
    assert negative["remaining_added_noise_min_eigenvalue"] > 0.1
    assert negative["commutator_deficit"] > 0.3
    assert negative["deficit_identity_error"] < ALGEBRAIC
    assert negative["gaussian_cp_min_eigenvalue_without_inflow"] == pytest.approx(
        -negative["expected_deficit_boundary"]/2, abs=ALGEBRAIC)
    assert negative["apparent_vacuum_variance_without_inflow"] < 0.4
    assert "analysis only" in negative["status"]


def test_control_is_strict_json_serializable_and_declares_scope_and_thresholds(control):
    decoded = json.loads(json.dumps(control, allow_nan=False))
    assert decoded["passed"] is True
    assert decoded["expected_controls_passed"] is True
    assert decoded["source"].endswith("reference/ballistic_linear_channel.py")
    assert decoded["frequency"] == {"omega": 0.8, "units": "rad/t", "fourier_convention": "exp(-i Omega t)"}
    assert decoded["units"]["v"] == "ell/t"
    assert decoded["schema_version"] == 1
    assert "finite-velocity" in decoded["scope"]
    assert "FWM squeezing" in decoded["exclusions"]
    assert decoded["precision_thresholds"]["algebraic_absolute"] == ALGEBRAIC


@pytest.mark.parametrize("field,value", [("v", 0), ("v", -1), ("c", 0),
                                         ("gamma", -1), ("omega", np.nan)])
def test_invalid_physical_parameters_rejected(field, value):
    with pytest.raises(ValueError):
        replace(Parameters(), **{field: value})


def test_invalid_lengths_and_thermal_states_rejected():
    with pytest.raises(ValueError):
        Segment(-1, 1)
    with pytest.raises(ValueError):
        Segment(1, complex(np.inf, 0))
    with pytest.raises(ValueError):
        atomic_source_covariance([0, -1])
    with pytest.raises(ValueError):
        optical_diagnostics(segment_channel(Segment(1, 1)), atomic_occupation=-0.1)
