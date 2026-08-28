"""Pump-only stationary-state and infinitesimal FWM response references."""

import numpy as np
import pytest

from gabes import core
from gabes.schemes import fwm


def _fixture():
    pump = fwm.rabi_freq(0.6, fwm.W_PUMP)
    atom = fwm.collisional_atom(394.15)
    delta = 2.0 * np.pi * np.array([-8e6])
    deff = 2.0 * np.pi * np.array([0.9e9])
    return pump, atom, delta, deff


def test_static_pump_state_matches_independent_svd_null_vector():
    pump, atom, delta, deff = _fixture()
    reference = fwm.pump_only_weak_response_reference(
        pump, pump, delta, deff, atom=atom)
    rho = reference.pump_state[0]
    L = core.build_liouvillian(
        fwm.pump_hamiltonian_at_Deff_zero(pump, pump), atom)
    L = L - deff[0] * atom.S_v

    _u, _s, vh = np.linalg.svd(L)
    null = vh.conj().T[:, -1]
    diagonal = np.arange(fwm.N_LEVELS) * (fwm.N_LEVELS + 1)
    null = null / np.sum(null[diagonal])
    np.testing.assert_allclose(rho.reshape(-1), null, rtol=2e-11, atol=2e-12)
    np.testing.assert_allclose(
        np.real(np.diag(rho)),
        [0.11574137, 0.84317396, 0.01515365, 0.02593102],
        rtol=4e-8, atol=4e-9)
    assert reference.diagnostics["max_pump_normalized_residual"] < 1e-12
    assert reference.diagnostics["max_pump_trace_error"] < 1e-12
    assert reference.diagnostics["max_pump_hermiticity_error"] < 1e-12
    assert reference.diagnostics["minimum_pump_eigenvalue"] > 0.0


def test_zero_pump_returns_thermal_reload_state_and_passive_response():
    _pump, atom, delta, _deff = _fixture()
    deff = 2.0 * np.pi * np.array([0.75e9, 0.9e9, 1.05e9])
    reference = fwm.pump_only_weak_response_reference(
        0.0, 0.0, delta, deff, atom=atom)
    target = np.diag([5.0 / 12.0, 7.0 / 12.0, 0.0, 0.0])
    np.testing.assert_allclose(
        reference.pump_state,
        np.broadcast_to(target, reference.pump_state.shape), atol=2e-15)
    np.testing.assert_allclose(reference.chi_cs, 0.0, atol=1e-30)
    np.testing.assert_allclose(reference.chi_sc, 0.0, atol=1e-30)
    assert np.all(np.imag(reference.chi_ss) <= 0.0)
    assert np.all(np.imag(reference.chi_cc) <= 0.0)


def test_pump_reference_default_atom_matches_production_dissipator():
    pump, _atom, delta, deff = _fixture()
    implicit = fwm.pump_only_weak_response_reference(
        pump, pump, delta, deff)
    explicit_atom = fwm.collisional_atom(fwm.T_CELL)
    explicit = fwm.pump_only_weak_response_reference(
        pump, pump, delta, deff, atom=explicit_atom)
    np.testing.assert_allclose(
        implicit.pump_state, explicit.pump_state, rtol=2e-13, atol=2e-15)
    np.testing.assert_allclose(
        implicit.chi_matrix, explicit.chi_matrix, rtol=2e-13, atol=2e-23)

    implicit_poles = fwm.pump_only_pole_residue_reference(
        pump, pump, delta[0], deff[0])
    explicit_poles = fwm.pump_only_pole_residue_reference(
        pump, pump, delta[0], deff[0], atom=explicit_atom)
    np.testing.assert_allclose(
        implicit_poles["response"], explicit_poles["response"],
        rtol=2e-13, atol=2e-23)


@pytest.mark.parametrize("n_f", [1, 3])
def test_static_pump_frame_is_unitarily_equal_to_minus_branch_floquet(n_f):
    pump, atom, _delta, _deff = _fixture()
    delta_axis = 2.0 * np.pi * np.array([-40e6, -8e6, 0.0, 25e6])
    deff = 2.0 * np.pi * np.array([0.75e9, 0.9e9, 1.05e9])
    static = fwm.pump_only_weak_response_reference(
        pump, pump, delta_axis[:1], deff, atom=atom).pump_state
    expected = fwm.pump_frame_to_seeded_harmonics(static, n_f=n_f, branch=-1)
    if n_f > 1:
        np.testing.assert_allclose(expected[:, :n_f - 1], 0.0, atol=0.0)
        np.testing.assert_allclose(expected[:, n_f + 2:], 0.0, atol=0.0)

    Cp, Cm = fwm.sideband_template(pump, pump, 0.0, branch=-1)
    for delta in delta_axis:
        L0 = core.build_liouvillian(
            fwm.static_hamiltonian_at_Deff_zero(
                pump, pump, 0.0, delta, branch=-1), atom)
        actual = core.floquet_solve_truncated(
            L0, Cp, Cm, fwm.seeded_sideband_beat(delta, -1), deff,
            atom.S_v, fwm.N_LEVELS, n_f, return_harmonics=True)
        np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=2e-13)


def test_plus_branch_reference_fails_closed_until_frame_is_corrected():
    pump, atom, delta, deff = _fixture()
    with pytest.raises(NotImplementedError, match="plus-branch"):
        fwm.pump_only_weak_response_reference(
            pump, pump, delta, deff, atom=atom, branch=+1)
    with pytest.raises(NotImplementedError, match="plus-branch"):
        fwm.pump_frame_to_seeded_harmonics(
            np.eye(fwm.N_LEVELS) / fwm.N_LEVELS, n_f=3, branch=+1)


def test_weak_nambu_response_golden_residual_and_frequency_separation():
    pump, atom, delta, deff = _fixture()
    analysis = 2.0 * np.pi * np.array([0.0, 0.1e6, 4.0e6])
    reference = fwm.pump_only_weak_response_reference(
        pump, pump, delta, deff, atom=atom,
        analysis_frequency_axis_rad_s=analysis)
    expected_relative = -fwm.OMEGA_HF + delta[None, :] + analysis[:, None]
    np.testing.assert_allclose(reference.relative_frequency_rad_s,
                               expected_relative, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(
        [reference.chi_ss[0, 0, 0], reference.chi_cs[0, 0, 0],
         reference.chi_sc[0, 0, 0], reference.chi_cc[0, 0, 0]],
        [-1.9081848375005e-11 - 2.4043614344959e-13j,
         -4.1419996461365e-11 + 1.3207594099432e-13j,
         -4.1419991825958e-11 - 8.1815140174397e-14j,
         -8.9363859101550e-12 + 7.0338967135415e-14j],
        rtol=2e-12, atol=2e-23)
    assert reference.diagnostics["max_response_normalized_residual"] < 1e-12
    assert reference.diagnostics["max_response_trace_error"] < 1e-12
    assert reference.provenance["analysis_frequency"].startswith("independent")


def test_finite_seed_converges_quadratically_to_pump_reference():
    pump, atom, _delta, _deff = _fixture()
    base_seed = fwm.rabi_freq(8e-6, fwm.W_PROBE)
    delta = 2.0 * np.pi * np.array([-30e6, -8e6, 0.0, 12e6, 40e6])
    deff = 2.0 * np.pi * np.array([0.75e9, 0.9e9, 1.05e9])
    reference = fwm.pump_only_weak_response_reference(
        pump, pump, delta, deff, atom=atom)
    target = (reference.chi_ss[0], reference.chi_cs[0],
              reference.chi_sc[0], reference.chi_cc[0])
    scale = max(float(np.max(np.abs(value))) for value in target)

    fractions = np.array([0.3, 0.1, 0.03])
    errors = []
    for fraction in fractions:
        finite = fwm.chi_matrix_table(
            pump, pump, base_seed * fraction, base_seed * fraction,
            delta, deff, -1, atom=atom, n_f=3)
        errors.append(max(float(np.max(np.abs(a - b)))
                          for a, b in zip(finite, target)) / scale)
    slopes = np.diff(np.log(errors)) / np.diff(np.log(fractions))
    assert np.all((slopes > 1.95) & (slopes < 2.05))
    assert errors[-1] < 1.4e-6

    finite_nf1 = fwm.chi_matrix_table(
        pump, pump, base_seed * fractions[-1], base_seed * fractions[-1],
        delta, deff, -1, atom=atom, n_f=1)
    nf1_error = max(float(np.max(np.abs(a - b)))
                    for a, b in zip(finite_nf1, target)) / scale
    assert nf1_error > 4e-3


def test_dc_response_uses_full_stationary_projection_and_matches_derivative():
    pump, atom, _delta, deff = _fixture()
    H = fwm.pump_hamiltonian_at_Deff_zero(pump, pump)
    L = core.build_liouvillian(H, atom) - deff[0] * atom.S_v
    rho = core.steady_state_from_liouvillian(L, fwm.N_LEVELS)
    dH = np.zeros_like(H)
    fwm._add_static_drive(dH, fwm.G1, 1.0)
    fwm._add_static_drive(dH, fwm.G2, 1.0)
    source = core.comm_super(dH) @ rho.reshape(-1)
    response = core.trace_zero_liouvillian_response(
        L, source, 0.0, fwm.N_LEVELS).reshape(fwm.N_LEVELS, fwm.N_LEVELS)

    epsilon = 1e-4 * pump
    states = []
    for shifted in (pump + epsilon, pump - epsilon):
        shifted_L = core.build_liouvillian(
            fwm.pump_hamiltonian_at_Deff_zero(shifted, shifted), atom)
        shifted_L = shifted_L - deff[0] * atom.S_v
        states.append(core.steady_state_from_liouvillian(
            shifted_L, fwm.N_LEVELS))
    finite_difference = (states[0] - states[1]) / (2.0 * epsilon)
    np.testing.assert_allclose(response, finite_difference, rtol=2e-7, atol=2e-18)
    assert abs(np.trace(response)) < 1e-14

    multiple_null = np.diag([0.0, -2.0, -3.0, 0.0]).astype(complex)
    coherence = core.trace_zero_liouvillian_response(
        multiple_null, np.array([0.0, 1.0, 0.0, 0.0]), 0.0, 2)
    np.testing.assert_allclose(coherence, [0.0, 0.5, 0.0, 0.0], atol=1e-15)
    with pytest.raises(np.linalg.LinAlgError, match="incompatible"):
        core.trace_zero_liouvillian_response(
            multiple_null, np.array([1.0, 0.0, 0.0, -1.0]), 0.0, 2)


def test_trace_zero_source_validation_is_independent_for_every_rhs_column():
    liouvillian = -np.eye(4, dtype=complex)
    sources = np.array([
        [0.0, 1e-12],
        [1.0, 0.0],
        [0.0, 0.0],
        [0.0, 0.0],
    ], dtype=complex)
    with pytest.raises(ValueError, match="trace-free"):
        core.trace_zero_liouvillian_response(
            liouvillian, sources, 1.0, 2)


def test_trace_zero_response_accepts_one_frequency_per_liouvillian_batch():
    liouvillian = np.broadcast_to(
        np.diag([0.0, -2.0, -3.0, -4.0]), (2, 4, 4)).astype(complex)
    sources = np.zeros((2, 4), dtype=complex)
    sources[:, 1] = 1.0
    frequencies = np.array([0.0, 1.0])
    batched = core.trace_zero_liouvillian_response(
        liouvillian, sources, frequencies, 2)
    separate = np.stack([
        core.trace_zero_liouvillian_response(
            liouvillian[index], sources[index], frequencies[index], 2)
        for index in range(2)
    ])
    np.testing.assert_allclose(batched, separate, rtol=1e-14, atol=1e-15)


def test_pole_sum_matches_direct_resolvent_and_excludes_stationary_visibility():
    pump, atom, delta, deff = _fixture()
    for analysis in 2.0 * np.pi * np.array([0.0, 0.1e6, 4.0e6]):
        direct = fwm.pump_only_weak_response_reference(
            pump, pump, delta, deff, atom=atom,
            analysis_frequency_axis_rad_s=[analysis])
        poles = fwm.pump_only_pole_residue_reference(
            pump, pump, delta[0], deff[0], atom=atom,
            analysis_frequency_rad_s=analysis)
        np.testing.assert_allclose(
            poles["response"], direct.chi_matrix[0, 0, 0],
            rtol=2e-12, atol=2e-23)
        assert poles["stationary_residue_max"] < 1e-12
        assert np.count_nonzero(poles["stationary_mask"]) == 1
        assert np.all(np.isfinite(
            poles["half_widths_rad_s"][poles["nonstationary_mask"]]))


def test_pole_sum_rejects_nonzero_stationary_residue_at_any_source_scale():
    liouvillian = np.diag([0.0, -1.0]).astype(complex)
    stationary_readout = np.array([1.0, 0.0])
    for source_scale in (1e-2, 1e-11):
        with pytest.raises(np.linalg.LinAlgError, match="nonzero residue"):
            core.liouvillian_pole_residue_response(
                liouvillian,
                np.array([source_scale, 0.0]),
                stationary_readout,
                0.0,
            )


def test_pole_zero_mode_threshold_does_not_erase_a_physical_slow_mode():
    liouvillian = np.diag([0.0, -1e-2, -1.0, -1e10]).astype(complex)
    source = np.array([0.0, 1.0, 0.0, 0.0])
    readout = np.array([0.0, 1.0, 0.0, 0.0])
    poles = core.liouvillian_pole_residue_response(
        liouvillian, source, readout, 0.0)
    np.testing.assert_allclose(poles["response"], [[100.0]], rtol=1e-14)
    np.testing.assert_allclose(np.max(poles["visibility"]), 100.0, rtol=1e-14)
    assert np.count_nonzero(poles["stationary_mask"]) == 1


@pytest.mark.parametrize("bad_argument", ["source", "readout"])
def test_pole_sum_rejects_nonfinite_contractions(bad_argument):
    liouvillian = np.diag([0.0, -1.0]).astype(complex)
    source = np.array([0.0, 1.0])
    readout = np.array([0.0, 1.0])
    if bad_argument == "source":
        source[1] = np.nan
    else:
        readout[1] = np.nan
    with pytest.raises(ValueError, match="finite"):
        core.liouvillian_pole_residue_response(
            liouvillian, source, readout, 1.0)


def test_block_response_solver_matches_independent_dense_floquet_reference():
    pump, atom, delta, deff = _fixture()
    L0 = core.build_liouvillian(
        fwm.static_hamiltonian_at_Deff_zero(
            pump, pump, 0.0, delta[0], branch=-1), atom)
    Cp, Cm = fwm.sideband_template(pump, pump, 0.0, branch=-1)
    drives, _readouts = fwm._pump_nambu_operators(-1)
    zeros = np.zeros_like(drives)
    continued = core.floquet_linear_response_truncated(
        L0, Cp, Cm, fwm.seeded_sideband_beat(delta[0], -1), deff,
        atom.S_v, fwm.N_LEVELS, drives, zeros, zeros, n_f=3)
    direct = core.floquet_linear_response_direct(
        L0, Cp, Cm, fwm.seeded_sideband_beat(delta[0], -1), deff,
        atom.S_v, fwm.N_LEVELS, drives, zeros, zeros, n_f=3)
    np.testing.assert_allclose(continued[0], direct[0], rtol=2e-11, atol=2e-13)
    np.testing.assert_allclose(continued[1], direct[1], rtol=2e-11, atol=2e-22)
