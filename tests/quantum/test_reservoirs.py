"""Independent operator-action and physicality checks for the S0 foundation."""

import numpy as np
import pytest

from gabes import atoms, core
from gabes.quantum.reservoirs import (
    CollapseChannel, ExplicitReservoirs, audit_generator, choi_matrix,
    thermal_reset_channels,
)
from gabes.schemes.fwm import collisional_atom, thermal_transit_reset_superoperator


def test_complex_jumps_and_hamiltonian_match_direct_operator_action():
    rng = np.random.default_rng(142)
    n = 3
    raw_h = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))
    h = (raw_h+raw_h.conj().T)/2
    jumps = [rng.normal(size=(n, n))+1j*rng.normal(size=(n, n)) for _ in range(2)]
    model = ExplicitReservoirs(n, tuple(
        CollapseChannel(str(i), j, "test", "random explicit GKSL fixture")
        for i, j in enumerate(jumps)))
    generator = model.generator(h)
    # A non-Hermitian input tests the entire complex-linear map, including
    # transpose/conjugation conventions that a diagonal state would not expose.
    rho = rng.normal(size=(n, n))+1j*rng.normal(size=(n, n))
    expected = -1j*(h@rho-rho@h)
    for jump in jumps:
        q = jump.conj().T@jump
        expected += jump@rho@jump.conj().T - (q@rho+rho@q)/2
    np.testing.assert_allclose((generator@rho.reshape(-1)).reshape(n, n), expected,
                               rtol=2e-14, atol=2e-14)
    assert audit_generator(generator).passed


def test_choi_indexing_matches_matrix_unit_definition_without_symmetrizing():
    rng = np.random.default_rng(115)
    matrix = rng.normal(size=(9, 9))+1j*rng.normal(size=(9, 9))
    expected = np.zeros((9, 9), complex)
    for i in range(3):
        for j in range(3):
            e = np.zeros((3, 3)); e[i, j] = 1
            expected += np.kron(e, (matrix@e.reshape(-1)).reshape(3, 3))
    np.testing.assert_array_equal(choi_matrix(matrix), expected)


def test_radiative_factory_has_exact_existing_assembly_and_total_decay():
    atom = atoms.double_lambda_rb85(gamma_gg=0)
    model = ExplicitReservoirs.from_atom(atom, source="GABES D1 radiative data")
    np.testing.assert_allclose(model.dissipator(), atom.lindblad, rtol=0, atol=0)
    assert audit_generator(model.dissipator()).passed
    assert {c.kind for c in model.channels} == {"spontaneous_emission"}
    for e in atom.excited:
        rho = np.zeros((4, 4)); rho[e, e] = 1
        derivative = (model.dissipator()@rho.reshape(-1)).reshape(4, 4)
        assert np.trace(derivative) == pytest.approx(0, abs=1e-8)
        assert derivative[e, e] == pytest.approx(-sum(
            rate for i, _j, rate in atom.decay if i == e))


def test_reset_jumps_match_analytic_replacement_and_legacy_reset_for_all_matrices():
    rate = 3.2e5
    p = np.array([5/12, 7/12, 0., 0.])
    model = ExplicitReservoirs(4, thermal_reset_channels(
        rate, p, source="declared thermal replacement fixture"))
    generator = model.dissipator()
    expected = thermal_transit_reset_superoperator(rate, p)
    np.testing.assert_allclose(generator, expected, rtol=2e-15, atol=2e-10)
    rng = np.random.default_rng(103)
    rho = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
    np.testing.assert_allclose((generator@rho.reshape(-1)).reshape(4, 4),
                               rate*(np.diag(p)*np.trace(rho)-rho),
                               rtol=3e-15, atol=3e-10)
    assert audit_generator(generator).passed
    assert thermal_reset_channels(0, p, source="zero-rate fixture") == ()


def test_default_counterexample_and_thermal_setting_have_separate_verdicts():
    default = atoms.double_lambda_rb85()
    bad = audit_generator(default.lindblad)
    assert bad.trace_relative_residual < 1e-12
    assert bad.hermiticity_relative_residual < 1e-12
    assert bad.minimum_conditional_choi_eigenvalue_s_inverse < -3e5
    assert not bad.passed
    thermal = collisional_atom(394.15)
    assert audit_generator(thermal.lindblad).passed
    # Combined CP is insufficient for importing unaccounted reservoirs.
    for atom in (default, thermal):
        with pytest.raises(ValueError, match="legacy dephasing"):
            ExplicitReservoirs.from_atom(atom, source="legacy settings")


def test_import_rejects_a_cp_but_untracked_generator_edit():
    atom = atoms.double_lambda_rb85(gamma_gg=0)
    atom.lindblad = atom.lindblad + thermal_transit_reset_superoperator(1e5)
    assert audit_generator(atom.lindblad).passed
    with pytest.raises(ValueError, match="untracked"):
        ExplicitReservoirs.from_atom(atom, source="radiative data only")


def test_each_generator_condition_is_enforced_and_zero_generator_is_valid():
    assert audit_generator(np.zeros((4, 4))).passed
    assert not audit_generator(np.eye(4)).passed  # trace not conserved
    hp_bad = audit_generator(1j*np.eye(4))
    assert hp_bad.hermiticity_relative_residual > 0.1
    assert not hp_bad.passed
    transpose = np.eye(4)[[0, 2, 1, 3]]
    non_cp = audit_generator(transpose-np.eye(4))
    assert non_cp.trace_relative_residual == 0
    assert non_cp.hermiticity_relative_residual == 0
    assert not non_cp.passed


def test_hamiltonian_does_not_mask_negative_dissipative_eigenvalue():
    atom = atoms.double_lambda_rb85()
    h = np.diag([0., 1e11, 2e11, -3e11])
    assert not audit_generator(atom.lindblad+core.comm_super(h)).passed
    assert audit_generator(core.comm_super(h)).passed


@pytest.mark.parametrize("rate,p", [(-1, [1, 0]), (np.nan, [1, 0]),
                                    (1, [0.2, 0.2]), (1, [1.1, -0.1])])
def test_invalid_reset_model_is_rejected(rate, p):
    with pytest.raises(ValueError):
        thermal_reset_channels(rate, p, source="fixture")


def test_explicit_reservoir_objects_do_not_alias_mutable_input():
    jump = np.array([[0., 1.], [0., 0.]])
    channel = CollapseChannel("decay", jump, "radiative", "fixture")
    jump[:] = 0
    assert channel.operator[0, 1] == 1
    with pytest.raises(ValueError):
        channel.operator.setflags(write=True)
