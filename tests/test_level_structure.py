"""Shared level-structure pieces: one D[L] and one transit/thermal reset."""
import numpy as np

from gabes import atoms, core, species, zeeman
from gabes.quantum import reservoirs


def _random_density(n, seed):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    rho = a @ a.conj().T
    return rho / np.trace(rho)


def test_dissipator_is_lindblad_form_in_row_major_vec():
    rng = np.random.default_rng(3)
    L = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    rho = _random_density(4, 4)
    LdL = L.conj().T @ L
    expected = L @ rho @ L.conj().T - 0.5 * (LdL @ rho + rho @ LdL)
    got = (core.dissipator(L) @ rho.reshape(-1)).reshape(4, 4)
    np.testing.assert_allclose(got, expected, atol=1e-13)


def test_with_reset_adds_thermal_replacement_and_keeps_builder_attributes():
    base = zeeman.zeeman_manifold(1, 2)
    p = np.array([0.5, 0.3, 0.2])
    rate = 2 * np.pi * 1e4
    reset = base.with_reset(rate, p, name="reset")
    rho = _random_density(base.n_levels, 5)
    rho_th = np.zeros_like(rho)
    rho_th[base.ground, base.ground] = p
    got = ((reset.lindblad - base.lindblad) @ rho.reshape(-1)).reshape(rho.shape)
    np.testing.assert_allclose(got, rate * (rho_th * np.trace(rho) - rho), atol=1e-9)
    assert reset.name == "reset" and reset.couplings is base.couplings
    np.testing.assert_array_equal(reset.m_ground, base.m_ground)
    assert reset.decay[:len(base.decay)] == base.decay


def test_with_reset_skips_zero_weight_targets():
    atom = species.build_manifold(species.RB87, "D1").atom
    reset = atom.with_reset(1.0, [0.0, 1.0])
    added = reset.decay[len(atom.decay):]
    assert {g for _s, g, _r in added} == {atom.ground[1]}
    assert len(added) == atom.n_levels


def test_reset_matches_explicit_quantum_reservoir_channels():
    atom = atoms.lambda3()
    rate = 3.7e5
    reset = atom.with_reset(rate, [0.25, 0.75])
    explicit = reservoirs.ExplicitReservoirs(3, reservoirs.thermal_reset_channels(
        rate, [0.25, 0.75, 0.0], source="test")).dissipator()
    np.testing.assert_allclose(reset.lindblad - atom.lindblad, explicit, atol=1e-9)
