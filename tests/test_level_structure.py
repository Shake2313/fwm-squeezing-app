"""Shared level-structure pieces: one D[L], one transit/thermal reset, and the
relations of every level builder to the signed dipoles of `species.level_structure`."""
import math

import numpy as np
import pytest

from gabes import atoms, constants, core, species, zeeman
from gabes.magneto_multilevel import build_multilevel_atom
from gabes.quantum import reservoirs

LINES = [(iso, line) for iso in ("Rb85", "Rb87", "Cs133") for line in ("D1", "D2")]


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


def _block(structure, fg, fe):
    g = [i for i, (f, _) in enumerate(structure.ground) if f == fg]
    e = [i for i, (f, _) in enumerate(structure.excited) if f == fe]
    return {q: structure.dipoles[q][np.ix_(g, e)] for q in (-1, 0, 1)}


def _block_sum(structure, fg, fe):
    return sum(float(np.sum(c ** 2)) for c in _block(structure, fg, fe).values())


@pytest.mark.parametrize("isotope, line", LINES)
def test_level_structure_dipoles_are_normalized_and_read_only(isotope, line):
    s = species.level_structure(isotope, line)
    assert s is species.level_structure(isotope, line)
    total = sum(c.T @ c for c in s.dipoles.values())
    np.testing.assert_allclose(total, np.eye(len(s.excited)), atol=1e-12)
    for array in (*s.dipoles.values(), s.ground_energy_mhz, s.excited_energy_mhz):
        assert not array.flags.writeable


@pytest.mark.parametrize("isotope, line", LINES[:4])
def test_multilevel_atom_shares_the_level_structure(isotope, line):
    s = species.level_structure(isotope, line)
    atom = build_multilevel_atom(isotope, line)
    assert (atom.ground, atom.excited) == (s.ground, s.excited)
    assert all(atom.dipoles[q] is s.dipoles[q] for q in (-1, 0, 1))
    assert atom.dipoles is not s.dipoles
    assert (atom.gamma, atom.frequency_hz, atom.dipole_sq, atom.mass) == (
        s.gamma, s.frequency_hz, s.dipole_sq, s.mass)


@pytest.mark.parametrize("isotope, line", LINES)
def test_f_resolved_manifold_branching_is_dipole_block_sums(isotope, line):
    # build_manifold keeps the analytic 6j strengths; the dipole block sums give
    # the same CG branching to rounding, and T·(2Je+1)/(2Jg+1) as block sums.
    iso = species.ISOTOPES[isotope]
    s = species.level_structure(isotope, line)
    man = species.build_manifold(iso, line)
    ratio = (2 * man.Je + 1) / (2 * man.Jg + 1)
    for fg in man.Fg:
        for fe in man.Fe:
            T = (2 * fg + 1) * species.line_strength(fg, fe, iso.I, man.Jg, man.Je)
            assert math.isclose(_block_sum(s, int(fg), int(fe)), ratio * T,
                                rel_tol=1e-14, abs_tol=1e-15)
    ng = len(man.Fg)
    for e, g, rate in man.atom.decay:
        fe, fg = int(man.Fe[e - ng]), int(man.Fg[g])
        branch = _block_sum(s, fg, fe) / (2 * fe + 1)
        assert math.isclose(rate, man.gamma * branch, rel_tol=1e-14)
    assert man.gamma == s.gamma


def test_double_lambda_is_the_rb85_d1_manifold_with_the_legacy_gamma():
    lumped = atoms.double_lambda_rb85().decay
    manifold = species.build_manifold(species.RB85, "D1").atom.decay
    assert [c[:2] for c in lumped] == [c[:2] for c in manifold]
    legacy = constants.GAMMA / (species.RB85.gamma_D1_mhz * species.MHZ)
    for (_e, _g, a), (_e2, _g2, b) in zip(lumped, manifold):
        assert math.isclose(a / b, legacy, rel_tol=1e-14)


def test_rb85_zeeman_manifold_is_the_dipoles_in_an_excited_f2_gauge():
    # Emission is √Γ·C_q; the drive is block-fitted to 3·C_F² per ground
    # sublevel, |c/C_q|² = (2Fg+1)/3, not the Wigner-Eckart amplitude. The
    # F'=2 states carry the opposite phase, which changes no observable.
    z = zeeman.rb85_d1_double_lambda_zeeman()
    s = species.level_structure("Rb85", "D1")
    ng = len(s.ground)
    assert [(int(f), int(m)) for f, m in zip(z.level_F, z.level_m)] == list(s.ground + s.excited)
    gauge = np.array([-1.0 if f == 2 else 1.0 for f, _ in s.excited])
    for qi, q in enumerate((-1, 0, 1)):
        expected = math.sqrt(constants.GAMMA) * s.dipoles[q] * gauge
        np.testing.assert_allclose(z.emission_ops[qi][:ng, ng:], expected,
                                   rtol=1e-14, atol=1e-15)
        assert len(z.couplings[q]) == np.count_nonzero(s.dipoles[q])
        for g, e, c in z.couplings[q]:
            fg = s.ground[g][0]
            d = s.dipoles[q][g, e - ng] * gauge[e - ng]
            assert math.isclose(c, math.sqrt((2 * fg + 1) / 3) * d, rel_tol=1e-14)


@pytest.mark.parametrize("fg, fe", [(1, 1), (1, 2), (2, 1), (2, 2)])
def test_single_pair_zeeman_is_a_normalized_dipole_block(fg, fe):
    zm = zeeman.zeeman_manifold(fg, fe)
    s = species.level_structure("Rb87", "D1")
    scale = math.sqrt(_block_sum(s, fg, fe) / (2 * fe + 1))
    ng = 2 * fg + 1
    block = _block(s, fg, fe)
    mine = [zm.emission_ops[qi][:ng, ng:] / math.sqrt(constants.GAMMA) for qi in range(3)]
    sign = np.sign(sum(np.sum(m.real * block[q]) for m, q in zip(mine, (-1, 0, 1))))
    for m, q in zip(mine, (-1, 0, 1)):  # one overall phase per closed pair
        np.testing.assert_allclose(m, sign * block[q] / scale, rtol=1e-14, atol=1e-15)
