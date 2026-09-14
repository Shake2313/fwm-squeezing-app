"""Microscopic optical theorem, normalization, ordered noise and propagation."""

from dataclasses import replace

import numpy as np
import pytest
from scipy.linalg import expm

from analysis.grand_challenge.reference.atomic_qrt import qrt_ordered_spectrum
from gabes import atoms, constants as c, observables
from gabes.fwm_quantum.field import reduced_field_pair, reduced_readout_operators
from gabes.fwm_quantum.model import reduced_pump_noise
from gabes.quantum.contracts import AnalysisFrequencyAxis, GeneratorFrequencyAxis, OpticalDetunings
from gabes.quantum.diffusion import stationary_atomic_noise
from gabes.quantum.reservoirs import CollapseChannel, ExplicitReservoirs, thermal_reset_channels
from gabes.quantum.traveling import (
    eliminate_atomic_noise, constant_segment, compose_segments, sum_independent_classes,
)
from gabes.schemes import fwm


@pytest.mark.parametrize("inverted", [False, True])
def test_two_level_passive_vacuum_and_inverted_quantum_limited_amplifier(inverted):
    gamma, resonance, g, linear_density = 2.4, 3.7, 0.7, 1.8
    lower = np.array([[0., 1.], [0., 0.]])
    jump = np.sqrt(gamma)*(lower.T if inverted else lower)
    atom = stationary_atomic_noise(np.diag([0., resonance]), ExplicitReservoirs(2, (
        CollapseChannel("reset", jump, "test", "analytic ground/excited reset"),)))
    axis = GeneratorFrequencyAxis(np.array([-2., 0., 3.7, 6.]), "two-level-lab")
    local = eliminate_atomic_noise(atom, [lower], [g], [1], axis,
                                   linear_density_m_inverse=linear_density, mode_labels=("signal",))
    denominator = gamma/2+1j*(resonance-axis.omega_rad_s)
    sign = 1 if inverted else -1
    expected_m = sign*linear_density*g*g/denominator
    expected_d = linear_density*g*g*gamma/np.abs(denominator)**2
    np.testing.assert_allclose(local.drift[:, 0, 0], expected_m, rtol=3e-15)
    np.testing.assert_allclose((local.noise_lesser if inverted else local.noise_greater)[:, 0, 0],
                               expected_d, rtol=3e-15)
    np.testing.assert_allclose(local.noise_greater if inverted else local.noise_lesser, 0, atol=2e-16)
    propagated = constant_segment(local, .61)
    power = np.exp(2*.61*expected_m.real)
    np.testing.assert_allclose(np.abs(propagated.transfer[:, 0, 0])**2, power, rtol=3e-15)
    np.testing.assert_allclose((propagated.noise_lesser if inverted else propagated.noise_greater)[:, 0, 0],
                               sign*(power-1), rtol=4e-15)
    cp, cm = propagated.vacuum_output()
    np.testing.assert_allclose(cp[:, 0, 0], power if inverted else np.ones_like(power), rtol=4e-15)
    np.testing.assert_allclose(cm[:, 0, 0], power-1 if inverted else np.zeros_like(power), atol=1e-15)


@pytest.fixture
def setup():
    pump = fwm.rabi_freq(.6, fwm.W_PUMP)
    optical = OpticalDetunings(2*np.pi*.9e9, -2*np.pi*8e6)
    axis = AnalysisFrequencyAxis.from_hz([-4e6, -1e6, 0, 1e6, 4e6])
    atom = reduced_pump_noise(pump, optical.one_photon_rad_s, transit_rate_s_inverse=c.GAMMA_GG)
    kwargs = dict(number_density_m3=1e18, uniform_area_m2=1.2e-7,
                  optical_omega_rad_s=[c.OMEGA_D1-c.OMEGA_HF, c.OMEGA_D1+c.OMEGA_HF],
                  effective_dipole_C_m=c.DIPOLE_D1/np.sqrt(12))
    return atom, optical, axis, kwargs


def test_response_and_photon_flux_transfer_match_existing_maxwell_reference(setup):
    atom, optical, axis, kwargs = setup
    pair = reduced_field_pair(atom, optical, axis, **kwargs)
    explicit = replace(atoms.double_lambda_rb85(gamma_gg=0), collapse_ops=tuple(
        channel.operator for channel in thermal_reset_channels(
            c.GAMMA_GG, [5/12, 7/12, 0, 0], source="same explicit parity fixture")))
    pump = fwm.rabi_freq(.6, fwm.W_PUMP)
    old = fwm.pump_only_weak_response_reference(
        pump, pump, [optical.two_photon_rad_s], [optical.one_photon_rad_s], atom=explicit,
        analysis_frequency_axis_rad_s=axis.omega_rad_s)
    wp, wc = kwargs["optical_omega_rad_s"]
    chi = [v[:, 0, 0] for v in (old.chi_ss, old.chi_sc, old.chi_cs, old.chi_cc)]
    _, _, field_t = observables.gain_from_chi(*chi, wp/c.C_LIGHT, wc/c.C_LIGHT,
                                           .0125, kwargs["number_density_m3"],
                                           dipole=c.DIPOLE_D1, line_strength=1/12)
    expected, _ = observables.canonical_transfer_from_field(
        field_t, wp, wc, kwargs["uniform_area_m2"], kwargs["uniform_area_m2"])
    actual = constant_segment(pair.main, .0125)
    np.testing.assert_allclose(actual.transfer, expected, rtol=2e-10, atol=1e-12)


def test_field_noise_matches_independent_atomic_qrt_projection(setup):
    atom, optical, axis, kwargs = setup
    pair = reduced_field_pair(atom, optical, axis, **kwargs)
    q = np.diag(observables.photon_flux_mode_matrix(*kwargs["optical_omega_rad_s"],
                                                  kwargs["uniform_area_m2"], kwargs["uniform_area_m2"]))
    g = kwargs["effective_dipole_C_m"]*q/(2*c.HBAR)
    w = np.einsum("jab,iba->ji", reduced_readout_operators(), atom.operators)
    c0 = -1j*np.diag([g[0], -g[1]])@w
    frequencies = pair.main.frequency_axis.omega_rad_s
    positive = qrt_ordered_spectrum(atom.generator, atom.stationary_state, atom.operators, frequencies)
    reversed_order = qrt_ordered_spectrum(atom.generator, atom.stationary_state, atom.operators,
                                         -frequencies).transpose(0, 2, 1)
    for source, actual in ((positive, pair.main.noise_greater), (reversed_order, pair.main.noise_lesser)):
        expected = kwargs["number_density_m3"]*kwargs["uniform_area_m2"]*(c0@source@c0.conj().T)
        error = np.linalg.norm(actual-expected, axis=(-2, -1))/np.linalg.norm(expected, axis=(-2, -1))
        assert error.max() < 2e-8


def test_density_dipole_and_common_area_normalization(setup):
    atom, optical, axis, kwargs = setup
    reference = reduced_field_pair(atom, optical, axis, **kwargs).main
    for key, factor, expected in (("number_density_m3", 2., 2.),
                                   ("uniform_area_m2", 3., 1.),
                                   ("effective_dipole_C_m", 2., 4.)):
        changed = reduced_field_pair(atom, optical, axis, **{**kwargs, key: kwargs[key]*factor}).main
        for name in ("drift", "noise_greater", "noise_lesser"):
            np.testing.assert_allclose(getattr(changed, name), expected*getattr(reference, name),
                                       rtol=3e-12, atol=1e-13)
    pieces = [reduced_field_pair(atom, optical, axis,
              **{**kwargs, "number_density_m3": kwargs["number_density_m3"]*fraction}).main
              for fraction in (.2, .3, .5)]
    combined = sum_independent_classes(pieces)
    for name in ("drift", "noise_greater", "noise_lesser"):
        np.testing.assert_allclose(getattr(combined, name), getattr(reference, name), rtol=3e-13)
    assert len(combined.reservoir_names) == 3*len(reference.reservoir_names)


def test_companion_and_both_noise_orderings_preserve_conjugation(setup):
    atom, optical, axis, kwargs = setup
    pair = reduced_field_pair(atom, optical, axis, phase_mismatch_rad_m=37., **kwargs)
    np.testing.assert_allclose(pair.companion.drift, pair.main.drift[::-1].conj(), rtol=3e-11, atol=1e-12)
    np.testing.assert_allclose(pair.companion.noise_greater, pair.main.noise_lesser[::-1].conj(), rtol=3e-11)
    np.testing.assert_allclose(pair.companion.noise_lesser, pair.main.noise_greater[::-1].conj(), rtol=3e-11)
    for local in (pair.main, pair.companion):
        assert local.audit()["passed"]
        segment = constant_segment(local, .0125)
        cp, cm = segment.vacuum_output()
        assert np.linalg.eigvalsh(cp).min() > -1e-12
        assert np.linalg.eigvalsh(cm).min() > -1e-12
        np.testing.assert_allclose(cp-cm, np.broadcast_to(np.diag(local.signs), cp.shape), atol=2e-12)


def test_constant_segment_matches_independent_quadrature_and_semigroup(setup):
    atom, optical, axis, kwargs = setup
    local = reduced_field_pair(atom, optical, axis, **kwargs).main
    whole = constant_segment(local, .0125)
    combined = compose_segments(constant_segment(local, .004), constant_segment(local, .0085))
    for name in ("transfer", "noise_greater", "noise_lesser"):
        np.testing.assert_allclose(getattr(whole, name), getattr(combined, name), rtol=5e-13, atol=2e-15)
    nodes, weights = np.polynomial.legendre.leggauss(32)
    for idx, drift in enumerate(local.drift):
        for source, noise in ((local.noise_greater, whole.noise_greater),
                              (local.noise_lesser, whole.noise_lesser)):
            expected = np.zeros((2, 2), complex)
            for node, weight in zip(nodes, weights):
                t = expm(drift*(node+1)*.0125/2)
                expected += weight*.0125/2*(t@source[idx]@t.conj().T)
            np.testing.assert_allclose(noise[idx], expected, rtol=5e-13, atol=1e-16)


def test_no_atom_zero_length_and_noncommuting_prescribed_segments(setup):
    atom, optical, axis, kwargs = setup
    local = reduced_field_pair(atom, optical, axis, **kwargs).main
    empty = reduced_field_pair(atom, optical, axis, **{**kwargs, "number_density_m3": 0.}).main
    for transfer in (constant_segment(local, 0.), constant_segment(empty, .0125)):
        np.testing.assert_array_equal(transfer.transfer, np.broadcast_to(np.eye(2), transfer.transfer.shape))
        np.testing.assert_array_equal(transfer.noise_greater, 0)
        np.testing.assert_array_equal(transfer.noise_lesser, 0)
    second = reduced_field_pair(atom, optical, axis, phase_mismatch_rad_m=1100., **kwargs).main
    a, b = constant_segment(local, .006), constant_segment(second, .007)
    ab, ba = compose_segments(a, b), compose_segments(b, a)
    assert np.linalg.norm(ab.transfer-ba.transfer) > .01
    assert np.linalg.norm(ab.noise_greater-ba.noise_greater) > 1e-4
    assert ab.audit()["passed"] and ba.audit()["passed"]


def test_bad_noise_coefficient_and_mismatched_coordinates_are_rejected(setup):
    atom, optical, axis, kwargs = setup
    local = reduced_field_pair(atom, optical, axis, **kwargs).main
    bad = replace(local, noise_greater_by_reservoir=local.noise_greater_by_reservoir/12,
                   noise_lesser_by_reservoir=local.noise_lesser_by_reservoir/12)
    assert not bad.audit()["passed"]
    with pytest.raises(ValueError, match="valid microscopic"):
        constant_segment(bad, .01)
    with pytest.raises(ValueError, match="identical"):
        sum_independent_classes([local, replace(local, mode_labels=("other", "conjugate"))])
    with pytest.raises(ValueError):
        reduced_field_pair(atom, optical, axis, **{**kwargs, "uniform_area_m2": -1})
    with pytest.raises(TypeError):
        reduced_field_pair(atom, optical, local.frequency_axis, **kwargs)
    with pytest.raises(ValueError):
        constant_segment(local, -.1)
    with pytest.raises(ValueError):
        local.drift.setflags(write=True)
