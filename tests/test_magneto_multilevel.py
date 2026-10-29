"""Physical limits of the optional multihyperfine coated-cell NMOR engine."""
import numpy as np
import pytest
from scipy.linalg import solve

from gabes import constants, species
from gabes.magneto_multilevel import (
    build_multilevel_atom, component_centroid_mhz, multilevel_response,
    velocity_quadrature, _optical_terms, _polarization, _ground_generator,
    _reservoir_states,
)


@pytest.mark.parametrize("isotope,line,ng,ne", [
    ("Rb85", "D1", 12, 12), ("Rb85", "D2", 12, 24),
    ("Rb87", "D1", 8, 8), ("Rb87", "D2", 8, 16),
])
def test_complete_dipole_emission_and_hyperfine_absorption(isotope, line, ng, ne):
    atom = build_multilevel_atom(isotope, line)
    iso = species.ISOTOPES[isotope]
    Je = iso.line(line)[0]
    assert len(atom.ground) == ng and len(atom.excited) == ne
    np.testing.assert_allclose(sum(C.T@C for C in atom.dipoles.values()),
                               np.eye(ne), atol=2e-15)
    for block in atom.ground_slices:
        F = atom.ground[block.start][0]
        expected, actual = [], []
        for Fe in species.f_values(iso.I, Je):
            columns = [i for i, (f, _) in enumerate(atom.excited) if f == Fe]
            actual.append(sum(np.sum(C[block][:, columns]**2) for C in atom.dipoles.values()))
            expected.append(species.line_strength(F, Fe, iso.I, iso.Jg, Je))
        np.testing.assert_allclose(np.array(actual)/sum(actual), expected, atol=1e-15)
    # Isotropic absorption per transverse polarization: D2 has twice D1 strength.
    strength = sum(np.sum(atom.dipoles[q]**2) for q in (-1, 1))/(2*ng)
    assert strength == pytest.approx((2*Je+1)/6, abs=2e-15)


def test_line_centroid_uses_strengths_and_declares_ground_component():
    atom = build_multilevel_atom("Rb85", "D1")
    upper = component_centroid_mhz(atom, 3)
    lower = component_centroid_mhz(atom, 2)
    assert lower-upper > 2900
    result = multilevel_response([0], isotope="Rb85", ground_f=3,
                                 intensity_mw_cm2=0, doppler=False)
    assert result["reference_centroid_mhz"] == upper
    assert "line-strength centroid" in result["detuning_reference"]


def test_optical_elimination_preserves_trace_and_hermiticity():
    atom = build_multilevel_atom("Rb87", "D2")
    L, _, _, _ = _optical_terms(atom, np.array([-100., 0., 200.]), 0.2, _polarization(13))
    scale = np.max(np.abs(L))
    np.testing.assert_allclose(np.einsum("i,vij->vj", atom.trace, L), 0, atol=scale*2e-14)
    rho = atom.thermal.astype(complex)
    evolved = L@rho
    for block in atom.reduced_slices:
        n = round(np.sqrt(block.stop-block.start))
        matrices = evolved[:, block].reshape(-1, n, n)
        np.testing.assert_allclose(matrices, matrices.conj().transpose(0, 2, 1), atol=scale*2e-14)


def test_shared_reservoir_schur_matches_independent_full_block_solve():
    atom = build_multilevel_atom("Rb87", "D1")
    optical, _, _, _ = _optical_terms(atom, np.array([-200., 100.]), 0.03, _polarization(0))
    light = optical + _ground_generator(atom, 1e-4, 3., 0.)
    dark = _ground_generator(atom, 1e-4, 1., 0.)
    w, gout, gin = np.array([0.3, 0.7]), 4000., 20.
    states, reservoir = _reservoir_states(light, dark, w, gout, gin,
                                         "complete", atom.trace, atom.thermal)
    M = atom.dimension
    A = np.zeros((3*M, 3*M), complex)
    for v in range(2):
        block = slice(v*M, (v+1)*M)
        A[block, block] = light[v]-gout*np.eye(M)
        A[block, 2*M:] = gout*np.eye(M)
        A[2*M:, block] = gin*w[v]*np.eye(M)
    A[2*M:, 2*M:] = dark-gin*np.eye(M)
    A[0] = 0
    A[0, 2*M:] = atom.trace
    rhs = np.zeros(3*M, complex); rhs[0] = 1
    direct = solve(A, rhs)
    np.testing.assert_allclose(states, direct[:2*M].reshape(2, M), atol=2e-12)
    np.testing.assert_allclose(reservoir, direct[2*M:], atol=2e-12)


def test_zero_return_retains_thermal_transit_relaxation():
    atom = build_multilevel_atom("Rb87", "D1")
    optical, _, _, _ = _optical_terms(atom, np.array([-400.]), .1, _polarization(0))
    light = optical+_ground_generator(atom, 0., 1., 0.)
    dark = _ground_generator(atom, 0., 1., 0.)
    states, thermal = _reservoir_states(light, dark, [1.], 1e5, 0., "complete", atom.trace, atom.thermal)
    residual = light[0]@states[0]+1e5*(atom.thermal-states[0])
    assert np.max(np.abs(residual)) < 1e-9
    np.testing.assert_allclose(thermal, atom.thermal)


def test_no_doppler_mixing_modes_are_identical_and_states_physical():
    args = dict(isotope="Rb87", line="D1", intensity_mw_cm2=.05,
                doppler=False, b_ut=[-1e-4, 0, 1e-4])
    shared = multilevel_response([-100., 150.], **args)
    separate = multilevel_response([-100., 150.], velocity_mixing="none", **args)
    np.testing.assert_allclose(shared["chi_probe"], separate["chi_probe"], rtol=2e-12)
    np.testing.assert_allclose(shared["rotation_rad"], separate["rotation_rad"], atol=2e-14)
    np.testing.assert_allclose(shared["rotation_rad"][0], -shared["rotation_rad"][-1], atol=2e-14)
    assert shared["diagnostics"]["min_ground_eigenvalue"] > -1e-12
    assert shared["diagnostics"]["trace_error"] < 1e-12


def test_zero_drive_matches_independent_linear_susceptibility_and_zero_density():
    atom = build_multilevel_atom("Rb85", "D1")
    detunings = np.array([-200., 0., 350.])
    density = 1e15
    result = multilevel_response(detunings, intensity_mw_cm2=0,
                                 density_m3=density, doppler=False)
    expected = np.zeros(3, complex)
    center = component_centroid_mhz(atom)
    for i, block in enumerate(atom.ground_slices):
        strengths = sum(np.sum(atom.dipoles[q][block]**2, axis=0) for q in (-1, 1))/2/len(atom.ground)
        delta = (center+detunings[:, None]-(atom.excited_energy_mhz-atom.ground_energy_mhz[i]))*species.MHZ
        expected += np.sum(strengths/(-delta-.5j*atom.gamma), axis=1)
    expected *= density*atom.dipole_sq/(constants.EPS_0*constants.HBAR)
    np.testing.assert_allclose(result["chi_probe"][0], expected, rtol=1e-12)
    assert np.all((result["transmission"] > 0) & (result["transmission"] <= 1))
    vacuum = multilevel_response([0], density_m3=0, doppler=False)
    assert vacuum["transmission"][0, 0] == 1
    assert vacuum["rotation_rad"][0, 0] == 0


def test_resonance_quadrature_refines_and_captures_maxwell_measure():
    atom = build_multilevel_atom("Rb85", "D1")
    x, w, sigma = velocity_quadrature(atom, component_centroid_mhz(atom), 293.15, order=6)
    assert len(x) > 100 and np.all(w > 0)
    assert w.sum() == pytest.approx(1)
    assert abs(w@x)/sigma < 2e-7
    assert (w@(x*x))/sigma**2 == pytest.approx(1, rel=2e-6)
    args = dict(b_ut=[-1e-5, 0., 1e-5], density_m3=4.5e15)
    coarse = multilevel_response([-600., 0., 600.], velocity_points=4, **args)
    fine = multilevel_response([-600., 0., 600.], velocity_points=8, **args)
    np.testing.assert_allclose(coarse["chi_probe"], fine["chi_probe"], rtol=.002, atol=1e-12)
    # A shared dark reservoir mixes pumping/probing classes; no forced optimum.
    np.testing.assert_allclose(coarse["rotation_rad"], fine["rotation_rad"], rtol=.01, atol=1e-11)
    assert abs(fine["rotation_rad"][0, 0]) > 5*abs(fine["rotation_rad"][0, 1])


@pytest.mark.parametrize("mixing", ["complete", "none"])
def test_analytic_coupled_derivative_matches_refined_field_differences(mixing):
    # These optical points were suppressed by the previous fixed B stencil.
    # An analytic response must retain them without interpolating missing data.
    detunings = [-100., 0., 300., 400.]
    args = dict(isotope="Rb85", line="D1", density_m3=4.5e15,
                wall_relaxation_hz=3., velocity_points=4,
                velocity_mixing=mixing)
    analytic = multilevel_response(detunings, b_ut=0., **args)
    step_ut = 1e-7
    finite = multilevel_response(detunings, b_ut=[-step_ut, step_ut],
                                 magnetic_derivative=False, **args)
    numerical = (finite["rotation_rad"][1]-finite["rotation_rad"][0])/(2e-6*step_ut)
    np.testing.assert_allclose(analytic["rotation_slope_rad_per_t"][0], numerical,
                               rtol=2e-5, atol=1e-2)
    assert np.all(analytic["slope_valid"])
    assert np.max(analytic["slope_relative_residual"]) < 1e-8
    assert "analytic" in analytic["slope_method"]


@pytest.mark.parametrize("gamma_out,gamma_in", [(0., 20.), (5e4, 0.)])
def test_analytic_derivative_in_transit_and_zero_return_limits(gamma_out, gamma_in):
    args = dict(isotope="Rb87", intensity_mw_cm2=.1, doppler=False,
                gamma_out_hz=gamma_out, gamma_in_hz=gamma_in,
                b_x_ut=.002, b_y_ut=.003)
    # The zero-return limit has a broad 50 kHz reset; a larger stencil avoids
    # cancellation of its very small rotation while remaining well resolved.
    bias, h = .001, (1e-4 if gamma_in == 0 else 1e-7)
    analytic = multilevel_response([-100.], b_ut=bias, **args)
    finite = multilevel_response([-100.], b_ut=[bias-h, bias+h],
                                 magnetic_derivative=False, **args)
    numerical = (finite["rotation_rad"][1, 0]-finite["rotation_rad"][0, 0])/(2e-6*h)
    assert analytic["rotation_slope_rad_per_t"][0, 0] == pytest.approx(numerical, rel=1e-5)
    assert analytic["slope_valid"][0, 0]


def test_transverse_field_preserves_hermiticity_and_circular_hanle_response():
    args = dict(isotope="Rb87", intensity_mw_cm2=.1, doppler=False,
                qwp_deg=45., density_m3=4.5e15,
                b_ut=[-.01, 0., .01])
    along = multilevel_response([0.], **args)
    transverse = multilevel_response([0.], b_x_ut=.001, b_y_ut=.002, **args)
    assert np.ptp(transverse["transmission"][:, 0]) > 1e-7
    assert np.max(np.abs(transverse["transmission"]-along["transmission"])) > 1e-7
    assert np.all(np.isnan(transverse["rotation_rad"]))
    assert np.all(np.isnan(transverse["rotation_slope_rad_per_t"]))
    assert np.all(np.isfinite(transverse["transmission_slope_per_t"]))
    assert np.all(transverse["slope_valid"])
    assert transverse["diagnostics"]["hermiticity_error"] < 1e-10
    assert transverse["diagnostics"]["min_ground_eigenvalue"] >= -1e-10
    bias, h = 0., 1e-7
    finite = multilevel_response([0.], b_ut=[bias-h, bias+h], b_x_ut=.001,
                                 b_y_ut=.002, magnetic_derivative=False,
                                 **{k:v for k,v in args.items() if k != "b_ut"})
    numerical = np.diff(finite["transmission"][:, 0])[0]/(2e-6*h)
    assert transverse["transmission_slope_per_t"][1, 0] == pytest.approx(numerical, rel=1e-4)


@pytest.mark.parametrize("qwp", [0., 90.])
def test_zero_field_linear_transmission_slope_is_exactly_zero_by_parity(qwp):
    result = multilevel_response([0., 300.], isotope="Rb85", doppler=False,
                                  qwp_deg=qwp, b_ut=[-1e-4, 0., 1e-4])
    np.testing.assert_array_equal(result["transmission_slope_per_t"][1], [0., 0.])
    np.testing.assert_array_equal(result["chi_probe_slope_per_t"][1].imag, [0., 0.])
    np.testing.assert_array_equal(result["transmission_slope_symmetry_zero"][1], [True, True])
    assert not np.any(result["transmission_slope_symmetry_zero"][[0, 2]])
    assert np.any(result["transmission_slope_per_t"][[0, 2]] != 0.)
    # A finite positive photon noise referred through exactly zero response
    # is unbounded field noise, not a spurious finite detectability claim.
    with np.errstate(divide="ignore"):
        field_asd = np.ones(2)/np.abs(result["transmission_slope_per_t"][1])
    assert np.all(np.isinf(field_asd))
    assert np.all(result["slope_valid"])


def test_transmission_parity_closure_does_not_zero_circular_response():
    result = multilevel_response([0.], isotope="Rb87", doppler=False,
                                  qwp_deg=45., b_ut=0., b_x_ut=.001)
    assert not result["transmission_slope_symmetry_zero"][0, 0]
    assert result["transmission_slope_per_t"][0, 0] != 0.


@pytest.mark.parametrize("changes", [
    {"intensity_mw_cm2": -1}, {"density_m3": -1}, {"velocity_points": 1},
    {"ground_relaxation_hz": 0}, {"velocity_mixing": "unknown"},
])
def test_invalid_settings_rejected(changes):
    with pytest.raises(ValueError):
        multilevel_response([0], doppler=False, **changes)
