"""Exact-model and shared-atom contracts for the magneto solver speed patch."""
import json

import numpy as np
import pytest

from gabes import magneto_multilevel as current
from analysis.magneto_speed_20261030.validation import (
    _json, array_check, branch_cases, compare_backend, compare_scheme_case,
    load_baselines, small_scheme_cases,
)


def test_frozen_scheme_is_bound_to_frozen_backend_only():
    import gabes
    before, scheme = load_baselines()
    assert scheme.magneto_multilevel is before
    assert before is not current
    assert gabes.magneto_multilevel is current


@pytest.mark.parametrize("name,kwargs", list(branch_cases().items()))
def test_fast_solver_matches_independent_frozen_complex_reference(name, kwargs):
    before, _ = load_baselines()
    reference = before.multilevel_response([-137.], **kwargs)
    result = current.multilevel_response([-137.], **kwargs)
    report = compare_backend(reference, result)
    assert report["passed"], name + ": " + json.dumps(_json(report), ensure_ascii=False)


@pytest.mark.parametrize("isotope,line", [("Rb85", "D1"), ("Rb85", "D2"),
                                         ("Rb87", "D1"), ("Rb87", "D2")])
def test_fast_calls_leave_shared_multilevel_atomic_arrays_unchanged(isotope, line):
    atom = current.build_multilevel_atom(isotope, line)
    arrays = {name: getattr(atom, name) for name in
              ("thermal", "trace", "ground_energy_mhz", "excited_energy_mhz")}
    arrays.update({f"dipole_{q}": value for q, value in atom.dipoles.items()})
    saved = {name: value.copy() for name, value in arrays.items()}
    current.multilevel_response([-137., 49.], isotope=isotope, line=line,
        intensity_mw_cm2=.08, density_m3=4.5e15, doppler=False,
        b_ut=[-.0002, 0., .0003], qwp_deg=13., b_x_ut=.0001, b_y_ut=-.0002)
    assert current.build_multilevel_atom(isotope, line) is atom
    for name, value in arrays.items():
        np.testing.assert_array_equal(value, saved[name], err_msg=f"mutated shared {isotope} {line} {name}")


def test_retained_complex_mode_matches_frozen_backend():
    before, _ = load_baselines()
    kwargs = dict(isotope="Rb87", line="D2", intensity_mw_cm2=.04,
                  doppler=False, b_ut=[0., .0004], b_x_ut=.0003, b_y_ut=.0002)
    reference = before.multilevel_response([-100.], **kwargs)
    result = current.multilevel_response([-100.], solver_basis="complex", **kwargs)
    assert compare_backend(reference, result)["passed"]


@pytest.mark.parametrize("name,params", list(small_scheme_cases().items()))
def test_physical_readout_snr_masks_and_best_detuning_match_frozen_scheme(name, params):
    report = compare_scheme_case(params)
    assert report["passed"], name + ": " + json.dumps(_json(report), ensure_ascii=False)


def test_comparison_uses_zero_floor_but_does_not_waive_availability_changes():
    assert array_check([0.], [1e-16], atol=1e-14)["passed"]
    report = array_check([None, 0.], [1., 0.], atol=1e-14)
    assert not report["passed"] and report["finite_mask_mismatches"] == 1


@pytest.mark.parametrize("isotope", ["Rb85", "Rb87"])
@pytest.mark.parametrize("even", [False, True])
def test_sparse_coordinates_match_independent_hermitian_basis_and_invariant_sector(isotope, even):
    from scipy.linalg import block_diag
    from gabes import core
    from gabes.magneto_coordinates import ground_coordinates

    atom = current.build_multilevel_atom(isotope, "D1")
    fs = tuple(atom.ground[block.start][0] for block in atom.ground_slices)
    blocks = []
    for f in fs:
        n = 2*f+1
        basis = core.hermitian_basis(n)
        if even:
            odd = ((np.arange(n)[:, None]-np.arange(n)[None, :]) % 2 != 0).ravel()
            basis = basis[:, np.all(basis[odd] == 0., axis=0)]
        blocks.append(basis)
    transform = block_diag(*blocks)
    coordinates = ground_coordinates(fs, even)
    expected_size = ({"Rb85": 38, "Rb87": 18}[isotope] if even else atom.dimension)
    assert coordinates.dimension == expected_size
    np.testing.assert_allclose(coordinates.expand(np.eye(expected_size)), transform.T, atol=3e-16)
    # Elliptical sigma input preserves parity; transverse fields require all
    # coherences. Compare the transform to the independent engine basis.
    optical, _, _, _ = current._optical_terms(atom, np.array([-137., 49.]), .08,
                                              current._polarization(13.))
    transverse = (0., 0.) if even else (.002, -.003)
    generator = optical+current._ground_generator(atom, .001, 3., 7., *transverse)
    expected = transform.conj().T @ generator @ transform
    scale = np.max(np.abs(generator))
    np.testing.assert_allclose(expected.imag, 0., atol=scale*3e-15)
    np.testing.assert_allclose(coordinates.generator(generator), expected.real,
                               rtol=2e-14, atol=scale*3e-15)
    projection = transform @ transform.conj().T
    leakage = (np.eye(atom.dimension)-projection) @ generator @ transform
    np.testing.assert_allclose(leakage, 0., atol=scale*3e-15)
