"""Exact optical assembly and bounded immutable angular-template reuse."""
from dataclasses import fields, is_dataclass, replace
from functools import lru_cache
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

from gabes import magneto_multilevel as current
from gabes import magneto_optical as optical


@lru_cache(maxsize=1)
def _reference():
    path = (Path(__file__).resolve().parents[1]/"analysis"/
            "magneto_optical_assembly_20261030"/"before_magneto_multilevel.py")
    name = "gabes._optical_assembly_reference"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


ATOMS = [("Rb85", "D1"), ("Rb85", "D2"), ("Rb87", "D1"), ("Rb87", "D2")]


def _assert_terms_equal(actual, reference):
    scale = max(float(np.max(np.abs(reference[0]))), 1.)
    np.testing.assert_allclose(actual[0], reference[0], rtol=3e-13, atol=scale*3e-14)
    for name, got, want in zip(("eliminated", "resolvents"), actual[1:3], reference[1:3]):
        assert len(got) == len(want), name
        for left, right in zip(got, want):
            np.testing.assert_array_equal(left, right, err_msg=name)
    np.testing.assert_array_equal(actual[3], reference[3])


@pytest.mark.parametrize("isotope,line", ATOMS)
@pytest.mark.parametrize("qwp", [0., 17.5, 45.])
@pytest.mark.parametrize("intensity", [0., .01, 1.])
def test_prepared_assembly_matches_frozen_full_jump_construction(isotope, line, qwp, intensity):
    atom = current.build_multilevel_atom(isotope, line)
    center = current.component_centroid_mhz(atom)
    # Both ground components and excited-hyperfine interference are present.
    lasers = center+np.array([-3500., -137., 0., 49., 3500.])
    pol = current._polarization(qwp)
    expected = _reference()._optical_terms(atom, lasers, intensity, pol)
    actual = current._optical_terms(atom, lasers, intensity, pol)
    _assert_terms_equal(actual, expected)


@pytest.mark.parametrize("isotope,line", ATOMS)
def test_phase_and_intensity_scaling_preserve_the_full_generator(isotope, line):
    atom = current.build_multilevel_atom(isotope, line)
    lasers = np.array([-2300., -10., 130.])
    pol = current._polarization(17.5)
    phased = {q: value*np.exp(.37j) for q, value in pol.items()}
    a = current._optical_terms(atom, lasers, .08, pol)
    b = current._optical_terms(atom, lasers, .16, phased)
    _assert_terms_equal(b, _reference()._optical_terms(atom, lasers, .16, phased))
    np.testing.assert_allclose(b[0], 2*a[0], rtol=3e-13, atol=np.max(abs(a[0]))*3e-14)


@pytest.mark.parametrize("isotope,line", ATOMS)
def test_generator_remains_trace_and_hermiticity_preserving(isotope, line):
    atom = current.build_multilevel_atom(isotope, line)
    generator, _, _, _ = current._optical_terms(
        atom, np.array([-137., 0., 130.]), .1, current._polarization(17.5))
    scale = np.max(abs(generator))
    np.testing.assert_allclose(atom.trace @ generator, 0., atol=scale*3e-14)
    rng = np.random.default_rng(314)
    rho = np.empty(atom.dimension, complex)
    for block in atom.reduced_slices:
        n = round(np.sqrt(block.stop-block.start))
        sample = rng.normal(size=(n, n))+1j*rng.normal(size=(n, n))
        rho[block] = (sample @ sample.conj().T).ravel()
    rho /= atom.trace @ rho
    evolved = generator @ rho
    for block in atom.reduced_slices:
        n = round(np.sqrt(block.stop-block.start))
        matrix = evolved[:, block].reshape(-1, n, n)
        np.testing.assert_allclose(matrix, matrix.conj().swapaxes(-1, -2), atol=scale*3e-14)


def _arrays(value):
    if isinstance(value, np.ndarray):
        yield value
    elif is_dataclass(value):
        for item in fields(value):
            yield from _arrays(getattr(value, item.name))
    elif isinstance(value, (tuple, list)):
        for item in value:
            yield from _arrays(item)


def test_templates_reuse_across_spectra_and_intensities_without_mutating_atoms():
    optical.coefficient_cache_clear()
    atom = current.build_multilevel_atom("Rb85", "D2")
    arrays = [atom.trace, atom.thermal, atom.ground_energy_mhz, atom.excited_energy_mhz,
              *atom.dipoles.values()]
    saved = [(a.copy(), a.flags.writeable) for a in arrays]
    pol = current._polarization(13.)
    prepared = optical.prepare_optical_coefficients(atom, pol)
    assert optical.prepare_optical_coefficients(atom, pol) is prepared
    # The entry retains the shared atom identity, without freezing its arrays.
    for a in _arrays(prepared.blocks):
        assert not a.flags.writeable
    current._optical_terms(atom, np.linspace(-300., 300., 180), .1, pol)
    current._optical_terms(atom, np.linspace(-900., 600., 97), 1., pol)
    info = optical.coefficient_cache_info()
    assert info["currsize"] == 1 and info["misses"] == 1
    assert info["hits"] >= 3
    for array, (before, writable) in zip(arrays, saved):
        np.testing.assert_array_equal(array, before)
        assert array.flags.writeable == writable


def test_atom_identity_and_content_changes_invalidate_templates():
    optical.coefficient_cache_clear()
    original = current.build_multilevel_atom("Rb87", "D1")
    atom = replace(original, dipoles={q: value.copy() for q, value in original.dipoles.items()})
    pol = current._polarization(17.5)
    old = optical.prepare_optical_coefficients(atom, pol)
    # Only pi emission changes here: the incident field still has sigma+/−.
    atom.dipoles[0] *= 1.01
    assert optical.prepare_optical_coefficients(atom, pol) is not old
    lasers = np.array([-2300., -137., 130.])
    _assert_terms_equal(current._optical_terms(atom, lasers, .08, pol),
                        _reference()._optical_terms(atom, lasers, .08, pol))
    new_atom = replace(atom, dipoles={q: value.copy() for q, value in atom.dipoles.items()})
    assert optical.prepare_optical_coefficients(new_atom, pol) is not old


def test_coefficient_cache_has_a_fixed_memory_bound():
    optical.coefficient_cache_clear()
    atom = current.build_multilevel_atom("Rb85", "D2")
    limit = optical.coefficient_cache_info()["maxsize"]
    for angle in range(limit+3):
        optical.prepare_optical_coefficients(atom, current._polarization(angle+.25))
    info = optical.coefficient_cache_info()
    assert info["currsize"] == limit
    assert 0 < info["coefficient_bytes"] < 64*1024**2
