"""Exact, reusable angular operators for weak-excitation magneto optics.

Only atom/polarization-dependent coefficients are cached. Resolvents, velocity
grids, intensity and assembled Liouvillians remain local to each calculation.
"""
from collections import OrderedDict
from dataclasses import dataclass
import hashlib
import threading

import numpy as np


_MAX_COEFFICIENTS = 8
_COEFFICIENTS = OrderedDict()
_CACHE_LOCK = threading.RLock()
_HITS = 0
_MISSES = 0


@dataclass(frozen=True)
class OpticalBlock:
    excited_indices: np.ndarray
    rows: np.ndarray
    columns: np.ndarray
    coefficients: np.ndarray


@dataclass(frozen=True)
class OpticalCoefficients:
    # Retaining the atom also prevents reuse of its Python id while cached.
    atom: object
    polarization: tuple
    blocks: tuple
    coefficient_bytes: int


def coefficient_cache_clear():
    """Clear static coefficients and counters without changing any atom."""
    global _HITS, _MISSES
    with _CACHE_LOCK:
        _COEFFICIENTS.clear()
        _HITS = _MISSES = 0


def coefficient_cache_info():
    with _CACHE_LOCK:
        return dict(hits=_HITS, misses=_MISSES, maxsize=_MAX_COEFFICIENTS,
                    currsize=len(_COEFFICIENTS),
                    coefficient_bytes=sum(value.coefficient_bytes for value in _COEFFICIENTS.values()))


def _atom_fingerprint(atom):
    """Detect in-place edits to the small mutable atomic inputs as well as ids."""
    digest = hashlib.blake2b(digest_size=16)
    digest.update(repr((atom.ground, atom.excited, float(atom.gamma),
                       tuple((s.start, s.stop, s.step) for s in atom.ground_slices),
                       tuple((s.start, s.stop, s.step) for s in atom.reduced_slices))).encode())
    for value in (atom.excited_energy_mhz, *(atom.dipoles[q] for q in (-1, 0, 1))):
        array = np.asarray(value)
        digest.update(repr((array.dtype.str, array.shape)).encode())
        digest.update(array.tobytes(order="C"))
    return digest.digest()


def _immutable(value, dtype=None):
    result = np.asarray(value, dtype=dtype)
    result.flags.writeable = False
    return result


def _prepare(atom, polarization, polarization_key):
    # The original effective-operator convention uses C_q.T here. Retain the
    # complex polarization and every coherent virtual-excited-F path.
    excitation = sum(polarization[q]*atom.dipoles[q].T for q in (-1, 1))
    unit_coupling = -0.5*excitation
    blocks = []
    for gs, rs in zip(atom.ground_slices, atom.reduced_slices):
        n = gs.stop-gs.start
        groups = [np.flatnonzero(atom.excited_energy_mhz == energy)
                  for energy in np.unique(atom.excited_energy_mhz)]
        groups = [indices for indices in groups
                  if np.any(unit_coupling[indices, gs] != 0)]
        count = len(groups)
        # Pair terms carry R_f R_h*, followed by Re(R_f) light-shift terms.
        coefficients = np.zeros((count*count+count, atom.dimension, n*n), complex)
        loss = np.zeros((count, count, n, n), complex)
        for target_gs, target_rs in zip(atom.ground_slices, atom.reduced_slices):
            target_n = target_gs.stop-target_gs.start
            for q in (-1, 0, 1):
                jumps = np.array([np.sqrt(atom.gamma)*atom.dipoles[q][target_gs][:, indices]
                                  @ unit_coupling[indices, gs] for indices in groups])
                if not count:
                    continue
                gain = np.einsum("fai,hbj->fhabij", jumps, jumps.conj())
                coefficients[:count*count, target_rs, :] += gain.reshape(
                    count*count, target_n*target_n, n*n)
                # The coefficient of R_f R_h* in K=ΣJ†J is J_h† J_f.
                # Full decay completeness makes f!=h cancel analytically, but
                # retaining this static form also preserves supplied dipoles
                # with roundoff or intentional changes to their decay metric.
                loss += np.einsum("hai,faj->fhij", jumps.conj(), jumps)
        eye = np.eye(n)
        for f in range(count):
            for h in range(count):
                K = loss[f, h]
                coefficients[f*count+h, rs, :] -= 0.5*(np.kron(K, eye)+np.kron(eye, K.T))
            V = unit_coupling[groups[f], gs]
            W = V.conj().T@V
            # H=-Ω² Re(R_f) W, hence -i[H,ρ]=+iΩ² Re(R_f)[W,ρ].
            coefficients[count*count+f, rs, :] = 1j*(np.kron(W, eye)-np.kron(eye, W.T))
        # Contract only the exact structural union, never threshold small
        # amplitudes. Source-F column blocks are disjoint in the full output.
        rows, columns = np.nonzero(np.any(coefficients != 0, axis=0))
        blocks.append(OpticalBlock(
            _immutable([indices[0] for indices in groups], int),
            _immutable(rows), _immutable(columns+rs.start),
            _immutable(coefficients[:, rows, columns])))
    blocks = tuple(blocks)
    memory = sum(array.nbytes for block in blocks for array in
                 (block.excited_indices, block.rows, block.columns, block.coefficients))
    return OpticalCoefficients(atom, polarization_key, blocks, memory)


def prepare_optical_coefficients(atom, polarization):
    """Return a bounded, immutable atom/polarization-specific angular cache.

    Identity guards builder cache clears; a content fingerprint also invalidates
    changes to dipoles, excited energies, decay rate or level/block ordering.
    No velocity-dependent tensor is stored in an entry.
    """
    global _HITS, _MISSES
    polarization_key = tuple(complex(polarization[q]) for q in (-1, 1))
    key = (id(atom), _atom_fingerprint(atom), polarization_key)
    with _CACHE_LOCK:
        cached = _COEFFICIENTS.get(key)
        if cached is not None and cached.atom is atom:
            _COEFFICIENTS.move_to_end(key)
            _HITS += 1
            return cached
        _MISSES += 1
    prepared = _prepare(atom, polarization, polarization_key)
    with _CACHE_LOCK:
        _COEFFICIENTS[key] = prepared
        _COEFFICIENTS.move_to_end(key)
        while len(_COEFFICIENTS) > _MAX_COEFFICIENTS:
            _COEFFICIENTS.popitem(last=False)
    return prepared


def assemble_optical_generator(atom, resolvents, polarization, omega_scale_squared):
    """Assemble the same full complex generator using fixed angular operators."""
    nv = len(resolvents[0])
    result = np.zeros((nv, atom.dimension, atom.dimension), complex)
    if omega_scale_squared == 0:
        return result
    prepared = prepare_optical_coefficients(atom, polarization)
    for block, resolvent in zip(prepared.blocks, resolvents):
        r = resolvent[:, block.excited_indices]
        count = r.shape[-1]
        if not count:
            continue
        pairs = (r[:, :, None]*r[:, None, :].conj()).reshape(nv, count*count)
        weights = np.concatenate((pairs, r.real), axis=1)
        result[:, block.rows, block.columns] = omega_scale_squared*(weights@block.coefficients)
    return result
