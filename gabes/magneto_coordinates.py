"""Exact sparse Hermitian coordinates for secular ground-F density blocks."""
from dataclasses import dataclass
from functools import lru_cache

import numpy as np


@dataclass(frozen=True)
class GroundCoordinates:
    full_dimension: int
    first: np.ndarray
    second: np.ndarray
    first_coefficient: np.ndarray
    second_coefficient: np.ndarray
    imaginary: np.ndarray
    row_scale: np.ndarray
    diagonal_indices: np.ndarray
    diagonal_coordinates: np.ndarray
    upper_indices: np.ndarray
    lower_indices: np.ndarray
    symmetric_coordinates: np.ndarray
    antisymmetric_coordinates: np.ndarray

    @property
    def dimension(self):
        return len(self.first)

    def generator(self, value):
        """U† L U, using the two nonzero entries of each Hermitian basis vector.

        The retained generators preserve Hermiticity, so the output is real.
        This sparse coordinate change avoids a dense cubic basis multiplication
        for every velocity. Its signs match core.hermitian_basis (row-major).
        """
        action = (value[..., self.first[:, None], self.first[None, :]]
                  * self.first_coefficient[None, :]
                  + value[..., self.first[:, None], self.second[None, :]]
                  * self.second_coefficient[None, :])
        result = action.real * self.row_scale[:, None]
        result[..., self.imaginary, :] = (action[..., self.imaginary, :].imag
                                         * self.row_scale[self.imaginary, None])
        return result

    def pack(self, state):
        """Coordinates of a Hermitian vector, including the thermal source."""
        value = state[..., self.first]
        result = value.real * self.row_scale
        result[..., self.imaginary] = (value[..., self.imaginary].imag
                                       * self.row_scale[self.imaginary])
        return result

    def expand(self, state):
        """Restore every original density entry for readout and validation."""
        if state is None:
            return None
        result = np.zeros(state.shape[:-1] + (self.full_dimension,), complex)
        result[..., self.diagonal_indices] = state[..., self.diagonal_coordinates]
        upper = (state[..., self.symmetric_coordinates]
                 + 1j*state[..., self.antisymmetric_coordinates])/np.sqrt(2.)
        result[..., self.upper_indices] = upper
        result[..., self.lower_indices] = upper.conj()
        return result


@lru_cache(maxsize=8)
def ground_coordinates(ground_f_values, even_coherences=False):
    """Blockwise Hermitian basis, optionally restricted to even m_i-m_j.

    A longitudinal field and σ± input preserve the even sector for any input
    ellipticity. A transverse field does not; callers must use the full basis.
    No angular dipoles or physical atom arrays are modified by this transform.
    """
    first, second, c1, c2, imaginary, scale = [], [], [], [], [], []
    diagonal, diagonal_c, upper, lower, symmetric, antisymmetric = [], [], [], [], [], []
    offset = 0
    for f in ground_f_values:
        n = 2*f+1
        for i in range(n):
            index = offset+i*n+i
            diagonal.append(index)
            diagonal_c.append(len(first))
            first.append(index); second.append(index)
            c1.append(1.); c2.append(0.); imaginary.append(False); scale.append(1.)
        for i in range(n):
            for j in range(i+1, n):
                if even_coherences and (i-j) % 2:
                    continue
                a, b = offset+i*n+j, offset+j*n+i
                upper.append(a); lower.append(b)
                symmetric.append(len(first)); antisymmetric.append(len(first)+1)
                first.extend((a, a)); second.extend((b, b))
                c1.extend((1/np.sqrt(2.), 1j/np.sqrt(2.)))
                c2.extend((1/np.sqrt(2.), -1j/np.sqrt(2.)))
                imaginary.extend((False, True)); scale.extend((np.sqrt(2.), np.sqrt(2.)))
        offset += n*n
    arrays = [np.asarray(x, dtype=dtype) for x, dtype in (
        (first, int), (second, int), (c1, complex), (c2, complex),
        (imaginary, bool), (scale, float), (diagonal, int), (diagonal_c, int),
        (upper, int), (lower, int), (symmetric, int), (antisymmetric, int))]
    for array in arrays:
        array.flags.writeable = False
    return GroundCoordinates(offset, *arrays)
