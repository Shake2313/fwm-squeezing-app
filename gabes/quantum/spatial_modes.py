"""Normalized transverse optical profiles and reciprocal local projection.

Profiles are physical collected modes, not atomic grating-index copies. The
local Markov construction requires stationary atomic centers; transport across
the aperture or between z slices needs its own correlated atomic noise model.
"""

from dataclasses import dataclass

import numpy as np

from .contracts import readonly_array
from .traveling import LocalNambuGenerator


@dataclass(frozen=True)
class TransverseModeGrid:
    positions_m: np.ndarray
    area_weights_m2: np.ndarray
    mode_values_m_inverse: np.ndarray
    source: str

    def __post_init__(self):
        positions = readonly_array(self.positions_m, real=True)
        weights = readonly_array(self.area_weights_m2, real=True)
        modes = readonly_array(self.mode_values_m_inverse)
        if (positions.ndim != 2 or positions.shape[1] != 2 or not len(positions)
                or weights.shape != (len(positions),) or np.any(weights <= 0)
                or modes.shape != (len(positions), 2)
                or not isinstance(self.source, str) or not self.source.strip()):
            raise ValueError('positions, positive area weights, two profiles and provenance required')
        norms = weights@abs(modes)**2
        if np.max(abs(norms-1)) > 1e-12:
            raise ValueError('each optical profile must satisfy integral |u|^2 dA=1; no silent normalization')
        for key, value in (('positions_m', positions), ('area_weights_m2', weights), ('mode_values_m_inverse', modes)):
            object.__setattr__(self, key, value)

    @classmethod
    def rectangle(cls, width_m, height_m, *, order_x, order_y=1):
        if not np.isfinite([width_m, height_m]).all() or width_m <= 0 or height_m <= 0:
            raise ValueError('positive finite aperture dimensions required')
        for order in (order_x, order_y):
            if isinstance(order, bool) or int(order) != order or order < 1:
                raise ValueError('positive integer spatial quadrature orders required')
        x, wx = np.polynomial.legendre.leggauss(int(order_x))
        y, wy = np.polynomial.legendre.leggauss(int(order_y))
        xx, yy = np.meshgrid(x*width_m/2, y*height_m/2, indexing='ij')
        weights = (wx[:, None]*wy[None, :]*width_m*height_m/4).ravel()
        return cls(np.column_stack([xx.ravel(), yy.ravel()]), weights,
            np.full((len(weights), 2), 1/np.sqrt(width_m*height_m), complex),
            'declared coextensive normalized top-hat probe/conjugate profiles; distinct optical carrier bands')


def project_local_generator(local, coefficients):
    """H drive uses eta*b, Maxwell readout eta*, noise eta* D eta.

    M and D use the SAME reciprocal profile factors. A diagonal E commutes
    with canonical J, hence E†(MJ+JM†+D>-D<)E=0 exactly. E need not be unitary;
    this is a coupling projection, not a lossy quantum channel on its own.
    """
    eta = readonly_array(coefficients)
    if not isinstance(local, LocalNambuGenerator) or eta.shape != (2,):
        raise TypeError('two optical profile coefficients and a full Nambu local generator required')
    if not np.array_equal(local.signs, [1., 1., -1., -1.]):
        raise ValueError('probe/conjugate annihilators then creators required')
    e = np.r_[eta, eta.conj()]
    def project(value):
        return e.conj()[:, None]*value*e[None, :]
    result = LocalNambuGenerator(local.frequency_axis, local.signs, local.mode_labels,
        project(local.drift), project(local.noise_greater_by_reservoir),
        project(local.noise_lesser_by_reservoir), local.reservoir_names)
    if not result.audit()['passed']:
        raise ValueError('reciprocal spatial projection failed canonical commutator/PSD')
    return result
