"""Direct ordered-Nambu bright-carrier current contraction (no quadrature adapter).

For coherent input sidebands only. Used to independently verify ordering,
carrier phase, opposite RF samples, detection loss and the factor of two.
"""

import numpy as np

from gabes.constants import ELEMENTARY_CHARGE


def direct_current_psd(main_greater, main_lesser, positive_indices, negative_indices,
                       carrier_amplitudes, transmissions, response, balance):
    eta = np.asarray(transmissions, float)
    amplitude = np.sqrt(eta)*np.asarray(carrier_amplitudes, complex)
    loss = np.diag(np.sqrt(eta))
    values = []
    for plus, minus, (hp, hc) in zip(positive_indices, negative_indices, response):
        greater = loss@main_greater[plus]@loss+np.diag([1-eta[0], 0])
        lesser = loss@main_lesser[minus]@loss+np.diag([0, 1-eta[1]])
        row_a = np.array([hp*amplitude[0].conjugate(), -balance*hc*amplitude[1]])
        row_b = np.array([hp*amplitude[0], -balance*hc*amplitude[1].conjugate()])
        value = row_a@greater@row_a.conj()+row_b@lesser.T@row_b.conj()
        values.append(2*ELEMENTARY_CHARGE**2*value.real)
    return np.array(values)
