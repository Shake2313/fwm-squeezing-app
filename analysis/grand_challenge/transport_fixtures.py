"""Declared moving two-level fixture; no Rb or experimental interpretation."""

import numpy as np

from gabes.quantum.reservoirs import CollapseChannel, ExplicitReservoirs, thermal_reset_channels
from gabes.quantum.transport import BallisticPath


def moving_fixture():
    path = BallisticPath([-.0003, 0., 0.], [150., 0., 100.], 4e-6,
        'synthetic ballistic ray through a spatially varying two-level drive')
    lower = np.array([[0., 1.], [0., 0.]])
    sx, sy, sz = lower+lower.T, -1j*(lower-lower.T), np.diag([1., -1.])
    rates = {'decay_s_inverse': 2*np.pi*.12e6, 'reset_s_inverse': 2*np.pi*.02e6,
             'detuning_rad_s': 2*np.pi*.17e6, 'peak_rabi_rad_s': 2*np.pi*.22e6}
    reservoirs = ExplicitReservoirs(2, (CollapseChannel('emission', np.sqrt(rates['decay_s_inverse'])*lower,
        'spontaneous_emission', 'declared toy decay'),)+thermal_reset_channels(rates['reset_s_inverse'], [.7, .3],
        source='declared internal bath, distinct from ballistic boundary renewal'))
    rho = np.array([[.65, .08+.04j], [.08-.04j, .35]])
    def envelope(r):
        return np.exp(-(r[0]/.00015)**2)*np.exp(1j*(9000*r[0]+5000*r[2]))
    def hamiltonian(age, r):
        drive = rates['peak_rabi_rad_s']/2*envelope(r)*lower.T
        return np.diag([0., rates['detuning_rad_s']])+drive+drive.conj().T
    def readout(age, r):
        first = np.exp(-(r[0]/.00018)**2)
        second = np.exp(-((r[0]-.00005)/.00016)**2)*np.exp(.4j*r[2]/.0004)
        return np.array([first*(.7*sz+.9*sx+.2j*sy)/np.sqrt(2)+.05*np.eye(2),
            second*(.15*sz+.2*sx+.9*sy)/np.sqrt(2)])
    def drives(age, r):
        value = .5*rates['peak_rabi_rad_s']*envelope(r)*lower.T
        return np.array([value+value.conj().T, .5*rates['detuning_rad_s']*sz])
    return {'path': path, 'hamiltonian': hamiltonian, 'readout': readout, 'drives': drives,
        'reservoirs': reservoirs, 'boundary_state': rho, 'parameters': rates,
        'scope': 'two-level prescribed spatial drive and readout; internal reset bath is not inferred transit'}
