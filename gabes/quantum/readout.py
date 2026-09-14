"""Bright-carrier linearized direct detection with explicit shot-noise units.

Uses four sidebands at positive RF and a separate coherent carrier solution.
Does not include quadratic fluctuation photocurrent, dark photons or RF-window
convolution. Coherent carrier flux approximates the mean detected photon flux.
"""

from dataclasses import dataclass

import numpy as np

from ..constants import ELEMENTARY_CHARGE
from .channels import GaussianChannel, canonical_commutator, covariance_uncertainty_minimum
from .contracts import AnalysisFrequencyAxis, readonly_array
from .sidebands import SIDEBAND_MODES
from .traveling import NambuTransfer


def coherent_carrier_output(main, analysis_axis, input_amplitudes_sqrt_flux):
    """Propagate DC carrier amplitudes [alpha_probe,alpha_conjugate], in s^-1/2.

    Select actual RF=0, never the first spectrum sample or generator-frame DC.
    Spontaneous mean photons are not included in this bright-carrier mean.
    """
    if not isinstance(main, NambuTransfer) or not isinstance(analysis_axis, AnalysisFrequencyAxis):
        raise TypeError("Nambu transfer and laboratory RF axis required")
    rf = analysis_axis.omega_rad_s
    if len(main.transfer) != len(rf) or main.mode_labels != ("probe", "conjugate") or not np.array_equal(main.signs, [1, -1]):
        raise ValueError("matching main probe/conjugate transfer required")
    offset = main.frequency_axis.omega_rad_s-rf
    if np.max(np.abs(offset-offset[0])) > 64*np.finfo(float).eps*max(1., abs(offset[0])):
        raise ValueError("RF axis does not match the main generator frame")
    zero = np.flatnonzero(rf == 0)
    if zero.size != 1 or not main.audit()["passed"]:
        raise ValueError("valid transfer with an explicit RF=0 carrier sample required")
    alpha = readonly_array(input_amplitudes_sqrt_flux)
    if alpha.shape != (2,):
        raise ValueError("two finite coherent input flux amplitudes required")
    b = main.transfer[int(zero[0])]@np.array([alpha[0], alpha[1].conjugate()])
    return readonly_array([b[0], b[1].conjugate()])


@dataclass(frozen=True)
class DetectorResponse:
    analysis_axis: AnalysisFrequencyAxis
    transmissions: np.ndarray              # optical power transmission times QE, two arms
    current_response: np.ndarray           # shape (frequency,2), calibrated A/A, complex
    balance: float                         # i_probe - balance*i_conjugate; declared, never optimized
    electronics_difference_psd_A2_Hz: np.ndarray # positive-frequency one-sided, after weighting
    source: str

    def __post_init__(self):
        if not isinstance(self.analysis_axis, AnalysisFrequencyAxis) or np.any(self.analysis_axis.omega_rad_s <= 0):
            raise ValueError("strictly positive laboratory RF axis required")
        eta = readonly_array(self.transmissions, real=True)
        h = readonly_array(self.current_response)
        electronics = readonly_array(self.electronics_difference_psd_A2_Hz, real=True)
        nf = len(self.analysis_axis.omega_rad_s)
        if eta.shape != (2,) or np.any((eta < 0) | (eta > 1)):
            raise ValueError("two intensity transmissions/QEs in [0,1] required")
        if h.shape != (nf, 2) or electronics.shape != (nf,) or np.any(electronics < 0):
            raise ValueError("frequency-resolved response and nonnegative electronics PSD required")
        g = float(self.balance)
        if not np.isfinite(g) or g < 0 or not isinstance(self.source, str) or not self.source.strip():
            raise ValueError("finite nonnegative balance and detector provenance required")
        for name, value in (("transmissions", eta), ("current_response", h),
                             ("electronics_difference_psd_A2_Hz", electronics), ("balance", g)):
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class IntensityDifferenceSpectrum:
    analysis_axis: AnalysisFrequencyAxis
    quantum_psd_A2_Hz: np.ndarray           # one-sided, positive f
    sql_psd_A2_Hz: np.ndarray               # same currents, balance and electronic response
    electronics_psd_A2_Hz: np.ndarray       # explicitly supplied independent contribution
    quantum_ratio: np.ndarray
    total_ratio: np.ndarray
    detected_carrier_flux_s_inverse: np.ndarray
    minimum_covariance_uncertainty_eigenvalue: float
    maximum_readout_commutator_absolute: float
    source: str

    def __post_init__(self):
        for name in ("quantum_psd_A2_Hz", "sql_psd_A2_Hz", "electronics_psd_A2_Hz",
                     "quantum_ratio", "total_ratio", "detected_carrier_flux_s_inverse"):
            object.__setattr__(self, name, readonly_array(getattr(self, name), real=True))

    @property
    def quantum_db(self):
        return 10*np.log10(self.quantum_ratio)

    @property
    def total_db(self):
        return 10*np.log10(self.total_ratio)


def intensity_difference_spectrum(sideband_covariances, output_carrier_amplitudes_sqrt_flux,
                                  detector, *, mode_labels, analysis_axis):
    """Linearized positive-RF current PSD and SQL, one-sided A^2/Hz.

    Mode order is (probe:+,conjugate:-,probe:-,conjugate:+), with (x,p) for each.
    Vacuum fluctuations enter both optical losses and quantum efficiency through
    a loss channel; coherent means get sqrt(eta) once. Electronic gain affects
    the optical PSD and SQL identically. Additive electronics is not subtracted.
    """
    if not isinstance(detector, DetectorResponse) or tuple(mode_labels) != SIDEBAND_MODES:
        raise ValueError("explicit detector response and ordered four-sideband labels required")
    if not isinstance(analysis_axis, AnalysisFrequencyAxis) or not np.array_equal(
            analysis_axis.omega_rad_s, detector.analysis_axis.omega_rad_s):
        raise ValueError("covariance and detector RF axes must match exactly")
    v = readonly_array(sideband_covariances, real=True)
    nf = len(detector.analysis_axis.omega_rad_s)
    if v.shape != (nf, 8, 8):
        raise ValueError("one real four-mode covariance per detector RF sample required")
    beta = readonly_array(output_carrier_amplitudes_sqrt_flux)
    if beta.shape != (2,):
        raise ValueError("two finite output carrier amplitudes required")
    eta = detector.transmissions
    detected_beta = np.sqrt(eta)*beta
    flux = np.abs(detected_beta)**2
    sql_two_sided = ELEMENTARY_CHARGE**2*(
        np.abs(detector.current_response[:, 0])**2*flux[0]
        +detector.balance**2*np.abs(detector.current_response[:, 1])**2*flux[1])
    if np.any(sql_two_sided <= 0) or not np.isfinite(sql_two_sided).all():
        raise ValueError("nonzero finite bright-carrier SQL required at every RF point")
    loss = GaussianChannel.vacuum_attenuator(SIDEBAND_MODES, eta[[0, 1, 0, 1]], source=detector.source)
    j = canonical_commutator(4)
    values, uncertainty, commutators = [], [], []
    for covariance, (hp, hc) in zip(v, detector.current_response):
        detected = loss.apply_covariance(covariance)
        uncertainty.append(covariance_uncertainty_minimum(detected))
        bp, bc = detected_beta
        # y(+Omega)=hp(bp* ap+ + bp ap-†)-g hc(bc* ac+ + bc ac-†).
        coefficients = np.array([hp*bp.conjugate(), -detector.balance*hc*bc,
                                 hp*bp, -detector.balance*hc*bc.conjugate()])
        row = np.empty(8, complex)
        row[::2] = coefficients/np.sqrt(2)
        row[1::2] = 1j*np.array([1, -1, -1, 1])*coefficients/np.sqrt(2)
        row *= ELEMENTARY_CHARGE
        value = row@detected@row.conj()
        commutator = float(abs(row@(1j*j)@row.conj()))
        if abs(value.imag) > 1e-10*max(abs(value.real), np.finfo(float).tiny) or value.real <= 0:
            raise ValueError("nonpositive or nonreal photocurrent PSD; no clipping applied")
        if commutator > 1e-10*np.linalg.norm(row)**2:
            raise ValueError("four-sideband photocurrent quadratures do not commute")
        values.append(2*value.real)  # S_two_sided_Hz(f)=S_omega(2*pi*f); positive one-sided doubles it.
        commutators.append(commutator)
    quantum = np.array(values)
    sql = 2*sql_two_sided
    return IntensityDifferenceSpectrum(
        detector.analysis_axis, quantum, sql, detector.electronics_difference_psd_A2_Hz,
        quantum/sql, (quantum+detector.electronics_difference_psd_A2_Hz)/sql,
        flux, float(min(uncertainty)), float(max(commutators)), detector.source)
