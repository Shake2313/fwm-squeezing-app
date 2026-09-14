"""Quadratic photocurrent of a stationary, phase-selected Gaussian beam pair.

Finite optical bins are physical top-hat collection bands with vacuum outside,
not a claim that uncomputed fluorescence vanishes. Frequencies here are optical
envelope offsets in Hz; photocurrent RF and generator-frame axes stay separate.
"""

from dataclasses import dataclass

import numpy as np

from ..constants import ELEMENTARY_CHARGE
from .contracts import AnalysisFrequencyAxis, readonly_array
from .readout import DetectorResponse, IntensityDifferenceSpectrum
from .sidebands import _validate_pair


@dataclass(frozen=True)
class PairedOpticalBins:
    """Piecewise constant n_p(f), n_c(f), m_pc(f)=<a_p(f)a_c(-f)>.

    No normal interbeam or same-beam anomalous moments are supplied. The band
    has equal bins, symmetric about zero, and independent vacuum outside.
    """

    centers_hz: np.ndarray
    width_hz: float
    occupations: np.ndarray
    anomalous_pc: np.ndarray
    source: str

    def __post_init__(self):
        f = readonly_array(self.centers_hz, real=True)
        n = readonly_array(self.occupations, real=True)
        m = readonly_array(self.anomalous_pc)
        width = float(self.width_hz)
        if (f.ndim != 1 or len(f) < 2 or n.shape != (len(f), 2) or m.shape != f.shape
                or not np.isfinite(width) or width <= 0 or not self.source.strip()):
            raise ValueError("finite symmetric optical bins and pair moments required")
        atol = 64*np.finfo(float).eps*max(width, float(np.max(abs(f))))
        if np.max(abs(f+f[::-1])) > atol or np.max(abs(np.diff(f)-width)) > atol:
            raise ValueError("optical bin centers must be symmetric with the declared uniform width")
        for key, value in (("centers_hz", f), ("occupations", n), ("anomalous_pc", m), ("width_hz", width)):
            object.__setattr__(self, key, value)
        if not self.audit()["passed"]:
            raise ValueError("unphysical paired spectral moments; no PSD repair applied")

    @property
    def half_width_hz(self):
        return float(self.centers_hz[-1]+self.width_hz/2)

    @property
    def spontaneous_flux_s_inverse(self):
        return readonly_array(self.width_hz*self.occupations.sum(axis=0), real=True)

    def audit(self):
        n, m = self.occupations, self.anomalous_pc
        matrices = np.zeros((2, len(n), 2, 2), complex)
        matrices[:, :, 0, 0] = n[:, 0]
        matrices[:, :, 1, 1] = n[::-1, 1]
        matrices[:, :, 0, 1], matrices[:, :, 1, 0] = m, m.conj()
        matrices[0, :, 0, 0] += 1
        matrices[1, :, 1, 1] += 1
        mins = np.linalg.eigvalsh(matrices)[..., 0]
        scale = np.maximum(1., np.linalg.norm(matrices, axis=(-2, -1)))
        return {"passed": bool(np.all(mins >= -1e-10*scale) and np.all(n >= 0)),
                "minimum_ordered_eigenvalue": float(mins.min()),
                "scope": "phase-selected stationary Gaussian two-mode closure"}

    def attenuate(self, transmissions):
        eta = readonly_array(transmissions, real=True)
        if eta.shape != (2,) or np.any(eta < 0) or np.any(eta > 1):
            raise ValueError("two intensity transmissions in [0,1] required")
        return PairedOpticalBins(self.centers_hz, self.width_hz, self.occupations*eta,
            self.anomalous_pc*np.sqrt(eta.prod()), self.source+"; vacuum attenuation")


def paired_optical_bins(main, companion, envelope_axis, *, bin_width_hz, source):
    """Extract ordered occupations without subtracting the vacuum identity."""
    _validate_pair(main, companion, envelope_axis)
    for a, b in ((main.transfer, companion.transfer[::-1].conj()),
                 (main.noise_greater, companion.noise_lesser[::-1].conj()),
                 (main.noise_lesser, companion.noise_greater[::-1].conj())):
        scale = np.maximum(np.linalg.norm(a, axis=(-2, -1)), np.finfo(float).tiny)
        if np.any(np.linalg.norm(a-b, axis=(-2, -1)) > 1e-8*scale):
            raise ValueError("reflected companion conjugation failed")
    greater, lesser = main.vacuum_output()
    n = np.column_stack([lesser[:, 0, 0].real, greater[::-1, 1, 1].real])
    return PairedOpticalBins(envelope_axis.frequency_hz, bin_width_hz, n, greater[:, 0, 1], source)


def box_autocorrelation(values, bin_width_hz, offsets_hz):
    """Integral x*(f)x(f+RF) df, exact for constant bins with zero padding.

    Fractional lags use overlap lengths, never circular FFT wraparound.
    Complex autocorrelations retain the relative interbeam delay phase.
    """
    x = readonly_array(values)
    offsets = readonly_array(offsets_hz, real=True)
    width = float(bin_width_hz)
    if x.ndim != 1 or offsets.ndim != 1 or not np.isfinite(width) or width <= 0:
        raise ValueError("one-dimensional bins/offsets and positive width required")
    def integral(k):
        if k >= len(x):
            return 0j
        return np.vdot(x[:len(x)-k], x[k:])
    result = []
    for offset in offsets:
        lag = abs(offset)/width
        k = int(np.floor(lag))
        fraction = lag-k
        value = width*((1-fraction)*integral(k)+fraction*integral(k+1))
        result.append(value if offset >= 0 else value.conjugate())
    return readonly_array(result)


def gaussian_count_covariance(bins, analysis_axis):
    """Two-sided photon-flux PSD [s^-1], including spontaneous shot noise.

    C_pp=Phi_p+int n_p(f)n_p(f+RF)df, likewise cc;
    C_pc=int m_pc*(f)m_pc(f+RF)df, C_cp=C_pc*. Wick factorization
    is an additional Gaussian assumption, not an atomic fourth-moment theorem.
    """
    if not isinstance(bins, PairedOpticalBins) or not isinstance(analysis_axis, AnalysisFrequencyAxis):
        raise TypeError("optical bins and a separate photocurrent RF axis required")
    f, n = analysis_axis.frequency_hz, bins.occupations
    normal = np.column_stack([box_autocorrelation(n[:, j], bins.width_hz, f).real for j in (0, 1)])
    anomalous = box_autocorrelation(bins.anomalous_pc, bins.width_hz, f)
    covariance = np.zeros((len(f), 2, 2), complex)
    covariance[:, 0, 0] = bins.spontaneous_flux_s_inverse[0]+normal[:, 0]
    covariance[:, 1, 1] = bins.spontaneous_flux_s_inverse[1]+normal[:, 1]
    covariance[:, 0, 1], covariance[:, 1, 0] = anomalous, anomalous.conj()
    minimum = np.linalg.eigvalsh(covariance)[:, 0]
    if np.any(minimum < -1e-10*np.maximum(np.linalg.norm(covariance, axis=(-2, -1)), 1.)):
        raise ValueError("quadratic photon-flux covariance is not positive")
    return readonly_array(covariance)


@dataclass(frozen=True)
class GaussianReadoutCorrection:
    analysis_axis: AnalysisFrequencyAxis
    spontaneous_detected_flux_s_inverse: np.ndarray
    quadratic_psd_A2_Hz: np.ndarray
    added_sql_psd_A2_Hz: np.ndarray
    count_covariance_s_inverse: np.ndarray
    optical_half_width_hz: float
    source: str

    def __post_init__(self):
        for name in ("spontaneous_detected_flux_s_inverse", "quadratic_psd_A2_Hz", "added_sql_psd_A2_Hz"):
            object.__setattr__(self, name, readonly_array(getattr(self, name), real=True))
        object.__setattr__(self, "count_covariance_s_inverse", readonly_array(self.count_covariance_s_inverse))
        nf = len(self.analysis_axis.omega_rad_s)
        if (self.spontaneous_detected_flux_s_inverse.shape != (2,)
                or self.quadratic_psd_A2_Hz.shape != (nf,) or self.added_sql_psd_A2_Hz.shape != (nf,)
                or self.count_covariance_s_inverse.shape != (nf, 2, 2)
                or not np.isfinite(self.optical_half_width_hz) or self.optical_half_width_hz <= 0):
            raise ValueError("consistent quadratic correction dimensions and positive optical band required")


def quadratic_readout_correction(bins, detector):
    if not isinstance(detector, DetectorResponse):
        raise TypeError("DetectorResponse required")
    detected = bins.attenuate(detector.transmissions)
    cov = gaussian_count_covariance(detected, detector.analysis_axis)
    weights = detector.current_response*np.array([1., -detector.balance])
    quadratic = 2*ELEMENTARY_CHARGE**2*np.einsum("fi,fij,fj->f", weights, cov, weights.conj()).real
    sql = 2*ELEMENTARY_CHARGE**2*(abs(weights)**2@detected.spontaneous_flux_s_inverse)
    return GaussianReadoutCorrection(detector.analysis_axis, detected.spontaneous_flux_s_inverse,
        quadratic, sql, cov, bins.half_width_hz, bins.source+"; Gaussian Wick quadratic photocurrent")


def corrected_gaussian_readout(bright, correction):
    """Keep the separately sampled carrier beat term; add quadratic PSD and SQL.

    Caller must use the same optical state, detector and collection band. The
    coherent carrier and both RF sidebands must lie inside that optical band.
    Linear/quadratic cross terms vanish for a displaced Gaussian state only.
    """
    if not isinstance(bright, IntensityDifferenceSpectrum) or not isinstance(correction, GaussianReadoutCorrection):
        raise TypeError("bright-carrier result and Gaussian correction required")
    if not np.array_equal(bright.analysis_axis.omega_rad_s, correction.analysis_axis.omega_rad_s):
        raise ValueError("bright and quadratic RF axes must match")
    if np.any(abs(bright.analysis_axis.frequency_hz) >= correction.optical_half_width_hz):
        raise ValueError("bright carrier beat sidebands must lie inside the declared optical collection band")
    quantum = bright.quantum_psd_A2_Hz+correction.quadratic_psd_A2_Hz
    sql = bright.sql_psd_A2_Hz+correction.added_sql_psd_A2_Hz
    if np.any(quantum <= 0) or np.any(sql <= 0):
        raise ValueError("positive full Gaussian PSD and SQL required")
    return {"quantum_psd_A2_Hz": readonly_array(quantum, real=True),
            "sql_psd_A2_Hz": readonly_array(sql, real=True),
            "quantum_ratio": readonly_array(quantum/sql, real=True),
            "total_ratio": readonly_array((quantum+bright.electronics_psd_A2_Hz)/sql, real=True),
            "quantum_db": readonly_array(10*np.log10(quantum/sql), real=True),
            "detected_total_flux_s_inverse": readonly_array(
                bright.detected_carrier_flux_s_inverse+correction.spontaneous_detected_flux_s_inverse, real=True)}
