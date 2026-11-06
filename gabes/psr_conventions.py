"""Classical polarization, measurement and literature conventions for PSR.

The propagation frame is (x,y,z), with the carrier along +z. Laboratory
optical detunings refer to the line-strength centroid of the upper ground-F
component. Spectral analyzer frequencies have their own, independent axis.
"""

from copy import deepcopy
import math

import numpy as np

from . import constants, magneto_multilevel
from .probe_saturation import gaussian_peak_intensity
from .quantum.contracts import AnalysisFrequencyAxis


def _unit_jones(jones_xy):
    field = np.asarray(jones_xy, dtype=complex)
    if field.shape != (2,) or not np.all(np.isfinite(field)):
        raise ValueError("Jones input requires two finite Cartesian components (Ex,Ey)")
    if not np.isclose(np.vdot(field, field).real, 1.0, rtol=0., atol=1e-12):
        raise ValueError("Jones polarization must have unit norm; intensity is a separate input")
    return field


def jones_from_angles(azimuth_rad=0., ellipticity_rad=0.):
    """Unit Jones carrier with independent azimuth ψ and ellipticity angle ε.

    At ψ=ε=0 the two tangent directions are δEy=δψ and δEy=i δε.
    Positive ellipticity gives positive S3=2 Im(Ex* Ey).
    """
    psi, epsilon = float(azimuth_rad), float(ellipticity_rad)
    if not math.isfinite(psi) or not -math.pi/4 <= epsilon <= math.pi/4:
        raise ValueError("finite azimuth and ellipticity in [-pi/4, pi/4] required")
    c, s = math.cos(psi), math.sin(psi)
    major, minor = math.cos(epsilon), 1j*math.sin(epsilon)
    return np.array([major*c - minor*s, major*s + minor*c])


def jones_to_circular(jones_xy):
    """Return unit polarization coefficients keyed by the signed dipole index q."""
    ex, ey = _unit_jones(jones_xy)
    return {1: (ex - 1j*ey)/math.sqrt(2), -1: (ex + 1j*ey)/math.sqrt(2)}


def circular_to_jones(plus, minus):
    """Transform circular field or polarization amplitudes to Cartesian order."""
    return np.array([plus + minus, 1j*(plus - minus)], dtype=complex)/math.sqrt(2)


def qwp_jones(qwp_deg):
    """Cartesian field from the magneto solver's quarter-wave-plate convention."""
    if not math.isfinite(qwp_deg):
        raise ValueError("finite QWP angle in degrees required")
    polarization = magneto_multilevel._polarization(qwp_deg)
    return circular_to_jones(polarization[1], polarization[-1])


def stokes_rates(jones_xy, derivative_xy_per_m):
    """Return ε, dψ/dz, dε/dz and −d ln S0/dz at the local carrier.

    The derivative is that of the propagating field, expressed in units of
    the local carrier amplitude per meter. It includes common attenuation.
    The azimuth/ellipticity chart is singular at circular polarization. Both
    angular rates there are NaN while the power derivative remains defined.
    """
    ex, ey = _unit_jones(jones_xy)
    derivative = np.asarray(derivative_xy_per_m, dtype=complex)
    if derivative.shape != (2,):
        raise ValueError("two Cartesian field derivatives per meter required")
    if not np.all(np.isfinite(derivative)):
        raise FloatingPointError("nonfinite Cartesian field derivative")
    dx, dy = derivative
    xx, yy, xy = abs(ex)**2, abs(ey)**2, ex.conjugate()*ey
    s0, s1, s2, s3 = xx + yy, xx - yy, 2*xy.real, 2*xy.imag
    dxx, dyy = 2*(ex.conjugate()*dx).real, 2*(ey.conjugate()*dy).real
    dxy = dx.conjugate()*ey + ex.conjugate()*dy
    ds1, ds2, ds3 = dxx - dyy, 2*dxy.real, 2*dxy.imag
    linear = math.hypot(s1, s2)
    epsilon = .5*math.atan2(s3, linear)
    alpha = -(dxx + dyy)/s0
    if linear <= 8*np.finfo(float).eps*s0:
        return epsilon, math.nan, math.nan, float(alpha)
    dpsi = .5*(s1*ds2 - s2*ds1)/linear**2
    dlinear = (s1*ds1 + s2*ds2)/linear
    depsilon = .5*(linear*ds3 - s3*dlinear)/(linear**2 + s3**2)
    return epsilon, float(dpsi), float(depsilon), float(alpha)


def jones_polarization_angles(jones_xy):
    """Return azimuth and ellipticity angles in radians for a unit Jones field.

    Azimuth is defined modulo pi on [-pi/2, pi/2]. Circular polarization
    has no azimuth and returns NaN there. Both angles are independent of
    the common optical phase. The ellipticity uses the shared Stokes rule.
    """
    jones = _unit_jones(jones_xy)
    epsilon = stokes_rates(jones, np.zeros(2, complex))[0]
    ex, ey = jones
    s0 = abs(ex)**2+abs(ey)**2
    s1 = abs(ex)**2-abs(ey)**2
    s2 = 2*(ex.conjugate()*ey).real
    azimuth = math.nan if math.hypot(s1, s2) <= 8*np.finfo(float).eps*s0 else .5*math.atan2(s2, s1)
    return float(azimuth), float(epsilon)


def atom_frame_detuning_mhz(detuning_mhz, *, velocity_m_s, optical_frequency_hz):
    """Laser detuning in one atom's frame, to first order in vz/c (Δ−kv)."""
    if not all(math.isfinite(v) for v in (detuning_mhz, velocity_m_s, optical_frequency_hz)):
        raise ValueError("finite laboratory detuning, velocity and optical frequency required")
    if optical_frequency_hz <= 0:
        raise ValueError("optical frequency must be positive in Hz")
    return float(detuning_mhz - velocity_m_s*optical_frequency_hz/(constants.C_LIGHT*1e6))


def gaussian_peak_intensity_mw_cm2(power_mw, *, intensity_radius_1e2_mm):
    """Peak intensity from matched local Gaussian power and 1/e² radius.

    For an unknown power or unestablished 1/e² radius pass None. In
    particular, a reported 'diameter' alone does not establish this radius.
    """
    if power_mw is None or intensity_radius_1e2_mm is None:
        return None
    return gaussian_peak_intensity(power_mw, intensity_radius_1e2_mm)/10.


def measurement_conventions():
    """Portable conventions shared by the scheme's model and benchmark exports."""
    return {
        "convention_id": "psr-measurement-v1",
        "optical_detuning": "laser minus upper ground-F line-strength centroid, MHz",
        "velocity": "positive along carrier +z in m/s; delta_atom = delta_lab - nu0*vz/c (matched units)",
        "drive": "local total intensity in mW/cm^2 and independent unit (Ex,Ey)",
        "rate": "reset rate in s^-1; Hamiltonian frequencies in rad/s",
        "polarization": "E+/-=(Ex -/+ i Ey)/sqrt(2); S3=2 Im(Ex* Ey); psi and epsilon in rad",
        "qwp": "shared QWP angle beta in degrees; near linear input psi~beta and epsilon~-beta in matched angle units",
        "input_mode": "bright carrier; orthogonal input vacuum for the conditional squeezing readout",
        "analysis_frequency": "independent f in Hz, Omega=2*pi*f in rad/s, inverse Fourier exp(-i*Omega*t)",
        "analyzer_bandwidth": "analysis filter width, instrument resolution bandwidth and video bandwidth retain their reported definitions in Hz",
        "phase": "readout theta in degrees relative to the carrier quadratures; orthogonal trace at theta+90 deg",
        "vacuum": "canonical vacuum covariance I/2, also used by channel added covariance; exported source_covariance and covariance use SQL=I (2V)",
        "noise_db": "10*log10(variance/SQL); positive reported squeezing depth maps to negative signed noise dB",
        "measurement_plane": "cell-output, direct detector, electronics-subtracted and loss-corrected observations retain distinct labels",
        "beam_geometry": "Gaussian peak intensity uses local power and 1/e^2 intensity radius; preserve unknown inputs as null",
    }


def benchmark_case(paper, experiment_id):
    """Select one canonical literature experiment with its evidence and planes.

    Original quantities, source locations, missing values and scan membership
    are retained. ``notes`` belongs to the selected experiment. ``paper_notes``
    and ``paper_ai_notes`` carry paper-wide interpretations and version cautions.
    Common RF units and signed noise-power dB are added for scalar quantities.
    """
    matches = [case for case in paper['experiments'] if case['experiment_id'] == experiment_id]
    if len(matches) != 1:
        raise ValueError(f"expected one experiment {experiment_id!r} in {paper['paper_id']!r}")
    result = deepcopy(matches[0])
    result['paper_id'] = paper['paper_id']
    result['paper_title'] = paper['title']
    result['paper_year'] = paper['year']
    result['publication'] = deepcopy(paper['publication'])
    result['paper_ai_notes'] = deepcopy(paper.get('ai_notes', []))
    if 'notes' in paper:
        result['paper_notes'] = deepcopy(paper['notes'])
    result['sources'] = deepcopy(paper['sources'])
    result['apparatus_group'] = result.get('apparatus_group', paper.get('apparatus_group'))
    for key, value in paper['scope'].items():
        result.setdefault(key, deepcopy(value))

    hz_per_unit = {'Hz': 1., 'kHz': 1e3, 'MHz': 1e6}
    scalar_kinds = {'reported', 'approximate', 'digitized', 'derived'}
    for item in result['conditions']:
        name, unit = item['quantity'], item['unit']
        if unit not in hz_per_unit or item['value_type'] not in scalar_kinds:
            continue
        if name == 'analysis_frequency':
            item['frequency_hz'] = item['value'] * hz_per_unit[unit]
            axis = AnalysisFrequencyAxis.from_hz([item['frequency_hz']])
            item['omega_rad_s'] = float(axis.omega_rad_s[0])
        elif name in {'analysis_bandwidth', 'resolution_bandwidth', 'video_bandwidth'}:
            item['bandwidth_hz'] = item['value'] * hz_per_unit[unit]

    depths = {'squeezing_measured', 'squeezing_before_sensor', 'squeezing_before_filter',
              'squeezing_electronic_noise_subtracted',
              'squeezing_loss_corrected', 'squeezing_detection_corrected'}
    excesses = {'anti_squeezing', 'anti_squeezing_before_filter',
                'anti_squeezing_reported_text', 'anti_squeezing_reported_caption'}
    for item in result['results']:
        if item['unit'] == 'dB' and item['value_type'] in scalar_kinds:
            if item['quantity'] in depths:
                item['noise_db_sql'] = -item['value']
            elif item['quantity'] in excesses:
                item['noise_db_sql'] = item['value']
    return result
