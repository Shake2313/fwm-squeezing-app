"""Narrowband, balanced homodyne spectra of the orthogonal PSR polarization."""
import numpy as np

from .psr_conventions import _unit_jones
from .psr_spectrum import PSR_NAMBU_MODES, PSR_RF_FRAME
from .quantum.channels import GaussianChannel
from .quantum.traveling import NambuTransfer


def _matrix_norm(values):
    return np.linalg.norm(values, axis=(-2, -1))


def _positive_sidebands(transfer):
    if not isinstance(transfer, NambuTransfer):
        raise TypeError("PSR homodyne requires a NambuTransfer")
    if (transfer.frequency_axis.frame != PSR_RF_FRAME or transfer.mode_labels != PSR_NAMBU_MODES
            or not np.array_equal(transfer.signs, [1, 1, -1, -1])):
        raise ValueError("PSR homodyne requires x, y, x-dagger, y-dagger in the pump rotating frame")
    omega = transfer.frequency_axis.omega_rad_s
    positive = np.flatnonzero(omega > 0)
    negative = np.flatnonzero(omega < 0)[::-1]
    if (not positive.size or positive.size != negative.size
            or not np.allclose(omega[positive], -omega[negative],
                               rtol=64*np.finfo(float).eps, atol=0)):
        raise ValueError("each nonzero analysis frequency requires a reflected RF partner")

    # a†(Ω) is a(−Ω)†. Reflection exchanges the annihilation/creation blocks
    # and the two orderings of the Langevin forces. Check all x and y ports
    # before selecting the measured polarization.
    swap = [2, 3, 0, 1]
    greater, lesser = transfer.noise_greater, transfer.noise_lesser
    for left, right in ((transfer.transfer, transfer.transfer), (greater, lesser), (lesser, greater)):
        actual = left[negative]
        expected = right[positive].conj()[:, swap][:, :, swap]
        scale = np.maximum(_matrix_norm(actual), _matrix_norm(expected))
        residual = _matrix_norm(actual - expected)
        # A vacuum ordering may vanish to roundoff even when its partner
        # ordering is finite. Use their joint size only for the roundoff floor.
        roundoff = 64*np.finfo(float).eps*(_matrix_norm(left[positive]) + _matrix_norm(right[negative]))
        if np.any(residual > 1e-8*scale + roundoff):
            raise ValueError("PSR reflected transfer or ordered noise fails sideband conjugation")
    return positive


def psr_homodyne_readout(transfer, *, detection_efficiency=1.0,
                         phase_degrees=None, carrier_phase_rad=0.0):
    """Return orthogonal Stokes/homodyne spectra for coherent x and vacuum y.

    For f>0 the complex photocurrent is proportional to
    exp(-iθ) a_y(+Ω) + exp(+iθ) a_y†(-Ω). θ is relative to the bright carrier.
    A stationary output has two sideband occupations and an anomalous moment
    m=<a_y(+Ω) a_y(-Ω)>. Those three quantities retain arbitrary sideband
    imbalance and complex squeezing angle. DC belongs to the mean carrier.

    The optical sideband covariance is canonical ([x,p]=i, vacuum I/2).
    The returned homodyne ellipse and phase variances use detected SQL=1.
    Equal optical/QE loss on both sidebands is applied once after the cell.
    The caller determines whether a bright detected carrier exists for a
    Stokes measurement. At zero efficiency this algebra returns vacuum.
    """
    positive = _positive_sidebands(transfer)
    efficiency, carrier_phase = float(detection_efficiency), float(carrier_phase_rad)
    if not np.isfinite(carrier_phase) or not 0 <= efficiency <= 1:
        raise ValueError("finite carrier phase and external efficiency in [0,1] required")
    phase = np.linspace(0., 180., 361) if phase_degrees is None else np.asarray(phase_degrees, dtype=float)
    if phase.ndim != 1 or not phase.size or not np.isfinite(phase).all():
        raise ValueError("phase_degrees must be a nonempty finite vector")

    greater, lesser = transfer.vacuum_output()  # includes the transfer commutator/PSD audit
    ordered_sym = (greater[positive] + lesser[positive])/2
    upper = ordered_sym[:, 1, 1].real
    lower = ordered_sym[:, 3, 3].real
    pair = (ordered_sym[:, 1, 3] + ordered_sym[:, 3, 1].conj())/2
    # Interleaved (x+,p+,x−,p−). The four entries coupling the sidebands are
    # fixed by stationarity and by the complex anomalous moment.
    source_v = np.zeros((len(positive), 4, 4))
    source_v[:, 0, 0] = source_v[:, 1, 1] = upper
    source_v[:, 2, 2] = source_v[:, 3, 3] = lower
    source_v[:, 0, 2] = source_v[:, 2, 0] = pair.real
    source_v[:, 1, 3] = source_v[:, 3, 1] = -pair.real
    source_v[:, 0, 3] = source_v[:, 3, 0] = pair.imag
    source_v[:, 1, 2] = source_v[:, 2, 1] = pair.imag

    detector = GaussianChannel.vacuum_attenuator(
        ("y:+", "y:-"), [efficiency, efficiency], source="PSR downstream optical and detector efficiency")
    detected_v = np.stack([detector.apply_covariance(v) for v in source_v])
    # A normalized cosine temporal mode is (a+ + a−)/sqrt(2). Multiplication
    # of its canonical covariance by 2 for SQL units cancels both sqrt(2)s.
    c, s = np.cos(carrier_phase), np.sin(carrier_phase)
    projection = np.array([[c, s, c, s], [-s, c, -s, c]])
    source_sql = projection @ source_v @ projection.T
    detected_sql = projection @ detected_v @ projection.T
    eigenvalues, eigenvectors = np.linalg.eigh(detected_sql)
    if not np.isfinite(eigenvalues).all() or np.any(eigenvalues <= 0):
        raise ValueError("nonpositive or nonfinite PSR homodyne spectral variance")
    minimum, maximum = eigenvalues.T
    if np.any(minimum <= 64*np.finfo(float).eps*maximum):
        raise ValueError("PSR minimum detected variance is unresolved at float64 precision")
    angle_defined = maximum - minimum > 64*np.finfo(float).eps*maximum
    min_angle = np.degrees(np.arctan2(eigenvectors[:, 1, 0], eigenvectors[:, 0, 0])) % 180
    min_angle[~angle_defined] = np.nan
    theta = np.deg2rad(phase)
    cos_theta, sin_theta = np.cos(theta), np.sin(theta)
    rows_a = np.column_stack((cos_theta, sin_theta))
    rows_b = np.column_stack((-sin_theta, cos_theta))
    noise_a = np.einsum("pi,fij,pj->fp", rows_a, detected_sql, rows_a)
    noise_b = np.einsum("pi,fij,pj->fp", rows_b, detected_sql, rows_b)
    if np.any(noise_a <= 0) or np.any(noise_b <= 0):
        raise ValueError("homodyne phase projection has unrepresentable variance")
    return {
        "frequencies_hz": transfer.frequency_axis.omega_rad_s[positive]/(2*np.pi),
        "phase_degrees": phase.copy(), "carrier_phase_rad": carrier_phase,
        "detection_efficiency": efficiency,
        "source_sideband_covariance": source_v, "sideband_covariance": detected_v,
        "source_covariance": source_sql, "covariance": detected_sql,
        "minimum_variance": minimum, "maximum_variance": maximum,
        "minimum_db": 10*np.log10(minimum), "maximum_db": 10*np.log10(maximum),
        "angle_defined": angle_defined, "min_angle_deg": min_angle,
        "noise_a": noise_a, "noise_b": noise_b,
        "noise_a_db": 10*np.log10(noise_a), "noise_b_db": 10*np.log10(noise_b),
        "convention": (
            "Narrowband positive RF, orthogonal y polarization referenced to the bright x carrier; "
            "sideband covariance in [x+,p+,x-,p-] uses vacuum I/2; "
            "homodyne covariance and phase variances use SQL=I, dB=10 log10(variance/SQL)."),
    }


def psr_output_homodyne_readout(transfer, output_jones_xy, *, detection_efficiency=1.,
                                phase_degrees=None, carrier_phase_rad=0.):
    """Read the polarization orthogonal to an arbitrary output Jones carrier.

    ``output_jones_xy`` and the transfer must use the same Cartesian optical
    phase frame. A separately removed common propagation phase must not be
    passed again. ``carrier_phase_rad`` is only an additional LO phase offset.

    The x component fixes the polarization chart when its magnitude exceeds
    1e-12; at the x pole the y component fixes it. Remove that component's phase
    gamma, choose v=(-u_y*,u_x*) for the canonical carrier u, and project onto
    (u,v). The homodyne LO phase includes gamma. Thus multiplying both the
    carrier and optical output by the same phase leaves the whole phase scan
    unchanged, not merely its eigenvalues. The chart changes at the x pole;
    phase angles there describe the declared y-reference chart.

    Both optical inputs have coherent-state vacuum fluctuations. Only output
    coordinates are rotated. The existing two-sideband SQL and detector-loss
    calculation is then reused without an additional vacuum reservoir.
    """
    _positive_sidebands(transfer)
    jones = _unit_jones(output_jones_xy).copy()
    pivot = 0 if abs(jones[0]) > 1e-12 else 1
    phase = float(np.angle(jones[pivot]))
    canonical = jones*np.exp(-1j*phase)
    orthogonal = np.array([-canonical[1].conjugate(), canonical[0].conjugate()])
    analyzer = np.stack((canonical.conj(), orthogonal.conj()))
    rotation = np.zeros((4, 4), complex)
    rotation[:2, :2] = analyzer
    rotation[2:, 2:] = analyzer.conj()
    adjoint = rotation.conj().T
    projected = NambuTransfer(
        transfer.frequency_axis, transfer.signs, transfer.mode_labels,
        rotation @ transfer.transfer,
        rotation @ transfer.noise_greater @ adjoint,
        rotation @ transfer.noise_lesser @ adjoint)
    result = psr_homodyne_readout(
        projected, detection_efficiency=detection_efficiency,
        phase_degrees=phase_degrees, carrier_phase_rad=phase+float(carrier_phase_rad))
    result.update(output_jones_xy=jones, canonical_output_jones_xy=canonical,
                  orthogonal_jones_xy=orthogonal, output_jones_phase_rad=phase,
                  lo_offset_rad=float(carrier_phase_rad),
                  polarization_phase_chart="x" if pivot == 0 else "y",
                  output_polarization_analyzer=analyzer)
    result["convention"] = (
        "Narrowband positive RF in the polarization orthogonal to the output carrier; "
        "x-first canonical Jones chart (y at the x pole), referenced to the same carrier phase; "
        "sideband covariance in [x+,p+,x-,p-] uses vacuum I/2; "
        "homodyne covariance and phase variances use SQL=I, dB=10 log10(variance/SQL).")
    return result
