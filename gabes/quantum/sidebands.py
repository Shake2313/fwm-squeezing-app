"""Phase-selected two-beam sidebands and normalized spectral wavepacket channels.

Positive RF uses four independent optical modes. No cross-sector correlations
outside the declared two-mode closure are invented. DC is reserved for carriers.
"""

from dataclasses import dataclass

import numpy as np

from .channels import GaussianChannel
from .contracts import AnalysisFrequencyAxis, readonly_array
from .traveling import NambuTransfer


SIDEBAND_MODES = ("probe:+", "conjugate:-", "probe:-", "conjugate:+")


def _validate_pair(main, companion, analysis_axis):
    if (not isinstance(main, NambuTransfer) or not isinstance(companion, NambuTransfer)
            or not isinstance(analysis_axis, AnalysisFrequencyAxis)):
        raise TypeError("main/companion Nambu transfers and a laboratory RF axis required")
    rf = analysis_axis.omega_rad_s
    if (main.mode_labels != ("probe", "conjugate") or companion.mode_labels != main.mode_labels
            or not np.array_equal(main.signs, [1, -1])
            or not np.array_equal(companion.signs, [-1, 1])):
        raise ValueError("ordered probe/conjugate main and companion sectors required")
    if len(main.transfer) != len(rf) or len(companion.transfer) != len(rf):
        raise ValueError("RF and generator-frequency dimensions must agree")
    offsets = [t.frequency_axis.omega_rad_s-rf for t in (main, companion)]
    tol = 64*np.finfo(float).eps*max(1., *(np.max(np.abs(v)) for v in offsets))
    if (np.max(np.abs(offsets[0]-offsets[0][0])) > tol
            or np.max(np.abs(offsets[1]-offsets[1][0])) > tol
            or abs(offsets[0][0]+offsets[1][0]) > tol):
        raise ValueError("main and companion must carry opposite constant optical beat offsets")
    if not main.audit()["passed"] or not companion.audit()["passed"]:
        raise ValueError("invalid Nambu transfer")


def _mirror_index(rf, index):
    target = -rf[index]
    found = np.flatnonzero(np.abs(rf-target) <= 64*np.finfo(float).eps*max(1., abs(target)))
    if found.size != 1:
        raise ValueError("every sideband sample needs exactly one reflected RF sample")
    return int(found[0])


def _pair_channel(main, companion, index, mirror, modes):
    for left, right in ((main.transfer[index], companion.transfer[mirror].conj()),
                        (main.noise_greater[index], companion.noise_lesser[mirror].conj()),
                        (main.noise_lesser[index], companion.noise_greater[mirror].conj())):
        scale = max(np.linalg.norm(left), np.linalg.norm(right), np.finfo(float).tiny)
        if np.linalg.norm(left-right) > 1e-8*scale:
            raise ValueError("companion conjugation failed; no sideband covariance constructed")
    t = main.transfer[index]
    u = np.diag([t[0, 0], t[1, 1].conjugate()])
    v = np.array([[0., t[0, 1]], [t[1, 0].conjugate(), 0.]])
    x = np.zeros((4, 4))
    x[::2, ::2], x[::2, 1::2] = (u+v).real, -(u-v).imag
    x[1::2, ::2], x[1::2, 1::2] = (u+v).imag, (u-v).real
    ordered_sym = (main.noise_greater[index]+main.noise_lesser[index])/2
    if np.linalg.norm(ordered_sym-ordered_sym.conj().T) > 1e-10*max(1., np.linalg.norm(ordered_sym)):
        raise ValueError("added ordered noise is not Hermitian")
    # Symmetric quadrature covariance: take the anticommutator, not a PSD repair.
    m = (ordered_sym[0, 1]+ordered_sym[1, 0].conjugate())/2
    y = np.diag(np.repeat(np.diag(ordered_sym).real, 2))
    cross = np.array([[m.real, m.imag], [m.imag, -m.real]])
    y[:2, 2:], y[2:, :2] = cross, cross.T
    result = GaussianChannel(x, y, modes, modes,
                             "declared paired-sideband sampling of microscopic Nambu transfer")
    if not result.audit().passed:
        raise ValueError("paired sideband channel violates complete positivity")
    return result


@dataclass(frozen=True)
class SidebandChannels:
    analysis_axis: AnalysisFrequencyAxis
    channels: tuple[GaussianChannel, ...]

    def __post_init__(self):
        if not isinstance(self.analysis_axis, AnalysisFrequencyAxis):
            raise TypeError("laboratory RF axis required")
        channels = tuple(self.channels)
        if (np.any(self.analysis_axis.omega_rad_s <= 0) or len(channels) != len(self.analysis_axis.omega_rad_s)
                or any(c.input_modes != SIDEBAND_MODES or c.output_modes != SIDEBAND_MODES
                       or not c.audit().passed for c in channels)):
            raise ValueError("one valid four-sideband channel per strictly positive RF required")
        object.__setattr__(self, "channels", channels)

    def vacuum_covariances(self):
        return readonly_array([c.apply_covariance(np.eye(8)/2) for c in self.channels], real=True)


def sideband_channels(main, companion, analysis_axis):
    """Four-mode channels for positive RF samples and their reflected partners.

    Each point is the narrowband limit. Use average_spectral_channels with a
    normalized quadrature for a finite wavepacket and retain vacuum mode leakage.
    """
    _validate_pair(main, companion, analysis_axis)
    rf = analysis_axis.omega_rad_s
    positives = np.flatnonzero(rf > 0)
    if not positives.size:
        raise ValueError("strictly positive RF samples required; DC is a carrier, not four modes")
    channels = []
    for index in positives:
        mirror = _mirror_index(rf, index)
        a = _pair_channel(main, companion, index, mirror, SIDEBAND_MODES[:2])
        b = _pair_channel(main, companion, mirror, index, SIDEBAND_MODES[2:])
        x, y = np.zeros((8, 8)), np.zeros((8, 8))
        x[:4, :4], x[4:, 4:] = a.transfer, b.transfer
        y[:4, :4], y[4:, 4:] = a.added_covariance, b.added_covariance
        channels.append(GaussianChannel(x, y, SIDEBAND_MODES, SIDEBAND_MODES,
                                        "four optical sidebands; cross-sector moments absent by declared closure"))
    return SidebandChannels(AnalysisFrequencyAxis(rf[positives]), channels)


@dataclass(frozen=True)
class TopHatBand:
    """Normalized flat spectral mode: integral |h(omega)|^2 d omega/(2 pi)=1.

    On the positive band, h=1/sqrt(bandwidth_hz). Its negative partner is the
    reflected conjugate filter. Gauss-Legendre probabilities approximate df/B.
    """
    center_hz: float
    bandwidth_hz: float
    order: int

    def __post_init__(self):
        center, width = float(self.center_hz), float(self.bandwidth_hz)
        if not np.isfinite([center, width]).all() or not 0 < width < 2*center:
            raise ValueError("positive band width with support strictly above DC required")
        if isinstance(self.order, bool) or int(self.order) != self.order or self.order < 2:
            raise ValueError("integer quadrature order >=2 required")
        object.__setattr__(self, "center_hz", center)
        object.__setattr__(self, "bandwidth_hz", width)
        object.__setattr__(self, "order", int(self.order))

    def quadrature(self):
        nodes, weights = np.polynomial.legendre.leggauss(self.order)
        return (AnalysisFrequencyAxis.from_hz(self.center_hz+self.bandwidth_hz*nodes/2),
                readonly_array(weights/2, real=True))


def average_spectral_channels(channels, probabilities, *, source):
    """Matched real flat spectral modes; all orthogonal input modes are vacuum.

    Xbar=<X>, Yeff=<Y>+0.5< (X-Xbar)(X-Xbar)^T >. The last term accounts for
    coupling from unobserved spectral modes; simply averaging X and Y loses it.
    This is a quadrature approximation to the declared filter integral.
    """
    channels = tuple(channels)
    p = readonly_array(probabilities, real=True)
    if (not channels or p.shape != (len(channels),) or np.any(p < 0)
            or abs(p.sum()-1) > 1e-12):
        raise ValueError("normalized nonnegative spectral probabilities required")
    first = channels[0]
    if any(c.input_modes != first.input_modes or c.output_modes != first.output_modes
           or not c.audit().passed for c in channels):
        raise ValueError("valid matching spectral channel coordinates required")
    xs = np.stack([c.transfer for c in channels])
    x = np.einsum("f,fij->ij", p, xs)
    ys = np.stack([c.added_covariance for c in channels])
    differences = xs-x
    y = np.einsum("f,fij->ij", p, ys)+.5*np.einsum("f,fij,fkj->ik", p, differences, differences)
    result = GaussianChannel(x, y, first.input_modes, first.output_modes, source)
    if not result.audit().passed:
        raise ValueError("spectral-mode projection violates complete positivity")
    return result
