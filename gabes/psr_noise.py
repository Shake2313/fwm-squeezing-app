"""Effective Gaussian polarization-squeezing readout for an atomic PSR slope.

The classical self-rotation slope supplies a uniform shear. Vacuum absorption
and optional isotropic excess diffusion form an explicit effective noise model;
they do not derive the driven atom's microscopic Langevin noise.
"""

import numpy as np
from scipy.special import gammainc

from .quantum.channels import (
    GaussianChannel, compose_channels, covariance_uncertainty_minimum,
)


_MODE = ("orthogonal_polarization",)
_MODEL_ID = "effective_psr_uniform_shear_distributed_loss_v1"


def _distributed_channel(shear, optical_depth, excess_noise_per_od):
    """Exact unit-length channel; canonical vacuum covariance is I/2."""
    k, a = float(shear), float(optical_depth)
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        transfer = np.exp(-a/2)*np.array([[1.0, k], [0.0, 1.0]])
        if a == 0.0:
            added = np.zeros((2, 2))
        else:
            # H_j = a * integral_0^1 exp(-a*t) t^j dt. For small a the
            # convergent Taylor series avoids subtracting nearly equal terms.
            if a < 1e-3:
                moments = np.zeros(3)
                term = 1.0
                for order in range(8):
                    moments += term/(np.arange(3)+order+1)
                    term *= -a/(order+1)
                h0 = a*moments[0]
                kh1 = k*a*moments[1]
                k2h2 = (k*np.sqrt(a*moments[2]))**2
            else:
                # Scaling by k/a before squaring remains finite at large OD,
                # including a shear that grows proportionally to cell length.
                ratio = k/a
                h0 = gammainc(1, a)
                kh1 = ratio*gammainc(2, a)
                k2h2 = 2*ratio*ratio*gammainc(3, a)
            added = .5*(1+excess_noise_per_od)*np.array(
                [[h0+k2h2, kh1], [kh1, h0]])
    return GaussianChannel(
        transfer, added, _MODE, _MODE,
        "uniform PSR shear plus vacuum absorption and assumed isotropic excess diffusion")


def psr_gaussian_readout(
        shear, optical_depth, *, excess_noise_per_od=0.0,
        detection_efficiency=1.0, phase_degrees=None):
    """Return SQL-normalized polarization noise and orthogonal phase traces.

    Numeric model inputs broadcast. Scalar inputs return a 2x2 covariance and
    one-dimensional traces; an input shape S returns covariances S+(2,2),
    extrema S and phase traces S+(number_of_phases,). Arrays are fresh.

    At normalized position z in [0,1], drift A=[[-OD/2,shear],[0,-OD/2]] and
    diffusion D=OD*(1+excess_noise_per_od)*I/2 define the canonical covariance
    equation dV/dz=A V+V A.T+D, with input vacuum V=I/2. The shear is gL and
    optical_depth is alpha L. The exact channel integrates both distributed
    reservoirs before an external vacuum-loss detector channel is applied once.

    If ellipticity epsilon and polarization angle phi are both small angles,
    shear=d(phi_out)/d(epsilon_in), with x_out proportional to phi_out and
    p_in proportional to epsilon_in. Input X-P noise is uncorrelated vacuum.
    Positive classical rotation slope is a shear x -> x+shear*p. Zero excess
    noise is an optimistic omission of atomic excess noise, not its prediction.
    """
    values = np.broadcast_arrays(*[
        np.asarray(value, dtype=float) for value in (
            shear, optical_depth, excess_noise_per_od, detection_efficiency)])
    k, od, excess, eta = values
    if not k.size or not all(np.all(np.isfinite(value)) for value in values):
        raise ValueError("nonempty finite PSR model inputs required")
    if np.any(od < 0) or np.any(excess < 0):
        raise ValueError("optical_depth and excess_noise_per_od must be nonnegative")
    if np.any((eta < 0) | (eta > 1)):
        raise ValueError("detection_efficiency must be in [0,1]")
    phase = (np.linspace(0.0, 180.0, 361) if phase_degrees is None
             else np.asarray(phase_degrees, dtype=float))
    if phase.ndim != 1 or not phase.size or not np.all(np.isfinite(phase)):
        raise ValueError("phase_degrees must be a nonempty finite vector")

    shape = k.shape
    source_covariance = np.empty(shape+(2, 2))
    covariance = np.empty_like(source_covariance)
    source_transfer = np.empty_like(source_covariance)
    source_added = np.empty_like(source_covariance)
    detected_transfer = np.empty_like(source_covariance)
    detected_added = np.empty_like(source_covariance)
    cp_minimum, uncertainty_minimum = np.empty(shape), np.empty(shape)
    audit_passed = np.empty(shape, dtype=bool)
    for index in np.ndindex(shape):
        try:
            source = _distributed_channel(k[index], od[index], excess[index])
        except FloatingPointError as error:
            raise ValueError("PSR channel exceeds floating-point range") from error
        detector = GaussianChannel.vacuum_attenuator(
            _MODE, [float(eta[index])], source="external detection efficiency applied once")
        combined = compose_channels(source, detector)
        source_v = source.apply_covariance(np.eye(2)/2)
        detected_v = combined.apply_covariance(np.eye(2)/2)
        audit = combined.audit()
        source_covariance[index], covariance[index] = 2*source_v, 2*detected_v
        source_transfer[index], source_added[index] = source.transfer, source.added_covariance
        detected_transfer[index], detected_added[index] = combined.transfer, combined.added_covariance
        cp_minimum[index] = audit.minimum_cp_eigenvalue
        uncertainty_minimum[index] = covariance_uncertainty_minimum(detected_v)
        audit_passed[index] = audit.passed

    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    if not np.all(np.isfinite(eigenvalues)) or np.any(eigenvalues <= 0):
        raise ValueError("PSR readout has nonpositive or unrepresentable variance; no repair applied")
    minimum, maximum = eigenvalues[..., 0], eigenvalues[..., 1]
    min_angle = np.degrees(np.arctan2(eigenvectors[..., 1, 0], eigenvectors[..., 0, 0])) % 180
    theta = np.deg2rad(phase)
    rows_a = np.column_stack((np.cos(theta), np.sin(theta)))
    rows_b = np.column_stack((-np.sin(theta), np.cos(theta)))
    noise_a = np.einsum("ni,...ij,nj->...n", rows_a, covariance, rows_a)
    noise_b = np.einsum("ni,...ij,nj->...n", rows_b, covariance, rows_b)
    crossing = (minimum+maximum)/2
    return {
        "model_id": _MODEL_ID,
        "shear": k.copy(), "optical_depth": od.copy(),
        "excess_noise_per_od": excess.copy(), "detection_efficiency": eta.copy(),
        "phase_degrees": phase.copy(),
        "source_covariance": source_covariance, "covariance": covariance,
        "source_transfer": source_transfer, "source_added_covariance": source_added,
        "detected_transfer": detected_transfer, "detected_added_covariance": detected_added,
        "noise_a": noise_a, "noise_b": noise_b,
        "noise_a_db": 10*np.log10(noise_a), "noise_b_db": 10*np.log10(noise_b),
        "minimum_variance": minimum, "maximum_variance": maximum,
        "minimum_db": 10*np.log10(minimum), "maximum_db": 10*np.log10(maximum),
        "min_angle_deg": min_angle,
        "crossing_variance": crossing, "crossing_db": 10*np.log10(crossing),
        "crossing_phase_degrees": (min_angle[..., None]+[45.0, 135.0]) % 180,
        "channel_audit_passed": audit_passed,
        "channel_cp_minimum": cp_minimum, "uncertainty_minimum": uncertainty_minimum,
        "channel_matrix_convention": "canonical [x,p]=i; channel added covariance uses vacuum I/2",
        "convention": (
            "returned source/detected covariances and phase noise use SQL=I; "
            "A=x*cos(theta)+p*sin(theta); B=A(theta+90 deg); "
            "theta relative to the bright carrier's amplitude/ellipticity axes; "
            "noise dB=10*log10(variance/SQL)"),
        "stokes_convention": (
            "bright x-polarized carrier with mean photon number N: "
            "delta S2=sqrt(2N)*x, delta S3=sqrt(2N)*p, "
            "[S2,S3] approximately 2iN; Var(S_theta)/N equals returned noise"),
        "model_formula": (
            "A=[[-OD/2,shear],[0,-OD/2]], D=OD*(1+excess_noise_per_od)*I/2; "
            "X=exp(A), Y=integral_0^1 exp(A*t)*D*exp(A.T*t) dt; "
            "Vout=X*(I/2)*X.T+Y; Vdet=eta*Vout+(1-eta)*I/2"),
        "assumptions": [
            "Uniform effective PSR shear and absorption across the cell; orthogonal input vacuum",
            "Atomic excess noise is an assumed nonnegative isotropic diffusion per absorption OD",
            "Zero excess noise optimistically omits atomic excess noise; it is not microscopic validation",
            "Bright-carrier Stokes linearization; detection phase does not change the atomic state",
        ],
        "scope": (
            "Effective Gaussian estimate from the classical PSR slope, distributed vacuum absorption "
            "and declared excess-noise assumption; not a microscopic atomic-noise or finite-RF "
            "spectrum prediction. Channel consistency does not validate the noise assumption."),
    }
