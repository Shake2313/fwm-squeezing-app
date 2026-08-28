"""
Maxwell velocity grid, Δ_eff axis construction, and Doppler averaging.

The Δ_eff trick (README "Speed"): the only velocity dependence of L₀ is the
excited-state diagonal shift, so all velocities share one L₀ template and differ
only by Δ_eff = Δ − k·v. R(δ, Δ_eff) is T- and Δ-independent. Ported verbatim;
mass / k_vec are now arguments (defaulting to the 85Rb values) so other isotopes
or geometries can reuse the routines.
"""
import numpy as np

from . import constants


def maxwell_legendre_grid(T, mass=constants.MASS_85RB, order=24,
                          cutoff_sigma=5.0):
    """Truncated 1-D Maxwell quadrature with explicit convergence controls.

    Gauss--Legendre nodes cover ``[-cutoff_sigma*sigma,+cutoff_sigma*sigma]``;
    the Gaussian density is folded into the weights and the retained interval is
    normalized to one.  Unlike :func:`velocity_grid`, ``order`` and cutoff are
    independent, which makes grid and tail-convergence audits unambiguous.
    """
    T = float(T)
    mass = float(mass)
    order = int(order)
    cutoff_sigma = float(cutoff_sigma)
    if not np.isfinite(T) or T <= 0.0:
        raise ValueError("T must be finite and positive")
    if not np.isfinite(mass) or mass <= 0.0:
        raise ValueError("mass must be finite and positive")
    if order < 2:
        raise ValueError("order must be at least 2")
    if not np.isfinite(cutoff_sigma) or cutoff_sigma <= 0.0:
        raise ValueError("cutoff_sigma must be finite and positive")

    sigma = np.sqrt(constants.KB * T / mass)
    nodes, legendre_weights = np.polynomial.legendre.leggauss(order)
    limit = cutoff_sigma * sigma
    velocity = limit * nodes
    pdf = np.exp(-0.5 * (velocity / sigma)**2) / (
        np.sqrt(2.0 * np.pi) * sigma)
    weights = legendre_weights * limit * pdf
    weights /= np.sum(weights)
    return velocity, weights


def noncollinear_raman_shift_rad_s(
        vx_m_s, vz_m_s, pump_k_rad_m, probe_k_rad_m, crossing_angle_rad):
    """Velocity shift of the atomic two-photon detuning in angular units.

    The pump is along ``+z`` and the probe lies in the ``x-z`` plane.  The
    returned quantity is added to the laboratory two-photon detuning:

    ``(k_pump-k_probe*cos(theta))*v_z-k_probe*sin(theta)*v_x``.

    This function does not accept or return the laboratory optical beat; that
    independent frequency must remain invariant under a velocity-grid change.
    """
    vx = np.asarray(vx_m_s, dtype=float)
    vz = np.asarray(vz_m_s, dtype=float)
    k_pump = float(pump_k_rad_m)
    k_probe = np.asarray(probe_k_rad_m, dtype=float)
    theta = float(crossing_angle_rad)
    if not (np.isfinite(vx).all() and np.isfinite(vz).all()
            and np.isfinite(k_pump) and np.isfinite(k_probe).all()
            and np.isfinite(theta)):
        raise ValueError("velocities, wave numbers, and angle must be finite")
    return ((k_pump - k_probe * np.cos(theta)) * vz
            - k_probe * np.sin(theta) * vx)


def noncollinear_atomic_detunings_rad_s(
        one_photon_lab_rad_s, two_photon_lab_rad_s, vx_m_s, vz_m_s,
        pump_k_rad_m, probe_k_rad_m, crossing_angle_rad):
    """Return ``(Delta_eff, delta_eff)`` without modifying any lab frequency."""
    one_photon = np.asarray(one_photon_lab_rad_s, dtype=float)
    two_photon = np.asarray(two_photon_lab_rad_s, dtype=float)
    vz = np.asarray(vz_m_s, dtype=float)
    if not (np.isfinite(one_photon).all() and np.isfinite(two_photon).all()):
        raise ValueError("laboratory detunings must be finite")
    delta_eff = two_photon + noncollinear_raman_shift_rad_s(
        vx_m_s, vz, pump_k_rad_m, probe_k_rad_m, crossing_angle_rad)
    Delta_eff = one_photon - float(pump_k_rad_m) * vz
    return Delta_eff, delta_eff


def noncollinear_raman_rms_budget(
        T, pump_k_rad_m, probe_k_rad_m, crossing_angle_rad,
        mass=constants.MASS_85RB):
    """Analytic full-Maxwell rms Raman-Doppler budget in rad/s and Hz."""
    T = float(T)
    mass = float(mass)
    k_pump = float(pump_k_rad_m)
    k_probe = float(probe_k_rad_m)
    theta = float(crossing_angle_rad)
    if not all(np.isfinite(value) for value in (
            T, mass, k_pump, k_probe, theta)):
        raise ValueError("Raman-Doppler inputs must be finite")
    if T <= 0.0 or mass <= 0.0:
        raise ValueError("T and mass must be positive")
    sigma_v = np.sqrt(constants.KB * T / mass)
    axial = abs(k_pump - k_probe * np.cos(theta)) * sigma_v
    transverse = abs(k_probe * np.sin(theta)) * sigma_v
    total = float(np.hypot(axial, transverse))
    return {
        "sigma_velocity_m_s": float(sigma_v),
        "axial_rms_rad_s": float(axial),
        "transverse_rms_rad_s": float(transverse),
        "total_rms_rad_s": total,
        "axial_rms_hz": float(axial / (2.0 * np.pi)),
        "transverse_rms_hz": float(transverse / (2.0 * np.pi)),
        "total_rms_hz": float(total / (2.0 * np.pi)),
    }


def velocity_grid(T, mass=constants.MASS_85RB, dv=1.0, cutoff_sigma=3.0):
    """1-D Maxwell velocity classes and normalised weights at temperature T."""
    sigma = np.sqrt(constants.KB * T / mass)
    v_limit = np.ceil(cutoff_sigma * sigma / dv) * dv
    v = np.arange(-v_limit, v_limit + 0.5 * dv, dv)
    pdf = np.exp(-v**2 / (2 * sigma**2)) / (np.sqrt(2 * np.pi) * sigma)
    w = pdf * dv
    return v, w / w.sum()


def build_Delta_eff_axis(Delta_min, Delta_max, v_grid, k_vec=constants.K_VEC):
    """Δ_eff sample axis covering Δ ∈ [Delta_min, Delta_max] over all velocities."""
    dv = v_grid[1] - v_grid[0]
    step = k_vec * dv
    lo = Delta_min - k_vec * v_grid.max()
    hi = Delta_max - k_vec * v_grid.min()
    n = int(np.ceil((hi - lo) / step)) + 1
    return np.linspace(lo, hi, n)


def doppler_average(chi_table, Delta_eff_axis, Delta, v_grid, weights,
                    k_vec=constants.K_VEC):
    """
    Σ_v weights(v) · χ_table[δ, Δ_eff = Δ − k·v] for every δ row.
    Linear interpolation along the Δ_eff axis.
    """
    idx_lo, frac = interpolation_weights(Delta_eff_axis, Delta, v_grid, k_vec)
    return apply_doppler_average(chi_table, idx_lo, frac, weights)


def interpolation_weights(Delta_eff_axis, Delta, v_grid, k_vec=constants.K_VEC):
    """
    Lower indices and interpolation fractions for Delta_eff = Delta - k*v.

    Current Delta_eff axes are uniformly sampled, so direct index arithmetic
    avoids rebuilding an arange axis and calling np.interp for every average.
    """
    n_de = Delta_eff_axis.size
    step = (Delta_eff_axis[-1] - Delta_eff_axis[0]) / (n_de - 1)
    deff_v = Delta - k_vec * v_grid
    idx_float = (deff_v - Delta_eff_axis[0]) / step
    idx_float = np.clip(idx_float, 0.0, n_de - 1)
    idx_lo = np.clip(np.floor(idx_float).astype(int), 0, n_de - 2)
    frac = idx_float - idx_lo
    return idx_lo, frac


def apply_doppler_average(chi_table, idx_lo, frac, weights):
    """Apply precomputed Doppler interpolation weights to every row of a table."""
    frac = frac.astype(chi_table.dtype, copy=False)
    lo_part = chi_table[:, idx_lo]
    hi_part = chi_table[:, idx_lo + 1]
    interp = lo_part * (1 - frac)[None, :] + hi_part * frac[None, :]
    return interp @ weights


def doppler_average_1d(chi_axis, Delta_eff_axis, Delta_axis, v_grid, weights,
                       k_vec=constants.K_VEC):
    """
    Doppler-average one chi(Delta_eff) table over many detuning samples at once.

    This is the scan-independent-H path used by OD / bare Voigt calculations:
    one fine Delta_eff table, many detunings, same velocity grid.
    """
    chi_axis = np.asarray(chi_axis)
    idx_lo, frac = interpolation_weights(
        Delta_eff_axis, np.asarray(Delta_axis)[:, None], v_grid[None, :], k_vec)
    frac = frac.astype(chi_axis.dtype, copy=False)
    interp = chi_axis[idx_lo] * (1 - frac) + chi_axis[idx_lo + 1] * frac
    return interp @ weights
