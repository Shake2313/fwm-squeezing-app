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
    deff_v = Delta - k_vec * v_grid
    idx_float = np.interp(deff_v, Delta_eff_axis, np.arange(Delta_eff_axis.size))
    n_de = Delta_eff_axis.size
    idx_lo = np.clip(np.floor(idx_float).astype(int), 0, n_de - 2)
    frac = (idx_float - idx_lo).astype(chi_table.dtype)
    lo_part = chi_table[:, idx_lo]
    hi_part = chi_table[:, idx_lo + 1]
    interp = lo_part * (1 - frac)[None, :] + hi_part * frac[None, :]
    return interp @ weights
