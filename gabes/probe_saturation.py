"""Finite-power Gaussian-probe saturation helpers for alkali D-line OD/SAS.

The numerical model follows ``references/AutoOD-NatRbD2``:

* ``I0 = 2 P / (pi w^2)`` for a Gaussian 1/e^2 intensity waist;
* each hyperfine transition has a resonant two-level saturation parameter;
* the homogeneous width grows as ``gamma * sqrt(1 + s)`` while the original
  numerator is retained, so the integrated absorption falls as
  ``1 / sqrt(1 + s)``;
* ``dI/dz = -alpha(I) I`` is propagated self-consistently along the cell and
  power-averaged over the full Gaussian cross-section.

The Gaussian intensity and propagation layer is species independent.  The
OD/SAS scheme supplies each line's optical frequency, natural/effective width,
and transition-resolved weak-probe profiles.  Rb D2 lies inside the supplied
calculator's isotope-composition scope; D1 and Cs reuse the same closed-two-level
construction only as explicitly labelled generalized estimates.  The scheme
also documents the separate approximation used when a counter-propagating pump
is present.
"""

from __future__ import annotations

import numpy as np

from . import constants


SATURATION_LOG_STEP = 0.1
RADIAL_QUADRATURE_ORDER = 32
NUMERICAL_TRANSPARENT_OD = 1e-8


def gaussian_peak_intensity(power_mw, waist_mm):
    """Return on-axis Gaussian intensity [W/m^2] for a 1/e^2 waist."""
    power_mw = float(power_mw)
    waist_mm = float(waist_mm)
    if not np.isfinite(power_mw) or power_mw < 0.0:
        raise ValueError("Probe power must be finite and non-negative.")
    if not np.isfinite(waist_mm) or waist_mm <= 0.0:
        raise ValueError("Probe waist must be finite and positive.")
    waist_m = waist_mm * 1e-3
    return 2.0 * power_mw * 1e-3 / (np.pi * waist_m**2)


def closed_two_level_saturation_intensity(nu_hz, gamma_rad_s):
    """Return ``pi*h*c*Gamma/(3*lambda^3)`` in W/m^2.

    ``Gamma`` is the angular natural linewidth.  This is the closed-two-level
    reference used by the Rb-D2 AutoOD implementation; hyperfine line strengths
    are applied separately by the caller.  For open D1 transitions and species
    outside that reference, the caller must label the result as a generalized
    estimate rather than a validated multilevel saturation model.
    """
    nu_hz = float(nu_hz)
    gamma_rad_s = float(gamma_rad_s)
    if not np.isfinite(nu_hz) or nu_hz <= 0.0:
        raise ValueError("Optical frequency must be finite and positive.")
    if not np.isfinite(gamma_rad_s) or gamma_rad_s <= 0.0:
        raise ValueError("Natural linewidth must be finite and positive.")
    wavelength = constants.C_LIGHT / nu_hz
    planck = 2.0 * np.pi * constants.HBAR
    return (np.pi * planck * constants.C_LIGHT * gamma_rad_s
            / (3.0 * wavelength**3))


def peak_saturation_parameter(
        peak_intensity, saturation_intensity, gamma_natural, gamma_effective):
    """Resonant strongest-line saturation ``s0`` at the Gaussian beam centre."""
    peak_intensity = float(peak_intensity)
    saturation_intensity = float(saturation_intensity)
    gamma_natural = float(gamma_natural)
    gamma_effective = float(gamma_effective)
    values = (peak_intensity, saturation_intensity, gamma_natural, gamma_effective)
    if any(not np.isfinite(value) for value in values):
        raise ValueError("Saturation inputs must be finite.")
    if peak_intensity < 0.0:
        raise ValueError("Peak intensity must be non-negative.")
    if saturation_intensity <= 0.0 or gamma_natural <= 0.0 or gamma_effective <= 0.0:
        raise ValueError("Saturation intensity and linewidths must be positive.")
    return peak_intensity * gamma_natural / (
        saturation_intensity * gamma_effective)


def saturation_log_grid(max_saturation):
    """Grid in ``x = -log(I/I_peak)`` resolving saturation into the weak tail."""
    max_saturation = float(max_saturation)
    if not np.isfinite(max_saturation) or max_saturation < 0.0:
        raise ValueError("Maximum saturation must be finite and non-negative.")
    x_max = max(12.0, np.log1p(max_saturation) + 10.0)
    points = max(97, int(np.ceil(x_max / SATURATION_LOG_STEP)) + 1)
    return np.linspace(0.0, x_max, points)


def power_broadened_profile_table(
        weak_profile, spacing, hwhm, peak_saturation, intensity_log_grid):
    """Return one local absorption profile per log-intensity sample.

    A Lorentzian of HWHM ``hwhm`` broadened by saturation ``s`` has HWHM
    ``hwhm*sqrt(1+s)`` and area ``1/sqrt(1+s)`` relative to the weak profile.
    Convolution with a normalized Lorentzian whose HWHM is the difference of
    those widths performs that transformation exactly for a weak Voigt profile.
    Applying it to a pump-prepared transition is the OD/SAS scheme's explicit
    separable pump/probe approximation.
    """
    profile = np.asarray(weak_profile, dtype=float)
    x_grid = np.asarray(intensity_log_grid, dtype=float)
    spacing = float(spacing)
    hwhm = float(hwhm)
    peak_saturation = float(peak_saturation)
    if profile.ndim != 1 or profile.size < 2:
        raise ValueError("Weak profile must be a one-dimensional sampled curve.")
    if x_grid.ndim != 1 or x_grid.size < 1 or np.any(~np.isfinite(x_grid)):
        raise ValueError("Intensity log grid must be a finite one-dimensional array.")
    if not np.isfinite(spacing) or spacing <= 0.0:
        raise ValueError("Profile spacing must be finite and positive.")
    if not np.isfinite(hwhm) or hwhm <= 0.0:
        raise ValueError("Homogeneous HWHM must be finite and positive.")
    if not np.isfinite(peak_saturation) or peak_saturation < 0.0:
        raise ValueError("Peak saturation must be finite and non-negative.")
    if np.any(~np.isfinite(profile)):
        raise ValueError("Weak profile must be finite.")
    if not np.any(profile):
        return np.zeros((x_grid.size, profile.size), dtype=float)

    saturation = peak_saturation * np.exp(-x_grid)
    width_ratio = np.sqrt(1.0 + saturation)
    extra_hwhm = hwhm * (width_ratio - 1.0)
    if not np.any(extra_hwhm > 0.0):
        return np.broadcast_to(profile, (x_grid.size, profile.size)).copy()

    # A Lorentzian of HWHM ``a`` has Fourier multiplier exp(-2*pi*a*|f|).
    # Applying that multiplier avoids a sampled narrow kernel collapsing to a
    # grid delta when the added width is smaller than one scan interval.  Centre
    # the curve in a >=4x zero-padded domain so its long Lorentzian tail cannot
    # wrap around into the opposite scan edge.
    n = profile.size
    fft_size = 1 << (4 * n - 1).bit_length()
    start = (fft_size - n) // 2
    padded = np.zeros(fft_size, dtype=float)
    padded[start:start + n] = profile
    frequencies = np.fft.rfftfreq(fft_size, d=spacing)
    multiplier = np.exp(
        -2.0 * np.pi * extra_hwhm[:, None] * frequencies[None, :])
    broadened_padded = np.fft.irfft(
        np.fft.rfft(padded)[None, :] * multiplier,
        n=fft_size, axis=1)
    broadened = broadened_padded[:, start:start + n]
    return broadened / width_ratio[:, None]


def _interpolate_log_intensity_table(x_grid, table, x_values):
    """Interpolate a shared-x table for several intensity coordinates."""
    x_grid = np.asarray(x_grid, dtype=float)
    table = np.asarray(table, dtype=float)
    x_values = np.asarray(x_values, dtype=float)
    indices = np.searchsorted(x_grid, x_values, side="right") - 1
    indices = np.clip(indices, 0, x_grid.size - 2)
    low = x_grid[indices]
    high = x_grid[indices + 1]
    fraction = (x_values - low) / (high - low)
    return table[indices] + fraction[:, None] * (
        table[indices + 1] - table[indices])


def propagate_gaussian_spectrum(
        intensity_log_grid, alpha_table, alpha_weak, length_m,
        radial_order=RADIAL_QUADRATURE_ORDER):
    """Power-average nonlinear absorption over radius and cell length.

    For ``u = I_in(r)/I_peak = exp(-2r^2/w^2)``, incident power is uniform in
    ``u``.  Along a ray, ``x = -log(I/I_peak)`` obeys ``dz = dx/alpha(x)``;
    integrating and inverting this coordinate avoids unstable tiny longitudinal
    steps at high optical depth.  Returns transmission and
    ``effective_alpha = -log(transmission)/length``.
    """
    x_grid = np.asarray(intensity_log_grid, dtype=float)
    alpha_table = np.asarray(alpha_table, dtype=float)
    alpha_weak = np.asarray(alpha_weak, dtype=float)
    length_m = float(length_m)
    radial_order = int(radial_order)
    if length_m < 0.0 or not np.isfinite(length_m):
        raise ValueError("Propagation length must be finite and non-negative.")
    if x_grid.ndim != 1 or x_grid.size < 2 or np.any(np.diff(x_grid) <= 0.0):
        raise ValueError("Intensity grid must be strictly increasing.")
    if alpha_table.shape != (x_grid.size, alpha_weak.size):
        raise ValueError("Absorption table shape is inconsistent with its axes.")
    if np.any(~np.isfinite(alpha_table)) or np.any(alpha_table < 0.0):
        raise ValueError("Absorption table must be finite and non-negative.")
    if np.any(~np.isfinite(alpha_weak)) or np.any(alpha_weak < 0.0):
        raise ValueError("Weak-probe absorption must be finite and non-negative.")
    if radial_order < 1:
        raise ValueError("Radial quadrature order must be positive.")

    nodes, weights = np.polynomial.legendre.leggauss(radial_order)
    input_fractions = 0.5 * (nodes + 1.0)
    weights = 0.5 * weights
    weights /= weights.sum()
    x_input = -np.log(input_fractions)
    entrance_alpha = _interpolate_log_intensity_table(
        x_grid, alpha_table, x_input)
    entrance_effective_alpha = weights @ entrance_alpha
    max_local_alpha = np.maximum(np.max(alpha_table, axis=0), alpha_weak)
    if length_m == 0.0:
        return np.ones(alpha_weak.size), entrance_effective_alpha

    # This is also the continuous L -> 0 limit.  It avoids forcing a numerical
    # alpha floor above a genuinely tiny curve when every possible ray optical
    # depth is already negligible.
    small_optical_depth = (
        max_local_alpha * length_m <= NUMERICAL_TRANSPARENT_OD)
    if np.all(small_optical_depth):
        return (np.exp(-entrance_effective_alpha * length_m),
                entrance_effective_alpha)

    # Exact zeros can occur at a passivity-limited transparency point.  Using a
    # machine-tiny floor makes the cumulative distance O(1e153 m), at which
    # adding an ordinary cell length is lost to floating-point cancellation and
    # can invert to a spuriously huge absorption.  A floor corresponding to only
    # 1e-8 optical depth across this cell is physically negligible, yet keeps the
    # distance coordinate numerically resolvable.  The rigorous maximum-local-
    # alpha bound below removes even that numerical floor from transparent tails.
    alpha_floor = NUMERICAL_TRANSPARENT_OD / length_m
    inverse_alpha = 1.0 / np.maximum(alpha_table, alpha_floor)
    intervals = np.diff(x_grid)[:, None]
    distance_grid = np.vstack((
        np.zeros((1, alpha_weak.size), dtype=float),
        np.cumsum(0.5 * (inverse_alpha[:-1] + inverse_alpha[1:]) * intervals,
                  axis=0),
    ))
    input_distance = _interpolate_log_intensity_table(
        x_grid, distance_grid, x_input)
    target_distance = input_distance + length_m
    x_output = np.empty_like(target_distance)

    for frequency_index in range(alpha_weak.size):
        distance_curve = distance_grid[:, frequency_index]
        targets = target_distance[:, frequency_index]
        inside = targets <= distance_curve[-1]
        x_output[inside, frequency_index] = np.interp(
            targets[inside], distance_curve, x_grid)
        if np.any(~inside):
            tail_distance = targets[~inside] - distance_curve[-1]
            x_output[~inside, frequency_index] = (
                x_grid[-1] + tail_distance * alpha_weak[frequency_index])

    log_ray_transmission = -np.maximum(x_output - x_input[:, None], 0.0)
    log_weighted = np.log(weights)[:, None] + log_ray_transmission
    log_scale = np.max(log_weighted, axis=0)
    log_transmission = log_scale + np.log(
        np.sum(np.exp(log_weighted - log_scale[None, :]), axis=0))
    log_transmission = np.minimum(log_transmission, 0.0)
    log_transmission = np.maximum(
        log_transmission, -max_local_alpha * length_m)
    log_transmission[small_optical_depth] = (
        -entrance_effective_alpha[small_optical_depth] * length_m)
    transparent = (alpha_weak == 0.0) & np.all(alpha_table == 0.0, axis=0)
    log_transmission[transparent] = 0.0
    transmission = np.exp(log_transmission)
    effective_alpha = -log_transmission / length_m
    effective_alpha[transparent] = 0.0
    return transmission, effective_alpha
