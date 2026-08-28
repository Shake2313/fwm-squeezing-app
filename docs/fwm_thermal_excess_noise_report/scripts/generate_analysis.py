"""Reproduce the numerical figures and tables for the FWM excess-noise report.

The script deliberately separates four evidence layers:

1. exact counting-statistics identities for coherent admixture;
2. phenomenological Fano-factor sensitivity inferred from published points;
3. deterministic GABES/SABES mean-response and optical-power ledgers; and
4. measured quantities transcribed from Sim, Kim, and Moon (2025).

Nothing in this script constructs the microscopic atomic Langevin diffusion
matrix needed for a physical squeezing spectrum.  In particular, the GABES
frequency response below is a deterministic susceptibility/gain diagnostic,
not S_-(Omega).
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import subprocess
import sys
from dataclasses import replace
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPORT_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = REPORT_ROOT.parents[1]
DATA_DIR = REPORT_ROOT / "data"
FIGURE_DIR = REPORT_ROOT / "figures"

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from gabes import constants, doppler, hyperfine, observables  # noqa: E402
from gabes.schemes import fwm  # noqa: E402
from sabes.beamline import SetupSettings, build_source_chain  # noqa: E402


COLORS = {
    "navy": "#183B56",
    "teal": "#207F79",
    "amber": "#A96500",
    "crimson": "#A12B3B",
    "gray": "#65727A",
    "light": "#E7EEF2",
}


def configure_plots() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 10.5,
            "axes.labelsize": 9.5,
            "legend.fontsize": 8.2,
            "figure.dpi": 150,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.04,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.22,
            "grid.linewidth": 0.6,
        }
    )


def db(value: np.ndarray | float) -> np.ndarray | float:
    return 10.0 * np.log10(np.maximum(value, np.finfo(float).tiny))


def git_metadata() -> dict[str, object]:
    def run(*args: str) -> str:
        return subprocess.check_output(
            ["git", *args],
            cwd=REPO_ROOT,
            text=True,
            encoding="utf-8",
            stderr=subprocess.DEVNULL,
        ).strip()

    try:
        commit = run("rev-parse", "HEAD")
        status = run("status", "--porcelain")
        tracked_status = run("status", "--porcelain", "--untracked-files=no")
        dirty = bool(status)
        tracked_dirty = bool(tracked_status)
        status_entry_count = len(status.splitlines()) if status else 0
        try:
            script_relative = Path(__file__).resolve().relative_to(REPO_ROOT).as_posix()
            run("ls-files", "--error-unmatch", script_relative)
            analysis_script_tracked = True
        except (OSError, subprocess.CalledProcessError, ValueError):
            analysis_script_tracked = False
    except (OSError, subprocess.CalledProcessError):
        commit = "unavailable"
        dirty = None
        tracked_dirty = None
        status_entry_count = None
        analysis_script_tracked = None
    return {
        "commit": commit,
        "worktree_dirty_before_output_writes": dirty,
        "tracked_files_dirty_before_output_writes": tracked_dirty,
        "status_entry_count_before_output_writes": status_entry_count,
        "analysis_script_tracked": analysis_script_tracked,
    }


def source_provenance() -> dict[str, object]:
    """Hash a conservative snapshot of every local Python runtime source.

    The analysis imports only a subset directly, but those modules have their
    own local imports and may acquire new ones over time.  Hashing the complete
    ``gabes`` and ``sabes`` Python trees is an intentional superset of the
    runtime import closure and records uncommitted edits as well.
    """

    paths = [Path(__file__).resolve()]
    for package in ("gabes", "sabes"):
        paths.extend(sorted((REPO_ROOT / package).rglob("*.py")))

    files: dict[str, object] = {}
    for path in paths:
        relative = path.relative_to(REPO_ROOT).as_posix()
        data = path.read_bytes()
        files[relative] = {
            "sha256": hashlib.sha256(data).hexdigest(),
            "bytes": len(data),
        }
    return {
        "scope": (
            "analysis script plus every .py file under the local gabes and "
            "sabes package trees (a conservative superset of the import closure)"
        ),
        "file_count": len(files),
        "files": files,
    }


def noise_mixing_analysis() -> dict[str, object]:
    # Sim Fig. 6 slope ratios.  The first point is labelled <<0.5% in the paper.
    sim_input_power_ratio = np.array([0.0, 0.125, 0.25])
    # The paper does not report post-cell detected shot weights for these modes.
    # Equality below is therefore an explicit proxy assumption, not a conversion.
    sim_detected_shot_weight_ratio_proxy = sim_input_power_ratio.copy()
    sim_R = np.array([0.134, 0.589, 0.964])
    R0 = float(sim_R[0])

    contaminant_ratio = np.linspace(0.0, 1.0, 401)
    input_ratio_proxy_coefficients = np.full(sim_R.shape, np.nan)
    mask = sim_detected_shot_weight_ratio_proxy > 0.0
    input_ratio_proxy_coefficients[mask] = (
        sim_R[mask] * (1.0 + sim_detected_shot_weight_ratio_proxy[mask]) - R0
    ) / sim_detected_shot_weight_ratio_proxy[mask]
    input_ratio_proxy_coefficient_mean = float(
        np.nanmean(input_ratio_proxy_coefficients)
    )

    fano_values = [1.0, 2.0, input_ratio_proxy_coefficient_mean]
    curves = {
        f"F={fano:.6g}": (R0 + fano * contaminant_ratio)
        / (1.0 + contaminant_ratio)
        for fano in fano_values
    }
    threshold_fano = np.full_like(contaminant_ratio, np.nan)
    positive = contaminant_ratio > 0.0
    threshold_fano[positive] = 1.0 + (1.0 - R0) / contaminant_ratio[positive]

    # Symmetric covariance-matched post-loss and gain/loss algebraic checks.
    G = 15.0
    ideal_noise = 1.0 / (2.0 * G - 1.0)
    eta_reported_total = 0.865
    eta_multiplicative = (1.0 - 0.055) * (1.0 - 0.080)
    detected_reported = 1.0 - eta_reported_total + eta_reported_total * ideal_noise
    detected_multiplicative = (
        1.0 - eta_multiplicative + eta_multiplicative * ideal_noise
    )

    direct_R0 = float(10.0 ** (-7.8 / 10.0))
    selected_ratios = np.array([0.005, 0.125, 0.25, 0.50, 1.0])
    selected_coherent = (direct_R0 + selected_ratios) / (1.0 + selected_ratios)

    return {
        "baseline_slope_ratio": R0,
        "baseline_slope_db": float(db(R0)),
        "sim_input_carrier_to_seed_power_ratio": sim_input_power_ratio,
        "sim_detected_shot_weight_ratio_proxy": sim_detected_shot_weight_ratio_proxy,
        "sim_ratio_proxy_assumption": (
            "input carrier/seed power ratio is set equal to the post-cell detected "
            "contaminant/wanted shot-weight ratio; mode-dependent transfer is unknown"
        ),
        "sim_R": sim_R,
        "sim_db": db(sim_R),
        "input_ratio_proxy_coefficients": input_ratio_proxy_coefficients,
        "input_ratio_proxy_coefficient_mean": input_ratio_proxy_coefficient_mean,
        "contaminant_ratio": contaminant_ratio,
        "curves": curves,
        "threshold_fano": threshold_fano,
        "selected_direct_ratios": selected_ratios,
        "selected_direct_coherent_R": selected_coherent,
        "selected_direct_coherent_db": db(selected_coherent),
        "gain_loss_check": {
            "gain": G,
            "ideal_R": ideal_noise,
            "ideal_db": float(db(ideal_noise)),
            "eta_reported_total": eta_reported_total,
            "detected_R_reported_total": detected_reported,
            "detected_db_reported_total": float(db(detected_reported)),
            "eta_multiplicative": eta_multiplicative,
            "detected_R_multiplicative": detected_multiplicative,
            "detected_db_multiplicative": float(db(detected_multiplicative)),
        },
    }


def vapor_analysis() -> dict[str, object]:
    temperature_c = np.array([90.0, 100.0, 110.0, 121.0, 131.0, 140.0, 150.0])
    temperature_k = temperature_c + 273.15
    density = np.array([hyperfine.number_density(T) for T in temperature_k])
    ground_collision = np.array(
        [
            constants.ground_coherence_dephasing(T, N, floor=0.0)
            for T, N in zip(temperature_k, density)
        ]
    )
    added_optical_full = np.array(
        [hyperfine.self_broadened_gamma(N) - constants.GAMMA for N in density]
    )
    added_optical_coherence = 0.5 * added_optical_full

    sigma_v = np.sqrt(constants.KB * temperature_k / constants.MASS_85RB)
    optical_doppler_fwhm = (
        2.0
        * np.sqrt(2.0 * np.log(2.0))
        * constants.NU_D1_85RB
        * sigma_v
        / constants.C_LIGHT
    )

    # Thermal-light and hyperfine energy-scale checks at the literature point.
    T_lit = 121.0 + 273.15
    wavelength = constants.WAVELENGTH_D1_85RB
    hbar = constants.HBAR
    h_planck = 2.0 * np.pi * hbar
    optical_x = h_planck * constants.C_LIGHT / (
        wavelength * constants.KB * T_lit
    )
    optical_occupation = 1.0 / math.expm1(optical_x)
    hyperfine_x = h_planck * constants.NU_GROUND_HF / (constants.KB * T_lit)

    probe_ghz = np.atleast_1d(
        0.9 - constants.NU_HF / 1e9 - 8.0e-3
    )
    k_pump, k_probe, _k_conj = fwm.seeded_option_a_wavenumbers(0.9, probe_ghz)
    raman = doppler.noncollinear_raman_rms_budget(
        T_lit, k_pump, float(k_probe[0]), math.radians(0.32)
    )

    sigma_lit = math.sqrt(constants.KB * T_lit / constants.MASS_85RB)
    mean_transverse_speed = math.sqrt(math.pi / 2.0) * sigma_lit
    transit = {}
    for waist_um in (330.0, 530.0):
        tau = waist_um * 1e-6 / mean_transverse_speed
        transit[str(int(waist_um))] = {
            "tau_s": tau,
            "corner_hz_1_over_2pi_tau": 1.0 / (2.0 * math.pi * tau),
        }

    return {
        "temperature_c": temperature_c,
        "temperature_k": temperature_k,
        "density_m3": density,
        "density_cm3": density / 1e6,
        "ground_collision_rad_s": ground_collision,
        "ground_collision_hz": ground_collision / (2.0 * np.pi),
        "inherited_transit_floor_hz": constants.GAMMA_GG / (2.0 * np.pi),
        "added_optical_full_hz": added_optical_full / (2.0 * np.pi),
        "added_optical_coherence_hwhm_hz": added_optical_coherence
        / (2.0 * np.pi),
        "sigma_velocity_m_s": sigma_v,
        "one_photon_doppler_fwhm_hz": optical_doppler_fwhm,
        "literature_point": {
            "temperature_K": T_lit,
            "optical_planck_exponent": optical_x,
            "optical_thermal_occupation": optical_occupation,
            "hyperfine_hnu_over_kT": hyperfine_x,
            "hyperfine_boltzmann_factor": math.exp(-hyperfine_x),
            "raman_doppler_rms_hz": raman["total_rms_hz"],
            "raman_doppler_axial_rms_hz": raman["axial_rms_hz"],
            "raman_doppler_transverse_rms_hz": raman["transverse_rms_hz"],
            "mean_transverse_speed_m_s": mean_transverse_speed,
            "transit_estimates": transit,
        },
    }


def eom_analysis() -> dict[str, object]:
    settings = SetupSettings()
    rf_dbm = np.linspace(6.0, 25.0, 77)
    wanted_unfiltered = []
    carrier_unfiltered = []
    wanted_cell = []
    carrier_cell = []
    other_cell = []
    purity_cell = []

    for rf in rf_dbm:
        chain = build_source_chain(replace(settings, eom_rf_dbm=float(rf)))
        eom_beam = next(
            stage.beam
            for stage in chain.stages
            if stage.path == "seed" and stage.name == "EOM sidebands"
        )
        wanted = eom_beam.power_at(chain.seed_offset_hz)
        carrier = eom_beam.power_at(0.0)
        wanted_unfiltered.append(wanted)
        carrier_unfiltered.append(carrier)
        wanted_cell.append(chain.wanted_seed_sideband_power_w)
        carrier_cell.append(chain.eom_residual_carrier_power_w)
        other_cell.append(chain.eom_other_sidebands_power_w)
        purity_cell.append(chain.seed_purity)

    wanted_unfiltered = np.asarray(wanted_unfiltered)
    carrier_unfiltered = np.asarray(carrier_unfiltered)
    wanted_cell = np.asarray(wanted_cell)
    carrier_cell = np.asarray(carrier_cell)
    other_cell = np.asarray(other_cell)
    purity_cell = np.asarray(purity_cell)

    default_chain = build_source_chain(replace(settings, eom_rf_dbm=18.0))
    default_eom = next(
        stage.beam
        for stage in default_chain.stages
        if stage.path == "seed" and stage.name == "EOM sidebands"
    )
    default_wanted_unfiltered = default_eom.power_at(default_chain.seed_offset_hz)

    return {
        "rf_dbm": rf_dbm,
        "wanted_unfiltered_w": wanted_unfiltered,
        "carrier_unfiltered_w": carrier_unfiltered,
        "carrier_ratio_unfiltered": carrier_unfiltered / wanted_unfiltered,
        "wanted_cell_w": wanted_cell,
        "carrier_cell_w": carrier_cell,
        "other_cell_w": other_cell,
        "carrier_ratio_cell": carrier_cell / wanted_cell,
        "other_ratio_cell": other_cell / wanted_cell,
        "purity_cell": purity_cell,
        "default_18_dbm": {
            "wanted_unfiltered_w": default_wanted_unfiltered,
            "carrier_ratio_unfiltered": default_eom.power_at(0.0)
            / default_wanted_unfiltered,
            "wanted_cell_w": default_chain.wanted_seed_sideband_power_w,
            "carrier_cell_w": default_chain.eom_residual_carrier_power_w,
            "other_cell_w": default_chain.eom_other_sidebands_power_w,
            "carrier_ratio_cell": default_chain.carrier_ratio,
            "other_ratio_cell": default_chain.eom_other_sidebands_to_wanted_ratio,
            "purity_cell": default_chain.seed_purity,
            "warnings": list(default_chain.warnings),
        },
        "status": (
            "deterministic Bessel-sideband and etalon power ledger only; "
            "no ECDL phase noise, RF noise, atomic transfer, or squeezing noise"
        ),
    }


def frequency_response_analysis() -> dict[str, object]:
    T = 121.0 + 273.15
    one_photon_ghz = 0.9
    two_photon_mhz = -8.0
    angle_deg = 0.32
    frequency_mhz = np.linspace(0.05, 4.50, 91)
    probe_ghz = np.atleast_1d(
        one_photon_ghz - constants.NU_HF / 1e9 + two_photon_mhz * 1e-3
    )
    k_pump, k_probe, k_conj = fwm.seeded_option_a_wavenumbers(
        one_photon_ghz, probe_ghz
    )
    pump_rabi = fwm.rabi_freq(0.6, fwm.W_PUMP)
    def calculate_reference(order: int):
        return fwm.pump_only_weak_response_noncollinear_reference(
            pump_rabi,
            pump_rabi,
            [2.0 * np.pi * two_photon_mhz * 1e6],
            2.0 * np.pi * one_photon_ghz * 1e9,
            T=T,
            pump_k_rad_m=k_pump,
            probe_k_axis_rad_m=k_probe,
            crossing_angle_rad=math.radians(angle_deg),
            analysis_frequency_axis_rad_s=2.0 * np.pi * frequency_mhz * 1e6,
            quadrature_order=order,
            cutoff_sigma=5.0,
        )

    quadrature_orders = (18, 24, 32)
    references = {order: calculate_reference(order) for order in quadrature_orders}
    reference = references[32]
    chi_reference = reference.chi_sc[:, 0]
    peak_chi = float(np.max(np.abs(chi_reference)))
    convergence = {}
    for order, candidate in references.items():
        chi = candidate.chi_sc[:, 0]
        difference = chi - chi_reference
        convergence[str(order)] = {
            "max_complex_error_relative_to_order32_peak": float(
                np.max(np.abs(difference)) / peak_chi
            ),
            "l2_complex_relative_error_to_order32": float(
                np.linalg.norm(difference) / np.linalg.norm(chi_reference)
            ),
            "max_abs_curve_error_relative_to_order32_peak": float(
                np.max(np.abs(np.abs(chi) - np.abs(chi_reference))) / peak_chi
            ),
        }

    density = hyperfine.number_density(T)
    line_strength = fwm.SEEDED_REFERENCE_RESIDUAL * fwm.physical_coupling_norm(-1)
    mismatch = fwm.seeded_phase_mismatch_z(
        one_photon_ghz, probe_ghz, angle_deg=angle_deg
    )
    gain_probe, gain_conjugate, _transfer = observables.gain_from_chi(
        reference.chi_ss[:, 0],
        reference.chi_sc[:, 0],
        reference.chi_cs[:, 0],
        reference.chi_cc[:, 0],
        k_probe,
        k_conj,
        fwm.L_CELL,
        density,
        line_strength=line_strength,
        delta_k_z=mismatch,
    )
    coupling = np.abs(reference.chi_sc[:, 0])

    return {
        "frequency_mhz": frequency_mhz,
        "chi_sc_abs_s": coupling,
        "chi_sc_normalized": coupling / coupling[0],
        "mean_field_probe_gain": gain_probe,
        "mean_field_conjugate_gain": gain_conjugate,
        "probe_gain_normalized": gain_probe / gain_probe[0],
        "conjugate_gain_normalized": gain_conjugate / gain_conjugate[0],
        "pump_rabi_2pi_ghz": pump_rabi / (2.0 * np.pi * 1e9),
        "raman_doppler_rms_mhz": float(
            reference.diagnostics["analytic_raman_rms_rad_s"][0]
            / (2.0 * np.pi * 1e6)
        ),
        "quadrature_raman_rms_mhz": float(
            reference.diagnostics["quadrature_raman_rms_rad_s"][0]
            / (2.0 * np.pi * 1e6)
        ),
        "diagnostics": {
            key: value
            for key, value in reference.diagnostics.items()
            if key
            not in {"analytic_raman_rms_rad_s", "quadrature_raman_rms_rad_s"}
        },
        "quadrature_convergence": {
            "orders_per_axis": list(quadrature_orders),
            "reference_order_per_axis": 32,
            "error_metrics": convergence,
        },
        "operating_point": {
            "temperature_C": 121.0,
            "one_photon_detuning_GHz": one_photon_ghz,
            "two_photon_detuning_MHz": two_photon_mhz,
            "crossing_angle_deg": angle_deg,
            "pump_power_mW": 600.0,
            "quadrature_order_per_axis": 32,
            "cutoff_sigma": 5.0,
        },
        "status": (
            "deterministic pump-only weak susceptibility and mean-field transfer; "
            "not a quantum-noise or squeezing spectrum"
        ),
    }


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_data(
    noise: dict[str, object],
    vapor: dict[str, object],
    eom: dict[str, object],
    response: dict[str, object],
) -> None:
    ratio = np.asarray(noise["contaminant_ratio"])
    curve_keys = list(noise["curves"])
    noise_rows = []
    for index, value in enumerate(ratio):
        row = {
            "contaminant_to_wanted_shot_weight_ratio": float(value),
            "threshold_fano_for_sql_crossing": float(noise["threshold_fano"][index])
            if value > 0.0
            else "",
        }
        for key in curve_keys:
            safe = key.replace("=", "_").replace(".", "p")
            row[f"R_{safe}"] = float(noise["curves"][key][index])
            row[f"dB_{safe}"] = float(db(noise["curves"][key][index]))
        noise_rows.append(row)
    noise_fields = list(noise_rows[0])
    write_csv(DATA_DIR / "noise_mixing.csv", noise_fields, noise_rows)

    vapor_rows = []
    for index, temp in enumerate(vapor["temperature_c"]):
        vapor_rows.append(
            {
                "temperature_C": float(temp),
                "density_m-3": float(vapor["density_m3"][index]),
                "density_cm-3": float(vapor["density_cm3"][index]),
                "ground_collision_Hz": float(vapor["ground_collision_hz"][index]),
                "inherited_transit_floor_Hz": float(
                    vapor["inherited_transit_floor_hz"]
                ),
                "added_optical_fullwidth_Hz": float(
                    vapor["added_optical_full_hz"][index]
                ),
                "added_optical_coherence_HWHM_Hz": float(
                    vapor["added_optical_coherence_hwhm_hz"][index]
                ),
                "sigma_velocity_m_s": float(vapor["sigma_velocity_m_s"][index]),
                "one_photon_Doppler_FWHM_Hz": float(
                    vapor["one_photon_doppler_fwhm_hz"][index]
                ),
            }
        )
    write_csv(DATA_DIR / "vapor_scaling.csv", list(vapor_rows[0]), vapor_rows)

    eom_rows = []
    for index, rf in enumerate(eom["rf_dbm"]):
        eom_rows.append(
            {
                "rf_dBm": float(rf),
                "wanted_unfiltered_W": float(eom["wanted_unfiltered_w"][index]),
                "carrier_to_wanted_unfiltered": float(
                    eom["carrier_ratio_unfiltered"][index]
                ),
                "wanted_at_cell_W": float(eom["wanted_cell_w"][index]),
                "carrier_at_cell_W": float(eom["carrier_cell_w"][index]),
                "other_sidebands_at_cell_W": float(eom["other_cell_w"][index]),
                "carrier_to_wanted_at_cell": float(eom["carrier_ratio_cell"][index]),
                "other_to_wanted_at_cell": float(eom["other_ratio_cell"][index]),
                "cell_spectral_purity": float(eom["purity_cell"][index]),
            }
        )
    write_csv(DATA_DIR / "eom_sweep.csv", list(eom_rows[0]), eom_rows)

    response_rows = []
    for index, frequency in enumerate(response["frequency_mhz"]):
        response_rows.append(
            {
                "analysis_frequency_MHz": float(frequency),
                "abs_chi_sc_s": float(response["chi_sc_abs_s"][index]),
                "abs_chi_sc_normalized": float(response["chi_sc_normalized"][index]),
                "mean_field_probe_gain": float(response["mean_field_probe_gain"][index]),
                "mean_field_conjugate_gain": float(
                    response["mean_field_conjugate_gain"][index]
                ),
                "probe_gain_normalized": float(response["probe_gain_normalized"][index]),
                "conjugate_gain_normalized": float(
                    response["conjugate_gain_normalized"][index]
                ),
            }
        )
    write_csv(
        DATA_DIR / "frequency_response.csv", list(response_rows[0]), response_rows
    )


def plot_noise_mixing(noise: dict[str, object]) -> None:
    ratio = 100.0 * np.asarray(noise["contaminant_ratio"])
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 3.05))

    styles = [
        ("F=1", COLORS["teal"], "-"),
        ("F=2", COLORS["amber"], "--"),
        (
            f"F={noise['input_ratio_proxy_coefficient_mean']:.6g}",
            COLORS["crimson"],
            "-",
        ),
    ]
    for key, color, linestyle in styles:
        axes[0].plot(
            ratio,
            db(noise["curves"][key]),
            color=color,
            linestyle=linestyle,
            linewidth=1.8,
            label=key,
        )
    axes[0].scatter(
        100.0 * np.asarray(noise["sim_detected_shot_weight_ratio_proxy"]),
        noise["sim_db"],
        s=32,
        marker="o",
        facecolor="white",
        edgecolor="black",
        linewidth=1.0,
        zorder=5,
        label="Sim slopes (input-ratio proxy)",
    )
    axes[0].axhline(0.0, color="black", linewidth=0.9)
    axes[0].set(
        xlabel="Contaminant / wanted shot weight [%]",
        ylabel="Normalized difference noise [dB]",
        xlim=(0.0, 100.0),
        ylim=(-9.5, 5.0),
        title="Coherent admixture cannot cross SQL",
    )
    axes[0].legend(loc="lower right", frameon=False)

    positive = ratio > 0.0
    axes[1].plot(
        ratio[positive],
        np.asarray(noise["threshold_fano"])[positive],
        color=COLORS["navy"],
        linewidth=1.9,
        label=r"$F_{\rm threshold}$",
    )
    valid = np.isfinite(noise["input_ratio_proxy_coefficients"])
    axes[1].scatter(
        100.0 * np.asarray(noise["sim_detected_shot_weight_ratio_proxy"])[valid],
        np.asarray(noise["input_ratio_proxy_coefficients"])[valid],
        color=COLORS["crimson"],
        s=38,
        zorder=4,
        label=r"$C_{proxy}$ from Sim slopes (not Fano)",
    )
    axes[1].axhline(1.0, color="black", linewidth=0.9, linestyle=":")
    axes[1].set(
        xlabel="Contaminant / wanted shot weight [%]",
        ylabel="Required contaminant Fano factor",
        xlim=(0.0, 100.0),
        ylim=(0.8, 10.0),
        title="Condition for super-SQL difference noise",
    )
    axes[1].legend(frameon=False)
    fig.tight_layout(w_pad=2.0)
    fig.savefig(FIGURE_DIR / "noise_mixing.pdf")
    plt.close(fig)


def plot_vapor_scaling(vapor: dict[str, object]) -> None:
    temp = np.asarray(vapor["temperature_c"])
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 3.05))
    axes[0].semilogy(
        temp,
        vapor["density_cm3"],
        "o-",
        color=COLORS["navy"],
        linewidth=1.8,
        markersize=4,
    )
    axes[0].axvline(121.0, color=COLORS["crimson"], linestyle="--", linewidth=1.0)
    axes[0].set(
        xlabel="Cell temperature [C]",
        ylabel=r"Pure-$^{85}$Rb density [cm$^{-3}$]",
        title="GABES CRC vapor-density model",
    )

    axes[1].semilogy(
        temp,
        np.asarray(vapor["ground_collision_hz"]) / 1e3,
        "o-",
        color=COLORS["teal"],
        linewidth=1.8,
        markersize=4,
        label="Rb-Rb ground collision / 2pi [kHz]",
    )
    axes[1].axhline(
        vapor["inherited_transit_floor_hz"] / 1e3,
        color=COLORS["gray"],
        linestyle=":",
        linewidth=1.5,
        label="inherited transit floor [kHz]",
    )
    axes[1].semilogy(
        temp,
        np.asarray(vapor["added_optical_coherence_hwhm_hz"]) / 1e3,
        "s-",
        color=COLORS["amber"],
        linewidth=1.7,
        markersize=3.8,
        label="optical coherence self-width [kHz]",
    )
    axes[1].axvline(121.0, color=COLORS["crimson"], linestyle="--", linewidth=1.0)
    axes[1].set(
        xlabel="Cell temperature [C]",
        ylabel="Rate / 2pi [kHz]",
        title="Implemented damping scales (not diffusion noise)",
    )
    axes[1].legend(frameon=False, loc="upper left")
    fig.tight_layout(w_pad=2.0)
    fig.savefig(FIGURE_DIR / "vapor_scaling.pdf")
    plt.close(fig)


def plot_eom(eom: dict[str, object]) -> None:
    rf = np.asarray(eom["rf_dbm"])
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 3.05))
    axes[0].semilogy(
        rf,
        eom["carrier_ratio_unfiltered"],
        color=COLORS["amber"],
        linewidth=1.8,
        label="carrier / wanted after EOM",
    )
    axes[0].semilogy(
        rf,
        eom["carrier_ratio_cell"],
        color=COLORS["teal"],
        linewidth=1.8,
        label="carrier / wanted at cell",
    )
    axes[0].semilogy(
        rf,
        eom["other_ratio_cell"],
        color=COLORS["navy"],
        linewidth=1.5,
        linestyle="--",
        label="other sidebands / wanted at cell",
    )
    axes[0].axvline(18.0, color=COLORS["crimson"], linestyle=":", linewidth=1.2)
    axes[0].set(
        xlabel="EOM RF drive [dBm]",
        ylabel="Power ratio",
        title="SABES deterministic spectral-power ledger",
    )
    axes[0].legend(frameon=False, loc="best")

    axes[1].plot(
        rf,
        1e6 * np.asarray(eom["wanted_cell_w"]),
        color=COLORS["navy"],
        linewidth=1.8,
        label="wanted seed",
    )
    axes[1].plot(
        rf,
        1e6 * np.asarray(eom["carrier_cell_w"]),
        color=COLORS["amber"],
        linewidth=1.6,
        label="residual carrier",
    )
    axes[1].plot(
        rf,
        1e6 * np.asarray(eom["other_cell_w"]),
        color=COLORS["teal"],
        linewidth=1.6,
        label="other sidebands",
    )
    axes[1].axvline(18.0, color=COLORS["crimson"], linestyle=":", linewidth=1.2)
    axes[1].set(
        xlabel="EOM RF drive [dBm]",
        ylabel="Cell-plane power [uW]",
        title="Etalon-filtered cell input",
    )
    axes[1].legend(frameon=False)
    fig.tight_layout(w_pad=2.0)
    fig.savefig(FIGURE_DIR / "eom_chain.pdf")
    plt.close(fig)


def plot_frequency_response(response: dict[str, object]) -> None:
    frequency = np.asarray(response["frequency_mhz"])
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 3.05))

    axes[0].plot(
        frequency,
        response["chi_sc_normalized"],
        color=COLORS["navy"],
        linewidth=1.9,
    )
    axes[1].plot(
        frequency,
        response["probe_gain_normalized"],
        color=COLORS["teal"],
        linewidth=1.8,
        label="probe gain / value at 0.05 MHz",
    )
    axes[1].plot(
        frequency,
        response["conjugate_gain_normalized"],
        color=COLORS["amber"],
        linewidth=1.5,
        linestyle="--",
        label="conjugate gain / value at 0.05 MHz",
    )
    for axis in axes:
        axis.axvline(
            response["raman_doppler_rms_mhz"],
            color=COLORS["gray"],
            linestyle=":",
            linewidth=1.1,
            label="Raman Doppler rms" if axis is axes[1] else None,
        )
        axis.axvline(2.0, color=COLORS["crimson"], linestyle="--", linewidth=1.0)
        axis.axvline(3.5, color="black", linestyle=":", linewidth=0.9)
        axis.set_xlabel("Analysis frequency [MHz]")
    axes[0].set(
        ylabel=r"$|\bar\chi_{sc}|/|\bar\chi_{sc}(0.05\,\mathrm{MHz})|$",
        title="2-D Maxwell-averaged atomic response",
    )
    axes[1].set(
        ylabel="Relative mean-field power gain",
        title="GABES transfer diagnostic, not squeezing",
    )
    axes[1].legend(frameon=False, loc="upper left")
    fig.tight_layout(w_pad=2.0)
    fig.savefig(FIGURE_DIR / "frequency_response.pdf")
    plt.close(fig)


def json_ready(value: object) -> object:
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, np.generic):
        return json_ready(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    return value


def main() -> None:
    # Capture input state before this run rewrites any generated artifact.
    repository_state = git_metadata()
    local_source_snapshot = source_provenance()
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    configure_plots()

    noise = noise_mixing_analysis()
    vapor = vapor_analysis()
    eom = eom_analysis()
    response = frequency_response_analysis()

    write_data(noise, vapor, eom, response)
    plot_noise_mixing(noise)
    plot_vapor_scaling(vapor)
    plot_eom(eom)
    plot_frequency_response(response)

    payload = {
        "generated_on": date.today().isoformat(),
        "repository": repository_state,
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "local_source_snapshot": local_source_snapshot,
        "claim_gate": {
            "physical_squeezing_prediction": False,
            "reason": (
                "GABES lacks a microscopic frequency-dependent atomic Langevin "
                "diffusion matrix and measured detector transfer functions"
            ),
            "allowed_claims": [
                "exact coherent-counting and symmetric covariance-matched post-loss identities",
                "phenomenological Fano-factor sensitivity under an explicit input-ratio proxy",
                "deterministic GABES atomic response and damping scales",
                "deterministic SABES EOM/filter power ledger",
            ],
        },
        "noise_mixing": noise,
        "vapor": vapor,
        "eom": eom,
        "frequency_response": response,
    }
    with (DATA_DIR / "analysis_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(
            json_ready(payload),
            handle,
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        )
        handle.write("\n")

    print(f"Wrote report analysis artifacts under {REPORT_ROOT}")
    print(
        "Claim gate: deterministic response and phenomenological sensitivity only; "
        "no physical squeezing spectrum."
    )


if __name__ == "__main__":
    main()
