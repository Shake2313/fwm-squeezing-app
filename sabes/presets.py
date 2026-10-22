"""Calculation profiles, separate from the dated apparatus inventory.

The lab reference contains measured/confirmed downstream inputs. Solved source
angles and effective EOM RF remain model inputs, never measured mount angles or
generator-to-port calibrations.
"""
from dataclasses import replace

from .beamline import SetupSettings, nominal_settings
from .calibration import default_calibration
from .detection import DetectionSettings

PAPER_PROFILE = "Sim 2025 model reference"
LAB_PROFILE = "Lab reference · 2026-10-16"
PROFILES = (LAB_PROFILE, PAPER_PROFILE)

LAB_PUMP_W = 0.294804
LAB_SEED_DETECTOR_W = 35.0 * (20.0 / 228.5) * 1e-6
LAB_CELL_OFF_TRANSMISSION = 0.8046
LAB_AFTER_CELL_TRANSMISSION = 0.9709
LAB_SEED_INPUT_PROXY_W = (LAB_SEED_DETECTOR_W
                          / LAB_CELL_OFF_TRANSMISSION
                          / LAB_AFTER_CELL_TRANSMISSION)


def profile_calibration(profile=PAPER_PROFILE):
    if profile not in PROFILES:
        raise ValueError(f"unknown SABES profile: {profile}")
    calibration = default_calibration()
    if profile == LAB_PROFILE:
        calibration = calibration.with_value(
            "t_noise_eater", 0.928, provenance="lab",
            note="Historical 2026-01-07 NEL03A/M transmission; retained as a "
                 "throughput estimate for the confirmed inserted/ON, "
                 "unoptimized device. Not a new 10/16 measurement.")
    return calibration


def power_targets(profile=PAPER_PROFILE):
    """Targets for inverse *model* source knobs, not measured optic angles."""
    if profile == LAB_PROFILE:
        return LAB_PUMP_W, 7e-3, LAB_SEED_INPUT_PROXY_W
    if profile == PAPER_PROFILE:
        return 0.600, 7e-3, 8e-6
    raise ValueError(f"unknown SABES profile: {profile}")


def profile_settings(profile=PAPER_PROFILE, calibration=None):
    calibration = calibration or profile_calibration(profile)
    if profile == PAPER_PROFILE:
        return SetupSettings(), DetectionSettings()
    if profile != LAB_PROFILE:
        raise ValueError(f"unknown SABES profile: {profile}")
    # The 11.8 dBm generator record does not establish EOM-port power. Leave
    # the legacy effective-port model setting unchanged and expose both planes.
    base = replace(
        SetupSettings(), cell_temp_c=113.0, opd_ghz=1.0,
        eom_frequency_hz=3.039e9,
        pump_waist_override_um=530.0, seed_waist_override_um=330.0)
    pump, eom, seed = power_targets(profile)
    settings = nominal_settings(calibration, pump_w=pump, eom_input_w=eom,
                                seed_w=seed, base=base)
    detection = DetectionSettings(
        post_cell_transmission_override=LAB_AFTER_CELL_TRANSMISSION,
        qe_pct_override=90.5)
    return settings, detection
