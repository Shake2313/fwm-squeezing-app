"""Explicit power-gain reference planes for the seeded atomic-medium map.

The Maxwell map contains atomic propagation; it does not contain cell-window
or external optical losses.  Pump-on/pump-off output ratios therefore require
the pump-off *atomic* transmission separately from external throughputs.  A
measured total-cell transmission cannot identify that atomic transmission.
"""

import numpy as np


def _power_array(value, name):
    value = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(value)) or np.any(value < 0.0):
        raise ValueError(f"{name} must be finite and non-negative")
    return value


def _transmission(value, name, *, positive=False):
    value = _power_array(value, name)
    if np.any(value > 1.0) or (positive and np.any(value <= 0.0)):
        interval = "(0, 1]" if positive else "[0, 1]"
        raise ValueError(f"{name} must be a power transmission in {interval}")
    return value


def pump_off_referenced_gains(
        G_s, G_c, *, pump_off_atomic_transmission=None,
        probe_output_transmission=1.0, conjugate_output_transmission=1.0,
        reference_output_transmission=1.0):
    """Refer atomic-medium power gains to a measured pump-off seed output.

    ``G_s`` and ``G_c`` are output *powers* at the vapor exit divided by seed
    power at the vapor entrance.  The output transmissions describe each path
    from that exit to its power-measurement plane, including output windows,
    collection and subsequent optics.  The reference path measures pump-off
    seed power.  On/off inputs must have the same vapor-entry seed power.

    With atomic pump-off transmission A0 and output throughputs tau_j,

        G_j,ref = G_j * tau_j / (A0 * tau_ref).

    Common output throughputs cancel.  The entrance window is already upstream
    of the declared seed input plane and cancels in an on/off ratio.  Electronic
    detector gains/responsivities must be calibrated before supplying powers.
    Photon-flux gains must be converted to powers before calling this helper.

    A0=None returns unavailable ratios rather than guessing A0 from a total-cell
    transmission.  The seed-power factor converts measured pump-off reference
    power to the required vapor-entry seed power.  Scalar/array inputs broadcast.
    This conversion neither changes the atomic response nor validates it.
    """
    probe = _power_array(G_s, "G_s")
    conjugate = _power_array(G_c, "G_c")
    tau_probe = _transmission(
        probe_output_transmission, "probe_output_transmission")
    tau_conjugate = _transmission(
        conjugate_output_transmission, "conjugate_output_transmission")
    tau_reference = _transmission(
        reference_output_transmission, "reference_output_transmission",
        positive=True)
    common = (probe, conjugate, tau_probe, tau_conjugate, tau_reference)
    provenance = {
        "gain_quantity": "power ratio, not photon-flux ratio",
        "source_input_plane": "vapor entrance after entrance window",
        "source_output_plane": "vapor exit before output window",
        "reference_plane": "measured pump-off seed output",
        "formula": "G_j_ref=G_j*tau_j/(A0*tau_ref)",
        "output_transmissions_include": (
            "output windows, collection and subsequent optical throughput"),
        "entrance_window": "cancels for equal vapor-entry on/off seed power",
        "total_cell_transmission_identifies_A0": False,
        "detector_response": "powers require prior responsivity calibration",
        "atomic_response_modified": False,
    }
    if pump_off_atomic_transmission is None:
        # Still reject incompatible supplied shapes even when A0 is unknown.
        np.broadcast_arrays(*common)
        return {
            "status": "UNAVAILABLE",
            "reason": "pump-off atomic transmission not supplied",
            "probe_power_gain": None,
            "conjugate_power_gain": None,
            "medium_input_seed_per_reference_power": None,
            "pump_off_atomic_transmission": None,
            "provenance": provenance,
        }

    atomic = _transmission(
        pump_off_atomic_transmission, "pump_off_atomic_transmission",
        positive=True)
    probe, conjugate, tau_probe, tau_conjugate, tau_reference, atomic = (
        np.broadcast_arrays(*common, atomic))
    reference_per_input = atomic * tau_reference
    return {
        "status": "CONDITIONAL_ON_SUPPLIED_TRANSMISSIONS",
        "reason": None,
        "probe_power_gain": probe * tau_probe / reference_per_input,
        "conjugate_power_gain": conjugate * tau_conjugate / reference_per_input,
        "medium_input_seed_per_reference_power": 1.0 / reference_per_input,
        "pump_off_atomic_transmission": atomic,
        "provenance": provenance,
    }
