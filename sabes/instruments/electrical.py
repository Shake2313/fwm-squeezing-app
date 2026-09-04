"""Photodiode, oscilloscope, and spectrum-analyser models.

Scan traces use the steady-state sweep. Time traces sample an input PSD. The
spectrum-analyser trace combines detector noise and supplied source PSDs; it is
synthetic because the atomic model has no frequency-dependent covariance.
"""
from dataclasses import dataclass

import numpy as np

from .base import (COMPUTED, HEAD_OPTICAL, HEAD_PHOTOCURRENT,
                   PhotocurrentSignal, Quantity, Reading, SYNTHESISED,
                   Trace, shot_noise_a_per_rthz)


@dataclass(frozen=True)
class Photodiode:
    """Beam in, photocurrent out. The head the electrical instruments need."""
    head = HEAD_OPTICAL
    name: str = "Photodiode"
    responsivity_a_per_w: float = 0.58
    transimpedance_v_per_a: float = 1e5
    bandwidth_hz: float = 4.0e6
    nep_w_per_rthz: float = 3.3e-12
    saturation_w: float = 9.0e-3

    def convert(self, beam) -> PhotocurrentSignal:
        power = beam.total_power_w
        current = self.responsivity_a_per_w * power
        return PhotocurrentSignal(
            optical_power_w=power,
            responsivity_a_per_w=self.responsivity_a_per_w,
            dc_current_a=current,
            shot_noise_a_per_rthz=shot_noise_a_per_rthz(
                self.responsivity_a_per_w, power),
            electronic_noise_a_per_rthz=(self.responsivity_a_per_w
                                         * self.nep_w_per_rthz),
            transimpedance_v_per_a=self.transimpedance_v_per_a,
            bandwidth_hz=self.bandwidth_hz,
            label=self.name,
        )

    def measure(self, beam) -> Reading:
        signal = self.convert(beam)
        warnings = []
        if signal.optical_power_w > self.saturation_w:
            warnings.append(
                f"{signal.optical_power_w * 1e3:.1f} mW exceeds the "
                f"{self.saturation_w * 1e3:.0f} mW CW saturation of the amplifier")
        quantities = (
            Quantity("Optical power", signal.optical_power_w * 1e6, "µW"),
            Quantity("Photocurrent", signal.dc_current_a * 1e6, "µA"),
            Quantity("DC output", signal.dc_voltage_v, "V"),
            Quantity("Shot noise", signal.shot_noise_a_per_rthz * 1e12,
                     "pA/√Hz"),
            Quantity("Amplifier noise",
                     signal.electronic_noise_a_per_rthz * 1e12, "pA/√Hz"),
            Quantity("Clearance", signal.clearance_db, "dB",
                     "Shot noise above the amplifier floor; detector headroom "
                     "only, independent of any atomic squeezing claim."),
        )
        return Reading(self.name, quantities, None, tuple(warnings), COMPUTED)


def _detector_response(frequency_hz, bandwidth_hz):
    """Single-pole rolloff, in power. The only frequency shape we can defend."""
    return 1.0 / (1.0 + (np.asarray(frequency_hz, float) / bandwidth_hz) ** 2)


@dataclass(frozen=True)
class Oscilloscope:
    """Two modes, and they are not equally trustworthy -- see the module note."""
    head = HEAD_PHOTOCURRENT
    name: str = "Oscilloscope"

    def scan(self, signal, x, optical_power_w, *, x_label="Probe detuning",
             x_unit="GHz") -> Reading:
        """A swept-laser trace: what a scope actually shows in this experiment.

        `optical_power_w` is the power reaching the diode at each sweep point,
        which the Bloch solve supplies. Nothing here is reconstructed.
        """
        x = np.asarray(x, dtype=float)
        power = np.asarray(optical_power_w, dtype=float)
        volts = power * signal.responsivity_a_per_w * signal.transimpedance_v_per_a
        quantities = (
            Quantity("Sweep span", float(x.max() - x.min()), x_unit),
            Quantity("Peak", float(volts.max()), "V"),
            Quantity("Trough", float(volts.min()), "V"),
            Quantity("Contrast", float(volts.max() - volts.min()), "V"),
        )
        trace = Trace(x=x, series={"signal": volts}, x_label=x_label,
                      x_unit=x_unit, y_label="Detector output", y_unit="V")
        return Reading(self.name, quantities, trace, (), COMPUTED,
                       note="Swept-laser trace, computed from the Bloch solve.")

    def timeseries(self, signal, *, duration_s=2e-5, sample_rate_hz=5e7,
                   excess_a_per_rthz=0.0, seed=0) -> Reading:
        """Noise samples drawn from the PSD the model already knows.

        Looks like a scope and is statistically right, but it re-expresses the
        noise budget rather than deriving anything, so it is SYNTHESISED.
        """
        count = max(int(duration_s * sample_rate_hz), 16)
        rng = np.random.default_rng(seed)
        # White within the detector bandwidth: sigma = ASD * sqrt(BW).
        density = signal.total_noise_a_per_rthz(excess_a_per_rthz)
        sigma_a = density * np.sqrt(min(signal.bandwidth_hz,
                                        sample_rate_hz / 2.0))
        noise = rng.normal(0.0, sigma_a, count)
        volts = (signal.dc_current_a + noise) * signal.transimpedance_v_per_a
        time = np.arange(count) / sample_rate_hz

        quantities = (
            Quantity("DC level", signal.dc_voltage_v, "V"),
            Quantity("RMS noise", float(np.std(volts)) * 1e3, "mV"),
            Quantity("Noise density", density * 1e12, "pA/√Hz"),
            Quantity("Bandwidth", signal.bandwidth_hz / 1e6, "MHz"),
        )
        trace = Trace(x=time * 1e6, series={"signal": volts},
                      x_label="Time", x_unit="µs",
                      y_label="Detector output", y_unit="V")
        return Reading(
            self.name, quantities, trace, (), SYNTHESISED,
            note="Samples drawn from the modelled noise PSD — statistically "
                 "correct, but it says nothing the PSD did not already.")


@dataclass(frozen=True)
class SpectrumAnalyzer:
    """Detector-headroom view against RF frequency.

    Shot-noise and electronics traces apply to any photodiode signal. Twin-beam
    diagnostic traces are added only when ``gain_referred_noise_db`` is supplied.
    No atomic covariance is present, so this is not a squeezing spectrum.
    """
    head = HEAD_PHOTOCURRENT
    name: str = "Spectrum analyser"
    start_hz: float = 1.0e5
    stop_hz: float = 5.0e6
    points: int = 401
    resolution_bandwidth_hz: float = 3.0e4
    #: Optional illustrative low-frequency transfer corner; never evidence for
    #: atomic noise. Zero by default.
    technical_corner_hz: float = 0.0

    def analyze(self, signal, gain_referred_noise_db=None, *, total_power_w=None,
                pump_leakage_dbm=None, pump_rin_db_per_hz=None,
                eom_noise=None) -> Reading:
        """Build a detector spectrum from the supplied powers and PSDs.

        ``total_power_w`` sets balanced shot noise. ``pump_leakage_dbm`` is the
        observed detector power; pump RIN is added only when supplied.
        ``eom_noise`` is applied at the detector and is flat because its input is
        a single RIN value rather than a measured spectrum.
        """
        frequency = np.linspace(self.start_hz, self.stop_hz, self.points)
        response = _detector_response(frequency, signal.bandwidth_hz)

        power = total_power_w if total_power_w is not None else signal.optical_power_w
        shot_a = shot_noise_a_per_rthz(signal.responsivity_a_per_w, power)
        if eom_noise is not None:
            # Same-condition SQL includes the unpaired EOM modes reaching the
            # probe diode as well as the wanted twin beams.
            shot_a = np.sqrt(max(eom_noise.sql_psd_a2_per_hz, 0.0))
        bandwidth = self.resolution_bandwidth_hz
        gain = signal.transimpedance_v_per_a
        clearance = (20.0 * np.log10(shot_a / signal.electronic_noise_a_per_rthz)
                     if signal.electronic_noise_a_per_rthz > 0 else float("inf"))

        # Observed pump leakage: shot noise is calculable from the power. The
        # classical RIN term remains absent unless the caller explicitly supplies
        # a measured/predeclared spectral density.
        leak_shot_a = leak_rin_a = 0.0
        if pump_leakage_dbm is not None and np.isfinite(pump_leakage_dbm):
            leak_w = 1e-3 * 10.0 ** (float(pump_leakage_dbm) / 10.0)
            leak_current = signal.responsivity_a_per_w * leak_w
            leak_shot_a = shot_noise_a_per_rthz(signal.responsivity_a_per_w,
                                                leak_w)
            if pump_rin_db_per_hz is not None:
                leak_rin_a = leak_current * np.sqrt(
                    10.0 ** (float(pump_rin_db_per_hz) / 10.0))
        leak_a = float(np.hypot(leak_shot_a, leak_rin_a))

        # Noise POWER into 50 ohm, which is what an analyser displays.
        def to_dbm(current_asd):
            volts = current_asd * gain * np.sqrt(bandwidth)
            return 10.0 * np.log10(np.maximum(volts ** 2 / 50.0, 1e-30) / 1e-3)

        floor_a = float(np.hypot(signal.electronic_noise_a_per_rthz, leak_a))
        electronic = to_dbm(floor_a * np.sqrt(response))
        snl = to_dbm(shot_a * np.sqrt(response))

        quantities = [
            Quantity("Shot-noise level", float(snl[0]), "dBm"),
            Quantity("Noise floor", float(electronic[0]), "dBm",
                     "Detector electronics plus leaked-pump noise when supplied."),
            Quantity("Clearance", clearance, "dB"),
            Quantity("Resolution bandwidth", bandwidth / 1e3, "kHz"),
        ]
        if pump_leakage_dbm is not None:
            quantities.extend((
                Quantity("Pump leakage", float(pump_leakage_dbm), "dBm",
                         "Observed at the detector, not derived from geometry."),
                Quantity("Pump RIN", float(pump_rin_db_per_hz)
                         if pump_rin_db_per_hz is not None else float("nan"),
                         "dBc/Hz", "Applied only when supplied."),
            ))
        warnings = []
        series = {"shot-noise level": snl, "noise floor": electronic}
        if gain_referred_noise_db is not None:
            diagnostic = np.full_like(frequency,
                                      float(gain_referred_noise_db))
            if self.technical_corner_hz > 0:
                spoil = 1.0 / (1.0
                               + (self.technical_corner_hz / frequency) ** 2)
                diagnostic *= spoil
            diagnostic_level = snl + diagnostic
            shot_only_offset = diagnostic.copy()
            rin_loaded_offset = diagnostic.copy()
            if eom_noise is not None and eom_noise.sql_psd_a2_per_hz > 0.0:
                baseline = 10.0 ** (diagnostic / 10.0)
                shot_only = (
                    baseline * eom_noise.twin_shot_psd_a2_per_hz
                    + eom_noise.unwanted_shot_psd_a2_per_hz
                ) / eom_noise.sql_psd_a2_per_hz
                rin_loaded = (shot_only
                              + eom_noise.classical_rin_psd_a2_per_hz
                              / eom_noise.sql_psd_a2_per_hz)
                shot_only_offset = 10.0 * np.log10(
                    np.maximum(shot_only, 1e-300))
                rin_loaded_offset = 10.0 * np.log10(
                    np.maximum(rin_loaded, 1e-300))
            shot_only_level = snl + shot_only_offset
            rin_loaded_level = snl + rin_loaded_offset
            observed = 10.0 * np.log10(
                10 ** (rin_loaded_level / 10.0)
                + 10 ** (electronic / 10.0))

            quantities.insert(1, Quantity(
                "Gain-diagnostic level", float(diagnostic_level[0]), "dBm",
                "Algebraic overlay only; physical squeezing unavailable."))
            quantities.extend((
                Quantity("Diagnostic after electronics",
                         float(observed[0] - snl[0]), "dB",
                         "Detector-headroom overlay; not an atomic noise "
                         "spectrum."),
            ))
            series.update({"gain diagnostic": diagnostic_level,
                           "diagnostic + electronics": observed})

            if eom_noise is not None:
                quantities.extend((
                    Quantity("EOM unwanted-mode RIN", eom_noise.rin_db_per_hz,
                             "dBc/Hz"),
                    Quantity("EOM fractional intensity RMS",
                             100.0 * eom_noise.fractional_intensity_rms, "%"),
                    Quantity("Unwanted EOM power",
                             eom_noise.unwanted_detector_power_w * 1e9, "nW"),
                    Quantity("Classical EOM RIN excess",
                             eom_noise.classical_rin_excess_sql, "×SQL"),
                    Quantity("RIN-loaded diagnostic",
                             float(rin_loaded_offset[0]), "dB"),
                ))
                series.update({
                    "EOM shot-noise-only": shot_only_level,
                    "EOM RIN-loaded": rin_loaded_level,
                })

            deficit = clearance - abs(float(rin_loaded_offset[0]))
            if deficit < 6.0:
                warnings.append(
                    f"the diagnostic overlay sits {deficit:.1f} dB above the "
                    "noise floor; electronics dominate below about 6 dB")
        trace = Trace(
            x=frequency / 1e6,
            series=series,
            x_label="RF frequency", x_unit="MHz",
            y_label="Noise power", y_unit="dBm",
        )
        return Reading(
            self.name, tuple(quantities), trace, tuple(warnings), SYNTHESISED,
            note=("MEAN_FIELD_DIAGNOSTIC / PHYSICAL_SQUEEZING_UNAVAILABLE. "
                  "The RF shape is not measured or modeled."
                  if gain_referred_noise_db is not None else
                  "Computed detector noise levels; no twin-beam gain trace "
                  "applies at this beam."))
