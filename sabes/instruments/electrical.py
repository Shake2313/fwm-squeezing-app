"""
The instruments that need a photodiode: scope and spectrum analyser.

A photodiode turns a `BeamState` into a `PhotocurrentSignal`; the scope and the
analyser take that, exactly as on the bench. Neither can be pointed at a bare
beam, which is the right restriction rather than an inconvenience.

Honesty, restated because this is where it bites
------------------------------------------------
The model is steady state. There is no stochastic time evolution and no quantum
Langevin noise, so:

  * the oscilloscope's SCAN mode is `COMPUTED` -- a swept-laser trace is what a
    scope genuinely shows here, and the sweep comes from the Bloch solve;
  * its TIME-SERIES mode is `SYNTHESISED` -- samples drawn from a PSD the model
    already knows. Statistically correct, and informationally empty;
  * the spectrum analyser is `SYNTHESISED` -- assembled from the noise budget
    (shot noise, squeezed level, amplifier floor), not simulated.

The squeezing is drawn flat in frequency because that is all the model supports.
A real trace degrades at low frequency from technical noise; `technical_corner_hz`
exposes that as an explicit knob and defaults to zero, so the shape stays honest
unless someone deliberately puts a measured number in.
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
                     "Shot noise above the amplifier floor. A squeezed trace "
                     "sits |S| dB below shot noise, so this has to exceed |S| "
                     "with margin."),
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
    """Intensity-difference noise power against RF frequency.

    This is the instrument that measures squeezing, so it is the one whose
    honesty matters most. The three traces are the shot-noise level, the squeezed
    level `SNL + S_dB`, and the amplifier floor, each rolled off by the
    detector's single-pole response. It reproduces the SHAPE Sim et al. report;
    it is not a quantum-Langevin simulation, and it is labelled SYNTHESISED.
    """
    head = HEAD_PHOTOCURRENT
    name: str = "Spectrum analyser"
    start_hz: float = 1.0e5
    stop_hz: float = 5.0e6
    points: int = 401
    resolution_bandwidth_hz: float = 3.0e4
    #: Low-frequency technical-noise corner. Zero by default: a real trace does
    #: degrade towards DC, but inventing that shape would dress a guess up as a
    #: measurement.
    technical_corner_hz: float = 0.0

    def analyze(self, signal, squeezing_db, *, total_power_w=None,
                pump_leakage_dbm=None, pump_rin_db_per_hz=-140.0) -> Reading:
        """`total_power_w` overrides the power setting the shot-noise level.

        A balanced measurement sees the shot noise of BOTH arms while the head is
        one diode, so the caller needs to be able to say so. The clearance is
        recomputed from the same power: quoting a shot-noise level from one
        number and a clearance from another would be quietly inconsistent.

        `pump_leakage_dbm` is the OBSERVED pump power reaching the detector. It
        is not derived from beam geometry -- a Gaussian tail claims 1e-9
        rejection, which is useless and contradicted -- and it enters twice: its
        own shot noise, and the classical intensity noise it carries at
        `pump_rin_db_per_hz`. The second is the one that matters, because leaked
        pump is not common-mode and so survives the balanced subtraction.
        """
        frequency = np.linspace(self.start_hz, self.stop_hz, self.points)
        response = _detector_response(frequency, signal.bandwidth_hz)

        power = total_power_w if total_power_w is not None else signal.optical_power_w
        shot_a = shot_noise_a_per_rthz(signal.responsivity_a_per_w, power)
        bandwidth = self.resolution_bandwidth_hz
        gain = signal.transimpedance_v_per_a
        clearance = (20.0 * np.log10(shot_a / signal.electronic_noise_a_per_rthz)
                     if signal.electronic_noise_a_per_rthz > 0 else float("inf"))

        # Observed pump leakage: shot noise of its own photocurrent, plus the
        # classical RIN it carries. Only the latter is large.
        leak_shot_a = leak_rin_a = 0.0
        if pump_leakage_dbm is not None and np.isfinite(pump_leakage_dbm):
            leak_w = 1e-3 * 10.0 ** (float(pump_leakage_dbm) / 10.0)
            leak_current = signal.responsivity_a_per_w * leak_w
            leak_shot_a = shot_noise_a_per_rthz(signal.responsivity_a_per_w,
                                                leak_w)
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

        squeezing = np.full_like(frequency, float(squeezing_db))
        if self.technical_corner_hz > 0:
            # Squeezing degrades towards DC; the corner is an explicit input.
            spoil = 1.0 / (1.0 + (self.technical_corner_hz / frequency) ** 2)
            squeezing = squeezing * spoil
        squeezed = snl + squeezing

        # The analyser sees the squeezed light and the floor together, so the
        # observed trace is their incoherent sum -- which is exactly why
        # clearance limits what a squeezed measurement can show.
        observed = 10.0 * np.log10(10 ** (squeezed / 10.0)
                                   + 10 ** (electronic / 10.0))

        quantities = (
            Quantity("Shot-noise level", float(snl[0]), "dBm"),
            Quantity("Squeezed level", float(squeezed[0]), "dBm"),
            Quantity("Noise floor", float(electronic[0]), "dBm",
                     "Amplifier noise and leaked-pump noise together."),
            Quantity("Pump leakage", float(pump_leakage_dbm)
                     if pump_leakage_dbm is not None else float("nan"), "dBm",
                     "Observed at the detector, not derived from geometry."),
            Quantity("Clearance", clearance, "dB"),
            Quantity("Observable squeezing", float(observed[0] - snl[0]), "dB",
                     "What survives once the noise floor is added in, before "
                     "any electronic-noise subtraction."),
            Quantity("Resolution bandwidth", bandwidth / 1e3, "kHz"),
        )
        warnings = []
        deficit = clearance - abs(float(squeezing_db))
        if deficit < 6.0:
            warnings.append(
                f"the squeezed trace sits {deficit:.1f} dB above the noise "
                f"floor; below about 6 dB the electronic-noise subtraction "
                f"dominates the reported number")

        trace = Trace(
            x=frequency / 1e6,
            series={"shot-noise level": snl, "squeezed": squeezed,
                    "observed": observed, "noise floor": electronic},
            x_label="RF frequency", x_unit="MHz",
            y_label="Noise power", y_unit="dBm",
        )
        return Reading(
            self.name, quantities, trace, tuple(warnings), SYNTHESISED,
            note="Assembled from the noise budget (shot, squeezed, amplifier, "
                 "observed pump leakage), not simulated. Reproduces the measured "
                 "shape; it is not a quantum-Langevin calculation.")
