"""EOM-sideband noise at the balanced detector.

The source/etalon chain determines how much power is present in every optical
frequency mode. This module applies an effective relative-intensity-noise PSD
(RIN) to the unwanted EOM modes.

For each post-etalon mode the small-signal field convention is
``E_n(t) = sqrt(P_n) [1 + epsilon_n(t)/2] exp(-i omega_n t)`` with
``PSD(epsilon_n) = RIN``.  Thus ``epsilon`` is the fractional-power fluctuation
and its field-amplitude fluctuation is one half as large.

This is not an atomic or microscopic EOM model. The wanted -1
sideband is the seeded twin-beam input; every other mode is treated as unpaired
light in the probe detector.  Its Poisson shot noise and a positive classical
``I^2 RIN`` term enter the balanced difference channel.
"""
from dataclasses import dataclass
from typing import Tuple

import math


ELEMENTARY_CHARGE_C = 1.602176634e-19
_MATCH_TOLERANCE_HZ = 1.0


def _db_to_linear(value_db):
    return 10.0 ** (float(value_db) / 10.0)


def _linear_to_db(value):
    return 10.0 * math.log10(max(float(value), 1.0e-300))


@dataclass(frozen=True)
class EOMModeNoise:
    """One unwanted EOM mode's detector-plane noise contribution."""
    order: int
    label: str
    offset_hz: float
    cell_power_w: float
    detector_power_w: float
    photocurrent_a: float
    shot_psd_a2_per_hz: float
    classical_rin_psd_a2_per_hz: float

    @property
    def thermal_like_psd_a2_per_hz(self):
        """Compatibility alias for :attr:`classical_rin_psd_a2_per_hz`."""
        return self.classical_rin_psd_a2_per_hz

    @property
    def fano_factor(self):
        """Mode noise divided by its Poisson shot-noise PSD."""
        if self.shot_psd_a2_per_hz <= 0.0:
            return 1.0
        return 1.0 + self.classical_rin_psd_a2_per_hz / self.shot_psd_a2_per_hz


@dataclass(frozen=True)
class EOMNoiseBudget:
    """Balanced-detector noise terms for unwanted post-etalon EOM modes.

    The normalization is the same-condition SQL: shot noise from the desired
    twins plus the unpaired unwanted light.  For a sub-SQL baseline ``R0``,

        R = (R0 S_twin + S_unwanted + S_RIN)
            / (S_twin + S_unwanted)

    Unwanted-mode shot noise moves a sub-SQL result toward one. Positive RIN
    raises it further.
    """
    modes: Tuple[EOMModeNoise, ...]
    desired_twin_detector_power_w: float
    unwanted_detector_power_w: float
    rin_db_per_hz: float
    rin_per_hz: float
    analysis_bandwidth_hz: float
    twin_shot_psd_a2_per_hz: float
    unwanted_shot_psd_a2_per_hz: float
    classical_rin_psd_a2_per_hz: float

    @property
    def thermal_like_psd_a2_per_hz(self):
        """Compatibility alias for :attr:`classical_rin_psd_a2_per_hz`."""
        return self.classical_rin_psd_a2_per_hz

    @property
    def sql_psd_a2_per_hz(self):
        return self.twin_shot_psd_a2_per_hz + self.unwanted_shot_psd_a2_per_hz

    @property
    def fractional_intensity_rms(self):
        """Declared unwanted-mode fractional-power RMS over analysis bandwidth."""
        return math.sqrt(self.rin_per_hz * self.analysis_bandwidth_hz)

    @property
    def fractional_field_amplitude_rms(self):
        """Small-signal field-amplitude RMS (half the fractional-power RMS)."""
        return 0.5 * self.fractional_intensity_rms

    @property
    def classical_rin_excess_sql(self):
        sql = self.sql_psd_a2_per_hz
        return self.classical_rin_psd_a2_per_hz / sql if sql > 0.0 else 0.0

    @property
    def thermal_like_excess_sql(self):
        """Compatibility alias for :attr:`classical_rin_excess_sql`."""
        return self.classical_rin_excess_sql

    def shot_noise_only_linear(self, baseline_db):
        """Poissonian unwanted-light penalty, normalized to same-condition SQL."""
        baseline = _db_to_linear(baseline_db)
        sql = self.sql_psd_a2_per_hz
        if sql <= 0.0:
            return baseline
        return ((baseline * self.twin_shot_psd_a2_per_hz
                 + self.unwanted_shot_psd_a2_per_hz) / sql)

    def rin_loaded_linear(self, baseline_db):
        """Shot-noise-only result plus classical EOM RIN excess."""
        sql = self.sql_psd_a2_per_hz
        if sql <= 0.0:
            return _db_to_linear(baseline_db)
        return (self.shot_noise_only_linear(baseline_db)
                + self.classical_rin_psd_a2_per_hz / sql)

    def shot_noise_only_db(self, baseline_db):
        return _linear_to_db(self.shot_noise_only_linear(baseline_db))

    def rin_loaded_db(self, baseline_db):
        return _linear_to_db(self.rin_loaded_linear(baseline_db))

    def rin_penalty_db(self, baseline_db):
        """Positive degradation relative to a sub-SQL baseline.

        Above SQL the same-condition normalization can move a ratio downward
        even though both unwanted shot noise and the classical RIN PSD
        are positive.  Calling that normalized shift a ``penalty`` would be
        misleading, so this convenience method is intentionally sub-SQL-only.
        """
        if float(baseline_db) > 0.0:
            raise ValueError("rin_penalty_db is defined only for sub-SQL baselines")
        return self.rin_loaded_db(baseline_db) - float(baseline_db)

    # Compatibility names retained for saved notebooks and downstream callers.
    def coherent_mixture_linear(self, baseline_db):
        return self.shot_noise_only_linear(baseline_db)

    def loaded_linear(self, baseline_db):
        return self.rin_loaded_linear(baseline_db)

    def coherent_mixture_db(self, baseline_db):
        return self.shot_noise_only_db(baseline_db)

    def loaded_db(self, baseline_db):
        return self.rin_loaded_db(baseline_db)

    def penalty_db(self, baseline_db):
        return self.rin_penalty_db(baseline_db)


def build_eom_noise_budget(chain, gains, optics_transmission,
                           responsivity_a_per_w, rin_db_per_hz,
                           analysis_bandwidth_hz, *, probe_transmission=None,
                           conjugate_transmission=None):
    """Build an independent-mode EOM noise budget at the photodiodes.

    One effective RIN is applied independently to each unwanted mode. The
    reduced atomic solver has no transfer law for those modes, so their powers
    pass through the probe arm without gain or a conjugate partner.
    """
    transmission = float(optics_transmission)
    probe_transmission = (transmission if probe_transmission is None else
                          float(probe_transmission))
    conjugate_transmission = (
        transmission if conjugate_transmission is None else
        float(conjugate_transmission))
    responsivity = float(responsivity_a_per_w)
    bandwidth = float(analysis_bandwidth_hz)
    if not 0.0 <= transmission <= 1.0:
        raise ValueError("optics_transmission must lie between zero and one")
    if not 0.0 <= probe_transmission <= 1.0:
        raise ValueError("probe_transmission must lie between zero and one")
    if not 0.0 <= conjugate_transmission <= 1.0:
        raise ValueError("conjugate_transmission must lie between zero and one")
    if responsivity <= 0.0 or bandwidth <= 0.0:
        raise ValueError("responsivity and analysis bandwidth must be positive")

    rin_db = float(rin_db_per_hz)
    if math.isnan(rin_db):
        raise ValueError("rin_db_per_hz cannot be NaN")
    rin = 0.0 if rin_db == -math.inf else _db_to_linear(rin_db)

    gain_s, gain_c = (max(float(value), 0.0) for value in gains)
    desired_detector_power = (
        chain.wanted_seed_sideband_power_w
        * (gain_s * probe_transmission + gain_c * conjugate_transmission))
    twin_current = responsivity * desired_detector_power
    twin_shot = 2.0 * ELEMENTARY_CHARGE_C * twin_current

    eom_frequency_hz = abs(chain.seed_offset_hz)
    modes = []
    for line in chain.seed.lines:
        if abs(line.offset_hz - chain.seed_offset_hz) <= _MATCH_TOLERANCE_HZ:
            continue
        detector_power = max(float(line.power_w), 0.0) * probe_transmission
        if detector_power <= 0.0:
            continue
        current = responsivity * detector_power
        shot_psd = 2.0 * ELEMENTARY_CHARGE_C * current
        classical_rin_psd = current ** 2 * rin
        order = (int(round(line.offset_hz / eom_frequency_hz))
                 if eom_frequency_hz > 0.0 else 0)
        modes.append(EOMModeNoise(
            order=order,
            label=line.label,
            offset_hz=line.offset_hz,
            cell_power_w=line.power_w,
            detector_power_w=detector_power,
            photocurrent_a=current,
            shot_psd_a2_per_hz=shot_psd,
            classical_rin_psd_a2_per_hz=classical_rin_psd,
        ))

    return EOMNoiseBudget(
        modes=tuple(modes),
        desired_twin_detector_power_w=desired_detector_power,
        unwanted_detector_power_w=sum(mode.detector_power_w for mode in modes),
        rin_db_per_hz=rin_db,
        rin_per_hz=rin,
        analysis_bandwidth_hz=bandwidth,
        twin_shot_psd_a2_per_hz=twin_shot,
        unwanted_shot_psd_a2_per_hz=sum(
            mode.shot_psd_a2_per_hz for mode in modes),
        classical_rin_psd_a2_per_hz=sum(
            mode.classical_rin_psd_a2_per_hz for mode in modes),
    )
