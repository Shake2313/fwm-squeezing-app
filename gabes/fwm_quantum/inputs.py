"""Power-to-readout adapter and evidence bound to the actual consumed SI values."""

from dataclasses import dataclass, fields

import numpy as np

from .. import constants as c
from ..quantum.contracts import InputAudit, OpticalDetunings, audit_independent_inputs
from .model import reduced_pump_noise
from .normalization import optical_carriers, reduced_dipoles
from .readout import reduced_seeded_readout


INPUT_UNITS = {
    "pump_power_W": "W", "pump_waist_m": "m", "number_density_m3": "m^-3",
    "uniform_area_m2": "m^2", "length_m": "m", "seed_power_W": "W",
    "one_photon_rad_s": "rad/s", "two_photon_rad_s": "rad/s",
    "transit_rate_s_inverse": "s^-1", "phase_mismatch_rad_m": "rad/m",
    "seed_phase_rad": "rad",
}


@dataclass(frozen=True)
class ReducedPowerInputs:
    pump_power_W: float
    pump_waist_m: float
    number_density_m3: float
    uniform_area_m2: float
    length_m: float
    seed_power_W: float
    one_photon_rad_s: float
    two_photon_rad_s: float
    transit_rate_s_inverse: float
    phase_mismatch_rad_m: float = 0.
    seed_phase_rad: float = 0.

    def __post_init__(self):
        for item in fields(self):
            value = float(getattr(self, item.name))
            if not np.isfinite(value):
                raise ValueError("finite SI inputs required")
            object.__setattr__(self, item.name, value)
        for key in ("pump_waist_m", "uniform_area_m2", "seed_power_W"):
            if getattr(self, key) <= 0:
                raise ValueError(f"{key} must be positive")
        for key in ("pump_power_W", "number_density_m3", "length_m", "transit_rate_s_inverse"):
            if getattr(self, key) < 0:
                raise ValueError(f"{key} must be nonnegative")

    @property
    def detunings(self):
        return OpticalDetunings(self.one_photon_rad_s, self.two_photon_rad_s)

    def consumed_inputs(self, detector, *, convention):
        """Scalar ledger including detector samples; units are exact SI labels.

        Density is supplied directly, without an implicit T->density law. Reset
        is supplied directly, without inferring collisions or beam transit from w.
        Zero phase/flat response are model assumptions requiring evidence too.
        """
        coupling = reduced_dipoles(convention)
        entries = {k: (getattr(self, k), unit) for k, unit in INPUT_UNITS.items()}
        entries.update({"natural_decay": (c.GAMMA, "s^-1"),
            "D1_centroid": (c.NU_D1_85RB, "Hz"),
            "ground_hyperfine_splitting": (c.NU_GROUND_HF, "Hz"),
            "excited_hyperfine_splitting": (c.NU_EXCITED_HF_D1, "Hz"),
            "reset_F2_population": (5/12, "1"), "reset_F3_population": (7/12, "1"),
            "detector_balance": (detector.balance, "1")})
        if convention == "legacy-reciprocal":
            entries["stored_D1_dipole"] = (c.DIPOLE_D1, "C m")
        for j, beam in enumerate(("probe", "conjugate")):
            entries[f"{beam}_transmission"] = (float(detector.transmissions[j]), "1")
        for k, freq in enumerate(detector.analysis_axis.frequency_hz):
            entries[f"rf_{k}"] = (float(freq), "Hz")
            for j, beam in enumerate(("probe", "conjugate")):
                response = detector.current_response[k, j]
                entries[f"{beam}_response_real_{k}"] = (float(response.real), "A/A")
                entries[f"{beam}_response_imag_{k}"] = (float(response.imag), "A/A")
            entries[f"electronics_psd_{k}"] = (float(detector.electronics_difference_psd_A2_Hz[k]), "A^2/Hz")
        return {"scalars": entries, "convention": convention,
            "derived": {"base_dipole_C_m": coupling.base_dipole_C_m,
                "dipole_matrix_C_m": coupling.matrix_C_m.tolist(),
                "pump_rabi_rad_s": coupling.pump_rabi_rad_s(self.pump_power_W, self.pump_waist_m),
                "optical_carriers_rad_s": optical_carriers(self.detunings).tolist()},
            "assumptions": [coupling.scope, "single velocity and uniform weak-field area",
                "Gaussian pump peak sampled uniformly; no transverse average",
                "Markov replacement to declared thermal populations",
                "bright seed, prescribed undepleted pump; no target-fit parameters"],
            "experimental_validation": False}


def audit_consumed_inputs(ledger, evidence, *, target_dataset_ids=()):
    """An evidence record for a different value or unit cannot approve this run.

    Exact equality is intentional after conversion to the declared SI units;
    numpy's absolute tolerance would falsely equate distinct atomic dipoles.
    Even a pass is a metadata audit, not proof of measurements/model validity.
    """
    evidence = tuple(evidence)
    base = audit_independent_inputs(evidence, required_ids=ledger["scalars"],
                                    target_dataset_ids=target_dataset_ids)
    reasons = list(base.reasons)
    by_id = {p.parameter_id: p for p in evidence}
    for key, (value, unit) in ledger["scalars"].items():
        if key not in by_id:
            continue
        record = by_id[key]
        if record.value != value:
            reasons.append(f"{key}: evidence value differs from consumed value")
        if record.unit != unit:
            reasons.append(f"{key}: evidence unit differs from consumed unit")
    return InputAudit(not reasons, tuple(reasons), "consumed scalar values and declared evidence only")


def power_normalized_readout(inputs, analysis_axis, detector, *, convention="uniform-zeeman-rms"):
    """Recompute rho, M, both D orderings and readout from one fixed dipole map."""
    if not isinstance(inputs, ReducedPowerInputs):
        raise TypeError("ReducedPowerInputs required")
    coupling = reduced_dipoles(convention)
    atom = reduced_pump_noise(coupling.pump_rabi_rad_s(inputs.pump_power_W, inputs.pump_waist_m),
        inputs.one_photon_rad_s, transit_rate_s_inverse=inputs.transit_rate_s_inverse,
        transition_scales=coupling.transition_scales)
    return reduced_seeded_readout(atom, inputs.detunings, analysis_axis, detector,
        number_density_m3=inputs.number_density_m3, uniform_area_m2=inputs.uniform_area_m2,
        optical_omega_rad_s=optical_carriers(inputs.detunings)[1:],
        effective_dipole_C_m=coupling.base_dipole_C_m, transition_scales=coupling.transition_scales,
        length_m=inputs.length_m, seed_power_W=inputs.seed_power_W,
        seed_phase_rad=inputs.seed_phase_rad, phase_mismatch_rad_m=inputs.phase_mismatch_rad_m)
