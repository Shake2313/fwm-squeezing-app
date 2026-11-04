"""Bounded spectroscopy parameter profiles, holding reported Omega_c fixed."""
from __future__ import annotations

import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
sys.path.insert(0, str(REPO))

from gabes import atoms
from gabes.schemes.rydberg import MHZ, RydbergEITScheme, _cascade_skeleton

POWERS = [.5, 2., 6., 10.]
GROUND_RATES = [.02, .04, .06, .08, .10, .15, .20, .25]
PROBE_ANCHORS = [1., 1.5, 2., 2.5, 3., 4.]
OPACITIES = [.001, .01, .1]


class Candidate(RydbergEITScheme):
    def __init__(self, variant, anchor, total_ground_rate):
        self.variant = variant
        self.anchor = anchor
        self.ground_rate = total_ground_rate
        super().__init__()

    def _scan(self, params):
        return np.linspace(-10., 10., 601) * MHZ

    def _probe_rabi(self, params):
        return self.anchor * np.sqrt(float(params["probe_power_uw"]) / 6.)

    def _transit_rate_mhz(self, params):
        return 0.

    def _atom(self, ground_deph, rf_deph):
        if self.variant == "legacy_pair":
            sk = _cascade_skeleton()
            return atoms.AtomModel(name="legacy_pair_snapshot", n_levels=4,
                labels=sk["labels"], ground=(0,), excited=(1,2,3),
                decay=sk["decay"], doppler_levels=sk["doppler_levels"],
                doppler_ratios=sk["doppler_ratios"],
                dephasing=((0,2,ground_deph),(2,0,ground_deph),
                           (0,3,ground_deph),(3,0,ground_deph),
                           (2,3,rf_deph),(3,2,rf_deph)))
        sk = _cascade_skeleton()
        reset = ground_deph if self.variant == "nested_cp_reset" else 0.
        # Reset-only ground coherence is one bounded candidate, not an
        # inference that intrinsic laser broadening is exactly zero.
        phase = max(ground_deph - reset, 0.)
        common = np.sqrt(2 * phase) * np.diag([0., 0., 1., 1.])
        relative = np.sqrt(2 * rf_deph) * np.diag([0., 0., 0., 1.])
        atom = atoms.AtomModel(name=self.variant, n_levels=4, labels=sk["labels"],
            ground=(0,), excited=(1,2,3), decay=sk["decay"], dephasing=(),
            doppler_levels=sk["doppler_levels"], doppler_ratios=sk["doppler_ratios"],
            collapse_ops=(common, relative))
        return atom.with_reset(reset, (1.,)) if reset else atom


def main():
    sources = [Path(__file__), REPO/"gabes/schemes/rydberg.py", REPO/"gabes/atoms.py"]
    before = {str(path.relative_to(REPO)).replace('\\','/'):hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}
    rows = []
    scheme = RydbergEITScheme()
    defaults = scheme.recommended_defaults(scheme.defaults())["EIT"]
    for variant in ("legacy_pair", "nested_cp", "nested_cp_reset"):
        for rate in GROUND_RATES:
            for anchor in PROBE_ANCHORS:
                candidate = Candidate(variant, anchor, rate)
                for power in POWERS:
                    params = dict(defaults, doppler="on", probe_power_uw=power,
                        rydberg_dephasing_mhz=rate, coupling_rabi_mhz=3.,
                        coupling_power_mw=30., temp_c=20., lo_rabi_mhz=0.)
                    raw = candidate.compute(params)
                    for opacity in OPACITIES:
                        scenario = dict(params, effective_opacity_participation=opacity)
                        x, curve, _ = candidate._transmission(raw["chi_bar"], raw, scenario)
                        width, contrast = candidate._eit_features(x, curve)
                        rows.append(dict(variant=variant, total_ground_coherence_rate_mhz=rate,
                            probe_rabi_anchor_mhz=anchor, probe_power_uw=power,
                            opacity_participation=opacity,
                            eit_fwhm_mhz=float(width), eit_contrast=float(contrast),
                            relative_voltage_amplitude=power*float(contrast),
                            transmission_at_zero=float(curve[len(curve)//2])))
        print(f"Completed {variant}", flush=True)
    after = {str(path.relative_to(REPO)).replace('\\','/'):hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}
    result = dict(origin="⟦AI⟧ Bounded static spectroscopy profiles, not a fitted final model",
        fixed=dict(coupling_rabi_mhz=3., coupling_power_mw=30., temperature_c=20.,
                   lo_rabi_mhz=0., wavelength_geometry="existing 780/481 counterpropagating"),
        variants=dict(legacy_pair="Existing pairwise non-GKSL phenomenology",
            nested_cp="Common optical phase on states2/3 plus relative RF phase on state3",
            nested_cp_reset="All fitted ground-coherence rate allocated to transit reset; common optical phase set0, relative RF phase retained. Reset changes populations and adds RF coherence loss; no doubled transit term."),
        scope="Parameter grid is a spectroscopy identification diagnostic. Gamma intrinsic/transit split is not identified at one temperature. Absolute photodetector voltage gain is not known.",
        rows=rows, source_sha256=before, source_sha256_after=after,
        changed_sources_during_run=[name for name in before if before[name]!=after[name]])
    ROOT.mkdir(parents=True, exist_ok=True)
    (ROOT/"linewidth_profiles.json").write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print("Rows:",len(rows))


if __name__ == "__main__":
    main()
