"""Gold-protocol static spectroscopy diagnostics; no production edits."""
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

from gabes.schemes.rydberg import MHZ, RydbergEITScheme
from gabes.schemes.rydberg import _cascade_skeleton
from gabes import atoms

POWERS = [.5, 1., 2., 4., 6., 8., 10.]


class FixedScan(RydbergEITScheme):
    def _scan(self, params):
        return np.linspace(-12., 12., 2401) * MHZ

    def _probe_rabi(self, params):
        return 2. * np.sqrt(float(params["probe_power_uw"]) / 6.)

    def _atom(self, ground_deph, rf_deph):
        sk = _cascade_skeleton()
        return atoms.AtomModel(name="legacy_pair_snapshot", n_levels=4,
            labels=sk["labels"], ground=(0,), excited=(1,2,3),
            decay=sk["decay"], doppler_levels=sk["doppler_levels"],
            doppler_ratios=sk["doppler_ratios"],
            dephasing=((0,2,ground_deph),(2,0,ground_deph),
                       (0,3,ground_deph),(3,0,ground_deph),
                       (2,3,rf_deph),(3,2,rf_deph)))


def main():
    scheme = FixedScan()
    defaults = json.loads((REPO / "analysis/rydberg_gold_map_20261103/resolved_parameters.json").read_text(encoding="utf-8"))["baseline_parameters"]
    rows, spectra = [], {}
    for mode in ("off", "on"):
        for lo in (0., 3.7):
            for power in POWERS:
                params = dict(defaults, doppler=mode, lo_rabi_mhz=lo,
                              view="AT electrometry", probe_power_uw=power,
                              coupling_power_mw=30., temp_c=20., mw_detuning_mhz=0.)
                raw = scheme.compute(params)
                x, transmission, _ = scheme._transmission(raw["chi_bar"], raw, params)
                _, uncomp, _ = scheme._transmission(raw["chi_bar_uncomp"], raw, params)
                key = f"{mode}_lo{lo:g}_pp{power:g}"
                spectra[key + "_x_mhz"] = x
                spectra[key + "_transmission"] = transmission
                spectra[key + "_transmission_uncomp"] = uncomp
                for state, curve in (("compensated", transmission), ("uncompensated", uncomp)):
                    if lo == 0:
                        width, contrast = scheme._eit_features(x, curve)
                    else:
                        width, contrast = float("nan"), float("nan")
                    maxima = scheme._transparency_maxima(x, curve, 6.)
                    row = dict(doppler=mode, lo_rabi_mhz=lo, probe_power_uw=power,
                        compensation=state, transmission_at_zero=float(curve[len(x)//2]),
                        local_floor_5mhz=float(np.min(curve[np.abs(x) <= 5.])),
                        eit_fwhm_mhz=float(width) if np.isfinite(width) else None,
                        eit_contrast=float(contrast) if np.isfinite(contrast) else None,
                        relative_voltage_amplitude=power*float(contrast) if np.isfinite(contrast) else None,
                        transparency_peaks=maxima, probe_rabi_mhz=float(raw["probe_rabi_mhz"]),
                        coupling_rabi_mhz=float(raw["coupling_rabi_mhz"]),
                        transit_dephasing_mhz=float(raw["transit_mhz"]))
                    rows.append(row)
    ROOT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(ROOT / "baseline_spectra.npz", **spectra)
    sources = [REPO / "gabes/schemes/rydberg.py", REPO / "gabes/affine_doppler.py",
               REPO / "gabes/atoms.py", Path(__file__)]
    result = dict(origin="⟦AI⟧ Static spectroscopy baseline, no fitted physical parameters",
        scope="LO-off EIT metrics kept separate from LO=3.7 MHz AT spectra; absolute detector voltage and optical detunings not supplied by this calculation",
        feature_definition="EIT centre minus minimum transmission in +/-5 MHz; linear half-height crossings; relative signal amplitude Pp*contrast",
        rows=rows, source_sha256={str(path.relative_to(REPO)).replace('\\','/'):hashlib.sha256(path.read_bytes()).hexdigest() for path in sources})
    (ROOT / "baseline_spectroscopy.json").write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    for row in rows:
        if row["lo_rabi_mhz"] == 0 and row["compensation"] == "compensated":
            print(json.dumps(row, ensure_ascii=False))


if __name__ == "__main__":
    main()
