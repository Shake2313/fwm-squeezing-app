"""Reproduce historical spectra with the explicit legacy atom snapshot."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
sys.path.insert(0, str(ROOT))
from baseline_spectroscopy import FixedScan


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    historical_path = ROOT / "baseline_spectroscopy.json"
    historical = json.loads(historical_path.read_text(encoding="utf-8"))
    legacy_archive = ROOT/"legacy_rydberg_v7.py"
    archive_matches = digest(legacy_archive) == historical["source_sha256"]["gabes/schemes/rydberg.py"]
    sources = [Path(__file__), ROOT/"baseline_spectroscopy.py", historical_path, legacy_archive,
        ROOT/"baseline_spectra.npz", REPO/"gabes/schemes/rydberg.py",
        REPO/"gabes/atoms.py", REPO/"gabes/core.py", REPO/"gabes/affine_doppler.py",
        REPO/"analysis/rydberg_gold_map_20261103/resolved_parameters.json"]
    before = {str(p.relative_to(REPO)).replace('\\','/'):digest(p) for p in sources}
    params = json.loads((REPO/"analysis/rydberg_gold_map_20261103/resolved_parameters.json").read_text(encoding="utf-8"))["baseline_parameters"]
    old_spectra = np.load(ROOT/"baseline_spectra.npz")
    scheme = FixedScan()
    cases, comparisons = {}, []
    spectrum_error = 0.
    for row in historical["rows"]:
        key = row["doppler"], row["lo_rabi_mhz"], row["probe_power_uw"]
        if key not in cases:
            mode, lo, power = key
            p = dict(params, doppler=mode, lo_rabi_mhz=lo, probe_power_uw=power,
                     view="AT electrometry", coupling_power_mw=30., temp_c=20., mw_detuning_mhz=0.)
            raw = scheme.compute(p)
            x, tc, _ = scheme._transmission(raw["chi_bar"], raw, p)
            _, tu, _ = scheme._transmission(raw["chi_bar_uncomp"], raw, p)
            cases[key] = x, tc, tu
            prefix = f"{mode}_lo{lo:g}_pp{power:g}"
            spectrum_error = max(spectrum_error,
                float(np.max(np.abs(tc-old_spectra[prefix+"_transmission"]))),
                float(np.max(np.abs(tu-old_spectra[prefix+"_transmission_uncomp"]))))
        x, tc, tu = cases[key]
        t = tc if row["compensation"] == "compensated" else tu
        width, contrast = scheme._eit_features(x,t) if row["lo_rabi_mhz"] == 0 else (None,None)
        now = dict(transmission_at_zero=float(t[len(t)//2]),
            local_floor_5mhz=float(np.min(t[np.abs(x)<=5.])),
            eit_fwhm_mhz=width, eit_contrast=contrast,
            relative_voltage_amplitude=row["probe_power_uw"]*contrast if contrast is not None else None)
        differences = {name:float(abs(value-row[name])) for name,value in now.items()
                       if value is not None and row[name] is not None}
        comparisons.append(dict(doppler=row["doppler"], lo_rabi_mhz=row["lo_rabi_mhz"],
            probe_power_uw=row["probe_power_uw"], compensation=row["compensation"],
            absolute_differences=differences))
    after = {str(p.relative_to(REPO)).replace('\\','/'):digest(p) for p in sources}
    changed = [name for name in before if before[name] != after[name]]
    maximum = max(value for row in comparisons for value in row["absolute_differences"].values())
    receipt = dict(origin="⟦AI⟧ Explicit legacy-physics reproduction after production model revision",
        status="PASS" if archive_matches and not changed and maximum<1e-10 and spectrum_error<1e-10 else "FAIL",
        historical_record_sha256=digest(historical_path),
        exact_legacy_engine_archive_sha256=digest(legacy_archive),
        exact_legacy_engine_matches_historical_source_hash=archive_matches,
        historical_source_version_note="Historical baseline JSON records the original execution source hashes, including its then-current helper. That old self-script version is no longer the current helper; historical self-hash equality is not claimed.",
        numerical_scope="Current helper explicitly freezes legacy pairwise atom, old probe anchor and old baseline parameters; production CP atom is not used.",
        checked_unique_atomic_cases=len(cases), checked_metric_rows=len(comparisons),
        maximum_absolute_metric_difference=maximum, maximum_absolute_spectrum_difference=spectrum_error,
        comparisons=comparisons, source_sha256_before=before, source_sha256_after=after,
        changed_sources_during_run=changed)
    (ROOT/"legacy_reproduction.json").write_text(json.dumps(receipt,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({key:receipt[key] for key in ["status","checked_unique_atomic_cases","checked_metric_rows","maximum_absolute_metric_difference","maximum_absolute_spectrum_difference","changed_sources_during_run"]}))


if __name__ == "__main__":
    main()
