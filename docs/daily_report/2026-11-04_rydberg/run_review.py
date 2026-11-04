"""Bounded, read-only physics review; outputs confined to this review directory."""
from __future__ import annotations
import hashlib
import dataclasses
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
from unittest.mock import patch
from datetime import datetime, timezone, timedelta

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
os.environ['MPLBACKEND'] = 'Agg'
os.environ['PYTHONIOENCODING'] = 'utf-8'
KST = timezone(timedelta(hours=9))

def stamp():
    return datetime.now(KST).isoformat()

def write(name, data):
    (OUT / name).write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    import numpy as np
    import scipy
    from threadpoolctl import threadpool_info
    from gabes import kernels, schemes
    from gabes.schemes.rydberg import RydbergEITScheme
    before = json.loads((OUT/'initial_sources.json').read_text('utf-8'))
    extra = ['tests/test_rydberg_model_physics.py', 'gabes/pole_doppler.py',
             'gabes/angular.py', 'gabes/report.py', 'gabes/__init__.py',
             'references/rydberg_experiments/reference_policy.json',
             'references/rydberg_experiments/papers/ju_2026_preprint.json',
             'references/rydberg_experiments/validation.json',
             'references/rydberg_experiments/validation_detail.json',
             'references/rydberg_experiments/REVIEW.md']
    before.update({p:digest(ROOT/p) for p in extra})
    snapshot = dict(started_at=stamp(), hashes=before,
        changes_since_initial=[p for p,h in before.items() if digest(ROOT/p)!=h],
        environment=dict(python=sys.version, executable=sys.executable, platform=platform.platform(),
            numpy=np.__version__, scipy=scipy.__version__, packages={p:importlib.metadata.version(p)
                for p in ['pytest','matplotlib','numba','threadpoolctl','streamlit']},
            threadpools=threadpool_info(), env={p:os.environ.get(p) for p in
                ['MPLBACKEND','NUMBA_DISABLE_JIT','NUMBA_CACHE_DIR','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']},
            numba_available=kernels.available()),
        registry=[s.name for s in schemes.all_schemes()])
    probe_only = '--probe-only' in sys.argv
    if not probe_only:
        write('snapshot.json',snapshot)
    else:
        before = json.loads((OUT/'snapshot.json').read_text('utf-8'))['hashes']
    command=[sys.executable,'-m','pytest','-q','tests/test_rydberg_model_physics.py',
        'tests/test_affine_doppler.py','tests/test_rydberg_eit.py','tests/test_rydberg_opacity.py',
        '--durations=8','--junitxml='+str(OUT/'pytest.xml')]
    if not probe_only:
        print('Snapshot saved; tests starting', flush=True)
        start=time.perf_counter()
        with (OUT/'pytest.txt').open('w',encoding='utf-8') as log:
            result=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        checks=dict(command=command,started_at=snapshot['started_at'],completed_at=stamp(),
                    returncode=result.returncode,wall_seconds=time.perf_counter()-start)
        write('checks.json',checks)
        print('Tests completed',result.returncode,flush=True)
    sc=RydbergEITScheme()
    cases=[]
    for view in ['EIT','AT electrometry']:
        for mode in ['off','on']:
            params=dict(sc.recommended_defaults(sc.defaults())[view],doppler=mode)
            raw=sc.compute(params)
            from gabes import rydberg_electrometry as em
            arms=[]
            original=em.balanced_detector_noise
            def detector_call(detector):
                arms.append([detector.signal.optical_power_w,
                    detector.reference.optical_power_w if detector.reference else 0.,
                    detector.reference_weight])
                return original(detector)
            with patch.object(em,'balanced_detector_noise',side_effect=detector_call):
                readout=sc._readout(raw,params)
            def native(value):
                if dataclasses.is_dataclass(value): return native(dataclasses.asdict(value))
                if isinstance(value,np.ndarray): return native(value.tolist())
                if isinstance(value,np.generic): return native(value.item())
                if isinstance(value,complex): return dict(real=value.real,imag=value.imag)
                if isinstance(value,dict): return {k:native(v) for k,v in value.items()}
                if isinstance(value,(list,tuple)): return [native(v) for v in value]
                if isinstance(value,float) and not np.isfinite(value): return None
                return value
            cases.append(dict(view=view,doppler=mode,params=params,
                raw_metadata={k:native(v) for k,v in raw.items() if not isinstance(v,np.ndarray)},
                scan_points=len(raw['scan']),
                transmission_min_max=[float(np.min(readout['T_trans'])),float(np.max(readout['T_trans']))],
                readout=native(readout['rydberg_readout']),
                receiver=dict(evaluations=len(arms),
                    reference_power_min_max_w=[min(a[1] for a in arms),max(a[1] for a in arms)] if arms else None,
                    max_abs_dc_imbalance_w=max(abs(a[0]-a[1]*a[2]) for a in arms) if arms else None),
                superhet_summary=native({k:v for k,v in (readout['superhet'] or {}).items()
                    if not isinstance(v,np.ndarray)})))
    write('probe.json',dict(origin='⟦AI⟧',completed_at=stamp(),cases=cases,
        scope='Current four default view/Doppler readouts; conditional software observations, no experimental validation or benchmark.'))
    after={p:digest(ROOT/p) for p in before}
    write('closure_sources.json',dict(completed_at=stamp(),hashes=after,
        changed=[p for p,h in before.items() if after[p]!=h]))
    print('Probe saved; source changes',[p for p,h in before.items() if after[p]!=h],flush=True)

if __name__=='__main__':
    main()
