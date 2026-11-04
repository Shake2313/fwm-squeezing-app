"""Read-only probe of the batch adapter's newly exposed optical-axis contract."""
import json
from pathlib import Path
import sys
import time
import subprocess
from datetime import datetime, timezone, timedelta
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from analysis.rydberg_cell_heating.adapters import StaticRydbergAdapter, ElectrometryAdapter
from gabes import rydberg_electrometry as em
from gabes.schemes.rydberg import MHZ

def stamp(): return datetime.now(timezone(timedelta(hours=9))).isoformat()
def write(name,data): (OUT/name).write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')

static=StaticRydbergAdapter()
cases=[]
for label,overrides in [
    ('default',{}),
    ('probe_scan_fixed_coupling_12',dict(coupling_detuning_mhz=12.)),
    ('coupling_scan_fixed_probe_12',dict(scan_axis='Coupling detuning',probe_detuning_mhz=12.)),
]:
    point=static.temperature_point(20.,'AT electrometry',overrides,include_spectrum=False)
    cfg=dict(enabled=True,mode='configured',status='ASSUMED',transition_dipole_c_m=2e-27,
             angular_factor=.5,field_amplitude_convention='peak',quantum_efficiency=.8,
             signal_probe_power_scale=1.,probe_detuning_search_points=5)
    configured=ElectrometryAdapter(static,cfg).evaluate(point)
    native=ElectrometryAdapter(static,dict(enabled=True,mode='scheme',status='ASSUMED')).evaluate(point)
    raw,params=point.raw,point.params
    atom=static.scheme._atom(sum(raw[k] for k in ('rydberg_dephasing_mhz','temperature_dephasing_mhz',
        'density_dephasing_mhz','transit_mhz'))*MHZ,raw['rf_dephasing_mhz']*MHZ)
    comparisons=[]
    for record in configured['detuning_scan']:
        coordinate=record['probe_detuning_mhz']
        h=static.scheme._hamiltonian_at_scan_detuning(raw,params,coordinate*MHZ)
        reference=em.weak_signal_response(atom,h,params['if_khz']*1e3*2*np.pi,
            signal_transition=(2,3)).real_observable_phasor_per_angular_rabi(static.scheme._absorption_quadrature_operator())
        actual=complex(record['atomic_phasor_per_angular_rabi_real'],record['atomic_phasor_per_angular_rabi_imag'])
        comparisons.append(dict(scan_coordinate_mhz=coordinate,actual=[actual.real,actual.imag],
            reference=[float(reference.real),float(reference.imag)],
            relative_error=float(abs(actual-reference)/max(abs(reference),1e-30))))
    coordinate=point.readout['superhet']['optimum_detuning_mhz']
    center=static.scheme._resonance_scan_mhz(params)
    cases.append(dict(label=label,params=params,configured_available=configured['available'],
        comparisons=comparisons,
        scheme_mode=dict(batch_optimum_probe_detuning_mhz=native['optimum_probe_detuning_mhz'],
            actual_optimum_probe_detuning_mhz=point.readout['superhet']['optimum_probe_detuning_mhz'],
            actual_optimum_coupling_detuning_mhz=point.readout['superhet']['optimum_coupling_detuning_mhz'],
            optimum_scan_coordinate_mhz=coordinate),
        resonance=dict(actual_scan_center_mhz=center,
            batch_transmission=point.result['metrics']['transmission_at_resonance']['value'],
            native_transmission=point.readout['rydberg_readout']['resonance_transmission'])))
write('integration_probe.json',dict(origin='⟦AI⟧',completed_at=stamp(),cases=cases,
    scope='Default versus 12 MHz fixed-counterpart optical scans; compare atomic phasors before detector scaling. Five search points are a bounded diagnostic, not converged sensitivity.'))
command=[sys.executable,'-m','pytest','-q','analysis/rydberg_cell_heating/test_adapter_integration_regressions.py',
    'analysis/rydberg_cell_heating/test_workflow.py::test_calibrated_finite_if_chain_can_emit_absolute_sensitivity',
    '--durations=5','--junitxml='+str(OUT/'integration_pytest.xml')]
started=stamp(); tick=time.perf_counter()
with (OUT/'integration_pytest.txt').open('w',encoding='utf-8') as log:
    r=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
write('integration_checks.json',dict(command=command,started_at=started,completed_at=stamp(),
    returncode=r.returncode,wall_seconds=time.perf_counter()-tick))
for case in cases:
    print(case['label'], 'worst relative error',max(p['relative_error'] for p in case['comparisons']),
          'scheme metadata',case['scheme_mode'],'resonance',case['resonance'])
print('Existing integration tests exit',r.returncode)
