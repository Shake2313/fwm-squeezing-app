"""Continuous Gaussian Rb: numerical refinements and independent raw QRT.

Each worker integrates the same declared physical equations. Process parallelism
is only an audit scheduling choice. Existing reports are never overwritten.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy

from gabes import core
from gabes.fwm_quantum.kinetic import CarrierGeometry
from gabes.fwm_quantum.smooth_transport import smooth_rb_problem, smooth_rb_wavepacket
from gabes.quantum.contracts import AnalysisFrequencyAxis
from gabes.quantum.transport import BallisticPath
from .normalization_audit import conditional_inputs
from .reference.smooth_qrt import smooth_qrt


ROOT=Path(__file__).resolve().parents[2]
REFINEMENTS=((3e-9,3e-12),(1e-9,1e-12),(3e-10,3e-13))
TOLERANCES={'successive_primary_relative':1e-3,'independent_QRT_relative':5e-6,
    'independent_QRT_refinement_relative':2e-6,'required_successive_primary_comparisons':2}
KEYS=('greater','lesser','mean_pulse','retarded_response','exit_state')
QKEYS=('greater','lesser','mean_pulse','exit_state')


def fixture():
    inputs=replace(conditional_inputs(),transit_rate_s_inverse=0.)
    geometry=CarrierGeometry.vacuum_beams(inputs,probe_angle_rad=.006,conjugate_angle_rad=-.005)
    path=BallisticPath([-1e-4,2e-5,0.],[150.,10.,100.],2e-6,'same declared path as rb_transport_report_v1')
    axis=AnalysisFrequencyAxis.from_hz([.1e6,1e6,4e6])
    return inputs,geometry,path,axis


def hashes():
    files=['gabes/'+x+'.py' for x in ('constants','core','atoms','hyperfine','species','observables','schemes/fwm')]
    files+=['gabes/quantum/'+x+'.py' for x in ('contracts','diffusion','reservoirs','transport','segmented_transport','smooth_transport')]
    files+=['gabes/fwm_quantum/'+x+'.py' for x in ('inputs','model','field','normalization','kinetic','transport','smooth_transport')]
    files+=['analysis/grand_challenge/'+x+'.py' for x in ('normalization_audit','smooth_transport_audit','reference/transport_qrt','reference/segmented_qrt','reference/smooth_qrt')]
    files+=['tests/quantum/test_smooth_transport.py','tests/quantum/test_smooth_qrt.py']
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in files}


def relative(a,b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),np.finfo(float).tiny))


def encode(x):
    if isinstance(x,np.ndarray):
        return {'real':x.real.tolist(),'imag':x.imag.tolist()} if np.iscomplexobj(x) else x.tolist()
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [encode(v) for v in x]
    if isinstance(x,np.generic):return x.item()
    return x


def worker(job):
    label,kind,rtol,atol=job
    before=hashes()
    inputs,geometry,path,axis=fixture()
    started=time.monotonic()
    with core.blas_single_thread():
        if kind=='primary':
            packet=smooth_rb_wavepacket(inputs,geometry,path,axis,rtol=rtol,atol=atol)
            keep=KEYS+('greater_by_source','lesser_by_source','audit','numerics','metadata','source_names')
        else:
            p=smooth_rb_problem(inputs,geometry,path,axis)
            packet=smooth_qrt(**{k:v for k,v in p.items() if k not in ('metadata','drives','drive_frequencies_rad_s')},rtol=rtol,atol=atol)
            keep=QKEYS+('evaluations','complex_dimension','rtol','atol')
    return {'label':label,'kind':kind,'rtol':rtol,'atol':atol,
        'elapsed_seconds':time.monotonic()-started,'source_stable_during_run':before==hashes(),
        'values':{key:packet[key] for key in keep if key in packet}}


def build_report(*,workers=3):
    before=hashes()
    jobs=[(f'primary_{j}','primary',r,a) for j,(r,a) in enumerate(REFINEMENTS)]
    jobs += [('QRT_coarse','reference',2e-10,2e-14),('QRT_fine','reference',2e-11,2e-15)]
    results={}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures={pool.submit(worker,job):job[0] for job in jobs}
        for future in as_completed(futures):
            row=future.result()
            results[row['label']]=row
            print(f"{row['label']} finished: {row['elapsed_seconds']:.2f}s",flush=True)
    primary=[results[f'primary_{j}'] for j in range(len(REFINEMENTS))]
    reference=results['QRT_fine']['values']
    refinements=[{key:relative(new['values'][key],old['values'][key]) for key in KEYS}
                 for old,new in zip(primary[:-1],primary[1:])]
    independent={key:relative(primary[-1]['values'][key],reference[key]) for key in QKEYS}
    qrefinement={key:relative(reference[key],results['QRT_coarse']['values'][key]) for key in QKEYS}
    old_path=ROOT/'docs/grand_challenge/rb_transport_report_v1.json'
    old=json.loads(old_path.read_text(encoding='utf-8'))
    frozen=[]
    for row in old['Gaussian_refinements']:
        values={key:np.array(row[key]['real'])+1j*np.array(row[key]['imag']) for key in KEYS}
        frozen.append({'segments':row['segments'],
            'relative_error_vs_continuous':{key:relative(values[key],primary[-1]['values'][key]) for key in KEYS}})
    converged=(len(refinements)==2 and all(max(r.values())<TOLERANCES['successive_primary_relative'] for r in refinements)
        and max(independent.values())<TOLERANCES['independent_QRT_relative']
        and max(qrefinement.values())<TOLERANCES['independent_QRT_refinement_relative'])
    stable=before==hashes() and all(row['source_stable_during_run'] for row in results.values())
    quantum=all(row['values']['audit']['passed'] for row in primary)
    inputs,geometry,path,axis=fixture()
    return {'schema':'gabes-continuous-rb-transport-v1','source_sha256':before,'source_stable_during_run':stable,
        'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
        'inputs_SI':asdict(inputs),'path':{'entry_position_m':path.entry_position_m.tolist(),
            'velocity_m_s':path.velocity_m_s.tolist(),'residence_time_s':path.residence_time_s},
        'wavevectors_rad_m':geometry.wavevectors_rad_m.tolist(),'analysis_frequencies_hz':axis.frequency_hz.tolist(),
        'tolerances':TOLERANCES,'primary_refinements':refinements,'independent_QRT_errors':independent,
        'independent_QRT_refinement':qrefinement,'cases':encode(results),
        'historical_frozen_report_sha256':hashlib.sha256(old_path.read_bytes()).hexdigest(),
        'frozen_envelope_error_vs_continuous':frozen,'quantum_controls_passed':quantum,
        'continuous_Gaussian_selected_path_converged':bool(converged),
        'all_declared_controls_passed':bool(converged and quantum and stable),
        'scope':{'same_pump_only_reduced_Rb_physics':True,'continuous_Gaussian_single_path':True,
            'finite_seed_saturation':False,'thermal_Rb_ensemble_converged':False,'nonlocal_Maxwell_field':False,
            'full_Zeeman_atom':False,'experimental_validation':False,'absolute_hot_vapor_squeezing':False},
        'limitations':['Refinement applies to this declared path and three RF samples, not all inputs or frequencies.',
            'ODE tolerances and raw-QRT agreement are observed numerical controls, not a rigorous global error bound.',
            'The continuous method resolves GHz oscillations and is expensive; it is not an Ultra scan replacement.',
            'Preserve source noise and validated observables in any subsequent acceleration.']}


def write_audit(output,plot,*,workers=3):
    output,plot=Path(output),Path(plot)
    if output.resolve()==plot.resolve() or output.exists() or plot.exists():
        raise FileExistsError('new distinct immutable report and plot paths required')
    if isinstance(workers,bool) or int(workers)!=workers or not 1<=workers<=8:
        raise ValueError('workers must be an integer from one to eight')
    report=build_report(workers=int(workers))
    output.parent.mkdir(parents=True,exist_ok=True); plot.parent.mkdir(parents=True,exist_ok=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.4),layout='constrained')
    for key in KEYS:
        axes[0].loglog([r['segments'] for r in report['frozen_envelope_error_vs_continuous']],
            [r['relative_error_vs_continuous'][key] for r in report['frozen_envelope_error_vs_continuous']],'o-',label=key)
        axes[1].semilogy([1,2],[r[key] for r in report['primary_refinements']],'o-',label=key)
    axes[0].set(xlabel='Frozen Gaussian subdivisions',ylabel='Relative difference from continuous solve',title='Same physics: envelope approximation error')
    axes[1].set(xlabel='Successive tolerance refinement',ylabel='Relative change',title='Continuous microscopic D propagation')
    axes[1].axhline(TOLERANCES['successive_primary_relative'],color='black',ls=':')
    axes[0].legend(fontsize=8); axes[1].legend(fontsize=8)
    fig.suptitle('Selected moving Rb path: '+('CONVERGED' if report['all_declared_controls_passed'] else 'UNCONVERGED')+'; no optical squeezing claim')
    fig.savefig(plot,dpi=160); plt.close(fig)
    with output.open('x',encoding='utf-8') as f:
        json.dump(report,f,ensure_ascii=False,indent=2,allow_nan=False); f.write('\n')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True); parser.add_argument('--plot',required=True)
    parser.add_argument('--workers',type=int,default=3)
    args=parser.parse_args()
    r=write_audit(args.output,args.plot,workers=args.workers)
    print(json.dumps({k:r[k] for k in ('quantum_controls_passed','continuous_Gaussian_selected_path_converged','all_declared_controls_passed')},indent=2))
    if not r['all_declared_controls_passed']:raise SystemExit(1)


if __name__=='__main__':main()
