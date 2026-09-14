"""Immutable moving reduced-Rb finite-transport audit, including failed gates."""

import argparse
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
from gabes.fwm_quantum.transport import reduced_ballistic_problem, reduced_ballistic_wavepacket, average_pump_only_entry_phase
from gabes.quantum.contracts import AnalysisFrequencyAxis
from gabes.quantum.segmented_transport import segmented_wavepacket
from gabes.quantum.transport import BallisticPath
from .normalization_audit import conditional_inputs
from .reference.segmented_qrt import segmented_qrt


ROOT = Path(__file__).resolve().parents[2]
TOLERANCES = {'independent_QRT_relative': 2e-9, 'phase_theorem_relative': 2e-9,
    'constant_segment_refinement_relative': 2e-9, 'smooth_envelope_relative': 1e-3,
    'required_consecutive_smooth_refinements': 2}
KEYS = ('greater', 'lesser', 'mean_pulse', 'retarded_response', 'exit_state')


def relative(a, b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b), np.finfo(float).tiny))


def encode(a):
    a = np.asarray(a)
    return {'real': a.real.tolist(), 'imag': a.imag.tolist()}


def source_hashes():
    files = ['gabes/'+s+'.py' for s in ('constants','core','atoms','hyperfine','species','observables','schemes/fwm')]
    files += ['gabes/quantum/'+s+'.py' for s in ('contracts','diffusion','reservoirs','transport','segmented_transport')]
    files += ['gabes/fwm_quantum/'+s+'.py' for s in ('inputs','model','field','normalization','kinetic','transport')]
    files += ['analysis/grand_challenge/'+s+'.py' for s in ('normalization_audit','rb_transport_audit','reference/segmented_qrt','reference/transport_qrt')]
    files += ['tests/quantum/test_segmented_transport.py','tests/quantum/test_segmented_qrt.py']
    return {f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files}


def pack(packet):
    return {**{key: encode(packet[key]) for key in KEYS},
        'greater_by_source': encode(packet['greater_by_source']), 'lesser_by_source': encode(packet['lesser_by_source']),
        'source_names': packet['source_names'], 'audit': packet['audit'], 'metadata': packet['metadata'],
        'matrix_exponentials': packet['matrix_exponentials'], 'identical_exponentials_reused': packet['identical_exponentials_reused']}


def qrt(problem):
    return segmented_qrt(**{key: value for key, value in problem.items()
                           if key not in ('metadata','drives','drive_frequencies_rad_s')})


def build_report():
    before = source_hashes()
    inputs = replace(conditional_inputs(), transit_rate_s_inverse=0.)
    geometry = CarrierGeometry.vacuum_beams(inputs, probe_angle_rad=.006, conjugate_angle_rad=-.005)
    path = BallisticPath([-1e-4,2e-5,0.], [150.,10.,100.], 2e-6, 'declared finite Rb path; not sampled thermal ensemble')
    axis = AnalysisFrequencyAxis.from_hz([.1e6,1e6,4e6])
    constant_path = replace(path, velocity_m_s=[0.,0.,100.])
    constant = []
    for segments in (1,4):
        p = reduced_ballistic_wavepacket(inputs,geometry,constant_path,axis,segments=segments)
        constant.append(p)
    cref = qrt(reduced_ballistic_problem(inputs,geometry,constant_path,axis,segments=1))
    constant_errors = {key:relative(constant[0][key], cref[key]) for key in ('greater','lesser','mean_pulse','exit_state')}
    subdivision_errors = {key:relative(constant[1][key],constant[0][key]) for key in KEYS}
    print('Constant physical Rb and independent QRT complete',flush=True)
    rows, packets, times, refined = [], [], [], []
    for segments in (8,16,32,64,128):
        started = time.monotonic()
        p = reduced_ballistic_wavepacket(inputs,geometry,path,axis,segments=segments)
        elapsed = time.monotonic()-started
        row = {'segments':segments,'elapsed_seconds':elapsed, **pack(p)}
        if packets:
            errors = {key:relative(p[key], packets[-1][key]) for key in KEYS}
            row['relative_change_from_previous'] = errors
            row['smooth_refinement_passed'] = all(x<TOLERANCES['smooth_envelope_relative'] for x in errors.values())
            refined.append(row['smooth_refinement_passed'])
        rows.append(row); packets.append(p); times.append(elapsed)
        print(f'Gaussian segments={segments}: {elapsed:.2f}s; '+str(row.get('relative_change_from_previous',{})),flush=True)
    problem = reduced_ballistic_problem(inputs,geometry,path,axis,segments=8)
    pref = qrt(problem)
    reference_errors = {key:relative(packets[0][key],pref[key]) for key in ('greater','lesser','mean_pulse','exit_state')}
    phi = .73
    pphi = reduced_ballistic_wavepacket(inputs,geometry,path,axis,segments=8,entry_phase_rad=phi)
    p = packets[0]
    charges = np.array([1,-1,-1,1])
    u = np.exp(1j*charges*phi)
    phase_errors = {key:relative(pphi[key],p[key]*u[:,None]*u.conj()[None,:]) for key in
                    ('greater','lesser','greater_by_source','lesser_by_source','retarded_response')}
    phase_errors['mean_pulse'] = relative(pphi['mean_pulse'],p['mean_pulse']*u)
    avg = average_pump_only_entry_phase(p)
    phase_row = {key:encode(value) for key,value in avg.items() if isinstance(value,np.ndarray)}
    # A single thermal beat phase preserves correlations (probe,conjugate†).
    wrong = avg['conditional_greater']*np.eye(4)
    phase_row['wrong_independent_port_phases_relative_change'] = relative(wrong,avg['conditional_greater'])
    phase_row['deleting_mean_outer_relative_change'] = relative(avg['conditional_greater'],avg['raw_greater_phase_average'])
    full = p['atomic_covariance_by_source'][-1,0].sum(axis=0)
    internal = p['atomic_covariance_by_source'][-1,0,1:].sum(axis=0)
    boundary = {'exit_atomic_commutator_relative_defect':relative(internal-internal.T,full-full.T),
        'wavepacket_greater_relative_change':relative(p['greater_by_source'][1:].sum(axis=0),p['greater']),
        'remaining_noise_minimum_eigenvalue':float(np.linalg.eigvalsh(p['greater_by_source'][1:].sum(axis=0)).min())}
    # Reinitialize fluctuations to each actual local density, retaining means
    # but erasing covariance between distinct segments. This is a negative
    # control, not an admissible transport approximation.
    diagonal = np.zeros_like(p['greater'])
    age = 0.
    for j, dt in enumerate(problem['durations_s']):
        local = segmented_wavepacket(problem['hamiltonians'][j:j+1],problem['reservoirs'],p['states'][j],
                    [dt],problem['readouts'][j:j+1],problem['frequencies_rad_s'])
        phase = np.exp(1j*problem['frequencies_rad_s']*age)
        diagonal += local['greater']*phase[:,:,None]*phase.conj()[:,None,:]
        age += dt
    negative = {'discard_cross_segment_covariance_relative_change':relative(diagonal,p['greater']),
        'remaining_noise_minimum_eigenvalue':float(np.linalg.eigvalsh(diagonal).min())}
    smooth_passed = len(refined)>=2 and all(refined[-2:])
    exact_passed = (max(constant_errors.values())<TOLERANCES['independent_QRT_relative']
        and max(reference_errors.values())<TOLERANCES['independent_QRT_relative']
        and max(subdivision_errors.values())<TOLERANCES['constant_segment_refinement_relative']
        and max(phase_errors.values())<TOLERANCES['phase_theorem_relative']
        and all(p['audit']['passed'] for p in constant+packets)
        and boundary['wavepacket_greater_relative_change']>1e-3
        and negative['discard_cross_segment_covariance_relative_change']>1e-3)
    after = source_hashes()
    return {'schema':'gabes-rb-finite-transport-v1','source_sha256':before,'source_stable_during_run':before==after,
        'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
        'inputs_SI':asdict(inputs),'path':{'entry_position_m':path.entry_position_m.tolist(),
            'velocity_m_s':path.velocity_m_s.tolist(),'residence_time_s':path.residence_time_s},
        'wavevectors_rad_m':geometry.wavevectors_rad_m.tolist(),'analysis_frequencies_hz':axis.frequency_hz.tolist(),
        'tolerances':TOLERANCES,'constant_reference_errors':constant_errors,'constant_subdivision_errors':subdivision_errors,
        'Gaussian_reference_errors_at_8_segments':reference_errors,'Gaussian_refinements':rows,
        'entry_phase_theorem_errors':phase_errors,'entry_phase_average_at_8_segments':phase_row,
        'inflow_omission':boundary,'cross_segment_omission':negative,
        'exact_segment_model_controls_passed':bool(exact_passed),
        'smooth_Gaussian_envelope_converged':bool(smooth_passed),
        'all_certification_gates_passed':bool(exact_passed and smooth_passed and before==after),
        'scope':{'moving_pump_only_RMS_Rb':True,'exact_piecewise_constant_noise_response':True,
            'finite_seed_saturation':False,'thermal_Rb_ensemble_evaluated':False,'nonlocal_Maxwell_closure':False,
            'absolute_experimental_squeezing':False},
        'status':'CONVERGED conditional characteristic' if smooth_passed else 'UNCONVERGED smooth Gaussian envelope; exact frozen-segment controls are separate'}


def write_audit(output, plot):
    output, plot = Path(output), Path(plot)
    if output.resolve()==plot.resolve() or output.exists() or plot.exists():
        raise FileExistsError('new distinct immutable output and plot paths required')
    with core.blas_single_thread():
        report = build_report()
    output.parent.mkdir(parents=True,exist_ok=True)
    plot.parent.mkdir(parents=True,exist_ok=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    rows = report['Gaussian_refinements']
    for key in KEYS:
        ax[0].loglog([r['segments'] for r in rows[1:]], [r['relative_change_from_previous'][key] for r in rows[1:]], 'o-',label=key)
    ax[0].axhline(TOLERANCES['smooth_envelope_relative'],color='black',ls=':',label='required tolerance')
    ax[0].set(xlabel='Gaussian envelope segments',ylabel='Relative change on doubling',title='Envelope refinement: all quantities must pass')
    ax[0].legend(fontsize=8)
    selected = rows[-1]
    greater = np.array(selected['greater']['real'])+1j*np.array(selected['greater']['imag'])
    for j,label in enumerate(('probe','conjugate','probe dagger','conjugate dagger')):
        ax[1].plot(report['analysis_frequencies_hz'],greater[:,j,j].real,'o-',label=label)
    ax[1].set(xlabel='Lab RF frequency (Hz)',ylabel='Single-atom wavepacket covariance (s^2)',title='128 segments: not a squeezing spectrum')
    ax[1].legend(fontsize=8)
    fig.suptitle('Moving reduced Rb, finite residence; '+report['status'],fontsize=11)
    fig.savefig(plot,dpi=160)
    plt.close(fig)
    with output.open('x',encoding='utf-8') as f:
        json.dump(report,f,ensure_ascii=False,indent=2,allow_nan=False)
        f.write('\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    parser.add_argument('--plot',required=True)
    args = parser.parse_args()
    report = write_audit(args.output,args.plot)
    print(json.dumps({k:report[k] for k in ('exact_segment_model_controls_passed','smooth_Gaussian_envelope_converged','source_stable_during_run','status')},indent=2))
    if not report['exact_segment_model_controls_passed'] or not report['source_stable_during_run']:
        raise SystemExit(1)


if __name__=='__main__':
    main()
