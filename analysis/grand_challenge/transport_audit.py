"""Immutable characteristic transport, boundary and cross-position noise audit.

python -m analysis.grand_challenge.transport_audit --output NEW.json --plot NEW.png
"""

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy

from gabes import core
from gabes.quantum.contracts import GeneratorFrequencyAxis
from gabes.quantum.diffusion import stationary_atomic_noise
from gabes.quantum.reservoirs import CollapseChannel, ExplicitReservoirs
from gabes.quantum.transport import BallisticPath, integrate_characteristic, poisson_beam_spectrum
from .transport_fixtures import moving_fixture
from .reference.transport_qrt import characteristic_qrt, mean_pulse
from .reference.ballistic_linear_channel import build_control


ROOT = Path(__file__).resolve().parents[2]
TOLERANCES = {'independent_QRT_relative':1e-9,'ODE_refinement_relative':1e-9,
    'mean_finite_difference_relative':1e-7,'sampled_double_integral_relative':5e-4,
    'constant_analytic_relative':1e-9}


def complex_array(a):
    a=np.asarray(a)
    return {'real':a.real.tolist(),'imag':a.imag.tolist()}


def relative(a,b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),np.finfo(float).tiny))


def hashes():
    files=['gabes/core.py']
    files+=['gabes/quantum/'+s+'.py' for s in ('contracts','diffusion','reservoirs','periodic','transport')]
    files+=['analysis/grand_challenge/'+s+'.py' for s in ('transport_audit','transport_fixtures',
        'reference/transport_qrt','reference/ballistic_linear_channel')]
    files+=['tests/quantum/test_transport.py','tests/quantum/test_ballistic_linear_channel.py']
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files}


def moving_control():
    f=moving_fixture()
    t=f['path'].residence_time_s
    started=time.monotonic()
    axis=GeneratorFrequencyAxis(2*np.pi*50000*np.arange(-10,11),'laboratory Fourier offset for age-stationary toy protocol')
    indices=[4,10,16]
    selected=GeneratorFrequencyAxis(axis.omega_rad_s[indices],axis.frame)
    a=integrate_characteristic(f['path'],f['hamiltonian'],f['reservoirs'],f['boundary_state'])
    packet=a.wavepacket(axis,f['readout'],drive_operators=f['drives'])
    print(f'Moving two-level wavepacket, 21 frequencies: {time.monotonic()-started:.2f}s',flush=True)
    refs=[characteristic_qrt(f['path'],f['hamiltonian'],f['reservoirs'],f['boundary_state'],f['readout'],
        selected.omega_rad_s,drives=f['drives'],rtol=rtol,atol=rtol/100,intervals=intervals)
        for rtol,intervals in ((2e-10,32),(2e-12,64))]
    keys=('greater','lesser','mean_pulse','retarded_response')
    errors={key:relative(packet[key][indices],refs[-1][key]) for key in keys}
    ref_changes={key:relative(refs[0][key],refs[1][key]) for key in keys}
    tightened=integrate_characteristic(f['path'],f['hamiltonian'],f['reservoirs'],f['boundary_state'],rtol=2e-12,atol=2e-14)
    refined=tightened.wavepacket(selected,f['readout'],drive_operators=f['drives'])
    ode_changes={key:relative(packet[key][indices],refined[key]) for key in keys}
    trajectory_errors=relative(a.states,tightened.states)
    fd=[]
    for step in (1e-4,5e-5):
        outputs=[]
        for sign in (-1,1):
            h=lambda age,r:f['hamiltonian'](age,r)+sign*step*f['drives'](age,r)[0]
            outputs.append(mean_pulse(f['path'],h,f['reservoirs'],f['boundary_state'],f['readout']))
        fd.append(relative((outputs[1]-outputs[0])/(2*step),packet['retarded_response'][10,:,0]))
    covariance=a.two_time_covariance(f['readout'])
    full=covariance.transpose(0,2,1,3).reshape(34,34)
    no_inflow=packet['greater']-packet['greater_by_source'][0]
    c_exact=a.covariance_by_source.sum(axis=0)
    missing=a.covariance_by_source[1:].sum(axis=0)
    defect=missing-missing.swapaxes(-1,-2)-c_exact+c_exact.swapaxes(-1,-2)
    boundary_defect=float(np.linalg.norm(defect)/np.linalg.norm(c_exact))
    beam=poisson_beam_spectrum(packet,2e6,source='declared Poisson entries at 2 million atoms/s; identical stationary spatial protocol')
    contributions=[float(np.trace(part[10]).real) for part in packet['greater_by_source']]
    grid_rows=[]
    for count in (17,33,65,129):
        b=a if count==17 else integrate_characteristic(f['path'],f['hamiltonian'],f['reservoirs'],f['boundary_state'],times_s=np.linspace(0,t,count))
        kernel=b.two_time_covariance(f['readout'])
        weight=np.full(count,t/(count-1));weight[[0,-1]]/=2
        full_integral=np.einsum('i,j,ijab->ab',weight,weight,kernel)
        diagonal=np.einsum('i,iiab->ab',weight**2,kernel)
        grid_rows.append({'samples':count,'full_double_integral_relative_error':relative(full_integral,packet['greater'][10]),
            'incorrect_diagonal_only_norm_ratio':float(np.linalg.norm(diagonal)/np.linalg.norm(packet['greater'][10])),
            'atomic_mean_exit_relative_change':relative(b.states[-1],a.states[-1])})
    # An interval partition changes neither the boundary ensemble nor the path.
    # State/propagator composition carries the SAME atom through every slice.
    split=integrate_characteristic(f['path'],f['hamiltonian'],f['reservoirs'],f['boundary_state'],times_s=np.linspace(0,t,33))
    split_packet=split.wavepacket(selected,f['readout'],drive_operators=f['drives'])
    partition={key:relative(packet[key][indices],split_packet[key]) for key in keys}
    causal=a.retarded_kernel(f['readout'],f['drives'])
    return {'path':{'entry_position_m':f['path'].entry_position_m.tolist(),'velocity_m_s':f['path'].velocity_m_s.tolist(),
        'residence_time_s':t,'source':f['path'].source},'parameters':f['parameters'],'boundary_state':complex_array(f['boundary_state']),
        'frequencies_Hz':(axis.omega_rad_s/(2*np.pi)).tolist(),'source_names':a.source_names,
        'mean_trajectory':complex_array(a.states),'sample_ages_s':a.times_s.tolist(),
        'sample_positions_m':a.path.position(a.times_s).tolist(),'sampled_two_time_greater':complex_array(covariance),
        'sampled_retarded_kernel':complex_array(causal),
        'greater':complex_array(packet['greater']),'lesser':complex_array(packet['lesser']),
        'greater_by_source':complex_array(packet['greater_by_source']),'lesser_by_source':complex_array(packet['lesser_by_source']),
        'mean_pulse':complex_array(packet['mean_pulse']),'retarded_response':complex_array(packet['retarded_response']),
        'atomic_audit':a.audit,'wavepacket_audit':packet['audit'],'joint_readout_minimum_eigenvalue':float(np.linalg.eigvalsh(full)[0]),
        'reference_frequencies_Hz':(selected.omega_rad_s/(2*np.pi)).tolist(),'independent_reference_relative_errors':errors,
        'reference_ODE_refinement_relative_changes':ref_changes,'joint_state_wavepacket_ODE_refinement':ode_changes,
        'state_trajectory_ODE_relative_change':trajectory_errors,'interval_partition_relative_changes':partition,
        'mean_response_finite_difference_steps':[1e-4,5e-5],'mean_response_finite_difference_relative_errors':fd,
        'sampled_kernel_quadrature':grid_rows,'missing_inflow_control':{'atomic_commutator_relative_defect':boundary_defect,
            'remaining_noise_minimum_eigenvalue':float(np.linalg.eigvalsh(no_inflow).min()),
            'spectrum_relative_change':relative(no_inflow,packet['greater']),
            'interpretation':'positive remaining noise does not justify discarding the input atomic state fluctuations'},
        'zero_frequency_source_trace_contributions_s2':contributions,
        'poisson_beam':{'arrival_rate_s_inverse':beam['arrival_rate_s_inverse'],'mean_occupancy':beam['mean_occupancy'],
            'greater':complex_array(beam['greater']),'lesser':complex_array(beam['lesser']),
            'number':complex_array(beam['poisson_number']),'internal':complex_array(beam['internal_greater']),
            'scope':beam['scope']},'elapsed_seconds':time.monotonic()-started,
        'passed':bool(a.audit['passed'] and packet['audit']['passed'] and np.linalg.eigvalsh(full)[0]>-1e-10
            and max(errors.values())<TOLERANCES['independent_QRT_relative']
            and max(ref_changes.values())<TOLERANCES['independent_QRT_relative']
            and max(ode_changes.values())<TOLERANCES['ODE_refinement_relative']
            and max(partition.values())<TOLERANCES['ODE_refinement_relative']
            and fd[-1]<TOLERANCES['mean_finite_difference_relative']
            and grid_rows[-1]['full_double_integral_relative_error']<TOLERANCES['sampled_double_integral_relative']
            and boundary_defect>.01 and np.linalg.eigvalsh(no_inflow).min()>0)}


def constant_limit_control():
    lower=np.array([[0.,1.],[0.,0.]])
    h=np.diag([0.,.7]);rho=np.diag([1.,0.])
    res=ExplicitReservoirs(2,[CollapseChannel('emission',lower,'decay','declared unit two-level decay')])
    axis=GeneratorFrequencyAxis([.2,1.1],'constant finite-residence Fourier frequency')
    atom=stationary_atomic_noise(h,res)
    c=np.einsum('ab,iba->i',lower,atom.operators)
    infinite=np.real(c@atom.spectrum(axis).ordered@c.conj())
    s=-.5+1j*(axis.omega_rad_s-.7)
    rows=[]
    for duration in (1.,3.,10.,30.,100.,300.):
        exact=2*((np.expm1(duration*s)-duration*s)/s**2).real
        row={'residence_time_s':duration,'fixed_length':1.,'velocity':1/duration,
            'analytic_single_atom_wavepacket_s2':exact.tolist(),'per_occupancy_spectrum':(exact/duration).tolist(),
            'stationary_limit_relative_error':relative(exact/duration,infinite)}
        if duration in (3.,30.):
            a=integrate_characteristic(BallisticPath([0,0,0],[0,0,1/duration],duration,'constant toy, fixed length'),lambda t,r:h,res,rho)
            p=a.wavepacket(axis,lambda t,r:lower[None])
            row['actual_analytic_relative_error']=relative(p['greater'][:,0,0],exact)
        rows.append(row)
    identity=integrate_characteristic(BallisticPath([0,0,0],[1,0,0],2.,'count-noise toy'),lambda t,r:np.zeros((2,2)),
        ExplicitReservoirs(2,()),np.diag([.6,.4]))
    ia=GeneratorFrequencyAxis([0.,.7],'identity finite residence')
    p=identity.wavepacket(ia,lambda t,r:np.eye(2)[None])
    beam=poisson_beam_spectrum(p,3.,source='declared Poisson identity pulses')
    pulse=np.array([2.,np.expm1(1.4j)/(.7j)])
    count_error=relative(beam['greater'][:,0,0],3*abs(pulse)**2)
    return {'scope':'constant stationary ground-state two-level atom, not full Rb or Maxwell transport limit',
        'decay_rate':1.,'transition_frequency':.7,'frequencies':axis.omega_rad_s.tolist(),
        'stationary_spectrum':infinite.tolist(),'residence_controls':rows,
        'identity_poisson_control':{'mean_pulse':complex_array(p['mean_pulse']),
            'internal_covariance_norm':float(np.linalg.norm(p['greater'])),'PSD':complex_array(beam['greater']),
            'analytic_relative_error':count_error},
        'DC_limit_caution':'nonzero mean-pulse Poisson transit peak may approach a delta distribution as v->0; no pointwise DC stationary claim',
        'passed':bool(max(r.get('actual_analytic_relative_error',0.) for r in rows)<TOLERANCES['constant_analytic_relative']
            and all(b['stationary_limit_relative_error']<a['stationary_limit_relative_error'] for a,b in zip(rows[:-1],rows[1:]))
            and count_error<TOLERANCES['constant_analytic_relative'] and np.linalg.norm(p['greater'])<1e-20)}


def build_report():
    before=hashes()
    with core.blas_single_thread():
        moving=moving_control()
        print('Constant residence/Poisson identity controls',flush=True)
        constant=constant_limit_control()
        print('Independent passive canonical optical transport bridge',flush=True)
        optical=build_control()
    predecessor='docs/grand_challenge/s1_spatial_field_report_v2.json'
    stable=before==hashes()
    return {'schema':'gabes.characteristic_transport_audit.v1','numerical_tolerances':TOLERANCES,
        'scope':'finite atomic characteristic response/noise and independent passive oscillator transport bridge',
        'absolute_hot_vapor_prediction':False,'moving_Rb_Maxwell_closure':False,'experimental_validation':False,
        'moving_atom':moving,'constant_residence':constant,'passive_optical_bridge':optical,
        'source_sha256':before,'source_stable_during_run':stable,
        'predecessor_evidence':{'path':predecessor,'sha256':hashlib.sha256((ROOT/predecessor).read_bytes()).hexdigest(),
            'scope':'immutable historical stationary-center evidence; not a claim of current source parity'},
        'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'BLAS_threads':1},
        'references':['https://arxiv.org/html/2608.15130v1','https://arxiv.org/abs/2301.11993'],
        'expected_controls_passed':bool(stable and moving['passed'] and constant['passed'] and optical['passed']),
        'remaining':['entry-phase and flux-weighted velocity/path ensemble integration',
            'moving periodic Rb with slowly varying optical envelopes',
            'self-consistent nonlocal Maxwell response/noise using the same inflow and internal reservoirs',
            'optical modes, full atom, pump depletion and independent experimental validation']}


def save_plot(report,path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    r=report['moving_atom'];freq=np.array(r['frequencies_Hz'])/1e6
    internal=np.array(r['poisson_beam']['internal']['real'])[:,0,0]
    number=np.array(r['poisson_beam']['number']['real'])[:,0,0]
    sources=np.array(r['greater_by_source']['real'])[:,:,0,0]*2e6
    fig,axes=plt.subplots(1,3,figsize=(13,3.8),layout='constrained')
    axes[0].plot(freq,internal+number,label='Poisson stream total')
    axes[0].plot(freq,internal,'--',label='Internal connected')
    axes[0].plot(freq,number,':',label='Poisson arrival counts')
    axes[0].plot(freq,sources[0],label='Internal inflow contribution',alpha=.7)
    axes[0].set(xlabel='Fourier frequency (MHz)',ylabel='Ordered polarization PSD (toy units s)',title='Moving two-level atom, 4 us residence')
    axes[0].legend(fontsize=7)
    k=r['sampled_two_time_greater'];k=np.array(k['real'])+1j*np.array(k['imag'])
    ages=np.array(r['sample_ages_s'])*1e6
    im=axes[1].imshow(abs(k[:,:,0,0]),origin='lower',extent=[ages[0],ages[-1],ages[0],ages[-1]],aspect='equal')
    axes[1].set(xlabel='Earlier / later age (us)',ylabel='Age (us)',title='Magnitude of cross-position covariance')
    fig.colorbar(im,ax=axes[1],shrink=.8)
    rows=report['constant_residence']['residence_controls']
    axes[2].loglog([x['residence_time_s'] for x in rows],[x['stationary_limit_relative_error'] for x in rows],'o-')
    axes[2].set(xlabel='Residence time (toy time units)',ylabel='Relative error to stationary spectrum',title='Fixed length, v -> 0; constant-atom limit')
    axes[0].grid(alpha=.2);axes[2].grid(alpha=.2)
    fig.suptitle('Boundary + internal reservoirs retain atomic memory; not a hot-Rb squeezing prediction',fontsize=12)
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('xb') as stream:
        fig.savefig(stream,format='png',dpi=160)
    plt.close(fig)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--plot',type=Path)
    args=parser.parse_args(argv)
    for path in (args.output,args.plot):
        if path is not None and path.exists():
            raise FileExistsError(path)
    report=build_report()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8') as stream:
        stream.write(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    if args.plot is not None:
        save_plot(report,args.plot)
    print(f'Wrote {args.output}; expected controls passed={report["expected_controls_passed"]}',flush=True)
    return 0 if report['expected_controls_passed'] else 1


if __name__=='__main__':
    raise SystemExit(main())
