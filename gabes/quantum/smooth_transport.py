"""Continuous affine-Hamiltonian characteristics with independent jump noise.

Sparse linear moment lifts are assembled once. Adaptive time integration
resolves internal and carrier dynamics; no staircase envelope or RWA is used.
"""

import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix

from .contracts import readonly_array
from .diffusion import _state, traceless_hermitian_basis
from .reservoirs import ExplicitReservoirs
from .segmented_transport import _ordered_audit
from .transport import _moments


def smooth_wavepacket(h0, h1, envelope, reservoirs, boundary_state, duration_s,
                       readouts, frequencies_rad_s, *, drives=None,
                       drive_frequencies_rad_s=None, rtol=1e-10, atol=1e-13,
                       max_step_s=None, sample_count=17):
    """H(t)=h0+envelope(t)*h1; fixed complex readout and drive columns.

    Y_j=int exp(i*w_j*t) O_j dt. Shapes: H(n,n), O(p,n,n), w(nf,p).
    The actual jump operators are constant. Finite atomic moments are neither
    canonical field ports nor a self-consistent Maxwell channel. Refinement
    of the returned values is required separately from quantum consistency.
    """
    if not isinstance(reservoirs, ExplicitReservoirs) or not callable(envelope):
        raise TypeError('explicit reservoirs and real envelope callback required')
    n = reservoirs.n_levels
    hs = readonly_array([h0,h1])
    if hs.shape != (2,n,n) or np.linalg.norm(hs-hs.conj().swapaxes(-1,-2)) > 1e-12*max(np.linalg.norm(hs),1e-300):
        raise ValueError('two finite Hermitian Hamiltonians required')
    rho0 = _state(boundary_state,n)
    duration = float(duration_s)
    ops, freq = readonly_array(readouts), readonly_array(frequencies_rad_s,real=True)
    if (not np.isfinite(duration) or duration<=0 or ops.ndim!=3 or ops.shape[1:]!=(n,n)
            or not len(ops) or freq.ndim!=2 or freq.shape[1]!=len(ops) or not len(freq)):
        raise ValueError('positive duration, fixed readouts and matching port frequencies required')
    if (not np.isfinite([rtol,atol]).all() or not 0<rtol<1 or not 0<atol<1
            or isinstance(sample_count,bool) or int(sample_count)!=sample_count or sample_count<2):
        raise ValueError('positive tolerances below one and integer sample_count >=2 required')
    maximum = duration/64 if max_step_s is None else float(max_step_s)
    if not np.isfinite(maximum) or maximum<=0:
        raise ValueError('positive finite max_step_s required')
    nf, ports = freq.shape
    if drives is None:
        if drive_frequencies_rad_s is not None:
            raise ValueError('drive frequencies require drive columns')
        v,wf = np.zeros((0,n,n),complex), np.zeros((nf,0))
    else:
        v,wf = readonly_array(drives), readonly_array(drive_frequencies_rad_s,real=True)
        if v.ndim!=3 or v.shape[1:]!=(n,n) or not len(v) or wf.shape!=(nf,len(v)):
            raise ValueError('matched drive columns and frequencies required')
    inputs = len(v)
    f = traceless_hermitian_basis(n)
    m,nr,ns = len(f),n*n,1+len(reservoirs.channels)
    basis = f.reshape(m,nr).T
    l0 = reservoirs.generator(hs[0])
    # L(H0+fH1)=L(H0)+f*(L(H1)-L(0)): include dissipation exactly ONCE.
    l1 = reservoirs.generator(hs[1])-reservoirs.dissipator()
    aa = np.array([basis.conj().T@l@basis for l in (l0,l1)])
    if np.linalg.norm(aa.imag)>1e-11*max(np.linalg.norm([l0,l1]),1e-300):
        raise ValueError('Hermitian atomic drift must be real')
    a0,a1 = aa.real
    c = np.einsum('jab,iba->ji',ops,f)
    atomic_end=nr+ns*m*m
    cross_size=m*ports
    window_size=2*cross_size+ports*ports
    window_shape=(2,ns,nf,window_size)
    cov_end=atomic_end+int(np.prod(window_shape))
    mean_end = cov_end+nf*ports
    response_dim = m+ports
    size = mean_end+nf*inputs*response_dim
    rows=[[],[]]; cols=[[],[]]; vals=[[],[]]
    def add(which,block,row,col):
        ii,jj=np.nonzero(block)
        rows[which].append(ii+row); cols[which].append(jj+col); vals[which].append(block[ii,jj])
    add(0,l0*duration,0,0); add(1,l1*duration,0,0)
    products=[]
    for channel in reservoirs.channels:
        comm=f@channel.operator-channel.operator@f
        products.append(np.einsum('iab,jbc->ijac',comm.conj().swapaxes(-1,-2),comm))
    products=np.array(products).reshape(ns-1,m,m,n,n)
    dmap=products.swapaxes(-1,-2).reshape(ns-1,m*m,nr)
    # Exact shared block: equal-time atomic C_r is independent of the RF
    # readout. For Hermitian F, the other ordering is precisely C_r.T.
    # Evolve it ONCE per physical source; do not restore 2*nf identical solves.
    for source in range(ns):
        start=nr+source*m*m
        for which,a in enumerate((a0,a1)):
            add(which,(np.kron(a,np.eye(m))+np.kron(np.eye(m),a))*duration,start,start)
        if source:
            add(0,dmap[source-1]*duration,start,0)
    transpose_columns=np.arange(m*m).reshape(m,m).T.ravel()
    # EXACT change of accumulator units Z_j -> eta_j Z_j, eta=1+|w_j|T.
    # This resolves small GHz-demodulated pulses under the ODE error norm;
    # undo eta and T below. It is not an optical efficiency or noise factor.
    eta=1+abs(freq)*duration
    for fi,w in enumerate(freq):
        read=c*eta[fi,:,None]
        iw=1j*np.diag(w*duration)
        x0=np.kron(a0*duration,np.eye(ports))+np.kron(np.eye(m),iw)
        x1=np.kron(a1*duration,np.eye(ports))
        y0=np.kron(np.eye(ports),a0*duration)-np.kron(iw,np.eye(m))
        y1=np.kron(np.eye(ports),a1*duration)
        z0=-np.kron(iw,np.eye(ports))+np.kron(np.eye(ports),iw)
        fx=np.kron(np.eye(m),read.conj())
        fy=np.kron(read,np.eye(m))
        for order in range(2):
            for source in range(ns):
                start=atomic_end+((order*ns+source)*nf+fi)*window_size
                add(0,x0,start,start); add(1,x1,start,start)
                add(0,y0,start+cross_size,start+cross_size); add(1,y1,start+cross_size,start+cross_size)
                add(0,z0,start+2*cross_size,start+2*cross_size)
                add(0,fx if order==0 else fx[:,transpose_columns],start,nr+source*m*m)
                add(0,fy if order==0 else fy[:,transpose_columns],start+cross_size,nr+source*m*m)
                add(0,np.kron(read,np.eye(ports)),start+2*cross_size,start)
                add(0,np.kron(np.eye(ports),read.conj()),start+2*cross_size,start+cross_size)
        start=cov_end+fi*ports
        add(0,ops.swapaxes(-1,-2).reshape(ports,nr)*eta[fi,:,None],start,0)
        add(0,-1j*np.diag(w*duration),start,start)
        for j in range(inputs):
            start=mean_end+(fi*inputs+j)*response_dim
            b=basis.conj().T@(-1j*(np.kron(v[j],np.eye(n))-np.kron(np.eye(n),v[j].T)))
            add(0,b*duration,start,0)
            add(0,(a0+1j*wf[fi,j]*np.eye(m))*duration,start,start)
            add(1,a1*duration,start,start)
            add(0,c*eta[fi,:,None],start+m,start)
            add(0,-1j*np.diag((w-wf[fi,j])*duration),start+m,start+m)
    matrices=[coo_matrix((np.concatenate(vals[q]),(np.concatenate(rows[q]),np.concatenate(cols[q]))),shape=(size,size)).tocsr() for q in range(2)]
    initial=np.zeros(size,complex); initial[:nr]=rho0.ravel()
    _,cin=_moments(rho0,f)
    initial[nr:nr+m*m]=cin.ravel()
    calls=0
    def rhs(u,y):
        nonlocal calls
        value=np.asarray(envelope(float(u)*duration))
        if value.ndim!=0 or np.iscomplexobj(value) or not np.isfinite(value):
            raise ValueError('envelope must return one finite real value')
        calls+=1
        return matrices[0]@y+float(value)*(matrices[1]@y)
    started=time.monotonic()
    times=np.linspace(0.,1.,int(sample_count))
    ode=solve_ivp(rhs,(0.,1.),initial,method='DOP853',rtol=rtol,atol=atol,
        max_step=maximum/duration,t_eval=times)
    if not ode.success:
        raise ValueError('continuous microscopic moment integration failed: '+ode.message)
    final=ode.y[:,-1]
    pc=final[atomic_end:cov_end].reshape(window_shape)[...,-ports*ports:].reshape(2,ns,nf,ports,ports)
    phase=np.exp(1j*freq*duration)/eta
    ordered=pc*phase[None,None,:,:,None]*phase.conj()[None,None,:,None,:]*duration**2
    mean=final[cov_end:mean_end].reshape(nf,ports)*phase*duration
    response=np.empty((nf,ports,inputs),complex)
    for fi in range(nf):
        for j in range(inputs):
            start=mean_end+(fi*inputs+j)*response_dim+m
            response[fi,:,j]=final[start:start+ports]*np.exp(1j*(freq[fi]-wf[fi,j])*duration)/eta[fi]*duration
    states=ode.y[:nr].T.reshape(-1,n,n)
    for state in states: _state(state,n)
    samples=ode.y[nr:atomic_end].T.reshape(len(times),ns,m,m)
    atomic=np.stack([samples,samples.swapaxes(-1,-2)],axis=1)
    exact=np.array([_moments(rho,f)[1] for rho in states])
    actual=atomic[:,0].sum(axis=1)
    scale=max(np.linalg.norm(exact),1e-300)
    floor=64*np.finfo(float).eps*duration**2*max(float(np.linalg.norm(ops,axis=(-2,-1)).max())**2,1e-300)
    audit={'atomic_covariance_relative_residual':float(np.linalg.norm(actual-exact)/scale),
        'atomic_commutator_relative_residual':float(np.linalg.norm(actual-actual.swapaxes(-1,-2)-exact+exact.swapaxes(-1,-2))/scale),
        'minimum_state_eigenvalue':float(np.linalg.eigvalsh(states).min()),
        'greater':_ordered_audit(ordered[0],floor),'lesser':_ordered_audit(ordered[1],floor)}
    audit['passed']=bool(audit['atomic_covariance_relative_residual']<2e-7 and
        audit['atomic_commutator_relative_residual']<2e-7 and audit['minimum_state_eigenvalue']>=-1e-10
        and audit['greater']['passed'] and audit['lesser']['passed'])
    if not audit['passed']:
        raise ValueError(f'continuous characteristic quantum audit failed: {audit}')
    return {'greater_by_source':readonly_array(ordered[0]),'lesser_by_source':readonly_array(ordered[1]),
        'greater':readonly_array(ordered[0].sum(axis=0)),'lesser':readonly_array(ordered[1].sum(axis=0)),
        'mean_pulse':readonly_array(mean),'retarded_response':readonly_array(response),
        'exit_state':readonly_array(states[-1]),'states':readonly_array(states),
        'atomic_covariance_by_source':readonly_array(atomic),'sample_ages_s':readonly_array(times*duration,real=True),
        'source_names':('atomic_inflow',)+tuple('jump:'+ch.name for ch in reservoirs.channels),
        'frequencies_rad_s':freq,'residence_time_s':duration,'audit':audit,
        'numerics':{'method':'DOP853 sparse affine microscopic moment lift','rtol':rtol,'atol':atol,
            'max_step_s':maximum,'ode_evaluations':calls,'elapsed_seconds':time.monotonic()-started,
            'complex_variables':size,'sparse_nonzeros':[x.nnz for x in matrices],
            'accumulator_scale':eta.tolist()},
        'scope':'continuous finite atomic covariance and response; numerical convergence is a separate gate'}
