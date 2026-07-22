"""FINAL analytic theory (frozen):
  carrier  : saturated rate equations (4x4)  +  B-block coherences via a
             3x3 Sylvester solve (pump-A mixing exact, given populations)
  response : exact first-order coherence sector given the carrier
             (x(3) + Y(9) linear block, closed-form entries)
All solves are <=12-dim linear blocks with explicit entries (Cramer-able);
no OBE integration, no Liouvillian.

Companion of docs: squeezing_analytic_reconstruction_v1.tex
Beat-corrected: Omega_beat = OMEGA_HF + branch*delta (seeded branch = -1).
At the beat-corrected optimum (Delta=-1.5 GHz, T=110 C, delta*=-280 MHz) this
leading-order truncation reproduces the optimum location and the phase structure
(Doppler-averaged dphase ~1e-2 rad) but its magnitude is ~0.5-0.8x the exact
susceptibility (chi_ss 0.52, chi_cs 0.59, chi_sc 0.72, chi_cc 0.77); the
exponentially amplified gain therefore needs the full compact block / Floquet
solve (ref_solver.py) for the quantitative value.
"""
import numpy as np

GAMMA = 2*np.pi*5.746e6
OMEGA_HF = 2*np.pi*3.035732439e9
OMEGA_EHF = 2*np.pi*361.58e6
CF2 = {(2,2):10/81,(2,3):35/81,(3,2):35/81,(3,3):28/81}
a1 = {2:np.sqrt(3*CF2[(2,2)]), 3:np.sqrt(3*CF2[(2,3)])}
a2 = {2:np.sqrt(3*CF2[(3,2)]), 3:np.sqrt(3*CF2[(3,3)])}
br1 = {2:10/45, 3:35/63}; br2 = {2:35/45, 3:28/63}
EHF = {2:OMEGA_EHF, 3:0.0}
IDX = {1:0, 2:1, 3:2}

def _HB(A, Deff):
    HB = np.zeros((3,3))
    for e in (2,3):
        HB[0,IDX[e]] = HB[IDX[e],0] = A[e]
        HB[IDX[e],IDX[e]] = -(Deff+EHF[e])
    return HB

def carrier(OpA, delta, Deff, Go, ggg, branch=-1):
    Ob = OMEGA_HF + branch*delta   # seeded_sideband_beat: nu_HF + branch*delta
    A = {e: OpA*a1[e]/2 for e in (2,3)}
    B = {e: OpA*a2[e]/2 for e in (2,3)}
    D1  = {e: (Deff+EHF[e])+1j*Go for e in (2,3)}
    D2p = {e: (Deff+EHF[e]+delta+Ob)+1j*Go for e in (2,3)}
    D2m = {e: (Deff+EHF[e]+delta-Ob)+1j*Go for e in (2,3)}
    # --- saturated rate equations ---
    WA = {e: 2*A[e]**2*Go/np.abs(D1[e])**2 for e in (2,3)}
    WB = {e: 2*B[e]**2*Go*(1/np.abs(D2p[e])**2+1/np.abs(D2m[e])**2) for e in (2,3)}
    R = np.zeros((4,4))
    for i,e in enumerate((2,3)):
        ie = 2+i
        R[ie,0]+=WA[e]; R[ie,ie]-=WA[e]; R[ie,1]+=WB[e]; R[ie,ie]-=WB[e]
        R[0,ie]+=WA[e]; R[0,0]-=WA[e];  R[1,ie]+=WB[e]; R[1,1]-=WB[e]
        R[0,ie]+=GAMMA*br1[e]; R[1,ie]+=GAMMA*br2[e]; R[ie,ie]-=GAMMA
    Ms = R.copy(); Ms[0,:]=1.0
    rhs4 = np.zeros(4); rhs4[0]=1.0
    p1,p2,pe2,pe3 = np.linalg.solve(Ms, rhs4)
    pe = {2:pe2, 3:pe3}
    # --- B-block off-diagonal coherences: [HB, D+O] - i g O = 0 ---
    HB = _HB(A, Deff)
    gO = np.zeros((3,3))
    for e in (2,3):
        gO[0,IDX[e]] = gO[IDX[e],0] = Go
        for ep in (2,3): gO[IDX[e],IDX[ep]] = GAMMA
    dvec = np.array([p1, pe2, pe3])
    MO = np.zeros((9,9), complex); rhsO = np.zeros(9, complex)
    def oi(r,b): return 3*r+b
    for r in range(3):
        for b in range(3):
            if r == b:
                MO[oi(r,b), oi(r,b)] = 1.0; rhsO[oi(r,b)] = 0.0
                continue
            MO[oi(r,b), oi(r,b)] += -1j*gO[r,b]
            for k in range(3):
                if k != b: MO[oi(r,b), oi(k,b)] += HB[r,k]
                if k != r: MO[oi(r,b), oi(r,k)] -= HB[k,b]
            rhsO[oi(r,b)] = -HB[r,b]*(dvec[b]-dvec[r])
    O = np.linalg.solve(MO, rhsO).reshape(3,3)
    rho_e1 = {e: O[IDX[e],0] for e in (2,3)}
    rho_1e = {e: O[0,IDX[e]] for e in (2,3)}
    r23 = O[1,2]
    rho_ee = {(2,2):pe2,(3,3):pe3,(2,3):r23,(3,2):np.conj(r23)}
    return dict(p1=p1,p2=p2,pe=pe,
                w1={e:p1-pe[e] for e in (2,3)}, w2={e:p2-pe[e] for e in (2,3)},
                rho_e1=rho_e1,rho_1e=rho_1e,rho_ee=rho_ee,
                A=A,B=B,D1=D1,D2p=D2p,D2m=D2m)

def response(OpA, Os, delta, Deff, Go, ggg, car=None, branch=-1):
    Ob = OMEGA_HF + branch*delta   # seeded_sideband_beat: nu_HF + branch*delta
    if car is None:
        car = carrier(OpA, delta, Deff, Go, ggg, branch)
    A,Bf = car["A"], car["B"]
    Se = {e: Os*a2[e]/2 for e in (2,3)}
    Ce = {e: Os*a1[e]/2 for e in (2,3)}
    HB = _HB(A, Deff)
    gx = np.array([ggg, Go, Go])
    gY = np.zeros((3,3))
    gY[0,0] = 0.0
    for e in (2,3):
        gY[IDX[e],0] = Go; gY[0,IDX[e]] = Go
        for ep in (2,3): gY[IDX[e],IDX[ep]] = GAMMA
    Bvec = np.zeros(3)
    for e in (2,3): Bvec[IDX[e]] = Bf[e]

    n = 12
    Mm = np.zeros((n,n), complex)
    def xi_(b): return b
    def yi(r,b): return 3 + 3*r + b
    for b in range(3):
        Mm[xi_(b), xi_(b)] += delta - 1j*gx[b]
        for bp in range(3):
            Mm[xi_(b), xi_(bp)] -= HB[bp,b]
        for e in (2,3):
            Mm[xi_(b), yi(IDX[e],b)] += Bvec[IDX[e]]
    for r in range(3):
        for b in range(3):
            Mm[yi(r,b), yi(r,b)] += -Ob - 1j*gY[r,b]
            for k in range(3):
                Mm[yi(r,b), yi(k,b)] += HB[r,k]
                Mm[yi(r,b), yi(r,k)] -= HB[k,b]
            if r != 0:
                Mm[yi(r,b), xi_(b)] += Bvec[r]

    p2 = car["p2"]; w2 = car["w2"]; rho_ee = car["rho_ee"]
    rho_e1 = car["rho_e1"]; rho_1e = car["rho_1e"]

    def solve(d_x, g_Y):
        rhs = np.zeros(n, complex)
        for b in range(3): rhs[xi_(b)] = -d_x[b]
        for r in range(3):
            for b in range(3): rhs[yi(r,b)] = -g_Y[r,b]
        U = np.linalg.solve(Mm, rhs)
        return U[:3], U[3:].reshape(3,3)

    d = np.zeros(3, complex)
    d[0] = sum(Se[e]*rho_e1[e] for e in (2,3))
    for ep in (2,3):
        d[IDX[ep]] = -Se[ep]*w2[ep] + sum(Se[e]*rho_ee[(e,ep)] for e in (2,3) if e!=ep)
    x, Y = solve(d, np.zeros((3,3), complex))
    chi_ss = sum(a2[e]*np.conj(x[IDX[e]]) for e in (2,3))/Os
    chi_cs = sum(a1[e]*Y[IDX[e],0] for e in (2,3))/Os

    d2 = np.zeros(3, complex)
    d2[0] = -sum(Ce[e]*Bf[e]*w2[e]/np.conj(car["D2p"][e]) for e in (2,3))
    g2Y = np.zeros((3,3), complex)
    for e in (2,3):
        r = IDX[e]
        g2Y[r,0] = Ce[e]*car["p1"] - sum(rho_ee[(e,ep)]*Ce[ep] for ep in (2,3))
        for bp in (2,3):
            g2Y[r,IDX[bp]] = Ce[e]*rho_1e[bp]
    g2Y[0,0] = -sum(rho_1e[e]*Ce[e] for e in (2,3))
    x2, Y2 = solve(d2, g2Y)
    chi_sc = sum(a2[e]*np.conj(x2[IDX[e]]) for e in (2,3))/Os
    chi_cc = sum(a1[e]*Y2[IDX[e],0] for e in (2,3))/Os
    return chi_ss, chi_cs, chi_sc, chi_cc
