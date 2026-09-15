"""Disposable Gold-point angular and transit diagnostics, not production models."""
import json
import math
from pathlib import Path
import time
import numpy as np
from gabes import observables, constants
from gabes.schemes import fwm


def snapshot(*, transit_khz=100., half_width_mhz=.1, n=3):
    original = observables._gain_matrix_from_chi
    captured = []
    def capture(*args, **kwargs):
        if not captured:
            captured.append((args, kwargs))
        return original(*args, **kwargs)
    center = fwm.branch_center_GHz(.9, -1) - .008
    observables._gain_matrix_from_chi = capture
    try:
        raw = fwm.compute_spectrum(
            .9, T=394.15, P_pump=.6, P_probe=8e-6,
            response_method=fwm.RESPONSE_POLE, phase_detail=fwm.PHASE_ULTRA,
            velocity_step=1., velocity_cutoff=4.,
            coarse_points=n, fine_points=0,
            transit_rate=2*math.pi*transit_khz*1e3,
            scan_min=center-half_width_mhz*1e-3,
            scan_max=center+half_width_mhz*1e-3)
    finally:
        observables._gain_matrix_from_chi = original
    return raw, captured[0]


def main():
    start = time.perf_counter()
    raw, (args, kwargs) = snapshot(half_width_mhz=30., n=401)
    ss, sc, cs, cc, ks, kc, density, dipole, ls = args
    weights = np.exp(-.5*(np.linspace(-30.,30.,401)/1.380173563663881)**2)
    weights /= np.sum(weights)
    offdiag = [np.dot(weights, sc), np.dot(weights, cs)]
    i = len(ss)//2
    gs, gc, _, _ = fwm._ultra_segmented_gain(
        ss[i:i+1], np.array([offdiag[0]]), np.array([offdiag[1]]), cc[i:i+1],
        ks[i:i+1],kc[i:i+1],.0125,density,ls,kwargs['delta_k_z'][i:i+1],
        np.ones(64),fwm._gaussian_overlap_profile(64,.0125,530e-6,330e-6,.32),.6,8e-6)
    gs, gc = observables.pump_depletion_saturation(gs, gc, .6,8e-6)
    baseline = float(raw['G_s'][i])
    record = {
        'unpatched_gold_gain':baseline,
        'naive_1d_curve_angular_convolution': {
            'angular_rms_mhz':1.380173563663881,
            'sc_amplitude_ratio':float(abs(offdiag[0]/sc[i])),
            'cs_amplitude_ratio':float(abs(offdiag[1]/cs[i])),
            'sc_phase_difference_deg':float(np.angle(offdiag[0]/sc[i])*180/math.pi),
            'cs_phase_difference_deg':float(np.angle(offdiag[1]/cs[i])*180/math.pi),
            'probe_gain':float(gs[0]),'conjugate_gain':float(gc[0]),
            'qualifier':'Diagnostic only: convolution of finite-seed one-dimensional lab-detuning curves changes both lab beat and Raman frequency; separate two-dimensional reference keeps lab beat fixed and is not equivalent. Do not adopt as validated angular model.',
        },
        'transit_rate_sensitivity':{},
    }
    for rate in (50.,100.,200.):
        point,_=snapshot(transit_khz=rate)
        record['transit_rate_sensitivity'][str(rate)+' kHz']={'G_s':float(point['G_s'][1]),'G_c':float(point['G_c'][1])}
    record['elapsed_seconds']=time.perf_counter()-start
    Path(__file__).with_name('sensitivity_results.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(record,indent=2))


if __name__ == '__main__':
    main()
