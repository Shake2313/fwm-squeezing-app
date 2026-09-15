"""Frozen off-diagonal closure compared with source-listed held-out points."""
import json
from pathlib import Path
import time
import numpy as np
from gabes import observables
from gabes.schemes import fwm

RESIDUAL = .7763997305


def evaluate(p, *, constant_factor=None):
    original = observables._gain_matrix_from_chi
    captured = []
    def capture(*args, **kwargs):
        if not captured:
            captured.append((args, kwargs))
        return original(*args, **kwargs)
    center = fwm.branch_center_GHz(p['opd'], -1) + p['tpd']*1e-3
    wp,ws = p['pump_waist_um']*1e-6,p['probe_waist_um']*1e-6
    pump,seed=p['pump_mw']*1e-3,p['probe_uw']*1e-6
    L=p['cell_mm']*1e-3
    observables._gain_matrix_from_chi = capture
    try:
        raw = fwm.compute_spectrum(
            p['opd'], T=p['temp_c']+273.15, P_pump=pump, P_probe=seed,
            w_pump=wp,w_probe=ws,L=L,pump_probe_angle_deg=p['seeded_angle_deg'],
            response_method=fwm.RESPONSE_POLE, phase_detail=fwm.PHASE_ULTRA,
            velocity_step=1., velocity_cutoff=4., coarse_points=3, fine_points=0,
            scan_min=center-.0001, scan_max=center+.0001)
    finally:
        observables._gain_matrix_from_chi = original
    args,kwargs=captured[0]
    ss,sc,cs,cc,ks,kc,density,dipole,ls=args
    factor=(RESIDUAL*wp*wp/(wp*wp+ws*ws) if constant_factor is None
            else float(constant_factor))
    gs,gc,_,_=fwm._ultra_segmented_gain(
        ss,sc*factor,cs*factor,cc,ks,kc,L,density,ls,kwargs['delta_k_z'],
        np.ones(64),fwm._gaussian_overlap_profile(64,L,wp,ws,p['seeded_angle_deg']),pump,seed)
    gs,gc=observables.pump_depletion_saturation(gs,gc,pump,seed)
    return {'input':p,'G_s_before':float(raw['G_s'][1]),'G_c_before':float(raw['G_c'][1]),
            'G_s_candidate':float(gs[1]),'G_c_candidate':float(gc[1]),
            'factor':factor,'floquet_status':raw['floquet_convergence']['status']}


def main():
    start=time.perf_counter()
    references=json.loads(Path(__file__).parents[1].joinpath('reference_points.json').read_text(encoding='utf-8'))
    data={'fixed_residual':RESIDUAL,'fit_target':15.5,'cases':{},'nearby_sweeps':{}}
    for case in [references['gold'],*references['held_out']]:
        row=evaluate(case['params'])
        row['literature_targets']=case['targets']
        data['cases'][case['id']]=row
    p=references['gold']['params']
    for key,values in {'pump_mw':[480,600,720],'temp_c':[116,121,126],
                       'pump_waist_um':[477,530,583],'probe_waist_um':[297,330,363],
                       'seeded_angle_deg':[.22,.32,.42],'opd':[.8,.9,1.0],
                       'tpd':[-10.,-8.,-6.]}.items():
        data['nearby_sweeps'][key]=[evaluate(dict(p,**{key:v})) for v in values]
    data['elapsed_seconds']=time.perf_counter()-start
    Path(__file__).with_name('heldout_results.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(data,indent=2))


if __name__ == '__main__':
    main()
