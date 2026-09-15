"""Disposable read-only production exploration: no physical model modifications."""
import json
import math
from pathlib import Path
import time
import numpy as np
from gabes.schemes import fwm
from gabes import observables, constants, doppler


def main():
    started = time.perf_counter()
    captured = []
    original = observables._gain_matrix_from_chi

    def capture(*args, **kwargs):
        if not captured:
            captured.append((args, kwargs))
        return original(*args, **kwargs)

    center = fwm.branch_center_GHz(0.9, -1) - .008
    observables._gain_matrix_from_chi = capture
    try:
        raw = fwm.compute_spectrum(
            .9, T=394.15, P_pump=.6, P_probe=8e-6,
            response_method=fwm.RESPONSE_POLE, phase_detail=fwm.PHASE_ULTRA,
            velocity_step=1., velocity_cutoff=4.,
            coarse_points=3, fine_points=0,
            scan_min=center-.0001, scan_max=center+.0001,
        )
    finally:
        observables._gain_matrix_from_chi = original
    args, kwargs = captured[0]
    ss, sc, cs, cc, ks, kc, density, dipole, ls = args
    assert len(ss) == 3
    profile = fwm._gaussian_overlap_profile(64, .0125, 530e-6, 330e-6, .32)

    def propagate(offdiag=1., all_chi=1.):
        gs, gc, _, _ = fwm._ultra_segmented_gain(
            ss, sc*offdiag, cs*offdiag, cc, ks, kc, .0125,
            density, ls*all_chi, kwargs['delta_k_z'],
            np.ones(64), profile, .6, 8e-6)
        gs, gc = observables.pump_depletion_saturation(gs, gc, .6, 8e-6)
        return float(gs[1]), float(gc[1])

    target = 15.5
    lo, hi = 0., 1.
    for _ in range(48):
        mid = (lo+hi)/2
        if propagate(offdiag=mid)[0] > target:
            hi = mid
        else:
            lo = mid
    fit = (lo+hi)/2
    kp, ks0, _ = fwm.seeded_option_a_wavenumbers(.9, np.array([center]))
    angular = doppler.noncollinear_raman_rms_budget(394.15, kp, ks0[0], math.radians(.32))
    atom = fwm.collisional_atom(394.15)
    data = {
        'operating_point': {'opd_GHz': .9, 'tpd_MHz': -8., 'T_K':394.15, 'pump_W':.6, 'seed_W':8e-6},
        'base_gain': [float(raw['G_s'][1]), float(raw['G_c'][1])],
        'propagator_parity': propagate(),
        'offdiag_only': {str(c): propagate(offdiag=c) for c in [0., .25, .5, .720626, 1.]},
        'whole_response': {str(c): propagate(all_chi=c) for c in [.4,.5,.6,.720626,1.]},
        'offdiag_fit': {'target':target,'coefficient':fit,'gain':propagate(offdiag=fit),
            'residual_after_transverse_factor':fit/(530**2/(530**2+330**2))},
        'transverse_factor':530**2/(530**2+330**2),
        'angular_budget':angular,
        'gamma_collision_over_2pi_hz': atom.ground_collision_dephasing_rate/(2*math.pi),
        'pump_rabi_over_2pi_MHz': fwm.rabi_freq(.6,530e-6)/(2*math.pi*1e6),
        'x_model': math.acosh(math.sqrt(raw['G_s'][1])),
        'x_exp': math.acosh(math.sqrt(target)),
        'raw_power_ratio':111/8,
        'x_raw_power_ratio':math.acosh(math.sqrt(111/8)),
        'exact_evaluation_settings': {'coarse_points':3,'fine_points':0,
            'scan_min_GHz':center-.0001,'scan_max_GHz':center+.0001,
            'response_method':fwm.RESPONSE_POLE,'velocity_step':1.,
            'velocity_cutoff':4., 'phase_detail':fwm.PHASE_ULTRA,
            'line_strength_residual':.74,'macroscopic_norm':1/12,
            'waist_pump_m':530e-6,'waist_probe_m':330e-6,'cell_length_m':.0125,
            'crossing_angle_deg':.32,'floquet_order':3,'gamma_transit_over_2pi_hz':1e5},
        'elapsed_seconds':time.perf_counter()-started,
        'model_boundary':'Offdiagonal mixing only; Gaussian mode participation assumes local nonlinear coupling proportional to pump intensity. Strong-pump correction is a Fermi closure, not controlled spatial integration. One residual calibrated; no separate Zeeman, polarization, angular, or Raman fits.',
    }
    data['x_ratio_canonical'] = data['x_exp']/data['x_model']
    data['x_ratio_raw_power'] = data['x_raw_power_ratio']/data['x_model']
    Path(__file__).with_name('candidate_results.json').write_text(json.dumps(data,indent=2)+'\n', encoding='utf-8')
    print(json.dumps(data, indent=2))


if __name__ == '__main__':
    main()
