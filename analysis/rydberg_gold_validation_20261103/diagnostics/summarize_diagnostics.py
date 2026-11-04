"""Summarize bounded model tests and retain scoped coarsening evidence."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]


def main():
    gaussian = json.loads((ROOT/'gaussian_trial.json').read_text(encoding='utf-8'))
    profiles = json.loads((ROOT/'linewidth_profiles.json').read_text(encoding='utf-8'))
    audit = json.loads((ROOT/'generator_audit.json').read_text(encoding='utf-8'))
    reproduction = json.loads((ROOT/'legacy_reproduction.json').read_text(encoding='utf-8'))
    source = list(csv.DictReader((ROOT.parent/'sources/fig2b_points.csv').open(encoding='utf-8-sig')))
    observed = {(float(row['probe_power_uw']),row['quantity']):float(row['value'])
                for row in source if row['condition']=='compensated'}
    by_key = {}
    for row in profiles['rows']:
        key = row['variant'],row['total_ground_coherence_rate_mhz'],row['probe_rabi_anchor_mhz'],row['opacity_participation']
        by_key.setdefault(key,[]).append(row)
    scores = []
    for key,rows in by_key.items():
        w=np.array([row['eit_fwhm_mhz']-observed[row['probe_power_uw'],'eit_fwhm'] for row in rows])
        amp=np.array([row['relative_voltage_amplitude'] for row in rows])
        data=np.array([observed[row['probe_power_uw'],'eit_peak'] for row in rows])
        gain=float(amp@data/(amp@amp))
        wr=float(np.sqrt(np.mean(w*w)))
        ar=float(np.sqrt(np.mean((gain*amp-data)**2)))
        scores.append(dict(variant=key[0],total_ground_coherence_rate_mhz=key[1],
            probe_rabi_anchor_mhz=key[2],opacity_participation=key[3],
            width_rmse_mhz=wr,amplitude_rmse_arb=ar,one_global_gain=gain,
            score=(wr/3.)**2+(ar/.735)**2))
    heldout=gaussian['heldout_rows']
    width_error=np.array([row['eit_fwhm_mhz']-row['measured_eit_fwhm_mhz'] for row in heldout])
    amp_error=np.array([row['predicted_amplitude_arb']-row['measured_amplitude_arb'] for row in heldout])
    convergence=gaussian['selected_hypothesis_convergence_rows']
    result=dict(origin='⟦AI⟧ Bounded diagnostics synthesis; physical validation remains incomplete',
        uniform_profile_best={variant:min((row for row in scores if row['variant']==variant),key=lambda row:row['score'])
            for variant in ('legacy_pair','nested_cp','nested_cp_reset')},
        profile_scope='Training observations at0.5,2,6,10uW; Omega_c=3MHz fixed; coarse bounded parameter profiles, not a continuous global fit',
        gaussian_selected=gaussian['selected_training_hypothesis'],
        gaussian_holdout_errors=dict(width_rmse_mhz=float(np.sqrt(np.mean(width_error**2))),
            maximum_absolute_width_error_mhz=float(abs(width_error).max()),
            amplitude_rmse_arb=float(np.sqrt(np.mean(amp_error**2))),
            maximum_absolute_amplitude_error_arb=float(abs(amp_error).max()),
            status='NOT_REPRODUCED_AT_FIGURE_RESOLUTION',
            interpretation='Readout resolution is not experimental standard uncertainty; no statistical rejection level assigned'),
        scoped_predictability_evidence=[dict(variable_or_numerical_configuration='Gaussian radial quadrature + optical scan discretization',
            coarse=dict(radial_nodes=12,optical_scan_step_mhz=1/30),
            reference=dict(radial_nodes=24,optical_scan_step_mhz=.01),
            maximum_absolute_width_difference_mhz=max(abs(row['width_difference_from_nodes12_step_1over30_mhz']) for row in convergence),
            maximum_relative_peak_amplitude_difference=max(abs(row['relative_amplitude_change_from_nodes12_step_1over30_mhz']) for row in convergence),
            powers_uw=[row['probe_power_uw'] for row in convergence],
            scope='Only selected Gaussian equal-waist, Omega_c,peak=3MHz hypothesis atT20,Pc30,LOoff, ground-phase rate0, probe peak anchor4MHz, opacity0.1',
            decision='Coarser joint numerical configuration is adequate relative to source graph resolution here; not a certificate for all radial profiles or parameter values',
            confounding='Both radial count and spectral step changed, so errors cannot be attributed to either independently')],
        legacy_generator=dict(conditional_choi_minimum_per_2pi_mhz=audit['conditional_choi_minimum_per_2pi_mhz'],
            minimum_sampled_stationary_eigenvalue=min(row['min_stationary_eigenvalue'] for row in audit['stationary_samples'])),
        legacy_reproduction_status=reproduction['status'],
        absolute_transmission_validation='Unavailable: source detector signal is arbitrary and absolute incident/dark calibration is not reported')
    paths=[Path(__file__),ROOT/'gaussian_trial.json',ROOT/'linewidth_profiles.json',
           ROOT/'generator_audit.json',ROOT/'legacy_reproduction.json',ROOT.parent/'sources/fig2b_points.csv']
    result['source_sha256']={str(p.relative_to(REPO)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (ROOT/'model_diagnostic_summary.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
