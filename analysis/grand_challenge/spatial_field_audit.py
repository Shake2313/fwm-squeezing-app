"""Immutable stationary-center finite-aperture field audit.

python -m analysis.grand_challenge.spatial_field_audit --output NEW.json --plot NEW.png
"""

import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy

from gabes import core, constants as c
from gabes.fwm_quantum.kinetic import CarrierGeometry, _pump_system
from gabes.fwm_quantum.normalization import optical_carriers
from gabes.fwm_quantum.spatial import spatial_velocity_atom
from gabes.fwm_quantum.spatial_cell import StationarySpatialMedium, stationary_spatial_cell
from gabes.quantum.contracts import AnalysisFrequencyAxis, GeneratorFrequencyAxis
from gabes.quantum.periodic import PeriodicAtomicNoise
from gabes.quantum.periodic_field import PeriodicFieldPorts, periodic_field_couplings
from gabes.quantum.readout import DetectorResponse
from gabes.quantum.spatial_modes import TransverseModeGrid
from gabes.quantum.traveling import LocalNambuGenerator
from .normalization_audit import conditional_inputs
from .reference.periodic_field import ordered_current_psd
from .reference.spatial_field import spatial_field_reference


ROOT = Path(__file__).resolve().parents[2]
TOLERANCES = {'local_M_D_relative': 1e-8, 'mean_reference_relative': 1e-9,
    'cell_gain_absolute': 1e-7, 'cell_spectrum_absolute_dB': 1e-6,
    'mean_tangent_relative': 1e-6, 'ordered_current_relative': 1e-10}


def source_manifest():
    # Explicit dependencies, including the existing atom/normalization code.
    # Unrelated concurrently developed downstream modules are not consumed.
    files = ['gabes/'+s+'.py' for s in ('core', 'kernels', 'constants', 'atoms', 'observables',
        'lineshape', 'hyperfine', 'species', 'zeeman', 'doppler', 'schemes/fwm', 'schemes/base')]
    files += ['gabes/quantum/'+s+'.py' for s in ('contracts', 'diffusion', 'reservoirs', 'periodic',
        'periodic_field', 'traveling', 'channels', 'readout', 'sidebands', 'spatial', 'spatial_modes')]
    files += ['gabes/fwm_quantum/'+s+'.py' for s in ('inputs', 'model', 'normalization', 'periodic',
        'periodic_cell', 'kinetic', 'spatial', 'spatial_cell', 'readout', 'field')]
    files += ['analysis/grand_challenge/'+s+'.py' for s in ('normalization_audit', 'spatial_field_audit',
        'reference/periodic_field', 'reference/periodic_qrt', 'reference/spatial_field',
        'reference/uncoupled_d1', 'reference/direct_readout')]
    files += ['tests/quantum/test_spatial_field.py']
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files}


def complex_array(value):
    value = np.asarray(value)
    return {'real': value.real.tolist(), 'imag': value.imag.tolist()}


def relative(a, b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b), np.finfo(float).tiny))


def medium(order, angles=(.005, -.004), **kwargs):
    inputs = conditional_inputs()
    geometry = CarrierGeometry.vacuum_beams(inputs, probe_angle_rad=angles[0], conjugate_angle_rad=angles[1])
    return StationarySpatialMedium(inputs, geometry,
        TransverseModeGrid.rectangle(.0004, .0003, order_x=order), **kwargs)


def comparison(a, b):
    return {name: relative(getattr(a, name), getattr(b, name))
            for name in ('drift', 'noise_greater', 'noise_lesser')}


def local_controls():
    original = medium(2)
    u = original.modes.mode_values_m_inverse.copy()
    u[:, 0] *= np.sqrt([1.3, .7])
    u[:, 1] *= np.sqrt([.8, 1.2])*np.exp(1j*np.array([.3, -.4]))
    a = StationarySpatialMedium(original.inputs, original.geometry, replace(original.modes, mode_values_m_inverse=u))
    z, beta = .006, np.array([5e6+1e6j, .5e6j])
    axis = GeneratorFrequencyAxis([0., 2*np.pi*1e6], 'local independent optical phase reference')
    actual = a.local_field(z, beta, axis)
    k0, kp, kc = a.geometry.wavevectors_rad_m
    refs = [spatial_field_reference(a.h0, a.reservoirs, a.beat_rad_s, a.ports.lowering_operators,
        a.ports.coupling_s_inverse_sqrt_flux, beta, a.modes.positions_m, a.modes.area_weights_m2,
        u, a.inputs.uniform_area_m2, np.array([kp-k0, kc-k0]), z, a.inputs.number_density_m3,
        axis.omega_rad_s, phase_samples=p) for p in (16, 24)]
    errors = {name: relative(getattr(actual, name), refs[-1][key]) for name, key in
              (('drift', 'drift'), ('noise_greater', 'greater'), ('noise_lesser', 'lesser'))}
    ref_change = {key: relative(refs[0][key], refs[1][key]) for key in ('greater', 'lesser')}
    mean_error = relative(a.mean_rate(z, beta), refs[-1]['mean_rate'])
    controls = [StationarySpatialMedium(a.inputs, a.geometry, a.modes, mean_order=m, response_order=r)
                .local_field(z, beta, axis) for m, r in ((5, 3), (5, 4))]
    refinement = [comparison(actual, controls[0]), comparison(controls[0], controls[1])]
    grids = [2, 4, 8, 12]
    grid_values = [medium(n).local_field(z, beta, axis) for n in grids]
    quadrature = [comparison(x, y) for x, y in zip(grid_values[:-1], grid_values[1:])]
    return {'z_m': z, 'carriers_sqrt_flux': complex_array(beta), 'frequencies_rad_s': axis.omega_rad_s.tolist(),
        'reference_profiles_m_inverse': complex_array(u), 'reference_positions_m': a.modes.positions_m.tolist(),
        'reference_area_weights_m2': a.modes.area_weights_m2.tolist(), 'raw_optical_phase_reference_errors': errors,
        'reference_phase_samples': [16, 24], 'reference_phase_changes': ref_change,
        'reference_mean_relative_error': mean_error, 'atomic_mean_response_orders': [[4, 3], [5, 3], [5, 4]],
        'adjacent_atomic_changes': refinement, 'uniform_profile_quadrature_orders_x': grids,
        'adjacent_quadrature_changes': quadrature, 'local_audit': actual.audit(),
        'passed': bool(max(errors.values()) < TOLERANCES['local_M_D_relative']
            and max(ref_change.values()) < TOLERANCES['local_M_D_relative']
            and mean_error < TOLERANCES['mean_reference_relative']
            and max(v for row in refinement for v in row.values()) < TOLERANCES['local_M_D_relative']
            and max(quadrature[-1].values()) < TOLERANCES['local_M_D_relative'] and actual.audit()['passed'])}


def frozen_moving_grating_control(*, mean_order=4, response_order=3):
    """Deliberately invalid analysis-only reduction; never return it as a field.

    Freezing the second phase of a moving torus drops omega_2*d_psi K from
    the dynamic commutator equation. PSD alone cannot validate this closure.
    """
    inputs = conditional_inputs()
    geometry = CarrierGeometry.vacuum_beams(inputs, probe_angle_rad=.005, conjugate_angle_rad=-.004)
    velocity = np.array([150., 0., 380.])
    beta = np.array([5e6, .5e6j])*np.sqrt(125.)
    model = spatial_velocity_atom(inputs, geometry, velocity, beta, mean_orders=(mean_order, mean_order))
    atom, phase = model.atom, .7

    def section(labels, data, harmonics):
        return np.array([np.einsum('q,q...->...', np.exp(-1j*labels[labels[:, 0] == h, 1]*phase),
            data[labels[:, 0] == h]) for h in harmonics])

    rho = section(atom.state_labels, atom.state_coefficients, range(-mean_order, mean_order+1))
    drift = section(atom.generator_labels, atom.drift_coefficients, [-1, 0, 1])
    diffusion = np.array([section(atom.state_labels, d, range(-mean_order, mean_order+1)) for d in atom.diffusion_by_reservoir])
    h0, reservoirs = _pump_system(inputs, inputs.one_photon_rad_s-geometry.wavevectors_rad_m[0]@velocity,
                                  'uniform-zeeman-rms')
    v = model.drive_operators_rad_s[0]+model.drive_operators_rad_s[1].conj().T*np.exp(1j*phase)
    invalid = PeriodicAtomicNoise(h0, v, reservoirs, model.convective_frequencies_rad_s[0], atom.operators,
        rho, drift, diffusion, {'scope': 'invalid frozen moving grating control'})
    ports = PeriodicFieldPorts(('probe', 'conjugate'), [1, -1], model.optical_lowering_operators,
                               model.coupling_s_inverse_sqrt_flux)
    lift, b, c0, signs, labels = periodic_field_couplings(invalid, ports, response_order=response_order)
    omega = 2*np.pi*1e6
    cr = np.linalg.solve((-1j*omega*np.eye(len(lift.drift))-lift.drift).T, c0.T).T
    candidate = LocalNambuGenerator(GeneratorFrequencyAxis([omega], 'invalid frozen grating control'), signs,
        labels, np.array([cr@b]), (cr@lift.greater_by_reservoir@cr.conj().T)[:, None],
        (cr@lift.lesser_by_reservoir@cr.conj().T)[:, None], lift.reservoir_names)
    audit = candidate.audit()
    rejected = False
    try:
        StationarySpatialMedium(inputs, geometry, TransverseModeGrid.rectangle(.0004, .0003, order_x=2),
                                velocity_m_s=velocity)
    except ValueError:
        rejected = True
    return {'velocity_m_s': velocity.tolist(), 'carriers_sqrt_flux': complex_array(beta),
        'grating_phase_rad': phase, 'convective_frequencies_rad_s': model.convective_frequencies_rad_s.tolist(),
        'response_order': response_order, 'mean_orders': [mean_order, mean_order], 'linear_density_m_inverse': 1.,
        'candidate_audit': audit, 'production_velocity_guard_rejected': rejected,
        'scope': 'negative control only; drops convective derivative; not evidence against every possible transport closure',
        'passed': bool(not audit['passed'] and audit['noise_greater_by_reservoir']['passed']
                       and audit['noise_lesser_by_reservoir']['passed'] and rejected)}


def summarize_cell(name, a, detector, *, rtol=2e-9):
    started = time.monotonic()
    result = stationary_spatial_cell(a, detector, propagation_rtol=rtol)
    transfer, spectrum = result['transfer'], result['spectrum']
    vac = np.diag([1., 1., 0., 0.])
    direct = []
    for j, w in enumerate(detector.analysis_axis.omega_rad_s):
        i = np.flatnonzero(transfer.frequency_axis.omega_rad_s == w)[0]
        greater = transfer.transfer[i]@vac@transfer.transfer[i].conj().T+transfer.noise_greater[i]
        direct.append(ordered_current_psd(greater, result['output_carriers_sqrt_flux'], detector.transmissions,
            detector.current_response[j], detector.balance, c.ELEMENTARY_CHARGE))
    readout_error = relative(np.array(direct), spectrum.quantum_psd_A2_Hz)
    row = {'name': name, 'area_quadrature_points': len(a.modes.positions_m), 'wavevectors_rad_m': a.geometry.wavevectors_rad_m.tolist(),
        'loop_wavevector_rad_m': (-a.geometry.mismatch_rad_m).tolist(),
        'probe_power_gain': result['probe_power_gain'], 'conjugate_power_gain': result['conjugate_power_gain'],
        'spectrum_dB': spectrum.quantum_db.tolist(), 'quantum_ratio': spectrum.quantum_ratio.tolist(),
        'quantum_psd_A2_Hz': spectrum.quantum_psd_A2_Hz.tolist(), 'direct_current_relative_error': readout_error,
        'output_carriers_sqrt_flux': complex_array(result['output_carriers_sqrt_flux']),
        'output_powers_W': result['output_powers_W'].tolist(),
        'generator_frequencies_rad_s': transfer.frequency_axis.omega_rad_s.tolist(),
        'transfer': complex_array(transfer.transfer), 'added_greater': complex_array(transfer.noise_greater),
        'added_lesser': complex_array(transfer.noise_lesser), 'sideband_covariances': result['sideband_covariances'].tolist(),
        'transfer_audit': transfer.audit(), 'minimum_channel_cp_eigenvalue': result['minimum_channel_cp_eigenvalue'],
        'minimum_covariance_uncertainty_eigenvalue': spectrum.minimum_covariance_uncertainty_eigenvalue,
        'propagation_rtol': rtol, 'quantum_propagation': result['quantum_propagation'],
        'mean_ode_evaluations': result['mean_ode_evaluations'], 'elapsed_seconds': time.monotonic()-started}
    row['passed'] = bool(transfer.audit()['passed'] and row['minimum_channel_cp_eigenvalue'] >= -1e-10
        and row['minimum_covariance_uncertainty_eigenvalue'] >= -1e-10
        and readout_error < TOLERANCES['ordered_current_relative'])
    print(f'{name}: {row["elapsed_seconds"]:.2f}s, Gp={row["probe_power_gain"]:.10f}, S(1MHz)={row["spectrum_dB"][1]:.10f}dB', flush=True)
    return row, result


def cell_change(a, b):
    return {'maximum_gain_absolute_change': max(abs(a[k]-b[k]) for k in ('probe_power_gain', 'conjugate_power_gain')),
        'maximum_spectrum_absolute_change_dB': float(np.max(abs(np.array(a['spectrum_dB'])-b['spectrum_dB'])))}


def change_passed(row):
    return (row['maximum_gain_absolute_change'] < TOLERANCES['cell_gain_absolute']
        and row['maximum_spectrum_absolute_change_dB'] < TOLERANCES['cell_spectrum_absolute_dB'])


def build_report():
    before = source_manifest()
    with core.blas_single_thread():
        report = _build_report()
    report['source_sha256'] = before
    report['source_stable_during_run'] = before == source_manifest()
    report['expected_controls_passed'] &= report['source_stable_during_run']
    return report


def _build_report():
    print('Independent raw-phase response and QRT, atomic and spatial refinement', flush=True)
    local = local_controls()
    print('Negative control: freezing a moving grating', flush=True)
    negative = frozen_moving_grating_control()
    negative_refined = frozen_moving_grating_control(mean_order=5, response_order=4)
    negative_change = abs(negative['candidate_audit']['maximum_commutator_relative_residual']
        -negative_refined['candidate_audit']['maximum_commutator_relative_residual'])
    inputs = conditional_inputs()
    detector = DetectorResponse(AnalysisFrequencyAxis.from_hz([.1e6, 1e6, 4e6]), [.85, .85],
        np.ones((3, 2)), 1., np.zeros(3), 'declared balanced conditional detector, no measured data')
    cases, selected = [], None
    for name, order, angles in (('collinear', 1, (0., 0.)), ('symmetric_5mrad', 2, (.005, -.005)),
                              ('unequal_Nx4', 4, (.005, -.004)), ('unequal_Nx8', 8, (.005, -.004)),
                              ('unequal_Nx12', 12, (.005, -.004))):
        row, result = summarize_cell(name, medium(order, angles), detector)
        cases.append(row)
        if order == 8:
            selected = result
    refined, _ = summarize_cell('unequal_Nx8_quantum_ODE_refined', medium(8), detector, rtol=2e-11)
    spatial_changes = [cell_change(a, b) for a, b in zip(cases[2:-1], cases[3:])]
    ode_change = cell_change(cases[3], refined)
    a = medium(8)
    base, scale = a.integrate_carriers()
    tighter, refined_scale = a.integrate_carriers(rtol=2e-12, atol=2e-14)
    z = np.linspace(0., inputs.length_m, 17)
    mean_change = relative(base.sol(z)*scale, tighter.sol(z)*refined_scale)
    initial = np.array([scale, 0.]); direction = np.array([.4+.2j, -.3+.7j])
    out = []
    for sign in (-1, 1):
        trajectory, normalization = a.integrate_carriers(initial_amplitudes=initial+sign*100*direction)
        out.append(trajectory.y[:, -1]*normalization)
    t = selected['transfer']; dc = np.flatnonzero(t.frequency_axis.omega_rad_s == 0)[0]
    tangent_error = relative((out[1]-out[0])/200, (t.transfer[dc]@np.r_[direction, direction.conj()])[:2])
    predecessor = 'docs/grand_challenge/s1_spatial_report.json'
    ledger = inputs.consumed_inputs(detector, convention='uniform-zeeman-rms')
    return {'schema': 'gabes.stationary_spatial_field_audit.v2',
        'scope': 'conditional stationary atomic centers, fixed classical pump, two fixed normalized transverse profiles',
        'absolute_hot_vapor_prediction': False, 'experimental_validation': False,
        'inputs': asdict(inputs), 'dipole_convention': 'uniform-zeeman-rms', 'velocity_m_s': [0., 0., 0.],
        'consumed_scalar_ledger': ledger['scalars'], 'derived_atomic_optical_inputs': ledger['derived'],
        'input_evidence': 'conditional declared values; no independent measurement evidence or target-data fits',
        'aperture_m': [.0004, .0003], 'profiles': 'identical normalized top-hats in a fixed rectangle; synthetic conditional mode choice',
        'detector_frequencies_Hz': (detector.analysis_axis.omega_rad_s/(2*np.pi)).tolist(),
        'detector_transmissions': detector.transmissions.tolist(), 'detector_balance': detector.balance,
        'mean_order': 4, 'response_order': 3, 'mean_ODE_rtol': 2e-10, 'mean_ODE_atol': 2e-12,
        'numerical_tolerances': TOLERANCES, 'independent_local_controls': local,
        'invalid_moving_grating_control': negative, 'cells': cases, 'refined_quantum_ODE_cell': refined,
        'refined_invalid_moving_grating_control': negative_refined,
        'invalid_grating_commutator_residual_absolute_refinement_change': negative_change,
        'cell_quadrature_orders_x': [4, 8, 12], 'adjacent_cell_quadrature_changes': spatial_changes,
        'quantum_ODE_change_at_Nx8': ode_change, 'mean_ODE_refinement_relative_trajectory_change': mean_change,
        'nonlinear_mean_map_DC_tangent_relative_error': tangent_error,
        'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
                        'blas_threads': 1, 'timing_scope': 'elapsed audit time, not a performance benchmark'},
        'predecessor_evidence': {'path': predecessor, 'sha256': hashlib.sha256((ROOT/predecessor).read_bytes()).hexdigest(),
                                 'scope': 'immutable historical atomic-only evidence'},
        'expected_controls_passed': bool(local['passed'] and negative['passed'] and negative_refined['passed']
            and negative_change < TOLERANCES['local_M_D_relative'] and all(row['passed'] for row in cases)
            and refined['passed'] and change_passed(spatial_changes[-1]) and change_passed(ode_change)
            and mean_change < TOLERANCES['mean_reference_relative'] and tangent_error < TOLERANCES['mean_tangent_relative']),
        'remaining': ['moving-atom transport and boundary/inter-slice noise', 'transverse walkoff, diffraction and additional optical modes',
            'full hyperfine/Zeeman and collisions', 'self-consistent pump depletion', 'measured inputs and held-out experimental spectra']}


def save_plot(report, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout='constrained')
    for index, label in ((0, 'Collinear'), (1, '+5 / -5 mrad'), (4, '+5 / -4 mrad, Nx=12')):
        axes[0].plot(np.array(report['detector_frequencies_Hz'])/1e6, report['cells'][index]['spectrum_dB'], 'o-', label=label)
    axes[0].set(xlabel='Analysis frequency (MHz)', ylabel='Intensity-difference noise (dB / SQL)',
                title='Three calculated frequencies; fixed 400 x 300 um aperture')
    axes[0].legend(fontsize=8)
    changes = report['adjacent_cell_quadrature_changes']
    axes[1].semilogy([0, 1], [r['maximum_spectrum_absolute_change_dB'] for r in changes], 'o-', label='Maximum change over three frequencies')
    axes[1].axhline(TOLERANCES['cell_spectrum_absolute_dB'], color='gray', linestyle='--', label='Declared tolerance')
    axes[1].set_xticks([0, 1], ['Nx 4 -> 8', 'Nx 8 -> 12'])
    axes[1].set(ylabel='Absolute change in S (dB)', title='Unequal-angle spatial quadrature refinement')
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle('Conditional stationary-center field: fixed pump and two transverse profiles', fontsize=12)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as stream:
        fig.savefig(stream, format='png', dpi=160)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plot', type=Path)
    args = parser.parse_args(argv)
    for path in (args.output, args.plot):
        if path is not None and path.exists():
            raise FileExistsError(path)
    report = build_report()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print(f'Wrote {args.output}; expected controls passed={report["expected_controls_passed"]}', flush=True)
    if args.plot is not None:
        save_plot(report, args.plot)
    return 0 if report['expected_controls_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
