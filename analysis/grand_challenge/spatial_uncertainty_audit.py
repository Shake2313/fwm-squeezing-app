"""Conditional joint gain/S-minus uncertainty for stationary spatial fields.

python -m analysis.grand_challenge.spatial_uncertainty_audit --output NEW.json --plot NEW.png
"""

import argparse
import copy
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy

from gabes import core
from gabes.fwm_quantum.inputs import ReducedPowerInputs
from gabes.fwm_quantum.kinetic import CarrierGeometry
from gabes.fwm_quantum.spatial_cell import StationarySpatialMedium, stationary_spatial_cell
from gabes.quantum.contracts import AnalysisFrequencyAxis, ParameterEvidence, readonly_array
from gabes.quantum.readout import DetectorResponse
from gabes.quantum.spatial_modes import TransverseModeGrid
from gabes.quantum.spectrum_analysis import SpectralTrace, analyze_squeezing_bands
from gabes.quantum.uncertainty import CorrelationEvidence, JointInputUncertainty, propagate_first_order


ROOT = Path(__file__).resolve().parents[2]
PARENT = ROOT/'docs/grand_challenge/s1_spatial_field_report_v2.json'
INPUT_IDS = ('pump_power_W', 'conjugate_angle_rad')
OUTPUT_IDS = ('probe_power_gain', 'conjugate_power_gain', 'R_0.1MHz', 'R_1MHz', 'R_4MHz')
SIGMA = np.array([.006, 1e-5])
COARSE_STEPS = np.array([.0006, 1e-6])
TOLERANCES = {'nominal_gain_absolute': 1e-7, 'nominal_spectrum_dB_absolute': 1e-6,
              'weighted_response_absolute': 1e-8, 'weighted_response_relative': 1e-3,
              'standard_uncertainty_absolute': 1e-8, 'standard_uncertainty_relative': 1e-3,
              'quantum_eigenvalue_absolute': 1e-10}


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    return value


def source_manifest():
    files = ['gabes/'+name+'.py' for name in (
        'constants', 'core', 'kernels', 'atoms', 'hyperfine', 'species', 'zeeman', 'doppler',
        'observables', 'lineshape', 'schemes/base', 'schemes/fwm')]
    files += ['gabes/quantum/'+name+'.py' for name in (
        'contracts', 'reservoirs', 'diffusion', 'channels', 'traveling', 'sidebands', 'readout',
        'periodic', 'periodic_field', 'spatial_modes', 'uncertainty', 'spectrum_analysis')]
    files += ['gabes/fwm_quantum/'+name+'.py' for name in (
        'inputs', 'normalization', 'model', 'field', 'readout', 'kinetic', 'periodic',
        'periodic_cell', 'spatial_cell')]
    files += ['analysis/grand_challenge/spatial_uncertainty_audit.py']
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in files}


def load_parent(path=PARENT):
    """Read historical evidence once; changed sources are reported, not concealed.

    Saved nominal outputs are comparisons only. Every current nominal and
    perturbed point is recomputed, including when the parent hashes still match.
    """
    path = Path(path)
    body = path.read_bytes()
    parent = json.loads(body)
    if (parent.get('schema') != 'gabes.stationary_spatial_field_audit.v2'
            or parent.get('expected_controls_passed') is not True
            or parent.get('source_stable_during_run') is not True):
        raise ValueError('passed immutable stationary spatial v2 parent required')
    if parent.get('velocity_m_s') != [0., 0., 0.]:
        raise ValueError('parent must describe stationary centers')
    if parent.get('detector_frequencies_Hz') != [1e5, 1e6, 4e6]:
        raise ValueError('parent must have the declared three RF frequencies')
    if parent.get('dipole_convention') != 'uniform-zeeman-rms':
        raise ValueError('parent dipole convention does not match this audit')
    inputs = ReducedPowerInputs(**parent['inputs'])
    if inputs.phase_mismatch_rad_m != 0:
        raise ValueError('explicit geometry excludes a separate scalar mismatch')
    rows = {row['name']: row for row in parent['cells']}
    for order in (8, 12):
        row = rows[f'unequal_Nx{order}']
        wavevectors = CarrierGeometry.vacuum_beams(inputs, probe_angle_rad=.005,
                                                   conjugate_angle_rad=-.004).wavevectors_rad_m
        if row.get('passed') is not True or not np.allclose(row['wavevectors_rad_m'], wavevectors, rtol=1e-13, atol=1e-9):
            raise ValueError('parent reference rows need passing +5/-4 mrad geometry')
    changed = []
    for name, digest in parent['source_sha256'].items():
        source = (ROOT/name).resolve()
        if not source.is_relative_to(ROOT.resolve()):
            raise ValueError('parent source paths must stay in the repository')
        current = hashlib.sha256(source.read_bytes()).hexdigest() if source.is_file() else None
        if current != digest:
            changed.append({'path': name, 'parent_sha256': digest, 'current_sha256': current})
    return parent, {'path': str(path.resolve()), 'sha256': hashlib.sha256(body).hexdigest(),
        'changed_parent_sources': changed, 'parent_code_matches_current_sources': not changed,
        'scope': 'historical fixture evidence; current nominal recomputed, not cached from this parent'}


def default_detector(parent):
    ledger = parent['consumed_scalar_ledger']
    frequencies = parent['detector_frequencies_Hz']
    response = [[complex(ledger[f'{arm}_response_real_{i}'][0], ledger[f'{arm}_response_imag_{i}'][0])
                 for arm in ('probe', 'conjugate')] for i in range(len(frequencies))]
    return DetectorResponse(AnalysisFrequencyAxis.from_hz(frequencies), parent['detector_transmissions'],
        response, parent['detector_balance'], [ledger[f'electronics_psd_{i}'][0] for i in range(len(frequencies))],
        'fixed conditional detector from parent ledger; no measured calibration claimed')


class SpatialUncertaintyEvaluator:
    """Rebuild actual geometry/atomic state for each point, cache within one run.

    An injected solver is for independent adapter tests. No results or cache
    entries are reused across source snapshots, grids or ODE tolerances.
    """

    def __init__(self, inputs, detector, *, probe_angle_rad=.005, aperture_m=(.0004, .0003),
                 mean_order=4, response_order=3, convention='uniform-zeeman-rms', solver=None):
        self.inputs, self.detector = inputs, detector
        self.probe_angle_rad, self.aperture_m = float(probe_angle_rad), tuple(aperture_m)
        self.mean_order, self.response_order, self.convention = mean_order, response_order, convention
        self.solver = stationary_spatial_cell if solver is None else solver
        self._cache = {}
        self.evaluations = []

    def medium(self, values, order_x):
        values = readonly_array(values, real=True)
        if values.shape != (2,):
            raise ValueError('ordered pump power and conjugate angle required')
        inputs = replace(self.inputs, pump_power_W=float(values[0]))
        geometry = CarrierGeometry.vacuum_beams(inputs, probe_angle_rad=self.probe_angle_rad,
                                               conjugate_angle_rad=float(values[1]))
        modes = TransverseModeGrid.rectangle(*self.aperture_m, order_x=order_x)
        return StationarySpatialMedium(inputs, geometry, modes, mean_order=self.mean_order,
            response_order=self.response_order, convention=self.convention)

    def evaluate(self, values, *, order_x, propagation_rtol=2e-9):
        # Validate even cached calls: a cache must not bypass input guards.
        medium = self.medium(values, order_x)
        tolerance = float(propagation_rtol)
        if not np.isfinite(tolerance) or tolerance <= 0:
            raise ValueError('positive finite propagation tolerance required')
        key = (int(order_x), tolerance, float(values[0]), float(values[1]))
        if key in self._cache:
            return copy.deepcopy(self._cache[key])
        started = time.monotonic()
        result = self.solver(medium, self.detector, propagation_rtol=tolerance)
        spectrum = result['spectrum']
        trace = SpectralTrace.from_intensity_difference(spectrum, contribution='quantum')
        if not np.array_equal(trace.analysis_axis.omega_rad_s, self.detector.analysis_axis.omega_rad_s):
            raise ValueError('solver returned a different RF axis')
        band = analyze_squeezing_bands(trace)
        transfer_audit = result['transfer'].audit()
        eigenvalues = np.array([result['minimum_channel_cp_eigenvalue'],
                               spectrum.minimum_covariance_uncertainty_eigenvalue], float)
        outputs = readonly_array(np.r_[result['probe_power_gain'], result['conjugate_power_gain'], trace.ratio], real=True)
        if (not transfer_audit['passed'] or not np.isfinite(eigenvalues).all()
                or eigenvalues.min() < -TOLERANCES['quantum_eigenvalue_absolute']
                or np.any(outputs < 0) or not band.complete_domain):
            raise ValueError('perturbed spatial solution failed physical/current-spectrum checks')
        band_record = _jsonable(asdict(band))
        for row, item in zip(band_record['bands'], band.bands):
            row.update(width_hz=item.width_hz, observed_span_hz=item.observed_span_hz)
        row = {'order_x': int(order_x), 'propagation_rtol': tolerance,
            'varied_inputs': dict(zip(INPUT_IDS, map(float, values))), 'outputs': outputs.tolist(),
            'wavevectors_rad_m': medium.geometry.wavevectors_rad_m.tolist(),
            'loop_wavevector_rad_m': (-medium.geometry.mismatch_rad_m).tolist(),
            'transfer_audit': transfer_audit, 'minimum_channel_cp_eigenvalue': float(eigenvalues[0]),
            'minimum_covariance_uncertainty_eigenvalue': float(eigenvalues[1]),
            'quantum_propagation': result['quantum_propagation'],
            'mean_ode_evaluations': result['mean_ode_evaluations'],
            'quantum_psd_A2_Hz': spectrum.quantum_psd_A2_Hz.tolist(),
            'sql_psd_A2_Hz': spectrum.sql_psd_A2_Hz.tolist(),
            'band_analysis': band_record, 'elapsed_seconds': time.monotonic()-started, 'passed': True}
        self._cache[key] = copy.deepcopy(row)
        self.evaluations.append(copy.deepcopy(row))
        print(f'Nx={order_x}, P={values[0]:.7g} W, theta_c={values[1]:.8g} rad: '
              f'{row["elapsed_seconds"]:.2f}s, Gp={outputs[0]:.10f}, R(1MHz)={outputs[3]:.10f}', flush=True)
        return row

    def callback(self, order_x, propagation_rtol=2e-9):
        return lambda values: np.array(self.evaluate(values, order_x=order_x,
                                         propagation_rtol=propagation_rtol)['outputs'])


def uncertainty_model(values=(.6, -.004), *, correlation=0.):
    parameters = tuple(ParameterEvidence(key, value, unit, 'assumed', 'declared spatial sensitivity fixture',
        sigma, 'user-independent diagnostic assumption; not apparatus measurements', (), (),
        'stationary fixed-pump two-profile calculation')
        for key, value, unit, sigma in zip(INPUT_IDS, values, ('W', 'rad'), SIGMA))
    evidence = CorrelationEvidence('assumed', 'explicit correlation sensitivity control',
        'predeclared correlation, never fitted to gain or squeezing', (), 'same two assumed input uncertainties')
    return JointInputUncertainty(parameters, [[1., correlation], [correlation, 1.]], evidence)


def run_linear_case(evaluator, model, *, order_x, steps, propagation_rtol=2e-9):
    return propagate_first_order(model, evaluator.callback(order_x, propagation_rtol),
        output_ids=OUTPUT_IDS, output_units=('1',)*len(OUTPUT_IDS), steps=steps)


def linear_summary(result):
    sigma = result.standard_uncertainties
    scale = sigma[:, None]*sigma[None, :]
    correlation = np.divide(result.covariance, scale, out=np.zeros_like(result.covariance), where=scale > 0)
    return {'nominal': result.nominal.tolist(), 'jacobian_per_input_unit': result.jacobian.tolist(),
        'steps_in_input_units': result.steps.tolist() if result.steps is not None else None,
        'weighted_response': (result.jacobian*result.model.standard_uncertainties).tolist(),
        'covariance': result.covariance.tolist(), 'standard_uncertainties': sigma.tolist(),
        'output_correlation': correlation.tolist(), 'zero_variance_output_indices': np.flatnonzero(sigma == 0).tolist()}


def compare_linear(coarse, refined):
    weighted_change = abs((coarse.jacobian-refined.jacobian)*refined.model.standard_uncertainties)
    weighted_reference = abs(refined.jacobian*refined.model.standard_uncertainties)
    sigma_change = abs(coarse.standard_uncertainties-refined.standard_uncertainties)
    weighted_budget = TOLERANCES['weighted_response_absolute']+TOLERANCES['weighted_response_relative']*weighted_reference
    sigma_budget = TOLERANCES['standard_uncertainty_absolute']+TOLERANCES['standard_uncertainty_relative']*refined.standard_uncertainties
    return {'weighted_response_absolute_change': weighted_change.tolist(),
        'standard_uncertainty_absolute_change': sigma_change.tolist(),
        'maximum_weighted_response_fraction_of_tolerance': float(np.max(weighted_change/weighted_budget)),
        'maximum_standard_uncertainty_fraction_of_tolerance': float(np.max(sigma_change/sigma_budget)),
        'covariance_relative_change': float(np.linalg.norm(coarse.covariance-refined.covariance)/max(np.linalg.norm(refined.covariance), np.finfo(float).tiny)),
        'passed': bool(np.all(weighted_change <= weighted_budget) and np.all(sigma_change <= sigma_budget))}


def nominal_parity(current, parent):
    gains = abs(current.nominal[:2]-[parent['probe_power_gain'], parent['conjugate_power_gain']])
    db = abs(10*np.log10(current.nominal[2:])-parent['spectrum_dB'])
    return {'maximum_gain_absolute_change': float(max(gains)),
        'maximum_spectrum_dB_absolute_change': float(max(db)),
        'passed': bool(max(gains) < TOLERANCES['nominal_gain_absolute']
                       and max(db) < TOLERANCES['nominal_spectrum_dB_absolute'])}


def build_report(parent_path=PARENT):
    before = source_manifest()
    parent, provenance = load_parent(parent_path)
    detector = default_detector(parent)
    inputs = ReducedPowerInputs(**parent['inputs'])
    evaluator = SpatialUncertaintyEvaluator(inputs, detector, aperture_m=parent['aperture_m'],
        mean_order=parent['mean_order'], response_order=parent['response_order'])
    model = uncertainty_model((inputs.pump_power_W, -.004))
    with core.blas_single_thread():
        coarse = run_linear_case(evaluator, model, order_x=8, steps=COARSE_STEPS)
        fine = run_linear_case(evaluator, model, order_x=8, steps=COARSE_STEPS/2)
        spatial = run_linear_case(evaluator, model, order_x=12, steps=COARSE_STEPS/2)
    step_change, quadrature_change = compare_linear(coarse, fine), compare_linear(fine, spatial)
    correlated = {}
    for rho in (-.6, 0., .6):
        changed = uncertainty_model(model.values, correlation=rho)
        result = propagate_first_order(changed, lambda x: spatial.nominal, jacobian=spatial.jacobian,
            output_ids=OUTPUT_IDS, output_units=('1',)*len(OUTPUT_IDS))
        correlated[str(rho)] = linear_summary(result)
    # Cross terms are signed: they are not an independent positive noise source.
    weighted = spatial.jacobian*model.standard_uncertainties
    contributions = {key: np.outer(weighted[:, i], weighted[:, i]).tolist() for i, key in enumerate(INPUT_IDS)}
    cross = .6*(np.outer(weighted[:, 0], weighted[:, 1])+np.outer(weighted[:, 1], weighted[:, 0]))
    parent_rows = {row['name']: row for row in parent['cells']}
    parity = {'Nx8': nominal_parity(fine, parent_rows['unequal_Nx8']),
              'Nx12': nominal_parity(spatial, parent_rows['unequal_Nx12'])}
    ledger = inputs.consumed_inputs(detector, convention='uniform-zeeman-rms')
    ledger['scalars'].update({'probe_angle_rad': (.005, 'rad'), 'conjugate_angle_rad': (-.004, 'rad'),
        'aperture_width_m': (evaluator.aperture_m[0], 'm'), 'aperture_height_m': (evaluator.aperture_m[1], 'm'),
        **{f'velocity_{axis}_m_s': (0., 'm/s') for axis in 'xyz'}})
    scope_audit = model.audit_evidence(required_ids=ledger['scalars'])
    stable = before == source_manifest()
    report = {'schema': 'gabes.stationary_spatial_uncertainty.v1',
        'scope': 'conditional joint five-output uncertainty, stationary centers, fixed pump and two fixed top-hat profiles',
        'absolute_hot_vapor_prediction': False, 'experimental_validation': False,
        'parent_artifact': provenance, 'source_sha256': before, 'source_stable_during_run': stable,
        'input_evidence': [asdict(item) for item in model.parameters],
        'correlation_evidence': asdict(model.correlation_evidence),
        'input_standard_uncertainties': SIGMA.tolist(), 'base_input_correlation': model.correlation.tolist(),
        'declared_input_evidence_audit': asdict(model.audit_evidence()),
        'all_consumed_input_evidence_audit': asdict(scope_audit),
        'consumed_input_ledger': _jsonable(ledger), 'output_ids': list(OUTPUT_IDS), 'output_units': ['1']*5,
        'rf_hz': detector.analysis_axis.frequency_hz.tolist(),
        'numerics': {'quadrature_orders_x': [8, 12], 'quadrature_order_y': 1,
            'mean_order': evaluator.mean_order, 'response_order': evaluator.response_order,
            'mean_ode_rtol': 2e-10, 'mean_ode_atol': 2e-12, 'quantum_propagation_rtol': 2e-9,
            'coarse_steps': COARSE_STEPS.tolist(), 'fine_steps': (COARSE_STEPS/2).tolist()},
        'tolerances': TOLERANCES, 'Nx8_coarse_steps': linear_summary(coarse),
        'Nx8_fine_steps': linear_summary(fine), 'Nx12_fine_steps': linear_summary(spatial),
        'step_refinement_at_Nx8': step_change, 'quadrature_refinement_at_fine_steps': quadrature_change,
        'nominal_parent_comparisons': parity, 'correlation_cases_using_same_J': correlated,
        'uncorrelated_covariance_contributions': contributions, 'cross_covariance_at_rho_plus_0_6': cross.tolist(),
        'current_evaluations': evaluator.evaluations, 'unique_spatial_cell_evaluations': len(evaluator.evaluations),
        'historical_ode_refinement': {'scope': 'parent nominal-only evidence; perturbed-point ODE sensitivity not re-refined here',
            'parent_quantum_ODE_change_at_Nx8': parent['quantum_ODE_change_at_Nx8'],
            'parent_mean_ODE_refinement_relative_trajectory_change': parent['mean_ODE_refinement_relative_trajectory_change']},
        'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
            'blas_threads': 1, 'timing_scope': 'elapsed audit time; not a benchmark'},
        'limitations': [
            'The two standard uncertainties and all tested correlations are assumed; no measured input or no-fit claim.',
            'Other consumed inputs, fixed top-hat profiles and detector calibration are conditioned upon; zero uncertainty is not assigned to their unknown errors.',
            'Central-difference and spatial-grid comparisons are deterministic errors, not added independent statistical covariance.',
            'Linear covariance is local; no Gaussian ensemble, nonlinear coverage interval, or large-input-excursion bound was computed.',
            'Three RF samples do not resolve intrinsic bandwidth or an RF minimum; censored bands retain width=None.',
            'Atomic harmonic and ODE refinement are inherited nominal evidence only, not jointly refined at every perturbed input.',
            'Moving atomic transport, inter-slice noise, walkoff, diffraction, extra collection modes and pump depletion remain outside this model.',
        ]}
    report['expected_controls_passed'] = bool(stable and step_change['passed'] and quadrature_change['passed']
        and all(row['passed'] for row in parity.values()) and all(row['passed'] for row in evaluator.evaluations)
        and not scope_audit.passed and all(band['width_hz'] is None for row in evaluator.evaluations
                                        for band in row['band_analysis']['bands']))
    return report


def save_plot(report, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    selected = report['Nx12_fine_steps']
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), layout='constrained')
    freq = np.array(report['rf_hz'])/1e6
    axes[0].errorbar(freq, selected['nominal'][2:], yerr=selected['standard_uncertainties'][2:],
        fmt='o', capsize=5, label='Assumed independent inputs: one standard uncertainty')
    axes[0].set(xlabel='RF frequency (MHz)', ylabel='Linear intensity noise / SQL',
        title='Stationary spatial field: three RF samples')
    axes[0].legend(fontsize=8)
    positions = np.arange(5)
    reference = np.array(selected['standard_uncertainties'])
    for i, (rho, row) in enumerate(report['correlation_cases_using_same_J'].items()):
        axes[1].bar(positions+(i-1)*.23, np.array(row['standard_uncertainties'])/reference,
                    width=.23, label=f'Input correlation {rho}')
    axes[1].set_xticks(positions, ['Gp', 'Gc', 'R(0.1)', 'R(1)', 'R(4)'])
    axes[1].set(ylabel='Standard uncertainty / independent-input value', title='Correlation effect at unchanged nominal inputs')
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2, axis='y')
    fig.suptitle('Conditional pump-power and conjugate-angle sensitivity | No experimental validation', fontsize=12)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as stream:
        fig.savefig(stream, format='png', dpi=160)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent', type=Path, default=PARENT)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--plot', type=Path)
    args = parser.parse_args(argv)
    if args.plot is not None and args.output.resolve() == args.plot.resolve():
        raise ValueError('report and plot require different paths')
    for path in (args.output, args.plot):
        if path is not None and path.exists():
            raise FileExistsError(path)
    report = build_report(args.parent)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    if args.plot is not None:
        save_plot(report, args.plot)
    print(f'Wrote {args.output}; controls={report["expected_controls_passed"]}', flush=True)
    return 0 if report['expected_controls_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
