"""Extend an immutable kinetic audit after its longitudinal convergence failure.

python -m analysis.grand_challenge.kinetic_refinement_audit --parent OLD.json --output NEW.json --plot NEW.png
"""

import argparse
import copy
import hashlib
import json
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np

from gabes import core
from gabes.fwm_quantum.inputs import ReducedPowerInputs
from gabes.fwm_quantum.kinetic import CarrierGeometry, VelocityQuadrature, thermal_pump_cell
from gabes.quantum.contracts import AnalysisFrequencyAxis
from gabes.quantum.readout import DetectorResponse
from .kinetic_audit import ROOT, TOLERANCES, _case, _comparison, _complex, save_plot


def _decode(value):
    return np.asarray(value['real'])+1j*np.asarray(value['imag'])


def build_report(parent_path):
    """Reuse only verified parent evidence; recalculate the failed axis and tail.

    Parent transverse and 6-to-7-sigma controls retain their exact old grids.
    Their provenance is explicit; they are not described as new-grid controls.
    """
    raw = parent_path.read_bytes()
    parent = json.loads(raw)
    if ('parent_artifact' in parent or
            parent['stage'] != 'S1 kinetic local mean/noise and conditional thermal pump-state cell'):
        raise ValueError('initial kinetic audit required')
    if parent['numerical_tolerances'] != TOLERANCES or not parent['local_controls']['passed']:
        raise ValueError('unchanged tolerances and passing local reference controls required')
    for name, expected in parent['source_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != expected:
            raise ValueError(f'parent source changed: {name}; regenerate the initial audit')
    report = copy.deepcopy(parent)
    inputs = ReducedPowerInputs(**{key: row['value'] for key, row in parent['base_input_ledger'].items()})
    geometry = CarrierGeometry(parent['synthetic_closed_geometry']['wavevectors_rad_m'], parent['synthetic_closed_geometry']['source'])
    detector = DetectorResponse(AnalysisFrequencyAxis.from_hz(parent['rf_hz']), parent['detector']['transmissions'],
        _decode(parent['detector']['current_response']), parent['detector']['balance'], np.zeros(len(parent['rf_hz'])), parent['detector']['source'])
    # The initial audit explicitly fixes electronics=0; no electronics fit is
    # imported from an output spectrum during refinement.
    stored = parent['selected_generator']
    gen = SimpleNamespace(drift=_decode(stored['drift_m_inverse']),
        noise_greater=_decode(stored['noise_greater_by_reservoir_m_inverse']).sum(axis=0),
        noise_lesser=_decode(stored['noise_lesser_by_reservoir_m_inverse']).sum(axis=0))
    selected_parent = parent['velocity_grid_cases'][parent['selected_case']]
    previous = {'local': {'generator': gen}, 'probe_power_gain': selected_parent['probe_power_gain'],
        'spectrum': SimpleNamespace(quantum_db=np.asarray(selected_parent['detected_S_minus_db']))}
    new_indices, new_results = [], []
    with core.blas_single_thread():
        for cut in (7, 8):
            started = time.monotonic()
            q = VelocityQuadrature.maxwell_xz(parent['kinetic_input_ledger']['temperature_K'],
                longitudinal_order=1536, transverse_order=48, cutoff_sigma=cut)
            result = thermal_pump_cell(inputs, geometry, q, detector)
            new_indices.append(len(report['velocity_grid_cases']))
            new_results.append(result)
            report['velocity_grid_cases'].append({'grid': {'longitudinal_order': 1536, 'transverse_order': 48, 'cutoff_sigma': cut},
                'omitted_Maxwell_probability': q.omitted_probability, 'seconds': time.monotonic()-started, **_case(result)})
            print(f'Completed refinement grid 1536 x 48, cutoff {cut} sigma', flush=True)
    report['initial_refinements'] = parent['refinements']
    report['refinements'] = [
        {'name': 'longitudinal', 'from_case': parent['selected_case'], 'to_case': new_indices[0],
            'evidence': 'new 1024-to-1536 comparison at Nx=48, cutoff=7', **_comparison(previous, new_results[0])},
        {**parent['refinements'][1], 'evidence': 'inherited Nx=32-to-48 comparison at Nz=1024, cutoff=7'},
        {**parent['refinements'][2], 'evidence': 'inherited cutoff=6-to-7 comparison at Nz=1024, Nx=48'},
        {'name': 'tail_7_to_8_sigma', 'from_case': new_indices[0], 'to_case': new_indices[1],
            'evidence': 'new cutoff=7-to-8 comparison at Nz=1536, Nx=48', **_comparison(*new_results)}]
    report['selected_case'] = new_indices[0]
    gen = new_results[0]['local']['generator']
    report['selected_generator'] = {'signed_lab_rf_rad_s': gen.frequency_axis.omega_rad_s.tolist(),
        'mode_labels': list(gen.mode_labels), 'signs': gen.signs.tolist(), 'drift_m_inverse': _complex(gen.drift),
        'reservoir_names': list(gen.reservoir_names),
        'noise_greater_by_reservoir_m_inverse': _complex(gen.noise_greater_by_reservoir),
        'noise_lesser_by_reservoir_m_inverse': _complex(gen.noise_lesser_by_reservoir)}
    report['coarse_to_selected_comparisons'] = parent['coarse_to_selected_comparisons']
    report['coarse_comparison_reference_case'] = parent['selected_case']
    report['parent_artifact'] = {'path': str(parent_path.resolve().relative_to(ROOT)).replace('\\', '/'),
        'sha256': hashlib.sha256(raw).hexdigest(), 'source_hashes_verified': len(parent['source_sha256']),
        'expected_controls_passed': parent['expected_controls_passed']}
    report['source_sha256'][str(Path(__file__).resolve().relative_to(ROOT)).replace('\\', '/')] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    report['expected_controls_passed'] = bool(parent['local_controls']['passed']
        and parent['vacuum_opposite_angle_geometry']['rejected_by_single_phase_solver']
        and all(row['passed'] for row in report['refinements'])
        and all(row['local_audit']['passed'] and row['global_audit']['passed']
            and row['minimum_channel_cp_eigenvalue'] >= -1e-10
            and row['minimum_detected_covariance_uncertainty_eigenvalue'] >= -1e-10 for row in report['velocity_grid_cases']))
    report['limits'].append('Refinement preserves the failed initial audit; transverse and 6-to-7-sigma controls remain at their recorded parent grids. No joint-grid or rigorous tail error bound is claimed.')
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--plot', type=Path)
    args = parser.parse_args(argv)
    for path in (args.output, args.plot):
        if path is not None and path.exists():
            raise FileExistsError(path)
    report = build_report(args.parent)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print(f'Wrote {args.output}')
    if args.plot is not None:
        save_plot(report, args.plot)
        print(f'Wrote {args.plot}')
    return 0 if report['expected_controls_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
