"""Immutable independent references for six new thermal Rb boundary paths.

This pilot uses the existing continuous Gaussian Rb model, with an explicitly
assumed square open column. It does not certify thermal ensemble convergence.
Run with --plan-only to validate sources and print identities without solving.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import platform

import numpy as np
import scipy
import threadpoolctl

from gabes import core
from gabes.constants import MASS_85RB
from gabes.fwm_quantum.smooth_transport import smooth_rb_problem
from gabes.quantum.inflow import maxwell_box_inflow
from . import adjoint_transport_audit as audit
from .reference.adjoint_transport import adjoint_wavepacket
from .smooth_transport_audit import ROOT, fixture


DRIVER = 'analysis/grand_challenge/rb_thermal_reference_jobs.py'
DEFAULT_OUTPUT = ROOT/'docs/grand_challenge/rb_thermal_reference_jobs_v1'
TEMPERATURE_K = 373.
POWER = 0
SEED = 11
WORKERS = 2
REFERENCE_TOLERANCES = ((1e-9, 1e-12), (3e-10, 3e-13))
METRICS = ('greater', 'lesser', 'greater_by_source', 'lesser_by_source',
           'mean_pulse', 'mean_outer', 'retarded_response', 'exit_state')
INFLOW_SOURCE = (
    'Conditional reduced-Rb open column: square transverse side sqrt(uniform_area_m2), '
    'length inputs.length_m, centered at origin; independently assumed 373 K; '
    'existing MASS_85RB and inputs.number_density_m3; no measured cell geometry')
SCOPE = {
    'kind': 'independent backward-observable references for six new boundary paths',
    'parent_use': 'physical fixture and source validation only; no selected-path solution reused',
    'geometry': 'assumed square open column; optical normalization area is not measured wall area',
    'boundary': 'fresh unpolarized inflow even where Gaussian pump remains nonzero; exterior history omitted',
    'pump': 'continuous constant-waist transverse Gaussian; no seed drive',
    'phase': 'one common laboratory beat; entry phase zero with geometric optical entry offsets',
    'density': 'boundary rates include density once; single-atom reference values are unweighted',
    'response': 'formal declared-drive atomic response; no spatial Maxwell response kernel',
    'trajectory': 'true first-exit chords; no residence clipping or geometry shrinkage',
    'thermal_ensemble_certified': False,
    'production_adapter_certified': False,
    'experimental_prediction': False,
}


def encode(value):
    """JSON encoding preserving all complex components and scalar precision."""
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {'real': value.real.tolist(), 'imag': value.imag.tolist()}
        return value.tolist()
    if isinstance(value, np.generic):
        return encode(value.item())
    if isinstance(value, complex):
        return {'real': value.real, 'imag': value.imag}
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def canonical(value):
    return json.dumps(encode(value), sort_keys=True, separators=(',', ':'),
                      ensure_ascii=False, allow_nan=False).encode('utf-8')


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def source_snapshot():
    """Verify the selected-path parent, then bind old physics plus this driver."""
    sources = audit.hashes()
    sources[DRIVER] = audit.file_hash(ROOT/DRIVER)
    return {'source_sha256': sources,
            'selected_path_parent': str(audit.PARENT.relative_to(ROOT)).replace('\\', '/'),
            'selected_path_parent_sha256': audit.file_hash(audit.PARENT)}


def require_sources(expected):
    current = source_snapshot()
    if current != expected:
        raise ValueError('source or selected-path parent changed; no result may be saved')
    return current


def environment():
    return {'python': platform.python_version(), 'numpy': np.__version__,
            'scipy': scipy.__version__, 'platform': platform.platform(),
            'threadpoolctl': threadpoolctl.__version__}


def physical_fixture():
    """Construct new boundary paths; the old selected path is discarded."""
    inputs, geometry, _selected_path, axis = fixture()
    side = np.sqrt(inputs.uniform_area_m2)
    lengths = np.array([side, side, inputs.length_m])
    inflow = maxwell_box_inflow(-lengths/2, lengths/2,
        temperature_K=TEMPERATURE_K, mass_kg=MASS_85RB,
        density_m3=inputs.number_density_m3, points_per_face_power=POWER,
        seed=SEED, source=INFLOW_SOURCE)
    return inputs, geometry, axis, inflow


def job_identity(path_index, reference_level, snapshot):
    """Return a fully declared job and solver kwargs, without running an ODE."""
    if type(path_index) is not int or not 0 <= path_index < 6:
        raise ValueError('path_index must select one of the six pilot paths')
    if type(reference_level) is not int or reference_level not in (0, 1):
        raise ValueError('reference_level must be 0 or 1')
    inputs, geometry, axis, inflow = physical_fixture()
    path = inflow.path(path_index)
    boundary = np.diag([5/12, 7/12, 0., 0.]).astype(complex)
    problem = smooth_rb_problem(inputs, geometry, path, axis,
        boundary_state=boundary, entry_phase_rad=0., pump_center_xy_m=(0., 0.))
    metadata = problem.pop('metadata')
    rtol, atol = REFERENCE_TOLERANCES[reference_level]
    settings = {'rtol': rtol, 'atol': atol, 'max_step_s': path.residence_time_s/64}
    reservoirs = problem['reservoirs']
    names = ('atomic_inflow',)+tuple('jump:'+channel.name for channel in reservoirs.channels)
    identity = encode({
        'schema': 'gabes-rb-thermal-reference-job-identity-v1',
        'path_index': path_index, 'reference_level': reference_level,
        'job_index': 2*path_index+reference_level,
        'inputs_SI': asdict(inputs),
        'box': {'lower_corner_m': inflow.lower_corner_m, 'upper_corner_m': inflow.upper_corner_m,
                'temperature_K': inflow.temperature_K, 'mass_kg': inflow.mass_kg,
                'density_m3': inflow.density_m3, 'source': inflow.source,
                'points_per_face_power': POWER, 'seed': SEED,
                'total_arrival_rate_s_inverse': inflow.total_arrival_rate_s_inverse,
                'equilibrium_atom_number': inflow.equilibrium_atom_number,
                'raw_mean_occupancy': inflow.mean_occupancy},
        'path': {'entry_position_m': path.entry_position_m,
                 'exit_position_m': path.entry_position_m+path.velocity_m_s*path.residence_time_s,
                 'velocity_m_s': path.velocity_m_s, 'residence_time_s': path.residence_time_s,
                 'face_index': int(inflow.face_index[path_index]),
                 'rate_s_inverse': float(inflow.rate_s_inverse[path_index]), 'source': path.source},
        'analysis_frequencies_hz': axis.frequency_hz,
        'wavevectors_rad_m': geometry.wavevectors_rad_m,
        'port_wavevectors_rad_m': np.concatenate((geometry.wavevectors_rad_m[1:]-geometry.wavevectors_rad_m[0],
                                                 geometry.wavevectors_rad_m[0]-geometry.wavevectors_rad_m[1:])),
        'port_frequencies_rad_s': problem['frequencies_rad_s'],
        'drive_frequencies_rad_s': problem['drive_frequencies_rad_s'],
        'physical_metadata': metadata,
        'atomic_problem': {'h0': problem['h0'], 'h1': problem['h1'],
            'boundary_state': problem['boundary_state'], 'readouts': problem['readouts'],
            'drives': problem['drives'], 'n_levels': reservoirs.n_levels,
            'envelope_formula': 'exp(-((entry_x+vx*age)^2+(entry_y+vy*age)^2)/pump_waist_m^2)',
            'pump_center_xy_m': [0., 0.],
            'channels': [asdict(channel) for channel in reservoirs.channels]},
        'source_names': names, 'numerical_settings': settings,
        'environment': environment(), 'source_snapshot': snapshot, 'scope': SCOPE,
    })
    return identity, dict(problem, **settings)


def job_filename(path_index, reference_level):
    return f'path{path_index:02d}_reference{reference_level}.json'


def decode_array(value):
    if not isinstance(value, dict) or set(value) != {'real', 'imag'}:
        raise ValueError('all raw metrics require explicit real and imaginary components')
    real, imag = np.asarray(value['real'], dtype=float), np.asarray(value['imag'], dtype=float)
    if real.shape != imag.shape or not np.all(np.isfinite(real)) or not np.all(np.isfinite(imag)):
        raise ValueError('raw metric components must be finite and shape-matched')
    return real+1j*imag


def validate_values(values, identity):
    if set(values) != set(METRICS):
        raise ValueError('all eight raw metrics required, with no substitutions')
    nf, ports = np.asarray(identity['port_frequencies_rad_s']).shape
    n = identity['atomic_problem']['n_levels']
    ns = len(identity['source_names'])
    shapes = {key: (nf, ports, ports) for key in METRICS}
    shapes.update(greater_by_source=(ns, nf, ports, ports),
                  lesser_by_source=(ns, nf, ports, ports),
                  mean_pulse=(nf, ports), exit_state=(n, n))
    arrays = {key: decode_array(values[key]) for key in METRICS}
    if any(arrays[key].shape != shape for key, shape in shapes.items()):
        raise ValueError('raw metric shape differs from the declared physical problem')
    # These are storage/algebra checks, not convergence evidence.
    for total, source in (('greater', 'greater_by_source'), ('lesser', 'lesser_by_source')):
        if not np.allclose(arrays[total], arrays[source].sum(axis=0), rtol=1e-12, atol=0):
            raise ValueError('stored source sum does not close')
    mean = arrays['mean_pulse']
    if not np.allclose(arrays['mean_outer'], mean[:, :, None]*mean[:, None, :].conj(),
                       rtol=1e-12, atol=0):
        raise ValueError('stored mean_outer differs from the complex mean pulse product')


def sealed(payload):
    record = encode(payload)
    if 'record_sha256' in record:
        raise ValueError('record hash must be computed from the complete unhashed payload')
    return dict(record, record_sha256=digest(record))


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError('duplicate JSON key in immutable record')
        result[key] = value
    return result


def read_sealed(path):
    record = json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=_unique_object)
    if not isinstance(record, dict) or 'record_sha256' not in record:
        raise ValueError('immutable record has no content hash')
    supplied = record['record_sha256']
    payload = {key: value for key, value in record.items() if key != 'record_sha256'}
    if supplied != digest(payload):
        raise ValueError('immutable record content hash mismatch')
    return record


def validate_job(path, identity, snapshot):
    require_sources(snapshot)
    record = read_sealed(path)
    return validate_record(record, identity, snapshot)


def validate_record(record, identity, snapshot):
    """Validate complete numerical storage before publishing or resuming it."""
    if record.get('record_sha256') != digest({
            key: value for key, value in record.items() if key != 'record_sha256'}):
        raise ValueError('reference record content hash mismatch')
    if (record.get('schema') != 'gabes-rb-thermal-reference-job-v1'
            or record.get('identity') != identity
            or record.get('full_identity_sha256') != digest(identity)
            or record.get('source_before') != snapshot or record.get('source_after') != snapshot
            or record.get('source_stable_during_run') is not True
            or record.get('source_names') != identity['source_names']
            or record.get('environment') != identity['environment']
            or record.get('scope') != SCOPE):
        raise ValueError('job identity/provenance mismatch')
    numerics = record.get('numerics', {})
    pools = record.get('runtime_blas_threads', [])
    if not pools or any(row.get('user_api') != 'blas' or row.get('num_threads') != 1 for row in pools):
        raise ValueError('Actual single-thread BLAS runtime evidence required')
    if any(numerics.get(key) != value for key, value in identity['numerical_settings'].items()):
        raise ValueError('stored numerical settings differ from the required reference level')
    for name in ('forward_density_evaluations', 'backward_evaluations',
                 'density_mesh_points', 'backward_complex_variables'):
        value = numerics.get(name)
        if type(value) is not int or value <= 0:
            raise ValueError('positive solver work counters required')
    elapsed = numerics.get('elapsed_seconds')
    if isinstance(elapsed, bool) or not isinstance(elapsed, (int, float)) or not np.isfinite(elapsed) or elapsed < 0:
        raise ValueError('finite nonnegative elapsed time required')
    if (record.get('frequencies_rad_s') != identity['port_frequencies_rad_s']
            or record.get('residence_time_s') != identity['path']['residence_time_s']):
        raise ValueError('stored frequency or residence differs from the new boundary path')
    validate_values(record['values'], identity)
    return record


def write_exclusive(path, record, snapshot):
    text = json.dumps(record, indent=2, ensure_ascii=False, allow_nan=False)+'\n'
    require_sources(snapshot)
    with Path(path).open('x', encoding='utf-8', newline='\n') as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())


def run_job(path_index, reference_level, output_dir, snapshot, expected_identity_hash):
    """A fresh independent solve, or a fully verified immutable resume."""
    before = require_sources(snapshot)
    identity, kwargs = job_identity(path_index, reference_level, snapshot)
    if digest(identity) != expected_identity_hash:
        raise ValueError('worker identity differs from the frozen plan')
    path = Path(output_dir)/job_filename(path_index, reference_level)
    if path.exists():
        record = validate_job(path, identity, snapshot)
        return {'file': path.name, 'record_sha256': record['record_sha256'], 'resumed': True}
    if core._threadpool_limits is None:
        raise RuntimeError('threadpoolctl is required to enforce one BLAS thread per reference job')
    with core.blas_single_thread():
        pools = [row for row in threadpoolctl.threadpool_info() if row['user_api'] == 'blas']
        if not pools or any(row['num_threads'] != 1 for row in pools):
            raise RuntimeError('Could not verify one thread for every active BLAS runtime')
        result = adjoint_wavepacket(**kwargs)
    after = require_sources(snapshot)
    values = encode({key: result[key] for key in METRICS})
    validate_values(values, identity)
    record = sealed({'schema': 'gabes-rb-thermal-reference-job-v1',
        'identity': identity, 'full_identity_sha256': digest(identity),
        'values': values, 'source_names': result['source_names'],
        'numerics': result['numerics'], 'frequencies_rad_s': result['frequencies_rad_s'],
        'residence_time_s': result['residence_time_s'], 'environment': environment(),
        'source_before': before, 'source_after': after,
        'source_stable_during_run': before == after, 'scope': SCOPE})
    record = sealed({**{key: value for key, value in record.items() if key != 'record_sha256'},
                     'runtime_blas_threads': pools})
    validate_record(record, identity, snapshot)
    write_exclusive(path, record, snapshot)
    validate_job(path, identity, snapshot)
    return {'file': path.name, 'record_sha256': record['record_sha256'], 'resumed': False,
            'elapsed_seconds': result['numerics']['elapsed_seconds']}


def prepare_plan():
    snapshot = source_snapshot()
    jobs = [job_identity(i, level, snapshot)[0] for i in range(6) for level in range(2)]
    require_sources(snapshot)
    return snapshot, jobs


def run(output_dir=DEFAULT_OUTPUT):
    snapshot, jobs = prepare_plan()
    output_dir = Path(output_dir).resolve()
    filenames = [job_filename(row['path_index'], row['reference_level']) for row in jobs]
    index_path = output_dir/'index.json'
    if output_dir.exists():
        unknown = {path.name for path in output_dir.iterdir()}-set(filenames)-{'index.json'}
        if unknown:
            raise ValueError(f'unexpected files in immutable output directory: {sorted(unknown)}')
    else:
        output_dir.mkdir(parents=True)
    pending = []
    for identity, filename in zip(jobs, filenames):
        if (output_dir/filename).exists():
            validate_job(output_dir/filename, identity, snapshot)
            print(f'verified resume: {filename}', flush=True)
        else:
            pending.append(identity)
    if index_path.exists() and pending:
        raise ValueError('existing completion index has missing job files')
    if index_path.exists():
        read_sealed(index_path)
    if pending:
        with ProcessPoolExecutor(max_workers=WORKERS) as pool:
            futures = [pool.submit(run_job, row['path_index'], row['reference_level'],
                str(output_dir), snapshot, digest(row)) for row in pending]
            for count, future in enumerate(as_completed(futures), start=1):
                result = future.result()
                print(json.dumps({'completed_this_run': count, 'scheduled_this_run': len(pending),
                                  **result}, allow_nan=False), flush=True)
    records = []
    for identity, filename in zip(jobs, filenames):
        path = output_dir/filename
        record = validate_job(path, identity, snapshot)
        records.append({'path_index': identity['path_index'], 'reference_level': identity['reference_level'],
            'file': filename, 'file_sha256': audit.file_hash(path),
            'record_sha256': record['record_sha256'], 'full_identity_sha256': record['full_identity_sha256']})
    after = require_sources(snapshot)
    index = sealed({'schema': 'gabes-rb-thermal-reference-job-index-v1',
        'path_count': 6, 'job_count': 12, 'all_jobs_complete': True,
        'workers': WORKERS, 'blas_threads_per_job': 1,
        'reference_tolerances': REFERENCE_TOLERANCES, 'default_max_step': 'true residence time / 64',
        'inputs_SI': jobs[0]['inputs_SI'], 'box': jobs[0]['box'],
        'analysis_frequencies_hz': jobs[0]['analysis_frequencies_hz'],
        'wavevectors_rad_m': jobs[0]['wavevectors_rad_m'],
        'paths': [row['path'] for row in jobs[::2]], 'jobs': records,
        'environment': environment(), 'source_before': snapshot, 'source_after': after,
        'source_stable_during_run': snapshot == after, 'scope': SCOPE})
    if index_path.exists():
        if read_sealed(index_path) != index:
            raise ValueError('existing completion index differs from the verified job manifest')
    else:
        write_exclusive(index_path, index, snapshot)
    print(json.dumps({'index': str(index_path), 'index_file_sha256': audit.file_hash(index_path),
                      'job_count': 12, 'thermal_ensemble_certified': False}), flush=True)
    return index


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--plan-only', action='store_true', help='print frozen job identities; no files or ODE solves')
    args = parser.parse_args()
    if args.plan_only:
        snapshot, jobs = prepare_plan()
        print(json.dumps({'source_snapshot': snapshot, 'workers': WORKERS,
            'jobs': [{'file': job_filename(row['path_index'], row['reference_level']),
                      'full_identity_sha256': digest(row), 'identity': row} for row in jobs]},
            indent=2, ensure_ascii=False, allow_nan=False))
    else:
        run(args.output_dir)


if __name__ == '__main__':
    main()
