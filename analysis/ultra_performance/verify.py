"""Revalidate the linked Ultra proposal and time the exact kernel replacement.

python -m analysis.ultra_performance.verify --output NEW.json
The existing pre-patch snapshots are immutable evidence, never new baselines.
"""

import argparse
import ast
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time

import numpy as np

from gabes import core, doppler, kernels, constants as c
from gabes.schemes import fwm


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
PROPOSAL = 'https://claude.ai/code/artifact/bf263a2e-0d04-4705-90b6-a9d14eca192a'


@contextmanager
def kernel_path(function):
    previous = kernels.floquet_chi_grid
    kernels.floquet_chi_grid = function
    try:
        yield
    finally:
        kernels.floquet_chi_grid = previous


def errors(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or not np.array_equal(np.isfinite(a), np.isfinite(b)):
        raise AssertionError('shape or finite-mask changed')
    finite = np.isfinite(a)
    diff = float(np.max(abs(a[finite]-b[finite]))) if finite.any() else 0.
    scale = float(np.max(abs(b[finite]))) if finite.any() else 0.
    return {'maximum_absolute': diff, 'maximum_relative_to_array_scale': diff/max(scale, 1e-300)}


def baseline_body_unchanged():
    old = subprocess.run(['git', 'show', 'c0c46f0:gabes/kernels.py'], cwd=ROOT,
        check=True, capture_output=True, encoding='utf-8').stdout
    new = (ROOT/'gabes/kernels.py').read_text(encoding='utf-8')
    a = next(n for n in ast.parse(old).body if isinstance(n, ast.FunctionDef) and n.name == 'floquet_chi_grid')
    b = next(n for n in ast.parse(new).body if isinstance(n, ast.FunctionDef) and n.name == '_floquet_chi_grid_general')
    a.name = b.name
    if ast.dump(a) != ast.dump(b):
        raise AssertionError('the general reference differs from the pre-patch kernel')
    return True


def timed(function):
    started = time.perf_counter()
    result = function()
    return time.perf_counter()-started, result


def array_comparisons(actual, expected):
    return {key: errors(actual[key], value) for key, value in expected.items()
            if isinstance(value, np.ndarray)}


def sampled_audit_counterexample():
    # A sampled comparison cannot certify an omitted scan row, irrespective of
    # how accurately either Floquet order is solved at the sampled rows.
    n = 401
    response = {key: np.full(n, 1e-12+2e-13j) for key in ('chi_ss', 'chi_sc', 'chi_cs', 'chi_cc')}
    high = np.broadcast_to(np.eye(2, dtype=complex), (n, 2, 2)).copy()
    low = high.copy(); low[205, 0, 0] *= 1.1*np.exp(1j*np.deg2rad(1.))
    sampled = np.arange(0, n, 10)
    def audit(ids):
        selected = {key: value[ids] for key, value in response.items()}
        return fwm.assess_floquet_scan_convergence(high_order=3, low_order=2,
            high_response=selected, low_response=selected, high_transfer=high[ids],
            low_transfer=low[ids], scan_axis=np.arange(n)[ids])['status']
    return {'sampled_status': audit(sampled), 'full_status': audit(np.arange(n)),
        'scope': 'logical counterexample to a full-scan certificate; not measured atom data'}


@contextmanager
def exact_legendre_nodes(order, cutoff):
    """Replace only the main FWM grid; leave shared absorption helpers intact."""
    previous = fwm.doppler

    class DiagnosticDoppler:
        velocities = None

        def __getattr__(self, name):
            return getattr(previous, name)

        def velocity_grid(self, temperature, mass=c.MASS_85RB, **kwargs):
            # compute_spectrum's first grid is its susceptibility grid. Later
            # diagnostic grids and other modules must retain their own rules.
            if self.velocities is not None:
                return previous.velocity_grid(temperature, mass=mass, **kwargs)
            self.velocities, weights = doppler.maxwell_legendre_grid(
                temperature, mass=mass, order=order, cutoff_sigma=cutoff)
            return self.velocities, weights

        def build_Delta_eff_axis(self, delta_min, delta_max, velocities, k_vec=c.K_VEC):
            if velocities is not self.velocities:
                return previous.build_Delta_eff_axis(delta_min, delta_max, velocities, k_vec)
            if delta_min != delta_max:
                raise ValueError('single-detuning diagnostic only')
            return np.sort(delta_min-k_vec*velocities)

        def interpolation_weights(self, axis, delta, velocities, k_vec=c.K_VEC):
            if velocities is not self.velocities:
                return previous.interpolation_weights(axis, delta, velocities, k_vec)
            ids = np.searchsorted(axis, delta-k_vec*velocities)
            np.testing.assert_array_equal(axis[ids], delta-k_vec*velocities)
            lo = np.minimum(ids, len(axis)-2)
            return lo, (ids-lo).astype(float)

    fwm.doppler = DiagnosticDoppler()
    try:
        yield
    finally:
        fwm.doppler = previous


def quadrature_checks():
    result = []
    for D, temperature, power in ((.9, 394.15, .6), (-2.1, 403.15, .6), (-.2, 373.15, .2)):
        center = fwm.branch_center_GHz(D, -1)
        # Include the requested operating point exactly; no readout interpolation.
        configs = []
        spectra = []
        for order, cutoff in ((128, 5), (256, 5), (256, 6)):
            with exact_legendre_nodes(order, cutoff):
                spec = fwm.compute_spectrum(D, T=temperature, P_pump=power, P_probe=8e-6,
                    coarse_points=81, fine_points=0, scan_min=center-.408, scan_max=center+.392,
                    phase_detail=fwm.PHASE_ULTRA, model_fidelity=fwm.FIDELITY_ULTRA)
            spectra.append(spec)
            configs.append({'order': order, 'cutoff_sigma': cutoff, 'G_s_at_minus8MHz': float(spec['G_s'][40]),
                'gain_referred_diagnostic_db_at_minus8MHz': float(spec['S_dB'][40])})
        comparisons = []
        for a, b in zip(spectra[:-1], spectra[1:]):
            comparisons.append({'gain': errors(a['G_s'], b['G_s']),
                'diagnostic_db': errors(a['S_dB'], b['S_dB']),
                'delta_db_at_minus8MHz': float(a['S_dB'][40]-b['S_dB'][40]),
                'full_scan_within_0p01_db': bool(np.max(abs(a['S_dB']-b['S_dB'])) < .01)})
        result.append({'D_GHz': D, 'temperature_K': temperature, 'pump_power_W': power,
            'configurations': configs, 'adjacent_comparisons': comparisons})
    return result


def build_report():
    if not kernels.available():
        raise RuntimeError('Numba is required for the timing comparison')
    import numba
    baseline_body_unchanged()
    parent = json.loads((HERE/'ultra_app_before.json').read_text(encoding='utf-8'))
    scheme = fwm.FWMScheme()
    params = parent['params']
    fast = kernels.floquet_chi_grid
    general = kernels._floquet_chi_grid_general
    times = {'general': [], 'exact_symmetry': []}
    with core.blas_single_thread():
        # Warm both compiled paths; compilation and cache loading are excluded.
        for function in (general, fast):
            with kernel_path(function):
                fwm.compute_spectrum(.9, coarse_points=5, fine_points=0, velocity_step=200., velocity_cutoff=1.)
        for sequence in (('general', 'exact_symmetry'), ('exact_symmetry', 'general')):
            for name in sequence:
                with kernel_path(general if name == 'general' else fast):
                    elapsed, spec = timed(lambda: scheme.compute(params))
                times[name].append(elapsed)
                if name == 'exact_symmetry':
                    selected = spec
                else:
                    reference = spec
                print(f'Actual app Ultra {name}: {elapsed:.4f} s', flush=True)
        full_errors = array_comparisons(selected, reference)
        with np.load(HERE/'ultra_app_before.npz') as frozen:
            frozen_errors = {key: errors(selected[key], frozen[key]) for key in frozen.files}
        stress = []
        for D, temperature, power, branch in ((.9, 394.15, .6, 1), (-1.5, 383.15, 2., -1),
            (-2.1, 403.15, .6, -1), (-.2, 373.15, .2, -1)):
            center = fwm.branch_center_GHz(D, branch)
            outputs = []
            for function in (general, fast):
                with kernel_path(function):
                    outputs.append(fwm.compute_spectrum(D, T=temperature, P_pump=power, P_probe=8e-6,
                        branch=branch, coarse_points=61, fine_points=0, velocity_step=10., velocity_cutoff=4.,
                        scan_min=center-.55, scan_max=center+.55, phase_detail=fwm.PHASE_ULTRA))
            stress.append({'D_GHz': D, 'T_K': temperature, 'pump_power_W': power, 'branch': branch,
                'arrays': array_comparisons(outputs[1], outputs[0]),
                'order_gate_equal': outputs[0]['floquet_convergence']['status'] == outputs[1]['floquet_convergence']['status']})
        quadrature = quadrature_checks()
    sampled = sampled_audit_counterexample()
    checked_keys = ('G_s', 'G_c', 'G_s_smallsignal', 'G_c_smallsignal', 'T_field_small_signal')
    all_errors = [full_errors, frozen_errors]+[row['arrays'] for row in stress]
    passed = (all(all(row[key]['maximum_relative_to_array_scale'] < 1e-9 for key in checked_keys)
        and row['S_dB']['maximum_absolute'] < 1e-7 for row in all_errors)
        and all(row['order_gate_equal'] for row in stress)
        and selected['floquet_convergence']['status'] == parent['floquet_convergence']['status']
        and selected['floquet_convergence']['full_scan_points'] == len(selected['probe_axis_GHz'])
        and sampled['sampled_status'] == 'CONVERGED' and sampled['full_status'] == 'UNCONVERGED')
    files = ['gabes/kernels.py', 'gabes/core.py', 'gabes/schemes/fwm.py', 'gabes/doppler.py',
        'gabes/observables.py', 'tests/test_floquet_symmetry.py', 'analysis/ultra_performance/verify.py']
    return {'proposal': PROPOSAL, 'scope': 'exact finite-Floquet performance patch; no new quantum/experimental squeezing claim',
        'expected_controls_passed': bool(passed), 'general_kernel_body_matches_commit_c0c46f0': True,
        'warm_seconds': times, 'median_speedup': float(np.median(times['general'])/np.median(times['exact_symmetry'])),
        'scan_points': len(selected['probe_axis_GHz']), 'velocity_points': selected['n_velocity'],
        'reported_floquet_order': selected['floquet_order'],
        'full_scan_order_comparison': selected['floquet_convergence']['comparison_order'],
        'full_scan_order_gate': selected['floquet_convergence']['status'],
        'compiled_general_vs_fast': full_errors, 'frozen_pre_patch_vs_fast': frozen_errors, 'stress_cases': stress,
        'sampled_audit_counterexample': sampled, 'plain_legendre_diagnostics': quadrature,
        'proposal_disposition': {
            'conjugate_symmetry_and_real_zero_block': 'implemented after exact input identities and finite-cutoff induction; general fallback and independent dense oracle retained',
            'sample_10_percent_of_order_audit': 'not adopted: cannot preserve the full-scan certificate; logical counterexample recorded',
            'reduce_nf_to_2': 'not adopted: a different truncated system; proposal itself reports changes beyond regression tolerance',
            '128_node_5sigma_velocity_default': 'not adopted: quadrature and cutoff change; raw 128/256 and 5/6-sigma diagnostics recorded, no universal error bound or validated composite panels',
            'commensurate_scan_delta_ladder': 'not adopted in this app patch: changes requested physical detunings by about 0.6 percent; exact reuse requires a separately declared scan grid',
            'cubic_operating_point_readout': 'separate accuracy change, not an identical computation or speed optimization; linear readout unchanged',
            'merge_probe_and_conjugate_seed_solves': 'rejected: different finite-drive states; misleading seed-amplitude-independence comment corrected'},
        'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'numba': numba.__version__,
            'numba_threads': numba.get_num_threads(), 'blas_threads_during_timing': 1},
        'source_sha256': {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in files},
        'pre_patch_artifact_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (HERE/'ultra_app_before.json', HERE/'ultra_app_before.npz')},
        'pre_patch_parameters': params}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = build_report()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    print(f'Wrote {args.output}')
    return 0 if report['expected_controls_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
