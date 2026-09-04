"""Independent finite-Floquet references and production convergence gates."""

import numpy as np
import pytest

from gabes import core, doppler, kernels
from gabes.schemes import fwm


def _atomic_fixture():
    pump = fwm.rabi_freq(0.6, fwm.W_PUMP)
    seed = fwm.rabi_freq(8e-6, fwm.W_PROBE)
    delta = 2.0 * np.pi * -8e6
    atom = fwm.collisional_atom(394.15)
    L0 = core.build_liouvillian(
        fwm.static_hamiltonian_at_Deff_zero(
            pump, pump, seed, delta, branch=-1),
        atom,
    )
    Cp, Cm = fwm.sideband_template(pump, pump, 0.0, branch=-1)
    deff = 2.0 * np.pi * np.array([0.82e9, 0.97e9])
    omega = fwm.seeded_sideband_beat(delta, branch=-1)
    return L0, Cp, Cm, omega, deff, atom.S_v


@pytest.mark.parametrize("n_f", [1, 2, 3])
def test_continued_fraction_matches_independent_dense_block(n_f):
    L0, Cp, Cm, omega, deff, S_v = _atomic_fixture()
    continued = core.floquet_solve_truncated(
        L0, Cp, Cm, omega, deff, S_v, fwm.N_LEVELS, n_f,
        return_harmonics=True,
    )
    direct = core.floquet_solve_direct(
        L0, Cp, Cm, omega, deff, S_v, fwm.N_LEVELS, n_f,
        return_harmonics=True,
    )
    assert np.max(np.abs(continued - direct)) < 1e-10

    diagnostics = core.floquet_solution_diagnostics(
        L0, Cp, Cm, omega, deff, S_v, continued)
    assert diagnostics["n_f"] == n_f
    assert diagnostics["max_normalized_residual"] < 1e-10
    assert diagnostics["max_trace_error"] < 1e-10
    assert diagnostics["max_nonzero_harmonic_trace"] < 1e-10
    assert diagnostics["max_rho0_hermiticity_error"] < 1e-10


def test_floquet_diagnostics_fail_closed_on_nonfinite_sideband_entry():
    L0, Cp, Cm, omega, deff, S_v = _atomic_fixture()
    harmonics = core.floquet_solve_truncated(
        L0, Cp, Cm, omega, deff, S_v, fwm.N_LEVELS, 3,
        return_harmonics=True,
    )
    harmonics[0, 0, 0, 1] = np.nan
    diagnostics = core.floquet_solution_diagnostics(
        L0, Cp, Cm, omega, deff, S_v, harmonics)
    assert diagnostics["finite"] is False
    assert diagnostics["nonfinite_entries"] == 1
    assert np.isinf(diagnostics["max_normalized_residual"])
    assert np.isinf(diagnostics["max_trace_error"])


@pytest.mark.skipif(not kernels.available(), reason="numba not available")
@pytest.mark.parametrize("n_f", [1, 2, 3])
def test_compiled_floquet_orders_match_numpy_reference(monkeypatch, n_f):
    pump = fwm.rabi_freq(0.6, fwm.W_PUMP)
    seed = fwm.rabi_freq(8e-6, fwm.W_PROBE)
    branch = -1
    D_GHz = 0.9
    Delta = 2.0 * np.pi * D_GHz * 1e9
    center = fwm.branch_center_GHz(D_GHz, branch)
    probe = np.linspace(center - 0.01, center + 0.01, 3)
    delta = fwm.two_photon_detuning_from_probe_scan(probe, D_GHz, branch)
    velocity, _ = doppler.velocity_grid(394.15, dv=100.0, cutoff_sigma=1.0)
    deff = doppler.build_Delta_eff_axis(Delta, Delta, velocity)
    atom = fwm.collisional_atom(394.15)

    fast = fwm.chi_matrix_table(
        pump, pump, seed, seed, delta, deff, branch, atom=atom, n_f=n_f)
    monkeypatch.setattr(fwm.kernels, "available", lambda: False)
    reference = fwm.chi_matrix_table(
        pump, pump, seed, seed, delta, deff, branch, atom=atom, n_f=n_f)
    for actual, expected in zip(fast, reference):
        scale = max(float(np.max(np.abs(expected))), 1e-300)
        assert np.max(np.abs(actual - expected)) / scale < 1e-10


def _identity_response(n=5):
    response = {
        name: np.full(n, 1e-12 + 2e-13j, dtype=complex)
        for name in ("chi_ss", "chi_cs", "chi_sc", "chi_cc")
    }
    transfer = np.broadcast_to(np.eye(2, dtype=complex), (n, 2, 2)).copy()
    return response, transfer


def test_scan_gate_rejects_a_hidden_phase_or_gain_change():
    response, transfer = _identity_response()
    passed = fwm.assess_floquet_scan_convergence(
        high_order=3,
        low_order=2,
        high_response=response,
        low_response={key: value.copy() for key, value in response.items()},
        high_transfer=transfer,
        low_transfer=transfer.copy(),
        scan_axis=np.arange(5.0),
    )
    assert passed["status"] == "CONVERGED"

    low_transfer = transfer.copy()
    low_transfer[2, 0, 0] = 1.1 * np.exp(1j * np.deg2rad(1.0))
    failed = fwm.assess_floquet_scan_convergence(
        high_order=3,
        low_order=2,
        high_response=response,
        low_response={key: value.copy() for key, value in response.items()},
        high_transfer=transfer,
        low_transfer=low_transfer,
        scan_axis=np.arange(5.0),
    )
    assert failed["status"] == "UNCONVERGED"
    assert not failed["gains"]["probe_power"]["passed"]
    assert not failed["wrapped_phase"]["passed"]

    low_response = {key: value.copy() for key, value in response.items()}
    low_response["chi_ss"][:] = 0.0
    response_failed = fwm.assess_floquet_scan_convergence(
        high_order=3,
        low_order=2,
        high_response=response,
        low_response=low_response,
        high_transfer=transfer,
        low_transfer=transfer.copy(),
        scan_axis=np.arange(5.0),
    )
    assert response_failed["status"] == "UNCONVERGED"
    assert not response_failed["response_components"]["chi_ss"]["passed"]


def test_scan_gate_cannot_be_certified_by_one_point_or_mismatched_response():
    response, transfer = _identity_response(n=1)
    with pytest.raises(ValueError, match="at least two scan points"):
        fwm.assess_floquet_scan_convergence(
            high_order=3, low_order=2,
            high_response=response, low_response=response,
            high_transfer=transfer, low_transfer=transfer,
            scan_axis=np.array([0.0]),
        )

    response, transfer = _identity_response(n=5)
    bad = dict(response)
    bad["chi_ss"] = bad["chi_ss"][:1]
    with pytest.raises(ValueError, match="full-scan shape"):
        fwm.assess_floquet_scan_convergence(
            high_order=3, low_order=2,
            high_response=response, low_response=bad,
            high_transfer=transfer, low_transfer=transfer,
            scan_axis=np.arange(5.0),
        )


def test_scan_gate_returns_unconverged_for_nonfinite_transfer():
    response, transfer = _identity_response(n=5)
    low_transfer = transfer.copy()
    low_transfer[:, 0, 0] = np.nan
    result = fwm.assess_floquet_scan_convergence(
        high_order=3, low_order=2,
        high_response=response, low_response=response,
        high_transfer=transfer, low_transfer=low_transfer,
        scan_axis=np.arange(5.0),
    )
    assert result["status"] == "UNCONVERGED"
    assert result["finite"] is False
    assert result["nonfinite_entries"] == 5
    assert result["probe_gain_optimum"]["passed"] is False
    assert result["probe_gain_optimum"]["high_index"] is None


def test_default_seeded_scan_uses_nf3_and_full_scan_gate():
    spectrum = fwm.compute_spectrum(
        0.9,
        T=394.15,
        coarse_points=9,
        fine_points=0,
        scan_min=-2.2,
        scan_max=-2.0,
        velocity_step=100.0,
        velocity_cutoff=1.0,
        phase_detail=fwm.PHASE_BALANCED,
    )
    convergence = spectrum["floquet_convergence"]
    assert spectrum["floquet_order"] == 3
    assert convergence["high_order"] == 3
    assert convergence["comparison_order"] == 2
    assert convergence["full_scan_points"] == spectrum["probe_axis_GHz"].size
    assert convergence["status"] == "CONVERGED"
    assert spectrum["atomic_solver_provenance"]["floquet_modes"] == tuple(range(-3, 4))
    scheme = fwm.FWMScheme()
    rendered = scheme.observables(
        spectrum, scheme.defaults(), include_figures=False)
    assert [metric["label"] for metric in rendered["metrics"]] == [
        "Squeezing indicator", "Seed gain G_s", "Conjugate gain G_c"
    ]
    assert not any(metric["label"] == "Floquet truncation"
                   for metric in rendered["metrics"])
    assert any(table["title"] == "Model diagnostics"
               for table in rendered["tables"])
    markdown = "\n".join(table["markdown"] for table in rendered["tables"])
    assert "| Floquet adjacent-order gate | CONVERGED vs N_F=2 |" in markdown


def test_unchecked_low_order_is_explicitly_historical():
    spectrum = fwm.compute_spectrum(
        0.9,
        T=394.15,
        coarse_points=5,
        fine_points=0,
        scan_min=-2.2,
        scan_max=-2.0,
        velocity_step=150.0,
        velocity_cutoff=1.0,
        phase_detail=fwm.PHASE_BALANCED,
        floquet_order=1,
        enforce_floquet_convergence=False,
    )
    assert spectrum["floquet_convergence"]["status"] == "NOT_EVALUATED"
    assert "FLOQUET_UNCONVERGED" in spectrum["claim_gate"]["badges"]
