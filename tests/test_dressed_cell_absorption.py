"""Independent positive ordered-QRT absorption/emission controls."""

import numpy as np
import pytest

from gabes import constants as c, core, hyperfine, observables
from gabes.dressed_cell_absorption import (
    local_dressed_absorption, production_collisional_reservoirs,
    selected_dressed_cell_absorption)
from gabes.schemes import fwm


def _local_settings(delta_mhz=23.7, velocity=0.0, seed_power=None):
    if seed_power is None:
        seed_power = 3.06e-6/(20.3063457/26)
    probe = fwm.branch_center_GHz(1., -1)+delta_mhz*1e-3
    _, ks, kc = fwm.seeded_option_a_wavenumbers(1., np.array([probe]))
    return dict(T=386.15, pump_rabi_rad_s=float(fwm.rabi_freq(.2948, 530e-6)),
        probe_rabi_rad_s=float(fwm.rabi_freq(seed_power, 330e-6)),
        one_photon_detuning_rad_s=2*np.pi*1e9-fwm.K_VEC*velocity,
        two_photon_detuning_rad_s=delta_mhz*2*np.pi*1e6,
        optical_omega_rad_s=c.C_LIGHT*np.array([ks[0], kc[0]]),
        uniform_area_m2=observables.gaussian_mode_area(330e-6),
        weak_line_strength=.74/12, mean_order=3, response_order=3,
        phase_samples=32)


def _rates(result):
    return np.array([result[arm][key] for arm in ("probe", "conjugate")
                     for key in ("alpha_abs_per_m", "alpha_em_per_m", "alpha_net_per_m")])


@pytest.mark.parametrize("temperature", (300.0, 386.15, 410.0))
def test_explicit_jumps_reproduce_whole_production_dissipator(temperature):
    reservoirs, ledger = production_collisional_reservoirs(temperature)
    production = fwm.collisional_atom(temperature)
    np.testing.assert_allclose(reservoirs.dissipator(), production.lindblad,
                               rtol=1e-13, atol=1e-7)
    assert ledger["dissipator_relative_error"] < 1e-13
    assert ledger["positive_dephasing_construction_margin_s_inverse"] >= 0


@pytest.mark.parametrize("delta,velocity", ((0.0, 200.0), (23.7, 0.0),
                                            (52.269382, 0.0)))
def test_connected_covariance_source_matches_independent_jump_diffusion(delta, velocity):
    kwargs = _local_settings(delta, velocity)
    with core.blas_single_thread():
        direct = local_dressed_absorption(**kwargs, method="connected-qrt")
        diffusion = local_dressed_absorption(**kwargs, method="ordered-diffusion")
    np.testing.assert_allclose(_rates(direct), _rates(diffusion), rtol=1e-9, atol=1e-8)
    assert diffusion["field_audit"]["passed"]
    assert direct["atomic_diagnostics"]["passed"]
    assert direct["atomic_diagnostics"]["source_covariance_harmonic_order"] == 6


def test_zero_seed_net_drift_matches_independent_susceptibility_response():
    kwargs = _local_settings(seed_power=0.)
    with core.blas_single_thread():
        result = local_dressed_absorption(**kwargs, method="connected-qrt")
        response = fwm.pump_only_weak_response_reference(
            kwargs["pump_rabi_rad_s"], kwargs["pump_rabi_rad_s"],
            [kwargs["two_photon_detuning_rad_s"]],
            [kwargs["one_photon_detuning_rad_s"]],
            atom=fwm.collisional_atom(kwargs["T"]))
    chi = response.chi_matrix[0, 0, 0]
    expected = observables._gain_matrix_from_chi(
        [chi[0, 0]], [chi[0, 1]], [chi[1, 0].conjugate()], [chi[1, 1].conjugate()],
        kwargs["optical_omega_rad_s"][0]/c.C_LIGHT,
        kwargs["optical_omega_rad_s"][1]/c.C_LIGHT,
        hyperfine.number_density(kwargs["T"]), c.DIPOLE_D1,
        kwargs["weak_line_strength"])[0]
    actual = np.array([result[arm]["alpha_net_per_m"] for arm in ("probe", "conjugate")])
    np.testing.assert_allclose(actual, -2*np.real(np.diag(expected)), rtol=1e-9)


def test_positive_absorption_survives_net_conjugate_gain():
    with core.blas_single_thread():
        result = local_dressed_absorption(**_local_settings(), method="connected-qrt")
    conjugate = result["conjugate"]
    assert conjugate["alpha_abs_per_m"] > 0.0
    assert conjugate["alpha_em_per_m"] > conjugate["alpha_abs_per_m"]
    assert conjugate["alpha_net_per_m"] < 0.0
    assert 0.0 < np.exp(-conjugate["alpha_abs_per_m"]*.0125) < 1.0
    assert np.exp(-conjugate["alpha_net_per_m"]*.0125) > 1.0


def test_finite_seed_harmonic_refinement_at_narrow_raman_pole():
    kwargs = _local_settings(52.269382)
    with core.blas_single_thread():
        coarse = local_dressed_absorption(**kwargs, method="connected-qrt")
        fine = local_dressed_absorption(**dict(kwargs, mean_order=4, response_order=4),
                                        method="connected-qrt")
    np.testing.assert_allclose(_rates(coarse), _rates(fine), rtol=1e-7)


def test_selected_maxwell_absorption_is_bounded_with_reported_quadrature_error():
    kwargs = dict(D_GHz=1., T=386.15, P_pump=.2948,
        P_probe=3.06e-6/(20.3063457/26), w_pump=530e-6, w_probe=330e-6,
        delta_mhz=23.7, L=.0125, weak_line_strength=.74/12)
    result = selected_dressed_cell_absorption(**kwargs)
    for arm in ("probe", "conjugate"):
        values = result[arm]
        assert 0 <= values["in_cell_absorption_transmission"] <= 1
        assert values["absorption_quadrature_relative_error"] <= result["relative_tolerance"]
        assert values["alpha_abs_per_m"]-values["alpha_em_per_m"] == pytest.approx(
            values["alpha_net_per_m"], abs=1e-10)
    assert result["quadrature_success"]
    assert result["maximum_floquet_relative_rate_difference"] < 1e-6
    # Mutating a readout copy must not corrupt cached physical coefficients.
    saved = result["conjugate"]["in_cell_absorption_transmission"]
    result["conjugate"]["in_cell_absorption_transmission"] = -99
    assert selected_dressed_cell_absorption(**kwargs)["conjugate"]["in_cell_absorption_transmission"] == saved


def test_pump_and_seed_off_is_the_jasperse_pure_loss_limit():
    result = selected_dressed_cell_absorption(D_GHz=1., T=386.15,
        P_pump=0., P_probe=0., w_pump=530e-6, w_probe=330e-6,
        delta_mhz=23.7, L=.0125, weak_line_strength=.74/12)
    for arm in ("probe", "conjugate"):
        assert result[arm]["alpha_em_per_m"] == pytest.approx(0., abs=1e-12)
        assert result[arm]["in_cell_absorption_transmission"] == pytest.approx(
            result[arm]["net_diagonal_transmission"], abs=1e-12)


@pytest.mark.parametrize("overrides", ({"branch": 1}, {"T": -1},
    {"relative_tolerance": 0.}, {"transit_rate": 0.}, {"L": -1.}))
def test_selected_invalid_or_unsupported_inputs_rejected(overrides):
    kwargs = dict(D_GHz=1., T=386.15, P_pump=.2948, P_probe=3e-6,
        w_pump=530e-6, w_probe=330e-6, delta_mhz=23.7, L=.0125,
        weak_line_strength=.74/12)
    with pytest.raises((ValueError, NotImplementedError)):
        selected_dressed_cell_absorption(**dict(kwargs, **overrides))
