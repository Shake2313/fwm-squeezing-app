"""Carrier centroids and passive-limit integration of the seeded mean map."""

import numpy as np
import pytest

from gabes import constants
from gabes.schemes import fwm


def _casimir_hyperfine_energy(f, splitting):
    """I=5/2,J=1/2 has A=hyperfine splitting/3; no quadrupole term."""
    nuclear, electronic = 2.5, 0.5
    k = f * (f + 1) - nuclear * (nuclear + 1) - electronic * (electronic + 1)
    return (splitting / 3) * k / 2


def _transition_origin():
    ground_f2 = _casimir_hyperfine_energy(2, constants.OMEGA_HF)
    excited_f3 = _casimir_hyperfine_energy(3, constants.OMEGA_EXCITED_HF)
    return excited_f3 - ground_f2


def test_optical_origin_matches_transition_between_independent_centroid_energies():
    for splitting in (constants.OMEGA_HF, constants.OMEGA_EXCITED_HF):
        lower = _casimir_hyperfine_energy(2, splitting)
        upper = _casimir_hyperfine_energy(3, splitting)
        assert upper - lower == pytest.approx(splitting, rel=2e-16)
        assert 5 * lower + 7 * upper == pytest.approx(0, abs=2e-5)
    actual, _, _ = fwm.seeded_option_a_wavenumbers(0, np.array([0.0]))
    expected = (constants.OMEGA_D1 + _transition_origin()) / constants.C_LIGHT
    assert actual == pytest.approx(expected, rel=0, abs=2e-9)
    legacy, _, _ = fwm.seeded_option_a_wavenumbers(
        0, np.array([0.0]), legacy_centroid=True)
    assert legacy == constants.OMEGA_D1 / constants.C_LIGHT
    assert (actual - legacy) * constants.C_LIGHT == pytest.approx(
        _transition_origin(), rel=3e-11)


@pytest.mark.parametrize("legacy", [False, True])
def test_carriers_preserve_energy_conservation_and_vacuum_geometric_phase(legacy):
    detuning = 1.0
    probes = np.array([-2.09, -2.035732439, -2.01])
    pump, seed, conjugate = fwm.seeded_option_a_wavenumbers(
        detuning, probes, legacy_centroid=legacy)
    np.testing.assert_allclose(seed + conjugate, 2 * pump, rtol=0, atol=4e-9)
    # With no crossing angle, the four-wave-mixing vacuum mismatch vanishes.
    collinear = fwm.seeded_phase_mismatch_z(
        detuning, probes, angle_deg=0, legacy_centroid=legacy)
    np.testing.assert_allclose(collinear, 0, rtol=0, atol=4e-9)
    angle = np.deg2rad(0.32)
    geometric = fwm.seeded_phase_mismatch_z(
        detuning, probes, angle_deg=0.32, legacy_centroid=legacy)
    expected = 4 * pump * np.sin(angle / 2)**2
    np.testing.assert_allclose(geometric, expected, rtol=0, atol=4e-9)


@pytest.fixture(scope="module")
def tiny_pump_off_spectra():
    # Three Maxwell points and three optical points exercise the real kernel.
    # These are diagnostic limits, not an ensemble convergence declaration.
    settings = dict(
        T=320.0, P_pump=0.0, P_probe=1e-12,
        phase_detail=fwm.PHASE_BALANCED, gain_closure_enabled=False,
        coarse_points=3, fine_points=0, scan_min=-3.04, scan_max=-3.03,
        velocity_cutoff=0.1, velocity_step=200.0,
        floquet_order=3, enforce_floquet_convergence=False)
    atomic_results = []
    original = fwm.chi_matrix_table

    def record_atomic_results(*args, **kwargs):
        result = original(*args, **kwargs)
        atomic_results.append(tuple(value.copy() for value in result))
        return result

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(fwm, "chi_matrix_table", record_atomic_results)
        modern = fwm.compute_spectrum(0.0, **settings)
        historical = fwm.compute_spectrum(
            0.0, legacy_mean_field_fixture=True, **settings)
    return modern, historical, atomic_results


def test_tiny_compute_preserves_real_pump_off_absorption_and_labels_old_floor(
        tiny_pump_off_spectra):
    modern, historical, _ = tiny_pump_off_spectra
    assert modern["n_velocity"] == 3
    assert len(modern["probe_axis_GHz"]) == 3
    assert np.all((modern["G_s"] > 0) & (modern["G_s"] < 1))
    np.testing.assert_array_equal(modern["G_s"], modern["G_s_smallsignal"])
    np.testing.assert_array_equal(modern["G_c"], 0)
    np.testing.assert_array_equal(historical["G_s"], 1)
    np.testing.assert_array_equal(historical["G_c"], 0)
    assert np.all(historical["G_s_smallsignal"] < 1)
    assert not modern["legacy_mean_field_fixture"]
    assert historical["legacy_mean_field_fixture"]
    assert "HISTORICAL_MEAN_FIELD_FIXTURE" in historical["claim_gate"]["badges"]


def test_fixture_flag_changes_carriers_and_readout_without_shifting_atomic_response(
        tiny_pump_off_spectra):
    modern, historical, atomic_results = tiny_pump_off_spectra
    assert len(atomic_results) == 2
    for present, old in zip(*atomic_results):
        np.testing.assert_array_equal(present, old)
    np.testing.assert_array_equal(modern["probe_axis_GHz"], historical["probe_axis_GHz"])
    for raw in (modern, historical):
        np.testing.assert_allclose(
            raw["omega_probe_rad_s"] / constants.C_LIGHT,
            raw["k_probe_propagation_per_m"], rtol=0, atol=2e-9)
        np.testing.assert_allclose(
            raw["omega_conjugate_rad_s"] / constants.C_LIGHT,
            raw["k_conjugate_propagation_per_m"], rtol=0, atol=2e-9)
    for frequency in ("omega_probe_rad_s", "omega_conjugate_rad_s"):
        np.testing.assert_allclose(
            modern[frequency] - historical[frequency],
            _transition_origin(), rtol=0, atol=1.0)


def test_raw_gain_reference_remains_unknown_with_explicit_medium_boundary_planes(
        tiny_pump_off_spectra):
    modern, historical, _ = tiny_pump_off_spectra
    for raw in (modern, historical):
        assert raw["input_seed_power_W"] == 1e-12
        assert raw["input_seed_plane"] == "atomic medium entrance, after entrance window"
        assert raw["output_gain_plane"] == "atomic medium exit, before exit window"
        reference = raw["gain_reference"]
        assert reference["status"] == "UNAVAILABLE"
        assert reference["probe_power_gain"] is None
        assert reference["conjugate_power_gain"] is None
        assert reference["pump_off_atomic_transmission"] is None
        assert not reference["provenance"]["total_cell_transmission_identifies_A0"]
        assert not reference["provenance"]["atomic_response_modified"]
