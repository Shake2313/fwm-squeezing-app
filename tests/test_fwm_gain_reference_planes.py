"""Independent optical power chains and passive limits for gain readout."""

import numpy as np
import pytest

from gabes import fwm_gain_closure, observables
from gabes.gain_reference import pump_off_referenced_gains


def test_reference_conversion_matches_independent_random_optical_power_chains():
    rng = np.random.default_rng(20261015)
    incident = rng.uniform(1e-6, 20e-6, 60)
    entrance_window = rng.uniform(0.8, 1.0, 60)
    atomic_off = rng.uniform(0.1, 1.0, 60)
    vapor_seed = incident * entrance_window
    probe_vapor_gain = rng.uniform(0.0, 100.0, 60)
    conjugate_vapor_gain = rng.uniform(0.0, 90.0, 60)
    window_s, window_c = rng.uniform(0.8, 1.0, (2, 60))
    collection_s, collection_c = rng.uniform(0.2, 1.0, (2, 60))
    downstream_s, downstream_c = rng.uniform(0.4, 1.0, (2, 60))
    reference_downstream = rng.uniform(0.4, 1.0, 60)

    # Independent end-to-end powers; entrance and exit windows are separate.
    probe_detector = (vapor_seed * probe_vapor_gain * window_s
                      * collection_s * downstream_s)
    conjugate_detector = (vapor_seed * conjugate_vapor_gain * window_c
                          * collection_c * downstream_c)
    reference_detector = (vapor_seed * atomic_off * window_s
                          * reference_downstream)

    converted = pump_off_referenced_gains(
        probe_vapor_gain, conjugate_vapor_gain,
        pump_off_atomic_transmission=atomic_off,
        probe_output_transmission=window_s * collection_s * downstream_s,
        conjugate_output_transmission=window_c * collection_c * downstream_c,
        reference_output_transmission=window_s * reference_downstream)
    np.testing.assert_allclose(
        converted["probe_power_gain"], probe_detector / reference_detector,
        rtol=6e-16)
    np.testing.assert_allclose(
        converted["conjugate_power_gain"], conjugate_detector / reference_detector,
        rtol=6e-16)
    np.testing.assert_allclose(
        reference_detector * converted["medium_input_seed_per_reference_power"],
        vapor_seed, rtol=6e-16)


def test_pump_off_absorption_has_unit_output_referenced_gain_with_common_optics():
    atomic = np.array([0.02, 0.4, 0.91, 1.0])
    common = np.array([0.9, 0.5, 0.75, 0.82])
    converted = pump_off_referenced_gains(
        atomic, 0.0, pump_off_atomic_transmission=atomic,
        probe_output_transmission=common,
        conjugate_output_transmission=common,
        reference_output_transmission=common)
    np.testing.assert_allclose(converted["probe_power_gain"], 1.0, atol=2e-16)
    np.testing.assert_array_equal(converted["conjugate_power_gain"], 0.0)


def test_total_cell_throughput_does_not_identify_pump_off_atomic_transmission():
    # Both complete cells transmit .72, yet their atomic transmissions differ.
    first_atomic, first_window = 0.8, 0.9
    second_atomic, second_window = 0.9, 0.8
    assert first_atomic * first_window == second_atomic * second_window
    first = pump_off_referenced_gains(
        12.0, 11.0, pump_off_atomic_transmission=first_atomic,
        probe_output_transmission=first_window,
        conjugate_output_transmission=first_window,
        reference_output_transmission=first_window)
    second = pump_off_referenced_gains(
        12.0, 11.0, pump_off_atomic_transmission=second_atomic,
        probe_output_transmission=second_window,
        conjugate_output_transmission=second_window,
        reference_output_transmission=second_window)
    assert first["probe_power_gain"] == pytest.approx(15.0)
    assert second["probe_power_gain"] == pytest.approx(13.333333333333334)
    unknown = pump_off_referenced_gains(12.0, 11.0)
    assert unknown["status"] == "UNAVAILABLE"
    assert unknown["probe_power_gain"] is None
    assert unknown["conjugate_power_gain"] is None
    assert unknown["medium_input_seed_per_reference_power"] is None
    assert not unknown["provenance"]["total_cell_transmission_identifies_A0"]


def test_passive_gap_bound_cannot_be_applied_to_pump_off_output_gains():
    # A passive attenuator after an ideal two-mode amplifier gives D_source=.792.
    gain = 5.0
    amplitude_map = np.sqrt(0.792) * np.array([
        [np.sqrt(gain), np.sqrt(gain - 1.0)],
        [np.sqrt(gain - 1.0), np.sqrt(gain)],
    ])
    J = np.diag([1.0, -1.0])
    # This is a physical passive-loss completion: the environment is one vacuum
    # annihilation channel per optical mode, with the correct Nambu sign.
    environment_commutator = 0.208 * J
    np.testing.assert_allclose(
        amplitude_map @ J @ amplitude_map.T + environment_commutator, J,
        atol=2e-15)
    source_s, source_c = np.abs(amplitude_map[:, 0])**2
    assert source_s - source_c == pytest.approx(0.792)
    converted = pump_off_referenced_gains(
        source_s, source_c, pump_off_atomic_transmission=0.6)
    assert (converted["probe_power_gain"] - converted["conjugate_power_gain"]
            == pytest.approx(1.32))


def test_reference_conversion_broadcasts_and_allows_blocked_output_arm():
    converted = pump_off_referenced_gains(
        np.array([[2.0], [3.0]]), np.array([0.5, 1.0, 1.5]),
        pump_off_atomic_transmission=np.array([0.5, 0.75, 1.0]),
        probe_output_transmission=0.0, conjugate_output_transmission=0.5,
        reference_output_transmission=0.5)
    assert converted["probe_power_gain"].shape == (2, 3)
    np.testing.assert_array_equal(converted["probe_power_gain"], 0.0)
    np.testing.assert_allclose(
        converted["conjugate_power_gain"], [[1.0, 4.0 / 3.0, 1.5]] * 2)
    np.testing.assert_allclose(
        converted["medium_input_seed_per_reference_power"],
        [[4.0, 8.0 / 3.0, 2.0]] * 2)


@pytest.mark.parametrize("overrides", [
    {"G_s": -0.1}, {"G_c": np.nan},
    {"pump_off_atomic_transmission": 0.0},
    {"pump_off_atomic_transmission": 1.01},
    {"reference_output_transmission": 0.0},
    {"probe_output_transmission": -0.1},
    {"conjugate_output_transmission": np.inf},
    {"probe_output_transmission": 1.01},
])
def test_reference_conversion_rejects_invalid_physical_inputs(overrides):
    inputs = dict(G_s=2.0, G_c=1.0, pump_off_atomic_transmission=0.8)
    inputs.update(overrides)
    with pytest.raises(ValueError):
        pump_off_referenced_gains(**inputs)


def test_depletion_cap_preserves_beer_lambert_absorption_even_with_zero_pump():
    passive = np.exp(-np.array([0.0, 0.1, 1.0, 7.0, 30.0]))
    for pump in (0.0, 0.2948, 0.6):
        capped_s, capped_c = observables.pump_depletion_saturation(
            passive, np.zeros_like(passive), pump, 3.92e-6)
        np.testing.assert_array_equal(capped_s, passive)
        np.testing.assert_array_equal(capped_c, 0.0)
    legacy, _ = observables.pump_depletion_saturation(
        passive, np.zeros_like(passive), 0.2948, 3.92e-6,
        legacy_probe_floor=True)
    np.testing.assert_array_equal(legacy, 1.0)


def test_corrected_cap_retains_positive_gain_and_pump_budget_regression():
    gs = np.array([0.02, 0.5, 1.0, 1.000001, 20.0, 1e9])
    gc = np.array([0.0, 0.1, 0.0, 0.000001, 19.0, 1e9 - 1.0])
    actual = observables.pump_depletion_saturation(gs, gc, 0.2948, 3.92e-6)
    legacy = observables.pump_depletion_saturation(
        gs, gc, 0.2948, 3.92e-6, legacy_probe_floor=True)
    np.testing.assert_array_equal(actual[0][gs >= 1.0], legacy[0][gs >= 1.0])
    np.testing.assert_array_equal(actual[1], legacy[1])
    assert np.all((actual[0] - 1.0 + actual[1]) * 3.92e-6 <= 0.2948)


def test_active_gain_fit_preserves_absorbed_probe_without_changing_fit_domain():
    calibration = fwm_gain_closure.laboratory_provenance(
        enabled=True, eligible=True, D_GHz=1.06, branch=-1)
    gs = np.array([0.01, 0.5, 1.0, 3.0, 9.0, 17.0])
    gc = np.maximum(gs - 1.0, 0.0)
    for pump in (0.0, 0.380):
        corrected = fwm_gain_closure.apply_power_gain(
            gs, gc, calibration=calibration, P_pump=pump, P_seed=3.7e-6)
        legacy = fwm_gain_closure.apply_power_gain(
            gs, gc, calibration=calibration, P_pump=pump, P_seed=3.7e-6,
            legacy_probe_floor=True)
        np.testing.assert_array_equal(corrected[0][:3], gs[:3])
        np.testing.assert_array_equal(legacy[0][:3], 1.0)
        np.testing.assert_array_equal(corrected[0][3:], legacy[0][3:])
        np.testing.assert_array_equal(corrected[1], legacy[1])


def test_calibration_distinguishes_simulation_seed_from_measured_reference():
    provenance = fwm_gain_closure.laboratory_provenance(
        enabled=True, eligible=True, D_GHz=1.06, branch=-1)
    conditions = provenance["calibration_conditions"]
    assert conditions["cell_input_seed_W"] is None
    assert conditions["assumed_simulation_seed_W"] == 3.7e-6
    assert conditions["measured_pump_off_reference_seed_W"] == 3.7e-6
    assert provenance["passive_cell_transmission"] is None
