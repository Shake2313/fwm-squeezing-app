"""Shared heated-cell layer: temperature state, cold-spot density, axial profile."""

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import constants, rydberg_experiment, species, vapor_cell  # noqa: E402
from gabes.schemes.rydberg import RydbergEITScheme  # noqa: E402


def test_rydberg_experiment_reexports_the_shared_cell_layer():
    for name in vapor_cell.__all__:
        if name in rydberg_experiment.__all__:
            assert getattr(rydberg_experiment, name) is getattr(vapor_cell, name), name
    moved = {"CellTemperatureState", "resolve_cell_temperature", "AxialCellProfile",
             "axial_cell_profile", "linear_axial_cell_profile", "SegmentedBeerLambert",
             "integrate_beer_lambert", "LinearTemperatureCalibration",
             "fit_effective_temperature_calibration"}
    assert moved <= set(rydberg_experiment.__all__) & set(vapor_cell.__all__)


@pytest.mark.parametrize("iso", [species.RB85, species.RB87, species.CS133],
                         ids=lambda iso: iso.name)
def test_cold_spot_limited_density_holds_the_reservoir_pressure(iso):
    cold_k = 311.4
    path_k = np.array([311.4, 320.0, 345.25, 372.0])
    n_cold = species.number_density(iso, cold_k)
    density = vapor_cell.cold_spot_limited_density(n_cold, cold_k, path_k)
    # Independent reference: the ideal gas along the beam sits at the reservoir's
    # saturated pressure, P = n k_B T for every path temperature.
    pressure = (species.rb_vapor_pressure_pa(cold_k) if iso.name.startswith("Rb")
                else n_cold * constants.KB * cold_k)
    np.testing.assert_allclose(density * constants.KB * path_k, pressure, rtol=1e-12)
    assert density[0] == pytest.approx(n_cold, rel=1e-15)
    assert np.all(np.diff(density) < 0)


@pytest.mark.parametrize("temps", [(45.0, 45.0, 45.0), (62.0, 58.0, 41.5), (90.0, 33.25, 21.0)])
def test_beam_path_density_matches_the_rydberg_scheme(temps):
    heater, effective, cold = temps
    state = vapor_cell.resolve_cell_temperature(
        heater, effective_vapor_temp_c=effective, cold_spot_temp_c=cold)
    # The Rydberg scheme still computes this inline (its file is being edited
    # elsewhere); this pins the two equal until it calls the shared method.
    assert state.beam_path_density_m3() == RydbergEITScheme._local_density_from_temperature_state(state)
    assert state.beam_path_density_m3(species.CS133) == vapor_cell.cold_spot_limited_density(
        species.number_density(species.CS133, cold + 273.15), cold + 273.15, effective + 273.15)


def test_axial_profile_uses_the_shared_density():
    temperature = np.array([52.0, 47.5, 44.0, 49.0, 61.0])
    profile = vapor_cell.axial_cell_profile(np.linspace(0.0, 0.075, 5), temperature,
                                            pressure_cold_spot_temp_c=39.0)
    expected = vapor_cell.cold_spot_limited_density(
        profile.cold_spot_density_m3, 39.0 + 273.15, temperature + 273.15)
    assert profile.density_m3.tobytes() == expected.tobytes()
