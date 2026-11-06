"""Planned Na/K records: present as data, not offered, consistent with their sources."""

import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import probe_saturation, species  # noqa: E402


def test_planned_elements_are_not_offered():
    planned = set(species.PLANNED_ISOTOPES)
    assert planned == {"Na23", "K39", "K40", "K41"}
    assert not planned & set(species.ISOTOPES)
    offered = {iso.name for entries in species.SPECIES.values() for iso, _ in entries}
    assert not planned & offered
    for iso in species.PLANNED_ISOTOPES.values():
        with pytest.raises(ValueError, match="planned, not yet supported"):
            species.number_density(iso, 330.0)


@pytest.mark.parametrize("name, splitting_mhz", [
    # Ground hyperfine splitting A(I+1/2): Steck Na Table/text 1771.626 MHz;
    # Tiecke Fig. 2 and Hanley Fig. 2 level diagrams for K.
    ("Na23", 1771.626), ("K39", 461.7), ("K40", 1285.8), ("K41", 254.0)])
def test_ground_hyperfine_splitting_matches_the_level_diagrams(name, splitting_mhz):
    iso = species.PLANNED_ISOTOPES[name]
    assert abs(iso.A_S) * (iso.I + 0.5) == pytest.approx(splitting_mhz, abs=0.05)


@pytest.mark.parametrize("name, line, isat_mw_cm2, rel", [
    ("Na23", "D2", 6.2600, 5e-4),   # Steck Table 7, cycling transition
    ("K39", "D2", 1.75, 5e-3),      # Tiecke Table 4
    ("K39", "D1", 1.71, 1e-2),      # Hanley et al. 2015, §4.1
])
def test_closed_two_level_saturation_intensity_matches_the_sources(name, line, isat_mw_cm2, rel):
    iso = species.PLANNED_ISOTOPES[name]
    _, nu, gamma_mhz, _, _ = iso.line(line)
    isat = probe_saturation.closed_two_level_saturation_intensity(nu, 2 * math.pi * gamma_mhz * 1e6)
    assert isat * 0.1 == pytest.approx(isat_mw_cm2, rel=rel)


@pytest.mark.parametrize("element, temperature_k, pressure, rel", [
    ("Na", 298.15, 2.38e-11, 0.05),   # Steck Table 2 (torr), ±5 % model accuracy
    ("K", 293.0, 1.3e-8, 0.02),       # Tiecke Table 2 (mbar)
])
def test_vapor_pressure_fits_reproduce_the_tabulated_room_temperature_value(
        element, temperature_k, pressure, rel):
    fit = species.PLANNED_VAPOR_PRESSURE[element]
    a, b = fit["solid"] if temperature_k < fit["melting_k"] else fit["liquid"]
    assert 10 ** (a - b / temperature_k) == pytest.approx(pressure, rel=rel)
