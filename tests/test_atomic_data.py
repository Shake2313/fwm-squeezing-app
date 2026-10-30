"""
Single-source guards for the atomic-data layer.

`species` holds the vapor pressure, densities and β; `hyperfine` is its ⁸⁵Rb D1
view; `angular` holds the Racah algebra. `constants` keeps legacy ⁸⁵Rb D1
values whose deviations from the species record are pinned here, so changing
either side is a visible decision (checklist
`atomic-data-legacy-value-unification`).

    python tests/test_atomic_data.py     # or: pytest tests/test_atomic_data.py
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import angular, atoms, constants, hyperfine, species, zeeman  # noqa: E402

TEMPS = np.linspace(250.0, 480.0, 47)


def test_hyperfine_is_a_bitwise_view_of_species():
    assert hyperfine.BETA_SELF == species.BETA_SELF
    assert hyperfine.vapor_pressure_pa is species.rb_vapor_pressure_pa
    for T in TEMPS:
        assert hyperfine.number_density(T) == species.number_density(species.RB85, T)
    N = np.array([0.0, 1e15, 3.3e17, 2e19, 7e19])
    np.testing.assert_array_equal(
        hyperfine.self_broadened_gamma(N),
        constants.GAMMA + species.self_broadening(species.RB85, N))


def test_data_layer_imports_without_level_schemes():
    probe = ("import sys, gabes.hyperfine; "
             "print(sorted(m for m in sys.modules if m.startswith('gabes')))")
    out = subprocess.run([sys.executable, "-c", probe], cwd=ROOT, check=True,
                         capture_output=True, text=True).stdout
    assert out.strip() == str(["gabes", "gabes.angular", "gabes.constants",
                               "gabes.hyperfine", "gabes.species"])


def test_angular_algebra_has_one_home():
    assert zeeman.clebsch_gordan is angular.clebsch_gordan
    assert species.wigner6j is angular.wigner6j
    assert species.wigner3j is angular.wigner3j


def test_steck_density_is_the_second_convention():
    assert atoms.rb85_density is species.rb85_steck_density
    ratio = [species.rb85_steck_density(T) / species.number_density(species.RB85, T)
             for T in np.linspace(313.15, 393.15, 9)]
    assert 0.50 < min(ratio) and max(ratio) < 0.60


@pytest.mark.parametrize("legacy, record, deviation", [
    (constants.GAMMA_2PI, species.RB85.gamma_D1_mhz * 1e6, -6.956521739e-4),
    (constants.MASS_85RB, species.RB85.mass, 4.651972889e-6),
    (constants.NU_EXCITED_HF_D1,
     3.0 * species.RB85.A_P12 * 1e6, -2.765632044e-6),
    (constants.NU_GROUND_HF, 3.0 * species.RB85.A_S * 1e6, 0.0),
    (constants.NU_D1_85RB, species.RB85.nu_D1, 0.0),
    (constants.RB85_ABUNDANCE, species.RB85.abundance, 0.0),
    (species.GJ_S12_LEGACY, species.RB87.gJ_S, -2.996507376e-8),
])
def test_legacy_constants_deviate_only_as_documented(legacy, record, deviation):
    assert legacy / record - 1.0 == pytest.approx(deviation, rel=1e-6, abs=1e-15)


@pytest.mark.parametrize("iso", [species.RB85, species.RB87, species.CS133])
def test_lande_gf_matches_the_j_half_closed_form(iso):
    I, gJ, gI = iso.I, iso.gJ_S, iso.g_I
    upper = species.lande_gf(gJ, gI, I, 0.5, I + 0.5)
    lower = species.lande_gf(gJ, gI, I, 0.5, I - 0.5)
    assert upper == pytest.approx((gJ + 2 * I * gI) / (2 * I + 1), rel=1e-14)
    assert lower == pytest.approx((-gJ + (2 * I + 2) * gI) / (2 * I + 1), rel=1e-14)
    assert species.lande_gf(gJ, gI, I, 0.5, 0) == 0.0


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
