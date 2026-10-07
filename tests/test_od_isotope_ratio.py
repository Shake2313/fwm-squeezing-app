"""Isotope-enriched Rb shares one vapor density and a stable frequency origin."""
from unittest.mock import patch

import numpy as np
import pytest

from gabes import constants, schemes, species
from gabes.schemes.sas import _probe_reference_matched


SAS = schemes.get("sas")


def _params(**overrides):
    params = SAS.defaults()
    params.update(pump_power_mw=0.0, line="D2", temp_c=25.0, scan_points=401)
    params.update(overrides)
    return params


def test_natural_composition_preserves_legacy_spectrum_exactly():
    legacy = _params()
    del legacy["rb_composition"]
    del legacy["rb85_fraction"]
    original = SAS.compute(legacy)
    natural = SAS.compute(_params(rb85_fraction=0.2))
    for key in ("scan", "alpha_unit", "chi_real_unit"):
        np.testing.assert_array_equal(natural[key], original[key])
    assert natural["rb85_fraction"] == species.RB85.abundance
    assert natural["isotope_fractions"] == {"Rb85": 0.7217, "Rb87": 0.2783}
    assert natural["frequency_reference"] == "Rb85"


def test_composition_controls_recompute_and_keep_existing_species_menu():
    specs = {spec.name: spec for spec in SAS.param_schema()}
    assert specs["species"].choices == tuple(species.SPECIES_ORDER)
    assert specs["rb_composition"].default == "Natural"
    assert specs["rb_composition"].choices == ("Natural", "Custom")
    assert specs["rb_composition"].recompute
    assert specs["rb85_fraction"].recompute
    assert specs["rb85_fraction"].visible_if == {
        "species": "Rb (natural)", "rb_composition": "Custom"}


@pytest.mark.parametrize("fraction", [-0.1, 1.1, np.nan, np.inf, -np.inf])
def test_invalid_custom_fractions_fail_before_solving(fraction):
    with patch.object(SAS, "_component_response", side_effect=AssertionError("solve")):
        with pytest.raises(ValueError, match="fraction"):
            SAS.compute(_params(rb_composition="Custom", rb85_fraction=fraction))


@pytest.mark.parametrize("fraction,expected", [(0.0, "Rb87"), (1.0, "Rb85")])
def test_zero_share_components_are_not_built(fraction, expected):
    with patch("gabes.schemes.sas.species.build_manifold",
               wraps=species.build_manifold) as build:
        raw = SAS.compute(_params(rb_composition="Custom", rb85_fraction=fraction))
    assert len(build.call_args_list) == 1
    assert build.call_args.args[0].name == expected
    assert raw["isotope_fractions"] == {expected: 1.0}
    assert all(species.ISOTOPES[expected].label in label for _, label in raw["markers"])


@pytest.mark.parametrize("line", ["D1", "D2"])
@pytest.mark.parametrize("fraction,pure_key", [(0.0, "⁸⁷Rb"), (1.0, "⁸⁵Rb")])
def test_mixture_endpoints_match_pure_isotopes_with_declared_centroid(line, fraction, pure_key):
    mixture = SAS.compute(_params(line=line, rb_composition="Custom", rb85_fraction=fraction))
    pure = SAS.compute(_params(line=line, species=pure_key))
    frequency_offset = (species.ISOTOPES[pure["frequency_reference"]].line(line)[1]
                        - species.RB85.line(line)[1])
    np.testing.assert_allclose(mixture["scan"], pure["scan"] + 2 * np.pi * frequency_offset,
                               rtol=1e-12, atol=1e-3)
    np.testing.assert_allclose(mixture["alpha_unit"], pure["alpha_unit"], rtol=1e-8)
    assert mixture["frequency_reference"] == "Rb85"
    assert mixture["k_vec"] == 2 * np.pi * species.RB85.line(line)[1] / constants.C_LIGHT


def test_custom_fraction_assigns_complementary_elemental_density_shares():
    params = _params(rb_composition="Custom", rb85_fraction=0.35)
    densities = {}
    component_response = SAS._component_response

    def capture(component, scan, temperature):
        densities[component["iso"].name] = component["N"]
        return component_response(component, scan, temperature)

    with patch.object(SAS, "_component_response", side_effect=capture):
        raw = SAS.compute(params)
    total = species.number_density(species.RB85, params["temp_c"] + 273.15)
    assert densities == pytest.approx({"Rb85": 0.35 * total, "Rb87": 0.65 * total})
    assert sum(densities.values()) == pytest.approx(total)
    assert raw["isotope_fractions"] == {"Rb85": 0.35, "Rb87": 0.65}


def test_dominant_isotope_change_does_not_move_spectral_features():
    low = SAS.compute(_params(rb_composition="Custom", rb85_fraction=0.49))
    high = SAS.compute(_params(rb_composition="Custom", rb85_fraction=0.51))
    assert low["frequency_reference"] == high["frequency_reference"] == "Rb85"
    assert low["markers"] == high["markers"]
    np.testing.assert_array_equal(low["scan"], high["scan"])
    x = low["scan"] / (2 * np.pi * 1e9)
    # The far positive-frequency 87Rb manifold lies outside the 85Rb groups.
    rb87_peak = int(np.argmax(low["alpha_unit"] * (x > 3.0)))
    rb85_peak = int(np.argmax(low["alpha_unit"] * ((x > 0.5) & (x < 2.0))))
    assert high["alpha_unit"][rb87_peak] < low["alpha_unit"][rb87_peak]
    assert high["alpha_unit"][rb85_peak] > low["alpha_unit"][rb85_peak]


@pytest.mark.parametrize("fraction", [0.0, 0.35, 0.7217, 1.0])
def test_finite_probe_custom_composition_is_passive_and_qualified(fraction):
    params = _params(rb_composition="Custom", rb85_fraction=fraction,
                     probe_power_uw=10.0)
    raw = SAS.compute(params)
    assert raw["probe_model"] == "finite_power_2level_gaussian_propagation"
    assert np.all(np.isfinite(raw["probe_alpha_table_unit"]))
    assert np.all(raw["probe_alpha_table_unit"] >= 0.0)
    assert _probe_reference_matched(raw, params) == (fraction in (0.0, 0.7217, 1.0))
    readout = SAS.headless_observables(raw, params)
    composition = next(metric for metric in readout["metrics"]
                       if metric["label"] == "Rb isotope composition")
    assert f"{fraction * 100:.2f}%" in composition["value"]
    if fraction == 0.35:
        model = next(table["markdown"] for table in readout["tables"]
                     if table["title"] == "Finite-probe model")
        assert "outside the existing natural-" in model
