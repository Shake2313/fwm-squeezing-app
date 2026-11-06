"""Condition comparison shares solver samples, units and resolution diagnostics."""
import json
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import pytest

from gabes import schemes
from gabes.core import blas_single_thread
from gabes.model_comparison import sample_values, same_grid
from gabes.schemes.sas import GENERIC, RB_COMPOSITION_CUSTOM


def _compare(**changes):
    scheme = schemes.get("sas")
    params = {**scheme.defaults(), **changes}
    solves = []

    def resolve(physical_params):
        solves.append(physical_params["pump_power_mw"])
        return scheme.compute({k: physical_params[k] for k in scheme.recompute_keys()})

    with blas_single_thread():
        current = resolve(params)
        view = scheme.headless_observables(current, params)
        solves.clear()
        before = plt.get_fignums()
        comparison = scheme.model_comparison(
            params, current_raw=current, current_view=view, resolve_raw=resolve)
        assert before == plt.get_fignums()
    json.dumps(comparison, allow_nan=False)
    return comparison, solves, view


def test_full_resolution_natural_rb_reuses_current_and_differences_at_every_point():
    model, solves, current_view = _compare(scan_points=4001, pump_power_mw=2.)
    assert solves == [0.]
    x = np.array(model["x"])
    np.testing.assert_array_equal(x, current_view["optical_readout"]["x"] * 1000.)
    assert len(x) == model["fixed"]["scan_points"] == 4001
    assert "pump_power_mw" not in model["fixed"]
    assert model["axis"]["reference"] == "Rb85 D1 line centroid"
    assert model["axis"]["overview_unit"] == "GHz"
    off = np.asarray(model["series"][0]["transmission"])
    for s in model["series"]:
        np.testing.assert_array_equal(s["delta"], np.asarray(s["transmission"]) - off)
    current = model["series"][-1]
    np.testing.assert_array_equal(current["transmission"],
                                  current_view["optical_readout"]["transmission"])
    np.testing.assert_array_equal(current["optical_depth"],
                                  current_view["optical_readout"]["optical_depth"])
    assert not current["effective"]
    assert current["resolution"] == current_view["optical_readout"]["resolution"]
    assert model["series"][0]["resolution"]["width_mhz"] is None
    assert [g["label"] for g in model["groups"]] == [
        "⁸⁷Rb F=2", "⁸⁵Rb F=3", "⁸⁵Rb F=2", "⁸⁷Rb F=1"]
    for g in model["groups"]:
        lo, hi = g["window_mhz"]
        actual = [m["x_mhz"] for m in g["markers"] if m["kind"] == "transition"]
        candidate = [m["x_mhz"] for m in g["markers"] if m["kind"] == "crossover_candidate"]
        assert len(actual) == 2 and candidate == [sum(actual)/2]
        assert x[0] <= lo < min(actual) < max(actual) < hi <= x[-1]


def test_mixed_rb_composition_isotope_absence_and_native_reference_do_not_conflate():
    model, solves, _ = _compare(line="D2", rb_composition=RB_COMPOSITION_CUSTOM,
                                rb85_fraction=0., scan_points=401, pump_power_mw=.7)
    assert solves == [0.]
    assert model["default_power_mw"] == .7
    assert len(model["groups"]) == 2 and all("⁸⁷Rb" in g["label"] for g in model["groups"])
    assert model["axis"]["reference"].startswith("Rb85 D2")
    for group in model["groups"]:
        lines = [m["x_mhz"] for m in group["markers"] if m["kind"] == "transition"]
        mids = [m["x_mhz"] for m in group["markers"] if m["kind"] == "crossover_candidate"]
        assert len(lines) == 3
        assert mids == sorted((a+b)/2 for a,b in combinations(lines, 2))


@pytest.mark.parametrize("atom,reference,count,unit", [
    ("¹³³Cs", "Cs133 D1 line centroid", 2, "GHz"),
    (GENERIC, "Generic model zero detuning", 0, "MHz"),
])
def test_other_manifolds_use_own_grid_and_no_invented_rb_groups(atom, reference, count, unit):
    model, solves, current = _compare(species=atom, scan_points=401,
                                      temp_c=65., pump_power_mw=0.)
    assert solves == []
    assert model["default_power_mw"] == 0
    assert [s["power_mw"] for s in model["series"]] == [0.]
    assert model["axis"]["overview_unit"] == unit
    assert model["axis"]["reference"] == reference
    assert len(model["groups"]) == count and model["fixed"]["temp_c"] == 65.
    factor = 1000 if unit == "GHz" else 1
    np.testing.assert_array_equal(model["x"], current["optical_readout"]["x"] * factor)


def test_axis_mismatch_and_missing_values_cannot_be_subtracted_as_valid_samples():
    with pytest.raises(ValueError, match="same sampled axis"):
        same_grid([0, 1, 2], [0, 1.01, 2])
    assert sample_values([0, np.nan, 1.0, np.inf]) == [0., None, 1., None]


def test_native_depth_survives_underflow_and_carries_finite_probe_scope():
    opaque, _, _ = _compare(scan_points=401, temp_c=170., cell_mm=150.)
    c = opaque["series"][-1]
    deep = int(np.argmax(c["optical_depth"]))
    assert c["transmission"][deep] == 0
    assert np.isfinite(c["optical_depth"][deep]) and c["optical_depth"][deep] > 745.

    finite, solves, current = _compare(scan_points=401, line="D2", probe_power_uw=.1,
                                       pump_power_mw=.83)
    assert solves == [0.]
    assert finite["varied"]["values"] == [0., .83]
    assert all(c["effective"] for c in finite["series"])
    np.testing.assert_array_equal(finite["series"][-1]["optical_depth"],
                                  current["optical_readout"]["optical_depth"])


def test_unused_scheme_does_not_resolve_an_optional_comparison():
    def unexpected(*args):
        pytest.fail("Unexpected solver call")
    other = schemes.get("fwm")
    assert other.model_comparison_parameter is None
    assert other.model_comparison({}, resolve_raw=unexpected) is None
