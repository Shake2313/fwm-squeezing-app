"""Lambda condition comparisons reuse native solves and existing numeric metrics."""
import json

import matplotlib.pyplot as plt
import numpy as np
import pytest

from gabes import schemes
from gabes.core import blas_single_thread
from gabes.schemes.lambda_comparison import _powers


def compare(mode, **changes):
    scheme = schemes.get("lambda")
    params = {**scheme.defaults(), **scheme.recommended_defaults({})[mode],
              "doppler": "off", **changes}
    resolved = []
    raws = {}
    with blas_single_thread():
        current_raw = scheme.compute(params)
        current_view = scheme.headless_observables(current_raw, params)

        def resolve(condition):
            power = condition["coupling_power_mw"]
            resolved.append(power)
            raws[power] = scheme.compute(condition)
            return raws[power]

        dataset = scheme.model_comparison(params, resolve_raw=resolve,
                                          current_raw=current_raw, current_view=current_view)
    raws[params["coupling_power_mw"]] = current_raw
    return scheme, params, dataset, resolved, raws


@pytest.mark.parametrize("mode", ["EIT", "AT", "CPT"])
def test_native_grids_and_existing_metric_values_are_preserved(mode, monkeypatch):
    monkeypatch.setattr(plt, "subplots", lambda *args, **kwargs: pytest.fail("Unexpected figure"))
    scheme, params, dataset, resolved, raws = compare(mode, cell_mm=1.0, temp_c=27.)
    assert len(resolved) == (2 if mode == "CPT" else 8)
    assert params["coupling_power_mw"] not in resolved
    assert dataset["varied"]["values"][0] == .25 and dataset["varied"]["values"][-1] == 4.
    assert dataset["fixed"]["cell_mm"] == 1.0
    assert dataset["fixed"]["coupling_diameter_mm"] == params["coupling_diameter_mm"]
    assert dataset["schema"] == "gabes-spectral-comparison-v1"
    for curve in dataset["series"]:
        power = curve["parameter"]["value"]
        view = scheme.headless_observables(raws[power], {**params, "coupling_power_mw": power})
        source = view["lambda_readout"]
        expected_x = source["x_mhz"]
        if mode == "CPT":
            expected_x = (expected_x - source["center_mhz"]) * 1000
        np.testing.assert_array_equal(curve["x"], expected_x)
        np.testing.assert_array_equal(curve["observables"]["transmission"], source["transmission"])
        np.testing.assert_array_equal(curve["observables"]["dispersion"], source["dispersion"])
        central = next(m["numeric_value"] for m in view["metrics"]
                       if m["label"] in ("Transmission at center", "Transmission at resonance"))
        assert next(m["value"] for m in curve["readouts"] if m["key"] == "central_transmission") == central
    assert not np.array_equal(dataset["series"][0]["x"], dataset["series"][-1]["x"])
    json.dumps(dataset, allow_nan=False)


def test_absorption_doublet_spacing_matches_existing_peak_extractor_and_input_reference():
    _, _, dataset, _, _ = compare("AT", coupling_rabi_mhz=15.)
    extracted = dataset["summary"]["series"][0]["observables"]["value"]
    reference = dataset["summary"]["series"][1]["observables"]["value"]
    np.testing.assert_allclose(extracted, reference, rtol=.03, atol=.1)
    _, _, weak, _, _ = compare("AT", coupling_rabi_mhz=.1)
    curve = next(c for c in weak["series"] if c["parameter"]["value"] == 1.)
    split = next(r for r in curve["readouts"] if r["key"] == "doublet_spacing")
    assert split["value"] is None and split["status"] == "unresolved"
    current_index = next(i for i,c in enumerate(weak["series"]) if c["parameter"]["value"] == 1.)
    assert weak["summary"]["series"][0]["observables"]["value"][current_index] is None


def test_cpt_two_photon_axis_subtracts_nonzero_coupling_detuning():
    _, _, data, _, _ = compare("CPT", coupling_detuning_mhz=2.7, coupling_rabi_mhz=.8)
    assert data["axis"]["unit"] == "kHz"
    assert "Δp − Δc" in data["axis"]["label"]
    for curve in data["series"]:
        assert abs(curve["x"][300]) < 1e-8
        width = next(r for r in curve["readouts"] if r["key"] == "window_width")
        assert width["unit"] == "kHz"
    assert data["summary"]["panels"][1]["unit"] == "kHz"


def test_zero_and_boundary_powers_preserve_current_and_finite_json():
    for power in (0., .01, 50.):
        specification = next(s for s in schemes.get("lambda").param_schema()
                             if s.name == "coupling_power_mw")
        levels = _powers(power, (specification.vmin, specification.vmax))
        assert power in levels and len(levels) == 3
        assert all(0 <= value <= 50 for value in levels)
    _, _, data, resolved, _ = compare("EIT", coupling_power_mw=0.)
    assert resolved == [.25, 1.]
    zero = next(c for c in data["series"] if c["parameter"]["value"] == 0)
    assert next(r["value"] for r in zero["readouts"] if r["key"] == "coupling_rabi") == 0
    json.dumps(data, allow_nan=False)
