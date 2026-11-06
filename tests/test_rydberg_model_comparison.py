"""Rydberg comparison preserves native spectra and separates fields from optical peaks."""
import json
import numpy as np
import pytest

from gabes.core import blas_single_thread
from gabes.schemes.rydberg import RydbergEITScheme


@pytest.fixture(scope="module")
def model():
    scheme = RydbergEITScheme()
    cache = {}
    actual = []

    def resolve(params):
        solve_params = {key:params[key] for key in scheme.recompute_keys()}
        key = tuple(sorted(solve_params.items()))
        if key not in cache:
            with blas_single_thread():
                cache[key] = scheme.compute(solve_params)
            actual.append(key)
        return cache[key]

    return scheme, resolve, actual


def test_eit_compensation_uses_one_solve_and_native_optical_readout(model):
    scheme, resolve, _ = model
    params = scheme.recommended_defaults(scheme.defaults())["EIT"]
    params["eit_compare"] = "Residual width"
    raw = resolve(params)

    def unexpected(_params):
        pytest.fail("Both compensation curves already exist in the current solve.")

    result = scheme.model_comparison(params, resolve_raw=unexpected, current_raw=raw)
    assert result["mode"] == "EIT"
    assert [s["parameter"]["value"] for s in result["series"]] == [0., 1.5]
    for curve, chi_key, expected_width in zip(
            result["series"], ("chi_bar", "chi_bar_uncomp"), (1.60349403255, 1.89749219860)):
        x, y, _ = scheme._transmission(raw[chi_key], raw, params)
        np.testing.assert_array_equal(curve["x"], x)
        np.testing.assert_array_equal(curve["observables"]["transmission"], y)
        assert curve["readouts"][0]["value"] == pytest.approx(expected_width, rel=1e-9)
        assert curve["readouts"][0]["status"] == "resolved"
    json.loads(json.dumps(result, allow_nan=False))


def test_rf_off_unresolved_and_resolved_are_distinct_without_interpolation(model):
    scheme, resolve, _ = model
    params = scheme.recommended_defaults(scheme.defaults())["AT electrometry"]
    current = resolve(params)
    calls = []

    def track(condition):
        calls.append(condition["lo_rabi_mhz"])
        return resolve(condition)

    with blas_single_thread():
        view = scheme.headless_observables(current, params)
        result = scheme.model_comparison(params, resolve_raw=track,
                                         current_raw=current, current_view=view)
    assert calls == [0., 1., 7.4]
    assert [s["parameter"]["value"] for s in result["series"]] == [0., 1., 3.7, 7.4]
    splits = [s["readouts"][0] for s in result["series"]]
    assert [s["status"] for s in splits] == ["not applicable", "unresolved", "resolved", "resolved"]
    assert splits[0]["value"] is None and splits[1]["value"] is None
    assert splits[2]["value"] == pytest.approx(3.1635, abs=1e-9)
    assert splits[3]["value"] == pytest.approx(6.882, abs=1e-9)
    curves = result["series"]
    assert len(curves[0]["x"]) == len(curves[3]["x"]) == 801
    assert curves[0]["x"][-1] != curves[3]["x"][-1]
    np.testing.assert_array_equal(curves[2]["x"], current["scan"]/(2*np.pi*1e6))
    assert result["source"]["interpolation"] == "none"
    assert curves[2]["readouts"][1]["value"] != curves[2]["readouts"][2]["value"]
    assert result["summary"]["series"][1]["observables"]["value"][:2] == [None, None]
    json.loads(json.dumps(result, allow_nan=False))


def test_peak_rms_reuses_atomic_cache_and_changes_field_only(model):
    scheme, resolve, actual = model
    params = scheme.recommended_defaults(scheme.defaults())["AT electrometry"]
    with blas_single_thread():
        rms = scheme.model_comparison(params, resolve_raw=resolve, current_raw=resolve(params))
        before = len(actual)
        peak_params = {**params, "rf_field_convention": "Peak"}
        peak = scheme.model_comparison(peak_params, resolve_raw=resolve, current_raw=resolve(peak_params))
    assert len(actual) == before
    for a, b in zip(rms["series"], peak["series"]):
        assert a["x"] == b["x"]
        assert a["observables"]["transmission"] == b["observables"]["transmission"]
        for name in ("if_response","if_phase","total_sensitivity","psn_sensitivity"):
            np.testing.assert_allclose(a["observables"][name],b["observables"][name],rtol=1e-12,atol=1e-12)
        assert a["readouts"][0] == b["readouts"][0]
        for index in (1, 2):
            expected = a["readouts"][index]["value"]
            assert b["readouts"][index]["value"] == (
                None if expected is None else pytest.approx(np.sqrt(2)*expected))
            assert b["readouts"][index]["unit"] == "mV/m Peak"


def test_nonzero_optical_detuning_axis_and_custom_rf_are_preserved(model):
    scheme, resolve, _ = model
    params = scheme.recommended_defaults(scheme.defaults())["AT electrometry"]
    params.update(scan_axis="Coupling detuning", probe_detuning_mhz=2., lo_rabi_mhz=5.25)
    with blas_single_thread():
        current = resolve(params)
        result = scheme.model_comparison(params, resolve_raw=resolve, current_raw=current)
    assert result["axis"]["label"] == "Coupling detuning"
    assert result["axis"]["reference"].endswith("-2 MHz")
    assert result["fixed"]["probe_detuning_mhz"] == 2.
    assert result["varied"]["values"] == [0., 1., 3.7, 5.25, 7.4]
    custom = next(s for s in result["series"] if s["parameter"]["value"] == 5.25)
    np.testing.assert_array_equal(custom["x"], current["scan"]/(2*np.pi*1e6))
    assert custom["readouts"][0]["status"] in ("unresolved", "resolved")
