"""New comparison displays retain the existing magneto and PSR physical readouts."""
import json

import numpy as np
import pytest

from gabes import schemes
from gabes.core import blas_single_thread
from gabes.psr import compute_psr_response
from gabes.schemes.magneto_unified import (
    FIELD_SCAN, OPTICAL_SCAN, ROTATION, TRANSMISSION, TRANSITION_COMPARISON,
)
from gabes.schemes.psr import CONDITIONAL, FULL_OBE, PHASE_VIEW, DETUNING_VIEW
from gabes.schemes.psr_rf import RF_MODEL, PHASE_VIEW as RF_PHASE_VIEW


def compare(scheme, params, raw=None):
    calls = []
    def solve(condition):
        calls.append(dict(condition))
        with blas_single_thread():
            return scheme.compute({key: condition[key] for key in scheme.recompute_keys()})
    if raw is None:
        raw = solve(params)
    view = scheme.headless_observables(raw, params)
    calls.clear()
    result = scheme.model_comparison(params, resolve_raw=solve, current_raw=raw, current_view=view)
    json.dumps(result, allow_nan=False)
    return result, calls, view


@pytest.fixture(scope="module")
def magneto_case():
    s=schemes.get("magneto")
    p={**s.defaults(), "scan_points":11, "doppler":"off", "b_max_ut":.01,
       "intensity_mw_cm2":.7, "b_offset_ut":.002, "residual_transverse_b_ut":.001}
    with blas_single_thread():
        r=s.compute(p)
    return s,p,r


def test_magneto_curves_and_photon_limit_are_existing_readouts(magneto_case):
    s,p,raw=magneto_case
    result,calls,view=compare(s,p,raw)
    assert [c["qwp_deg"] for c in calls] == [9., 18., 27., 36., 45.]
    assert all(c["intensity_mw_cm2"] == p["intensity_mw_cm2"] for c in calls)
    assert len(result["series"]) == 6
    assert result["varied"]["key"] == "qwp_deg"
    assert result["varied"]["values"] == [0., 9., 18., 27., 36., 45.]
    transmission_panel = next(panel for panel in result["panels"] if panel["key"] == "transmission")
    assert transmission_panel["stackable"] and transmission_panel["wide"]
    assert transmission_panel["stack_label"] == "QWP angle"
    assert transmission_panel["stack_scale_unit"] == "percentage points T"
    rotation_panel = next(panel for panel in result["panels"] if panel["key"] == "rotation")
    assert "defined only for linear input" in rotation_panel["note"]
    native=view["magneto_readout"]
    current=next(c for c in result["series"] if c["parameter"]["value"] == p["qwp_deg"])
    np.testing.assert_array_equal(current["x"],native["x"])
    np.testing.assert_array_equal(current["observables"]["transmission_pct"],
                                  np.asarray(native["transmission"])*100)
    np.testing.assert_array_equal(current["observables"]["rotation_mrad"],
                                  np.asarray(native["rotation_rad"])*1000)
    np.testing.assert_allclose(np.asarray(result["source"]["physical_b_ut"]),
                               np.asarray(result["source"]["command_b_ut"])+.002,rtol=0,atol=1e-18)
    assert result["markers"][1]["x"] == -.002
    assert "residual transverse" in " ".join(result["notes"])
    assert "photon shot-noise" in current["readouts"][-1]["reason"].lower()


def test_magneto_optical_axis_and_unsupported_rotation_use_native_scope():
    s=schemes.get("magneto")
    p={**s.defaults(),"scan_axis":OPTICAL_SCAN,"detuning_points":3,
       "detuning_max_mhz":10.,"doppler":"off","intensity_mw_cm2":.7,
       "qwp_deg":45.,"signal_type":"Transmission"}
    result,_,view=compare(s,p)
    assert result["axis"]["unit"] == "MHz"
    assert result["panels"][0]["key"] == "transmission"
    assert result["markers"] == []
    circular = next(c for c in result["series"] if c["parameter"]["value"] == 45.)
    assert all(v is None for v in circular["observables"]["rotation_mrad"])
    assert "unsupported" in circular["scope"]["rotation"].lower()
    assert circular["observable_status"]["rotation_mrad"] == [circular["scope"]["rotation"]]*len(circular["x"])
    assert all(v is not None for v in circular["observables"]["transmission_pct"])


def test_magneto_optical_optimum_reports_native_coordinate_and_search_scope():
    s=schemes.get("magneto")
    p={**s.defaults(),"scan_axis":OPTICAL_SCAN,"detuning_points":3,
       "detuning_max_mhz":10.,"doppler":"off","intensity_mw_cm2":.7,"b_offset_ut":.0002}
    result,_,view=compare(s,p)
    native=view["magneto_readout"]
    assert native["optimum_index"] is not None
    current=next(c for c in result["series"] if c["parameter"]["value"] == p["qwp_deg"])
    readings={r["key"]:r for r in current["readouts"]}
    assert readings["sensitivity_optical_detuning"]["value"] == native["optimum_detuning_mhz"]
    assert readings["sensitivity_optimum_status"]["value"] == native["optimum_status"]
    assert readings["photon_field_asd"]["value"] == native["optimum_field_asd_t_per_sqrt_hz"]*1e15
    assert "Best sampled" in readings["photon_field_asd"]["label"]
    assert current["scope"]["sensitivity_selection"]["at_scan_edge"] == native["optimum_at_scan_edge"]


@pytest.fixture(scope="module")
def conditional_case():
    s=schemes.get("psr")
    p={**s.defaults(),"calculation":CONDITIONAL,"view":PHASE_VIEW,"detuning_mhz":-900.,
       "doppler":False}
    raw=compute_psr_response([-1000.,-900.,-800.], doppler=False)
    return s,p,raw


@pytest.mark.parametrize("view_name",[PHASE_VIEW,DETUNING_VIEW])
def test_conditional_noise_assumptions_share_one_atomic_response(conditional_case,view_name):
    s,p,raw=conditional_case
    p={**p,"view":view_name}
    result,calls,view=compare(s,p,raw)
    assert calls == []
    assert result["varied"]["values"] == [0.,1.,2.]
    current=next(c for c in result["series"] if c["parameter"]["value"] == p["excess_noise_per_od"])
    native=view["psr_readout"]
    y=(native["phase_scan"]["noise_a_db"] if view_name==PHASE_VIEW
       else native["spectrum"]["minimum_db"])
    np.testing.assert_array_equal(current["observables"]["noise_a_db"],y)
    minimum=[c["readouts"][0]["value"] for c in result["series"]]
    assert minimum[0] < minimum[1] < minimum[2]
    assert result["axis"]["unit"] == ("°" if view_name==PHASE_VIEW else "MHz")
    assert result["panels"][0]["reference"] == 0


def test_full_obe_curves_preserve_all_sixteen_populations():
    s=schemes.get("psr")
    p={**s.defaults(),"calculation":FULL_OBE,"obe_intensity_mw_cm2":.1}
    result,calls,view=compare(s,p)
    assert [c["obe_intensity_mw_cm2"] for c in calls] == [.05,.2]
    assert len(result["axis"]["tick_labels"]) == 16
    for curve in result["series"]:
        assert len(curve["x"]) == 16
        assert sum(curve["observables"]["population"]) == pytest.approx(1,abs=1e-10)
        assert curve["readouts"][0]["unit"] == "m⁻¹"
        assert "jones_xy" in curve["scope"]["conditions"]


@pytest.mark.parametrize("phase",[False,True])
def test_rf_curves_reuse_native_spectrum_and_missing_phase(phase):
    s=schemes.get("psr")
    p={**s.defaults(),"calculation":RF_MODEL,"cell_mm":0.,
       "rf_view":RF_PHASE_VIEW if phase else "Analysis frequency"}
    result,calls,view=compare(s,p)
    assert calls == []
    assert result["series"][0]["readouts"][0]["value"] == pytest.approx(0,abs=1e-12)
    assert result["series"][0]["readouts"][-1]["value"] is None
    native=view["psr_readout"]["phase_scan" if phase else "spectrum"]
    key="noise_a_db" if phase else "minimum_db"
    np.testing.assert_allclose(result["series"][0]["observables"][key],native[key],rtol=0,atol=0)
    assert result["axis"]["unit"] == ("°" if phase else "MHz")


def test_unavailable_rf_is_gap_with_reason_instead_of_zero_db():
    s=schemes.get("psr")
    p={**s.defaults(),"calculation":RF_MODEL,"obe_intensity_mw_cm2":0.}
    result,calls,view=compare(s,p)
    curve=result["series"][0]
    assert calls == []
    assert not curve["scope"]["available"]
    assert curve["scope"]["reason"]
    assert all(v is None for v in curve["observables"]["minimum_db"])
    assert curve["readouts"][0]["value"] is None


@pytest.fixture(scope="module")
def conditional_screen_case():
    s = schemes.get("psr")
    p = {**s.defaults(), "calculation": CONDITIONAL, "view": DETUNING_VIEW,
         "detuning_mhz": -900., "doppler": False}
    with blas_single_thread():
        raw = compute_psr_response([-900., -412.5, 0., 412.5, 900.], doppler=False)
    return s, p, raw


def test_conditional_spectrum_retains_each_samples_native_applicability(conditional_screen_case):
    s, p, raw = conditional_screen_case
    result, calls, view = compare(s, p, raw)
    assert calls == []
    assert view["psr_readout"]["selected_point"]["approximation_valid"]
    expected = ["available", "available outside approximation screens", "available",
                "available outside approximation screens", "available"]
    for curve in result["series"]:
        assert curve["sample_status"] == expected
        assert curve["sample_reasons"] == raw["approximation_reasons"]
        native = s.headless_observables(
            raw, {**p, "excess_noise_per_od": curve["parameter"]["value"]})["psr_readout"]["spectrum"]
        np.testing.assert_array_equal(curve["observables"]["noise_a_db"], native["minimum_db"])
        np.testing.assert_array_equal(curve["observables"]["noise_b_db"], native["maximum_db"])
        assert curve["observables"]["noise_a_db"][1] is not None
        assert "uniform-pump screen" in " ".join(curve["sample_reasons"][1])


def test_conditional_phase_repeats_selected_points_native_applicability(conditional_screen_case):
    s, p, raw = conditional_screen_case
    p = {**p, "view": PHASE_VIEW, "detuning_mhz": -412.5}
    result, calls, view = compare(s, p, raw)
    assert calls == []
    point = view["psr_readout"]["selected_point"]
    assert point["available"] and not point["approximation_valid"]
    for curve in result["series"]:
        assert curve["sample_status"] == ["available outside approximation screens"] * len(curve["x"])
        assert curve["sample_reasons"] == [point["approximation_reasons"]] * len(curve["x"])
        assert all(value is not None for value in curve["observables"]["noise_a_db"])


def test_conditional_failed_atomic_response_keeps_gaps_and_sample_reasons(monkeypatch):
    s = schemes.get("psr")
    p = {**s.defaults(), "view": DETUNING_VIEW, "doppler": False, "detuning_mhz": -900.}
    def failed_solve(*args, **kwargs):
        raise np.linalg.LinAlgError("test failed solve")
    monkeypatch.setattr("gabes.psr.magneto_multilevel.multilevel_response", failed_solve)
    raw = compute_psr_response([-1000., -900., -800.], doppler=False)
    view = s.headless_observables(raw, p)
    result = s.model_comparison(p, resolve_raw=s.compute, current_raw=raw, current_view=view)
    for curve in result["series"]:
        assert curve["sample_status"] == ["unavailable"] * len(curve["x"])
        assert all("atomic solve failed: LinAlgError" in reasons for reasons in curve["sample_reasons"])
        assert all(value is None for value in curve["observables"]["noise_a_db"])
        assert all(value is None for value in curve["observables"]["noise_b_db"])



def test_magneto_angle_comparison_includes_arbitrary_current_angle_and_noise_toggle():
    s = schemes.get("magneto")
    p = {**s.defaults(), "qwp_deg": 17., "scan_points": 11, "doppler": "off", "show_shot_noise": False}
    result, calls, view = compare(s, p)
    assert result["varied"]["values"] == [0., 9., 17., 18., 27., 36., 45.]
    assert [c["qwp_deg"] for c in calls] == [0., 9., 18., 27., 36., 45.]
    assert [panel["key"] for panel in result["panels"]] == ["rotation", "transmission"]
    assert all(not any(item["key"].startswith("photon_") for item in curve["readouts"])
               for curve in result["series"])
    current = next(curve for curve in result["series"] if curve["parameter"]["value"] == 17.)
    assert current["observables"]["transmission_pct"] == pytest.approx(
        np.asarray(view["magneto_readout"]["transmission"])*100)
    assert "023416" in result["source"]["comparison_reference"]["doi"]


@pytest.mark.parametrize("axis", [FIELD_SCAN, OPTICAL_SCAN])
def test_magneto_transition_comparison_changes_native_branch_and_retains_source(axis):
    s = schemes.get("magneto")
    p = {**s.defaults(), "isotope": "Rb87", "scan_axis": axis, "comparison_mode": TRANSITION_COMPARISON,
         "scan_points": 11, "detuning_points": 3, "detuning_max_mhz": 10., "doppler": "off"}
    result, calls, view = compare(s, p)
    assert [condition["line"] for condition in calls] == ["D2"]
    assert result["varied"]["key"] == "line"
    assert result["varied"]["values"] == ["D1", "D2"]
    assert len(result["series"]) == 2
    for curve in result["series"]:
        condition = {**p, "line": curve["parameter"]["value"]}
        with blas_single_thread():
            raw = s.compute(condition)
        native = s.headless_observables(raw, condition)["magneto_readout"]
        assert curve["source"]["line"] == condition["line"]
        assert curve["source"]["reference_centroid_mhz"] == native["reference_centroid_mhz"]
        assert curve["x"] == native["x"]
        np.testing.assert_array_equal(curve["observables"]["transmission_pct"], np.asarray(native["transmission"])*100)
    assert any("complete hyperfine manifold" in note for note in result["notes"])
    assert not any(panel.get("stackable") for panel in result["panels"])
    assert [panel["key"] for panel in result["panels"]][-1] == "photon_field_asd"
    assert result["primary_panel"] == ("rotation" if axis == FIELD_SCAN else "transmission")


def test_magneto_legacy_current_view_regenerates_axis_selected_sensitivity(magneto_case):
    s, p, raw = magneto_case
    native_view = s.headless_observables(raw, p)
    stale_view = {**native_view, "magneto_readout": {**native_view["magneto_readout"], "signal_type": TRANSMISSION}}
    calls = []
    def solve(condition):
        calls.append(condition)
        with blas_single_thread():
            return s.compute(condition)
    result = s.model_comparison({**p, "signal_type": TRANSMISSION}, resolve_raw=solve,
                                current_raw=raw, current_view=stale_view)
    current = next(curve for curve in result["series"] if curve["parameter"]["value"] == p["qwp_deg"])
    native = native_view["magneto_readout"]["rotation_sensitivity"]["field_asd_t_per_sqrt_hz"]
    np.testing.assert_array_equal(current["observables"]["photon_field_asd_ft_per_sqrt_hz"], np.asarray(native)*1e15)
    assert result["fixed"]["signal_type"] == ROTATION
    assert len(calls) == 5
