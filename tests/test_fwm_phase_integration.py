"""FWM phase readouts reuse the selected gains without changing the atomic scan."""
from copy import deepcopy
import json

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from gabes.schemes import fwm


def phase_params(**changes):
    return dict(
        fwm.FWMScheme().defaults(),
        **dict({
            "tpd": -3.7,
            "seed_classical_rin_scaled": 0.0,
            "after_cell_eff_pct": 100.0,
            "qe_pct": 100.0,
            "jasperse_loss_pct": 0.0,
            "experimental_gain_calibration": False,
            "measured_gain_calibration": False,
            "squeezing_view": fwm.SQUEEZING_QUADRATURE_PHASE,
        }, **changes),
    )


def phase_raw(params=None):
    """Small physical lossless amplifier scan; no atomic solve in these tests."""
    params = phase_params() if params is None else params
    delta = np.array([-30.0, -10.0, 0.0, 10.0, 30.0])
    gs = np.array([2.0, 3.0, 4.0, 5.0, 6.0])
    gc = gs - 1.0
    transfer = np.empty((delta.size, 2, 2), dtype=complex)
    transfer[:, 0, 0] = transfer[:, 1, 1] = np.sqrt(gs)
    transfer[:, 0, 1] = transfer[:, 1, 0] = np.sqrt(gc)
    center = fwm.branch_center_GHz(params["opd"], -1)
    return {
        "D_GHz": params["opd"],
        "probe_axis_GHz": center + delta * 1e-3,
        "raman_center_minus_GHz": center,
        "raman_center_plus_GHz": fwm.branch_center_GHz(params["opd"], 1),
        "temperature_K": params["temp_c"] + 273.15,
        "input_pump_power_W": params["pump_mw"] * 1e-3,
        "input_seed_power_W": params["probe_uw"] * 1e-6,
        "cell_length_m": params["cell_mm"] * 1e-3,
        "w_pump_m": params["pump_waist_um"] * 1e-6,
        "w_probe_m": params["probe_waist_um"] * 1e-6,
        "pump_probe_angle_deg": params["seeded_angle_deg"],
        "branch": -1,
        "legacy_mean_field_fixture": False,
        "gain_closure": {"historical_fixture_requested": False},
        "T_canonical_small_signal": transfer,
        "G_s": gs,
        "G_c": gc,
        "gain_referred_noise_dB": np.full(delta.size, -5.0),
        "S_dB": np.full(delta.size, -5.0),
        "eta": 1.0,
        "N_atoms": 1e18,
        "sigma_v": 200.0,
        "n_velocity": 3,
        "Op_A_2pi_GHz": 0.5,
        "Os_2pi_MHz": 0.1,
        "phase_detail": fwm.PHASE_LEGACY,
        "physical_squeezing_dB": None,
        "claim_gate": {"physical_squeezing_prediction": False},
    }


def test_selector_changes_only_readout_and_preserves_atomic_compute_arguments(monkeypatch):
    scheme = fwm.FWMScheme()
    specs = {spec.name: spec for spec in scheme.param_schema()}
    selector = specs["squeezing_view"]
    assert selector.recompute is False
    assert "squeezing_view" not in scheme.recompute_keys()
    assert selector.visible_if["mode"] == fwm.MODE_SEEDED
    assert set(selector.choices) == {
        fwm.SQUEEZING_INTENSITY_DIFFERENCE, fwm.SQUEEZING_QUADRATURE_PHASE}
    assert scheme.defaults()["squeezing_view"] == fwm.SQUEEZING_INTENSITY_DIFFERENCE
    assert specs["noise_model"].visible_if["squeezing_view"] == fwm.SQUEEZING_INTENSITY_DIFFERENCE
    calls = []
    monkeypatch.setattr(fwm, "compute_spectrum", lambda *a, **kw: calls.append((a, kw)) or {})
    scheme.compute(phase_params(squeezing_view=fwm.SQUEEZING_INTENSITY_DIFFERENCE))
    scheme.compute(phase_params())
    assert calls[0] == calls[1]


def test_phase_switch_preserves_raw_and_the_saved_intensity_noise_model():
    scheme = fwm.FWMScheme()
    params = phase_params(noise_model=fwm.NOISE_PASSIVE_L1)
    raw = phase_raw(params)
    original = deepcopy(raw)
    id_params = dict(params, squeezing_view=fwm.SQUEEZING_INTENSITY_DIFFERENCE)
    before = scheme.headless_observables(raw, id_params)
    phase = scheme.headless_observables(raw, params)
    after = scheme.headless_observables(raw, id_params)
    assert phase["noise_model"] == fwm.NOISE_EXPERIMENTAL_JASPERSE
    assert phase["squeezing_view"] == fwm.SQUEEZING_QUADRATURE_PHASE
    assert phase["quadrature_phase_scan"]["available"]
    assert params["noise_model"] == fwm.NOISE_PASSIVE_L1
    assert after["noise_model"] == fwm.NOISE_PASSIVE_L1
    assert before["metrics"] == after["metrics"]
    np.testing.assert_array_equal(before["jasperse_readout"]["NF"], after["jasperse_readout"]["NF"])
    for key in ("G_s", "G_c", "T_canonical_small_signal", "S_dB", "probe_axis_GHz"):
        np.testing.assert_array_equal(raw[key], original[key])
    assert raw["claim_gate"] == original["claim_gate"]


def test_missing_view_in_older_saved_params_retains_intensity_difference():
    scheme = fwm.FWMScheme()
    params = phase_params(squeezing_view=fwm.SQUEEZING_INTENSITY_DIFFERENCE)
    raw = phase_raw(params)
    explicit = scheme.headless_observables(raw, params)
    legacy = dict(params)
    legacy.pop("squeezing_view")
    implicit = scheme.headless_observables(raw, legacy)
    assert implicit["squeezing_view"] == fwm.SQUEEZING_INTENSITY_DIFFERENCE
    assert explicit["metrics"] == implicit["metrics"]
    assert explicit["noise_operating_point"] == implicit["noise_operating_point"]
    assert implicit["metrics"][0]["label"] == "Squeezing [dB]"
    np.testing.assert_array_equal(explicit["jasperse_readout"]["NF"], implicit["jasperse_readout"]["NF"])


@pytest.mark.parametrize("calibrated", [False, True])
def test_phase_uses_the_selected_operating_gain_pair_including_calibration(calibrated):
    scheme = fwm.FWMScheme()
    params = phase_params(experimental_gain_calibration=calibrated)
    raw = phase_raw(params)
    view = scheme.headless_observables(raw, params)
    point = view["gain_operating_point"]
    if calibrated:
        assert view["experimental_gain_calibration"]["applied"]
        expected = view["experimental_gain_calibration"]["selected_point"]
    else:
        axis = (raw["probe_axis_GHz"] - raw["raman_center_minus_GHz"]) * 1e3
        expected = {key: np.interp(params["tpd"], axis, raw[key]) for key in ("G_s", "G_c")}
    for key in ("G_s", "G_c"):
        assert point[key] == pytest.approx(expected[key], rel=1e-12)
        assert np.ndim(point[key]) == 0  # LO phase does not create a gain scan.
    assert point["phase_independent"]
    assert point["two_photon_detuning_mhz"] == params["tpd"]
    assert view["quadrature_phase_scan"]["available"]
    if not calibrated:
        # The covariance is interpolated in linear units. Off-grid mixtures of
        # pure endpoint covariances need not themselves be pure squeezed states.
        endpoint_minima = (np.sqrt(raw["G_s"]) - np.sqrt(raw["G_c"])) ** 2
        expected_minimum = np.interp(params["tpd"], axis, endpoint_minima)
        assert view["quadrature_phase_scan"]["minimum_variance"] == pytest.approx(expected_minimum, rel=1e-11)
    else:
        # The selected calibrated point is inserted into the model's scan.
        assert "calibrated-gain approximation" in view["metrics"][0]["delta"]


def test_phase_plot_headless_metrics_and_strict_export_share_the_same_scan():
    import matplotlib.pyplot as plt
    from gabes_ui.export import build_export_payload
    from gabes_ui.liveplot import preview_data

    scheme = fwm.FWMScheme()
    params = phase_params()
    raw = phase_raw(params)
    headless = scheme.headless_observables(raw, params)
    view = scheme.observables(raw, params)
    try:
        scan = view["quadrature_phase_scan"]
        assert view["metrics"] == headless["metrics"]
        assert headless["figure"] is None
        assert len(view["figure"].axes) == 1
        axis = view["figure"].axes[0]
        assert "phase" in axis.get_xlabel().lower()
        assert "db" in axis.get_ylabel().lower()
        for index, key in enumerate(("noise_a_db", "noise_b_db")):
            np.testing.assert_array_equal(axis.lines[index].get_xdata(), scan["phase_degrees"])
            np.testing.assert_array_equal(axis.lines[index].get_ydata(), scan[key])
            np.testing.assert_array_equal(scan[key], headless["quadrature_phase_scan"][key])
        assert any(np.all(np.asarray(line.get_ydata()) == 0.0) for line in axis.lines[2:])
        assert scan["minimum_variance"] < 1.0 < scan["maximum_variance"]
        assert view["metrics"][0]["value"] == f"{scan['minimum_db']:.2f} dB"
        assert view["metrics"][1]["value"] == f"{scan['maximum_db']:.2f} dB"
        assert all(metric["tier"] == "hero" for metric in view["metrics"][:2])
        assert all(metric["tier"] == "secondary" for metric in view["metrics"][2:])
        np.testing.assert_allclose(scan["noise_a_db"], 10 * np.log10(scan["noise_a"]))
        np.testing.assert_allclose(scan["noise_b_db"], 10 * np.log10(scan["noise_b"]))
        table = view["tables"][0]["markdown"]
        assert "Common homodyne LO rotation" in table
        assert "Shot-noise normalization" in table
        assert "not a measured optical phase" in table
        assert "independent of LO phase" in table
        assert "shot-noise seed phase variance = 1" in table
        assert preview_data(view["figure"], "fwm", canonical_noise=True) == (None, [])
        payload = build_export_payload(
            scheme=scheme, raw=raw, params=params, view=view,
            plot_payloads=[], comparison_payload=None, cache_version=scheme.cache_version,
            readout_cache_version="phase-integration-test")
        decoded = json.loads(json.dumps(payload, allow_nan=False))
        exported = decoded["view"]["quadrature_phase_scan"]
        assert exported["available"]
        np.testing.assert_array_equal(exported["phase_degrees"], scan["phase_degrees"])
        np.testing.assert_array_equal(exported["noise_a_db"], scan["noise_a_db"])
        assert decoded["view"]["gain_operating_point"]["G_s"] == view["gain_operating_point"]["G_s"]
    finally:
        plt.close(view["figure"])


def test_changed_operating_detuning_updates_phase_envelope_and_gain_together():
    scheme = fwm.FWMScheme()
    params = phase_params()
    raw = phase_raw(params)
    first = scheme.headless_observables(raw, params)
    second = scheme.headless_observables(raw, dict(params, tpd=7.3))
    assert second["gain_operating_point"]["G_s"] > first["gain_operating_point"]["G_s"]
    assert second["quadrature_phase_scan"]["minimum_variance"] < first["quadrature_phase_scan"]["minimum_variance"]
    assert second["quadrature_phase_scan"]["maximum_variance"] > first["quadrature_phase_scan"]["maximum_variance"]


def test_phase_outside_uncalibrated_scan_is_unavailable_without_inventing_a_curve():
    scheme = fwm.FWMScheme()
    params = phase_params(tpd=40.0)
    view = scheme.headless_observables(phase_raw(params), params)
    scan = view["quadrature_phase_scan"]
    assert not scan["available"]
    assert scan["status"]
    assert all(metric["value"] == "unavailable" for metric in view["metrics"][:2])
    assert not np.any(np.isfinite(scan["noise_a_db"]))
    assert not np.any(np.isfinite(scan["noise_b_db"]))


@pytest.mark.parametrize("endpoint", [-30.0, 30.0])
def test_phase_accepts_physical_scan_endpoints_after_ghz_roundtrip(endpoint):
    scheme = fwm.FWMScheme()
    params = phase_params(tpd=endpoint)
    raw = phase_raw(params)
    view = scheme.headless_observables(raw, params)
    scan = view["quadrature_phase_scan"]
    assert scan["available"]
    index = 0 if endpoint < 0 else -1
    for key in ("G_s", "G_c"):
        assert view["gain_operating_point"][key] == pytest.approx(raw[key][index])
    ideal_minimum = (np.sqrt(raw["G_s"][index]) - np.sqrt(raw["G_c"][index])) ** 2
    assert scan["minimum_variance"] == pytest.approx(ideal_minimum, rel=1e-12)
    outside = scheme.headless_observables(raw, dict(params, tpd=endpoint + np.sign(endpoint) * 1e-6))
    assert not outside["quadrature_phase_scan"]["available"]
    assert np.isnan(outside["gain_operating_point"]["G_s"])


def test_missing_covariance_reports_unavailable_but_retains_operating_gain(monkeypatch):
    original_readout = fwm.experimental_jasperse_readout

    def readout_without_covariance(*args, **kwargs):
        readout = original_readout(*args, **kwargs)
        readout.pop("quadrature_covariance")
        return readout

    monkeypatch.setattr(fwm, "experimental_jasperse_readout", readout_without_covariance)
    params = phase_params()
    view = fwm.FWMScheme().headless_observables(phase_raw(params), params)
    scan = view["quadrature_phase_scan"]
    assert not scan["available"]
    assert scan["status"] == "EFFECTIVE_QUADRATURE_INPUT_UNAVAILABLE"
    assert all(metric["value"] == "unavailable" for metric in view["metrics"][:2])
    assert np.isfinite(view["gain_operating_point"]["G_s"])
    assert np.isfinite(view["gain_operating_point"]["G_c"])
    assert not np.any(np.isfinite(scan["noise_a_db"]))
