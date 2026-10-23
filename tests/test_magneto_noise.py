"""Photon counting normalization and integration of the magneto noise overlay."""
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from gabes import constants, observables, schemes
from gabes.magneto_noise import shot_noise_reference
from gabes.schemes.magneto import SIGNAL_NMOR, SIGNAL_TRANSMISSION
from gabes_ui.export import build_export_payload, collect_view_plot_payloads
from gabes_ui.readout import partition_readout, strip_html


def _reference(**changes):
    settings = dict(reference_transmission=0.5, frequency_hz=377e12,
                    intensity_mw_cm2=0.8, beam_diameter_mm=2.0,
                    efficiency=0.8, bandwidth_hz=1000.0)
    settings.update(changes)
    return shot_noise_reference(np.array([0.25, 0.5, 1.0]), **settings)


def test_noise_matches_independent_poisson_counts_and_one_sided_enbw():
    data = _reference()
    # An average over tau has one-sided ENBW = 1/(2 tau). Independent
    # Poisson counts in that interval fix both the factor of two and units.
    tau = 1 / (2 * 1000.0)
    energy = 2 * np.pi * constants.HBAR * 377e12
    power = 0.8e-3 * np.pi * 0.1 ** 2  # mW/cm² -> W, diameter 2 mm
    n_incident = 0.8 * power * tau / energy
    n_detected = n_incident * np.array([0.25, 0.5, 1.0])
    assert data["incident_power_w"] == pytest.approx(power)
    np.testing.assert_allclose(data["transmission_rms"], np.sqrt(n_detected) / n_incident)
    np.testing.assert_allclose(data["rotation_rms_rad"], 1 / (2 * np.sqrt(n_detected)))
    assert data["transmission_reference_rms"] == pytest.approx(
        np.sqrt(n_incident * 0.5) / n_incident)


@pytest.mark.parametrize("change, factor", [
    ({"bandwidth_hz": 4000}, 2),
    ({"efficiency": 0.2}, 2),
    ({"intensity_mw_cm2": 3.2}, 0.5),
    ({"beam_diameter_mm": 4.0}, 0.5),
])
def test_noise_scales_with_bandwidth_and_detected_power(change, factor):
    base, changed = _reference(), _reference(**change)
    for key in ("transmission_rms", "rotation_rms_rad"):
        np.testing.assert_allclose(changed[key], np.asarray(base[key]) * factor)
    for key in ("transmission_reference_trace", "rotation_reference_trace_rad"):
        center = 0.5 if key.startswith("transmission") else 0.0
        np.testing.assert_allclose(np.asarray(changed[key]) - center,
                                   (np.asarray(base[key]) - center) * factor)


def test_gaussian_reference_has_declared_rms_and_is_reproducible():
    data = shot_noise_reference(np.ones(100000), reference_transmission=1.0,
                                frequency_hz=377e12, intensity_mw_cm2=0.8)
    t_noise = np.asarray(data["transmission_reference_trace"]) - 1.0
    theta_noise = np.asarray(data["rotation_reference_trace_rad"])
    for noise, sigma in [(t_noise, data["transmission_reference_rms"]),
                         (theta_noise, data["rotation_rms_rad"][0])]:
        assert abs(np.mean(noise)) < sigma * 0.01
        assert np.std(noise) == pytest.approx(sigma, rel=0.01)
    assert _reference() == _reference()


def test_zero_light_and_unsupported_polarization_do_not_claim_zero_rotation_noise():
    dark = shot_noise_reference(np.zeros(3), reference_transmission=0,
                                frequency_hz=377e12, intensity_mw_cm2=0.8)
    assert dark["rotation_rms_rad"] == [None] * 3
    assert dark["rotation_reference_trace_rad"] == [None] * 3
    assert dark["transmission_reference_trace"] == [None] * 3
    assert "insufficient detected photons" in dark["rotation_status"]
    elliptical = _reference(linear_input=False)
    assert elliptical["rotation_reference_trace_rad"] == [None] * 3
    assert elliptical["rotation_status"] == "unsupported polarization"
    assert elliptical["transmission_status"] == "available"


@pytest.mark.parametrize("settings", [{"intensity_mw_cm2": 0.0}, {"efficiency": 0.0}])
def test_zero_incident_photon_budget_is_unavailable_without_division_by_zero(settings):
    with np.errstate(divide="raise", invalid="raise"):
        data = _reference(**settings)
    assert data["transmission_reference_rms"] is None
    for key in ("transmission_rms", "transmission_reference_trace",
                "rotation_rms_rad", "rotation_reference_trace_rad"):
        assert data[key] == [None] * 3
    assert "insufficient detected photons" in data["transmission_status"]


def test_opaque_scan_samples_are_gaps_without_overflow():
    with np.errstate(over="raise", divide="raise", invalid="raise"):
        data = shot_noise_reference(np.array([1e-320, 0.5, 1.0]),
                                    reference_transmission=0.5,
                                    frequency_hz=377e12, intensity_mw_cm2=0.8)
    assert data["rotation_rms_rad"][0] is None
    assert data["rotation_reference_trace_rad"][0] is None
    assert data["rotation_rms_rad"][1] > 0
    assert data["rotation_status"].startswith("partial:")


@pytest.fixture(scope="module")
def magneto_result():
    scheme = schemes.get("magneto")
    params = dict(scheme.defaults(), doppler="off", scan_points=21)
    return scheme, params, scheme.compute(params)


def _controlled_scan(magneto_result, transmission, rotation):
    """Three known detector signals, independent of the atomic line shape."""
    scheme, params, raw = magneto_result
    params = dict(params, b_offset_ut=0.25)
    transmission, rotation = np.asarray(transmission), np.asarray(rotation)
    length = params["cell_mm"] * 1e-3
    coupling = observables.chi_phys(1.0, raw["N_eff"], dipole=raw["dipole"],
                                     line_strength=params["line_strength"])
    command = np.array([-1.0, 0.0, 1.0])
    return scheme, params, dict(
        raw, b_ut=command, b_physical_ut=command + 0.25, b_offset_ut=0.25,
        chi_probe=-1j * np.log(transmission) / (length * raw["k_vec"] * coupling),
        chi_p=(rotation / (0.25 * raw["k_vec"] * length * coupling)).astype(complex),
        chi_m=np.zeros(3, dtype=complex))


def test_transmission_snr_is_contrast_over_background_noise(magneto_result):
    scheme, params, raw = _controlled_scan(magneto_result, [0.8, 0.7, 0.8], [0, 0, 0])
    view = scheme.headless_observables(raw, params)
    data = view["magneto_readout"]
    snr = data["shot_noise_snr"]
    # Symmetric physical wings fall at command -1 and +0.5: background 0.775.
    background = 0.775
    sigma = data["shot_noise"]["transmission_reference_rms"]
    assert snr["reference_level"] == pytest.approx(background)
    assert snr["peak_amplitude_ratio"] == pytest.approx(0.075 / sigma)
    assert snr["peak_signal_amplitude"] == pytest.approx(0.075)
    assert snr["peak_command_field_ut"] == 0.0
    assert snr["peak_physical_field_ut"] == 0.25
    # The large DC transmission must not be counted as resonance signal.
    assert snr["peak_amplitude_ratio"] < 0.7 / sigma / 5


def test_rotation_peak_snr_uses_noise_at_the_same_field(magneto_result):
    scheme, params, raw = _controlled_scan(magneto_result, [0.01, 0.5, 1.0], [0.004, 0, 0.003])
    params["signal_type"] = SIGNAL_NMOR
    view = scheme.headless_observables(raw, params)
    data = view["magneto_readout"]
    snr = data["shot_noise_snr"]
    # The largest angle is at -1 uT, but its photon flux is 100 times lower.
    # The correct largest SNR is at +1 uT, not at the largest rotation angle.
    energy = constants.HBAR * raw["k_vec"] * constants.C_LIGHT
    power = params["intensity_mw_cm2"] * 10 * np.pi * (params["beam_diameter_mm"] * 1e-3 / 2) ** 2
    sigma = np.sqrt(energy * params["noise_bandwidth_hz"] /
                    (2 * params["detection_efficiency"] * power))
    assert snr["peak_amplitude_ratio"] == pytest.approx(0.003 / sigma)
    assert snr["peak_command_field_ut"] == 1.0
    assert snr["peak_physical_field_ut"] == 1.25
    assert snr["peak_noise_rms"] == pytest.approx(sigma)
    assert snr["signal_unit"] == "rad"


@pytest.mark.parametrize("factor, visibility", [(0.0, "Below shot noise"),
                                                 (0.5, "Below shot noise"),
                                                 (8.0, "Above shot noise")])
def test_snr_labels_compare_amplitude_with_one_rms(magneto_result, factor, visibility):
    scheme, params, raw = _controlled_scan(magneto_result, [0.5] * 3, [0, 0, 0])
    params["signal_type"] = SIGNAL_NMOR
    baseline = scheme.headless_observables(raw, params)
    sigma = baseline["magneto_readout"]["shot_noise"]["rotation_rms_rad"][1]
    scheme, params, raw = _controlled_scan(magneto_result, [0.5] * 3, [0, sigma * factor, 0])
    params["signal_type"] = SIGNAL_NMOR
    view = scheme.headless_observables(raw, params)
    snr = view["magneto_readout"]["shot_noise_snr"]
    assert snr["peak_amplitude_ratio"] == pytest.approx(factor)
    assert snr["visibility"] == visibility
    metric = next(m for m in view["metrics"] if m["label"] == "Peak shot-noise SNR")
    assert metric["delta"] == visibility


@pytest.mark.parametrize("signal", [SIGNAL_TRANSMISSION, SIGNAL_NMOR])
def test_snr_replaces_command_zero_in_visible_strip_and_keeps_its_data(magneto_result, signal):
    scheme, params, raw = magneto_result
    view = scheme.headless_observables(raw, dict(params, signal_type=signal))
    heroes, secondary, statuses, overflow = partition_readout(view["metrics"], view["hero_count"])
    assert any(m["label"] == "Peak shot-noise SNR" for m in heroes + secondary)
    assert not any("command B=0" in m["label"] for m in view["metrics"])
    assert len(heroes + secondary) <= 5
    assert not overflow
    html = strip_html(heroes, secondary, statuses)
    assert "Peak shot-noise SNR" in html and "Above shot noise" in html
    assert "command B=0" not in html
    assert view["magneto_readout"]["command_zero_transmission"] is not None
    derived = next(t["markdown"] for t in view["tables"] if t["title"] == "Derived quantities")
    assert "Transmission at command B=0" in derived and "Rotation at command B=0" in derived


@pytest.mark.parametrize("signal", [SIGNAL_TRANSMISSION, SIGNAL_NMOR])
def test_no_detected_light_has_unavailable_snr(magneto_result, signal):
    scheme, params, raw = magneto_result
    coupling = observables.chi_phys(1.0, raw["N_eff"], dipole=raw["dipole"],
                                     line_strength=params["line_strength"])
    opaque = dict(raw, chi_probe=np.full(raw["b_ut"].shape, 1e6j / coupling))
    view = scheme.headless_observables(opaque, dict(params, signal_type=signal))
    snr = view["magneto_readout"]["shot_noise_snr"]
    assert snr["peak_amplitude_ratio"] is None
    assert snr["valid_samples"] == 0
    assert all(value is None for value in snr["amplitude_ratio"])
    metric = next(m for m in view["metrics"] if m["label"] == "Peak shot-noise SNR")
    assert metric["value"] == "unavailable" and metric["delta"] == "Unavailable"
    assert metric["kind"] == "status"


@pytest.mark.parametrize("signal", [SIGNAL_TRANSMISSION, SIGNAL_NMOR])
def test_zero_drive_preserves_atomic_off_state_and_withholds_noise_metrics(magneto_result, signal):
    scheme, params, _ = magneto_result
    params = dict(params, intensity_mw_cm2=0.0, signal_type=signal)
    raw = scheme.compute(params)
    np.testing.assert_array_equal(raw["chi_probe"], 0.0)
    with np.errstate(divide="raise", invalid="raise"):
        view = scheme.headless_observables(raw, params)
    assert view["magneto_readout"]["shot_noise_snr"]["peak_amplitude_ratio"] is None
    assert view["magneto_readout"]["shot_noise"]["incident_power_w"] == 0
    heroes, secondary, _, _ = partition_readout(view["metrics"], view["hero_count"])
    assert all(m["value"] != "unavailable" for m in heroes + secondary)


@pytest.mark.parametrize("signal", [SIGNAL_TRANSMISSION, SIGNAL_NMOR])
def test_overlay_uses_axis_units_and_exports_the_displayed_trace(magneto_result, signal):
    scheme, params, raw = magneto_result
    params = dict(params, signal_type=signal)
    view = scheme.observables(raw, params)
    try:
        axis = view["figure"].axes[0]
        noise = view["magneto_readout"]["shot_noise"]
        line = next(line for line in axis.lines if line.get_label() == "Shot noise (simulated)")
        rotation = signal == SIGNAL_NMOR
        key = "rotation_reference_trace_rad" if rotation else "transmission_reference_trace"
        np.testing.assert_allclose(line.get_xdata(), raw["b_ut"])
        np.testing.assert_allclose(line.get_ydata(), np.asarray(noise[key]) * (1e3 if rotation else 1))
        assert axis.get_ylabel().endswith("[mrad]") if rotation else axis.get_ylabel() == "Transmission"
        plots = collect_view_plot_payloads(view, view["figure"])
        exported = next(trace for trace in plots[0]["traces"]
                        if trace["figure_curve_label"] == "Shot noise (simulated)")
        np.testing.assert_allclose(exported["y_data"], line.get_ydata())
        assert exported["y_unit"] == ("mrad" if rotation else "1")
        payload = build_export_payload(
            scheme=scheme, raw=raw, params=params, view=view, plot_payloads=plots,
            comparison_payload=None, cache_version=scheme.cache_version, readout_cache_version="test")
        assert payload["view"]["magneto_readout"]["shot_noise"] == noise
        assert payload["view"]["magneto_readout"]["shot_noise_snr"] == view["magneto_readout"]["shot_noise_snr"]
    finally:
        plt.close(view["figure"])


def test_noise_controls_render_toggle_and_advanced_settings():
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string(
        "import streamlit as st\n"
        "from gabes import schemes\n"
        "from gabes_ui.controls import render_rail\n"
        "scheme = schemes.get('magneto')\n"
        "for name, value in scheme.defaults().items():\n"
        "    st.session_state.setdefault(f'magneto__{name}', value)\n"
        "names = {'show_shot_noise', 'beam_diameter_mm', 'detection_efficiency', 'noise_bandwidth_hz'}\n"
        "specs = [sp for sp in scheme.param_schema() if sp.name in names]\n"
        "params = render_rail(st, scheme, specs)\n"
        "st.session_state['noise_controls'] = params\n"
    ).run()
    assert not app.exception
    assert app.checkbox(key="magneto__show_shot_noise").value is True
    app.checkbox(key="magneto__show_shot_noise").uncheck().run()
    assert app.session_state["noise_controls"]["show_shot_noise"] is False
    app.button(key="magneto___advanced_toggle").click().run()
    assert not app.exception
    for name in ("beam_diameter_mm", "detection_efficiency", "noise_bandwidth_hz"):
        assert app.slider(key=f"magneto__{name}").value == schemes.get("magneto").defaults()[name]
    app.slider(key="magneto__noise_bandwidth_hz").set_value(4000.0).run()
    assert not app.exception
    assert app.session_state["noise_controls"]["noise_bandwidth_hz"] == 4000.0


def test_readout_controls_reuse_atomic_scan_and_toggle_preserves_signal(magneto_result):
    scheme, params, raw = magneto_result
    controls = {"show_shot_noise", "beam_diameter_mm", "detection_efficiency", "noise_bandwidth_hz"}
    assert controls.isdisjoint(scheme.recompute_keys())
    with patch.object(scheme, "compute", side_effect=AssertionError("unexpected atomic solve")):
        on = scheme.observables(raw, params)
        off = scheme.observables(raw, dict(params, show_shot_noise=False))
        wider = scheme.headless_observables(raw, dict(params, noise_bandwidth_hz=4000))
    try:
        np.testing.assert_array_equal(on["figure"].axes[0].lines[0].get_ydata(),
                                      off["figure"].axes[0].lines[0].get_ydata())
        assert not any(line.get_label() == "Shot noise (simulated)"
                       for line in off["figure"].axes[0].lines)
        assert wider["magneto_readout"]["shot_noise"]["transmission_reference_rms"] == pytest.approx(
            2 * on["magneto_readout"]["shot_noise"]["transmission_reference_rms"])
        assert wider["magneto_readout"]["shot_noise_snr"]["peak_amplitude_ratio"] == pytest.approx(
            0.5 * on["magneto_readout"]["shot_noise_snr"]["peak_amplitude_ratio"])
        assert on["magneto_readout"]["shot_noise_snr"] == off["magneto_readout"]["shot_noise_snr"]
    finally:
        plt.close(on["figure"])
        plt.close(off["figure"])


@pytest.mark.parametrize("qwp", [20.0, 45.0])
def test_rotation_overlay_is_withheld_for_elliptical_or_circular_input(magneto_result, qwp):
    scheme, params, raw = magneto_result
    view = scheme.observables(dict(raw, qwp_deg=qwp), dict(params, signal_type=SIGNAL_NMOR))
    try:
        assert view["magneto_readout"]["shot_noise"]["rotation_status"] == "unsupported polarization"
        assert view["magneto_readout"]["shot_noise_snr"]["peak_amplitude_ratio"] is None
        assert view["magneto_readout"]["shot_noise_snr"]["valid_samples"] == 0
        heroes, secondary, statuses, _ = partition_readout(view["metrics"], view["hero_count"])
        assert all(m["value"] != "unavailable" for m in heroes + secondary)
        assert any(m["label"] == "Peak shot-noise SNR" for m in statuses)
        assert not any(line.get_label() == "Shot noise (simulated)"
                       for line in view["figure"].axes[0].lines)
    finally:
        plt.close(view["figure"])
