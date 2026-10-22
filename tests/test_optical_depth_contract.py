"""Separate natural-log optical depth from the legacy decadic OD contract."""

from __future__ import annotations

import csv
from io import StringIO
import pickle
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from gabes import observables, probe_saturation, schemes, species  # noqa: E402
from gabes.schemes.sas import GENERIC  # noqa: E402
from gabes_ui.export import (  # noqa: E402
    _curve_csv_bytes, build_export_payload, collect_view_plot_payloads,
    to_json_primitive,
)

SAS = schemes.get("sas")


def test_optical_helpers_keep_both_log_conventions_and_legacy_values():
    # Independent known transmission: exp(-2), with alpha in inverse metres.
    alpha = np.array([0.0, 40.0, -20.0])
    length = 0.05
    tau = observables.optical_depth(alpha, length)
    a10 = observables.decadic_absorbance(alpha, length)
    np.testing.assert_array_equal(tau, [0.0, 2.0, -1.0])
    np.testing.assert_array_equal(a10, alpha * length / np.log(10.0))
    np.testing.assert_array_equal(observables.optical_density(alpha, length), a10)
    transmission = observables.transmission(alpha, length)
    np.testing.assert_allclose(tau, -np.log(transmission), atol=2e-16)
    np.testing.assert_allclose(a10, -np.log10(transmission), atol=2e-16)
    assert observables.optical_depth(40.0, length) == 2.0
    assert observables.optical_density(40.0, length) == 2.0 / np.log(10.0)


def test_log_readouts_remain_finite_when_transmission_underflows():
    alpha, length = np.array([0.0, 1000.0, 10000.0]), 1.0
    transmission = observables.transmission(alpha, length)
    assert transmission.tolist() == [1.0, 0.0, 0.0]
    tau = observables.optical_depth(alpha, length)
    a10 = observables.decadic_absorbance(alpha, length)
    assert np.all(np.isfinite(tau)) and np.all(np.isfinite(a10))
    np.testing.assert_array_equal(tau, alpha)
    np.testing.assert_array_equal(observables.optical_density(alpha, length), a10)


@pytest.fixture(scope="module", params=[
    (species.SPECIES_ORDER[1], 0.0, 0.0),
    (species.SPECIES_ORDER[1], 0.5, 0.0),
    (GENERIC, 0.0, 0.0),
    (GENERIC, 0.5, 0.0),
    ("Rb (natural)", 0.0, 40.0),
], ids=["species-OD", "species-SAS", "generic-OD", "generic-SAS", "finite-probe"])
def attenuation_case(request):
    medium, pump, probe = request.param
    params = SAS.defaults()
    params.update(species=medium, line="D2", pump_power_mw=pump,
                  probe_power_uw=probe, temp_c=30.0, scan_points=401,
                  cell_mm=50.0, line_strength=1.0)
    return SAS.compute(params), params


def test_live_length_and_line_strength_update_both_logs_without_raw_mutation(
        attenuation_case):
    raw, original_params = attenuation_case
    raw_before = pickle.dumps(raw, protocol=5)
    finite = raw.get("probe_alpha_table_unit") is not None
    for length_mm, strength in [(25.0, 0.6), (80.0, 1.3)]:
        params = dict(original_params, cell_mm=length_mm, line_strength=strength)
        view = SAS.headless_observables(raw, params)
        readout = view["optical_readout"]
        length = length_mm * 1e-3
        weak_alpha = raw["alpha_unit"] * strength
        if finite:
            expected_t, expected_alpha = probe_saturation.propagate_gaussian_spectrum(
                raw["probe_intensity_log_grid"],
                raw["probe_alpha_table_unit"] * strength, weak_alpha, length)
            # Collected nonlinear attenuation differs from the unsaturated one.
            peak = np.argmax(weak_alpha)
            assert expected_alpha[peak] < weak_alpha[peak]
        else:
            expected_alpha = weak_alpha
            expected_t = np.exp(-weak_alpha * length)
        expected_tau = expected_alpha * length
        np.testing.assert_array_equal(readout["optical_depth"], expected_tau)
        np.testing.assert_array_equal(readout["decadic_absorbance"],
                                      expected_tau / np.log(10.0))
        np.testing.assert_allclose(readout["transmission"], expected_t,
                                   rtol=1e-14, atol=0.0)
        np.testing.assert_allclose(readout["optical_depth"],
                                   -np.log(readout["transmission"]), atol=2e-16)
        assert readout["effective"] is finite
        assert readout["cell_length_m"] == length
        assert readout["line_strength"] == strength
        assert readout["unit"] == "dimensionless"
        assert readout["x_unit"] == ("MHz" if medium_is_generic(raw) else "GHz")
        expected_x = raw["scan"] / (2 * np.pi) / (
            1e6 if medium_is_generic(raw) else 1e9)
        np.testing.assert_array_equal(readout["x"], expected_x)
        if finite:
            assert "Collected Gaussian power" in readout["attenuation_scope"]
        else:
            # The added readout leaves the existing phase convention intact.
            phase = observables.single_pass_phase(
                raw["chi_real_unit"] * strength, raw["k_vec"], length)
            phase_metric = next(m for m in view["metrics"]
                                if m["label"] == "Peak phase shift")
            assert float(phase_metric["value"].split()[0]) == pytest.approx(
                np.max(np.abs(phase)) * (
                    1e3 if phase_metric["value"].endswith("mrad") else 1.0),
                abs=0.05)
    assert pickle.dumps(raw, protocol=5) == raw_before


def medium_is_generic(raw):
    return raw["mode"] == "generic"


def test_figures_headless_json_and_csv_share_explicit_quantity_contract(
        attenuation_case):
    raw, params = attenuation_case
    try:
        view = SAS.observables(raw, params)
        headless = SAS.headless_observables(raw, params)
        readout = view["optical_readout"]
        for quantity in ("transmission", "decadic_absorbance", "optical_depth"):
            np.testing.assert_array_equal(readout[quantity],
                                          headless["optical_readout"][quantity])
        assert view["comparison"] == headless["comparison"]
        assert view["comparison"]["axis_index"] == 0
        np.testing.assert_array_equal(view["figure"].axes[0].lines[0].get_ydata(),
                                      readout["transmission"])
        assert view["figure_views"][0]["label"] == "Transmission"
        assert view["figure_views"][1]["label"] == "Decadic absorbance A10"
        assert view["figure_views"][2]["label"] == "Optical depth tau"
        plots = collect_view_plot_payloads(view, view["figure"])
        quantity_traces = {
            trace["quantity_metadata"]["quantity"]: trace
            for plot in plots for trace in plot["traces"]
            if "quantity_metadata" in trace
        }
        for quantity, equation in (
                ("decadic_absorbance", "Decadic absorbance A10 = -log10(T)"),
                ("optical_depth", "Optical depth tau = -ln(T)")):
            trace = quantity_traces[quantity]
            assert trace["quantity_metadata"]["definition"] == (
                readout["conventions"][quantity])
            assert trace["quantity_metadata"]["effective"] == readout["effective"]
            assert trace["y_unit"] == "dimensionless"
            assert trace["x_unit"] == readout["x_unit"]
            rows = list(csv.reader(StringIO(_curve_csv_bytes(trace).decode("utf-8"))))
            assert rows[0][1] == f"{equation} [dimensionless]"
            assert len(rows[0]) == 2
            np.testing.assert_array_equal(
                np.asarray(rows[1:], dtype=float)[:, 1], readout[quantity])
        payload = build_export_payload(
            scheme=SAS, raw=raw, params=params, view=view, plot_payloads=plots,
            comparison_payload=None, cache_version=SAS.cache_version,
            readout_cache_version="test-optical-contract")
        assert payload["view"]["optical_readout"] == to_json_primitive(readout)
        assert payload["raw"] == to_json_primitive(raw)
        if not readout["effective"]:
            dispersion = next(item for item in view["figure_views"]
                              if item["label"] == "Dispersion")
            np.testing.assert_array_equal(
                dispersion["figure"].axes[0].lines[0].get_ydata(),
                (observables.refractive_index(
                    raw["chi_real_unit"] * params["line_strength"]) - 1.0) * 1e6)
    finally:
        plt.close("all")
