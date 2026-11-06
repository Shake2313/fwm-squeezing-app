"""CSV downloads must contain the native curves displayed in the model workspace."""
import csv
from io import BytesIO, StringIO
import zipfile

import numpy as np
import pytest

from gabes import schemes
from gabes.core import blas_single_thread
from gabes_ui.export import _build_plot_zip_bytes
from gabes_ui.model_compare_section import comparison_plot_payloads


@pytest.mark.parametrize("power", [0., .83])
def test_comparison_csv_exports_current_and_baseline_samples_without_marker_traces(power):
    scheme = schemes.get("sas")
    params = {**scheme.defaults(), "pump_power_mw": power, "scan_points": 401}
    with blas_single_thread():
        raw = scheme.compute({key: params[key] for key in scheme.recompute_keys()})
        view = scheme.headless_observables(raw, params)
        model = scheme.model_comparison(
            params, current_raw=raw, current_view=view,
            resolve_raw=lambda values: scheme.compute({key: values[key] for key in scheme.recompute_keys()}))
    plots = comparison_plot_payloads(model)
    assert [p["label"] for p in plots] == ["Transmission T", "Optical depth OD"]
    assert all(p["source"] == "model_comparison" for p in plots)
    assert all(len(p["traces"]) == len(model["series"]) for p in plots)
    archive = zipfile.ZipFile(BytesIO(_build_plot_zip_bytes(plots)))
    assert len(archive.namelist()) == 2 * len(model["series"])
    assert not any("line-" in name or "Dispersion" in name or "absorbance" in name
                   for name in archive.namelist())
    for plot, key in zip(plots, ("transmission", "optical_depth")):
        for trace, condition in zip(plot["traces"], model["series"]):
            assert trace["quantity_metadata"]["power_mw"] == condition["power_mw"]
            assert trace["quantity_metadata"]["frequency_reference"] == model["axis"]["reference"]
            np.testing.assert_array_equal(trace["x_data"], model["x"])
            np.testing.assert_array_equal(trace["y_data"], condition[key])
    # Parse actual ZIP bytes, including the native MHz axis and both full curves.
    for name, key in (("Transmission-T/Pump-off-0-mW.csv", "transmission"),
                      ("Optical-depth-OD/Pump-off-0-mW.csv", "optical_depth")):
        rows = list(csv.reader(StringIO(archive.read(name).decode("utf-8"))))
        assert rows[0][0] == "Relative frequency [MHz]"
        assert len(rows) == 402
        np.testing.assert_array_equal([float(row[0]) for row in rows[1:]], model["x"])
        np.testing.assert_array_equal([float(row[1]) for row in rows[1:]], model["series"][0][key])
