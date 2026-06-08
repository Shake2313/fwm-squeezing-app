"""
Reference checks for the Rydberg-EIT electrometry scheme.

The public UI shows the static spectrum only. Experimental sensitivity values
from arXiv:2606.04354 are kept as internal constants and tested here.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import schemes  # noqa: E402
from gabes.schemes.rydberg import RydbergEITScheme  # noqa: E402


def _metric_value(view, label):
    for metric in view["metrics"]:
        if metric["label"] == label:
            return float(metric["value"].split()[0])
    raise AssertionError(f"missing metric {label!r}: {view['metrics']}")


def test_reference_defaults_match_rydberg_eit_paper():
    sc = schemes.get("rydberg_eit")
    sets = sc.recommended_defaults(sc.defaults())
    ref = sets["AT electrometry"]
    assert ref["probe_power_uw"] == 6.0
    assert ref["coupling_power_mw"] == 30.0
    assert ref["beam_diameter_mm"] == 0.15
    assert ref["cell_mm"] == 50.0
    assert ref["coupling_rabi_mhz"] == 3.0
    assert ref["lo_rabi_mhz"] == 3.7
    assert ref["mw_frequency_ghz"] == 37.0
    assert ref["if_khz"] == 40.0


def test_reference_eit_linewidth_near_experiment():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["EIT"]
    view = sc.observables(sc.compute(params), params)
    linewidth = _metric_value(view, "EIT linewidth")
    plt.close(view["figure"])
    assert 1.3 <= linewidth <= 1.9


def test_microwave_at_splitting_tracks_lo_rabi():
    sc = schemes.get("rydberg_eit")
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    view = sc.observables(sc.compute(params), params)
    split = _metric_value(view, "RF AT splitting")
    plt.close(view["figure"])
    assert abs(split / params["lo_rabi_mhz"] - 1.0) < 0.03


def test_sensitivity_reference_constants_are_internal_only():
    sc = schemes.get("rydberg_eit")
    assert isinstance(sc, RydbergEITScheme)
    assert sc.REFERENCE_SENSITIVITY_NV_CM_SQRT_HZ == 12.5
    assert sc.REFERENCE_PSN_LIMIT_NV_CM_SQRT_HZ == 11.2
    params = sc.recommended_defaults(sc.defaults())["AT electrometry"]
    view = sc.observables(sc.compute(params), params)
    labels = " ".join(metric["label"].lower() for metric in view["metrics"])
    tables = " ".join(table["markdown"].lower() for table in view["tables"])
    plt.close(view["figure"])
    assert "sensitivity" not in labels
    assert "psn" not in labels
    assert "12.5" not in tables and "11.2" not in tables


if __name__ == "__main__":
    test_reference_defaults_match_rydberg_eit_paper()
    test_reference_eit_linewidth_near_experiment()
    test_microwave_at_splitting_tracks_lo_rabi()
    test_sensitivity_reference_constants_are_internal_only()
    print("Rydberg-EIT reference checks OK.")
