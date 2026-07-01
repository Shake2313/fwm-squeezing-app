"""
Smoke tests for figureless metric/table readout.

The normal render smoke test still proves each scheme can build a Matplotlib
figure. This file proves the complementary path: headless_observables() must not
call figure construction, while preserving metrics/tables for batch/report use.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import schemes  # noqa: E402
from gabes.schemes.absorption import ODScheme  # noqa: E402
from gabes.schemes.base import Scheme  # noqa: E402


def _fast(params):
    p = dict(params)
    if "doppler" in p:
        p["doppler"] = "off"
    if "resolution" in p:
        p["resolution"] = "Fast  (~3 s)"
    return p


def _scheme_cases():
    cases = list(schemes.all_schemes())
    cases.append(ODScheme())
    return cases


def test_current_schemes_support_headless_observables():
    unsupported = [
        scheme.name for scheme in _scheme_cases()
        if not getattr(scheme, "supports_headless_observables", False)
    ]
    assert unsupported == []


def test_headless_observables_do_not_build_figures():
    original_subplots = plt.subplots

    def _fail_subplots(*_args, **_kwargs):
        raise AssertionError("headless observables must not build figures")

    try:
        plt.close("all")
        plt.subplots = _fail_subplots
        for scheme in _scheme_cases():
            params = _fast(scheme.defaults())
            raw = scheme.compute(params)
            view = scheme.headless_observables(raw, params)
            assert view.get("figure") is None, f"{scheme.name}: built a figure"
            assert not view.get("figures", []), f"{scheme.name}: built extra figures"
            assert view.get("metrics"), f"{scheme.name}: no metrics"
            for metric in view["metrics"]:
                assert "label" in metric and "value" in metric, (
                    f"{scheme.name}: bad metric {metric}")
    finally:
        plt.subplots = original_subplots
        plt.close("all")


def test_base_headless_observables_passes_include_figures_false():
    class ProbeScheme(Scheme):
        supports_headless_observables = True

        def param_schema(self):
            return []

        def compute(self, params):
            return {}

        def observables(self, raw, params, include_figures=True):
            return dict(
                metrics=[dict(label="include_figures", value=str(include_figures))],
                figure="figure" if include_figures else None,
                tables=[],
            )

    view = ProbeScheme().headless_observables({}, {})
    assert view["figure"] is None
    assert view["metrics"][0]["value"] == "False"


if __name__ == "__main__":
    test_current_schemes_support_headless_observables()
    test_headless_observables_do_not_build_figures()
    test_base_headless_observables_passes_include_figures_false()
    print("All current schemes expose headless observables.")
