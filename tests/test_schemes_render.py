"""
Smoke test for the scheme/UI contract: every registered scheme must compute and
render (metrics + a matplotlib figure) without error. Catches rendering-layer
regressions independent of the physics tests.

    python tests/test_schemes_render.py   # or: pytest tests/test_schemes_render.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import schemes  # noqa: E402


def _fast(params):
    """Cheapen heavy knobs so the smoke test stays quick."""
    p = dict(params)
    if "doppler" in p:
        p["doppler"] = "off"
    if "resolution" in p:
        p["resolution"] = "Fast  (~4 s)"
    return p


def test_all_schemes_render():
    for scheme in schemes.all_schemes():
        params = _fast(scheme.defaults())
        raw = scheme.compute(params)
        view = scheme.observables(raw, params)
        assert view.get("figure") is not None, f"{scheme.name}: no figure"
        assert view.get("metrics"), f"{scheme.name}: no metrics"
        for m in view["metrics"]:
            assert "label" in m and "value" in m, f"{scheme.name}: bad metric {m}"
        plt.close(view["figure"])
        for _title, extra in view.get("figures", []):
            plt.close(extra)


if __name__ == "__main__":
    test_all_schemes_render()
    print("All schemes compute + render OK:",
          ", ".join(s.name for s in schemes.all_schemes()))
