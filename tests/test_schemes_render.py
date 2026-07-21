"""
Smoke test for the scheme/UI contract: every registered scheme must compute and
render (metrics + a matplotlib figure) without error. Catches rendering-layer
regressions independent of the physics tests.

    python tests/test_schemes_render.py   # or: pytest tests/test_schemes_render.py
"""
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import schemes  # noqa: E402
from gabes.schemes import fwm  # noqa: E402
from gabes.schemes.sas import GENERIC  # noqa: E402
from gabes.ui_metrics import partition_metrics  # noqa: E402


_UNAVAILABLE_VALUE = re.compile(
    r"(?:^|\b)(?:nan|n/a|[+\-]?inf(?:inity)?|unresolved|unconverged)(?:\b|$)",
    re.I,
)
_UNAVAILABLE_PLACEHOLDERS = {"—", "–", "-"}


def _fast(params):
    """Cheapen heavy knobs so the smoke test stays quick."""
    p = dict(params)
    if "doppler" in p:
        p["doppler"] = "off"
    if "resolution" in p:
        p["resolution"] = "Fast  (~4 s)"
    if "velocity_classes" in p:
        p["velocity_classes"] = 1
    if "scan_points" in p:
        p["scan_points"] = 51
    return p


def _assert_metric_contract(name, metrics, hero_count=2):
    assert metrics, f"{name}: no metrics"
    assert type(hero_count) is int and hero_count in (1, 2), (
        f"{name}: hero_count must be 1 or 2, got {hero_count!r}"
    )
    for metric in metrics:
        assert "label" in metric and "value" in metric, f"{name}: bad metric {metric}"

    explicit_heroes = [metric for metric in metrics if metric.get("tier") == "hero"]
    assert len(explicit_heroes) == hero_count, (
        f"{name}: expected {hero_count} explicitly prioritized hero metrics, "
        f"got {[metric.get('label') for metric in explicit_heroes]}"
    )
    heroes, _ribbon = partition_metrics(metrics, hero_count=hero_count)
    assert heroes == explicit_heroes, f"{name}: partition changed explicit hero order"

    for metric in heroes:
        value = str(metric["value"]).strip()
        assert value, f"{name}: empty hero value: {metric}"
        if value in _UNAVAILABLE_PLACEHOLDERS or _UNAVAILABLE_VALUE.search(value):
            assert metric.get("kind") == "status", (
                f"{name}: unavailable numeric result was promoted as a hero: {metric}"
            )


def test_metric_contract_rejects_placeholder_numeric_hero():
    with pytest.raises(AssertionError, match="unavailable numeric result"):
        _assert_metric_contract("placeholder", [
            {"label": "Width", "value": "—", "tier": "hero"},
            {"label": "Transmission", "value": "0.5", "tier": "hero"},
        ])


def _defaults_for(scheme, label):
    params = scheme.defaults()
    params.update(scheme.recommended_defaults(params)[label])
    return _fast(params)


def _internal_mode_cases():
    cases = []

    sas = schemes.get("sas")
    pump_off = _fast(sas.defaults())
    pump_off["pump_power_mw"] = 0.0
    cases.append(("SAS / pump-off OD", sas, pump_off))
    generic = _fast(sas.defaults())
    generic.update(species=GENERIC, transitions="single line", pump_power_mw=0.5)
    cases.append(("SAS / Generic pump-on", sas, generic))

    lam = schemes.get("lambda")
    for label in ("EIT", "AT", "CPT"):
        cases.append((f"Lambda / {label}", lam, _defaults_for(lam, label)))
    for label in ("EIT", "CPT"):
        opaque = _defaults_for(lam, label)
        opaque.update(temp_c=200.0, cell_mm=200.0, doppler="off")
        cases.append((f"Lambda / {label} opaque edge", lam, opaque))
    weak_at = _defaults_for(lam, "AT")
    weak_at.update(coupling_rabi_mhz=0.1, coupling_power_mw=0.01,
                   coupling_diameter_mm=5.0, doppler="off")
    cases.append(("Lambda / unresolved weak AT", lam, weak_at))

    rydberg = schemes.get("rydberg_eit")
    cases.append(("Rydberg / EIT", rydberg, _defaults_for(rydberg, "EIT")))
    at = _defaults_for(rydberg, "AT electrometry")
    cases.append(("Rydberg / resonant AT", rydberg, at))
    detuned_at = dict(at, mw_detuning_mhz=4.0)
    cases.append(("Rydberg / detuned AT", rydberg, detuned_at))
    lo_off = dict(at, lo_rabi_mhz=0.0)
    cases.append(("Rydberg / AT view with LO off", rydberg, lo_off))

    magneto = schemes.get("magneto")
    defaults = magneto.defaults()
    for index, params in enumerate(magneto.recommended_defaults(defaults).values()):
        cases.append((f"Magneto / regime {index + 1}", magneto,
                      _fast(dict(defaults, **params))))
    for offset in (2.0, 50.0):
        shifted = _fast(dict(defaults, b_offset_ut=offset))
        cases.append((f"Magneto / B offset {offset:g} uT", magneto, shifted))

    fwm_scheme = schemes.get("fwm")
    biphoton_seed = dict(fwm_scheme.defaults(), mode=fwm.MODE_BIPHOTON,
                         topology=fwm.TOPOLOGY_RB87_TELECOM,
                         biphoton_model=fwm.BIPHOTON_CALIBRATED)
    biphoton = fwm_scheme.recommended_defaults(biphoton_seed)[fwm.MODE_BIPHOTON]
    biphoton["biphoton_model"] = fwm.BIPHOTON_CALIBRATED
    cases.append(("FWM / calibrated Biphoton", fwm_scheme, biphoton))
    return cases


_INTERNAL_MODE_CASES = _internal_mode_cases()


def test_all_schemes_render():
    for scheme in schemes.all_schemes():
        params = _fast(scheme.defaults())
        raw = scheme.compute(params)
        view = scheme.observables(raw, params)
        assert view.get("figure") is not None, f"{scheme.name}: no figure"
        _assert_metric_contract(
            scheme.name, view.get("metrics"), view.get("hero_count", 2))
        figure_views = view.get("figure_views", [])
        if figure_views:
            assert figure_views[0]["figure"] is view["figure"], (
                f"{scheme.name}: first figure view must be the default figure"
            )
        figures_to_close = [view["figure"]]
        figures_to_close.extend(
            item["figure"] for item in figure_views if item.get("figure") is not None
        )
        closed = set()
        for figure in figures_to_close:
            if id(figure) not in closed:
                plt.close(figure)
                closed.add(id(figure))
        for _title, extra in view.get("figures", []):
            plt.close(extra)


@pytest.mark.parametrize(
    "case_name,scheme,params",
    _INTERNAL_MODE_CASES,
    ids=[case[0] for case in _INTERNAL_MODE_CASES],
)
def test_internal_modes_have_useful_heroes(case_name, scheme, params):
    raw = scheme.compute(params)
    view = scheme.headless_observables(raw, params)
    assert view.get("figure") is None, f"{case_name}: headless path built a figure"
    assert not view.get("figure_views"), f"{case_name}: headless path built carousel views"
    _assert_metric_contract(
        case_name, view.get("metrics"), view.get("hero_count", 2))


if __name__ == "__main__":
    test_all_schemes_render()
    print("All schemes compute + render OK:",
          ", ".join(s.name for s in schemes.all_schemes()))
