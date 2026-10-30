"""Frozen comparison modules must not contaminate later source inventories."""
from dataclasses import replace
import importlib
import sys
from types import ModuleType

import pytest

from analysis.magneto_speed_20261030 import validation as speed
from analysis.magneto_optical_assembly_20261030 import validation as optical


ALIASES = (
    "gabes._magneto_speed_before", "gabes.schemes._magneto_speed_before",
    "gabes._magneto_optical_before", "gabes.schemes._magneto_optical_before",
    "gabes._optical_assembly_reference",
)
FACTORIES = (
    (speed, "_BASELINES", "load_baselines", ALIASES[:2]),
    (optical, "_REFERENCE", "load_reference", ALIASES[2:4]),
)


@pytest.mark.parametrize("validator", [speed, optical], ids=["speed", "optical"])
@pytest.mark.parametrize("previous_kind", ["absent", "module", "none"])
@pytest.mark.parametrize("fail", [False, True], ids=["success", "failure"])
def test_loader_restores_module_registry_on_success_and_failure(
        validator, previous_kind, fail, tmp_path, monkeypatch):
    name = "gabes._reference_loader_isolation_test"
    source = tmp_path/"reference.py"
    body = ("from __future__ import annotations\n"
            "from dataclasses import dataclass\nimport sys\n"
            "registered_during_execution = sys.modules[__name__]\n"
            "@dataclass\nclass Value:\n    amount: int\n")
    if fail:
        body += "raise RuntimeError('fixture load failed')\n"
    source.write_text(body, encoding="utf-8")
    if previous_kind == "absent":
        monkeypatch.delitem(sys.modules, name, raising=False)
        previous = None
    else:
        previous = ModuleType(name) if previous_kind == "module" else None
        monkeypatch.setitem(sys.modules, name, previous)
    if fail:
        with pytest.raises(RuntimeError, match="fixture load failed"):
            validator._load(name, source)
    else:
        module = validator._load(name, source)
        assert module.registered_during_execution is module
        assert replace(module.Value(3), amount=4).amount == 4
    if previous_kind == "absent":
        assert name not in sys.modules
    else:
        assert name in sys.modules and sys.modules[name] is previous


@pytest.mark.parametrize("validator,cache,factory,names", FACTORIES,
                         ids=["speed", "optical"])
def test_frozen_factories_and_cache_hits_leave_only_live_package_bindings(
        validator, cache, factory, names, monkeypatch):
    import gabes
    from gabes.schemes import magneto_unified

    live = gabes.magneto_multilevel
    assert not any(name in sys.modules for name in ALIASES)
    monkeypatch.setattr(validator, cache, None)
    modules = getattr(validator, factory)()
    backend, scheme = modules
    assert scheme.magneto_multilevel is backend
    assert gabes.magneto_multilevel is live
    assert magneto_unified.magneto_multilevel is live
    assert not any(name in sys.modules for name in ALIASES)
    assert getattr(validator, factory)() is modules
    assert not any(name in sys.modules for name in ALIASES)
    # Dataclass construction remains usable after its temporary registration.
    atom = backend.build_multilevel_atom("Rb87", "D1")
    assert atom.dimension == 34


@pytest.mark.parametrize("validator,cache,factory,names", FACTORIES,
                         ids=["speed", "optical"])
def test_partial_factory_failure_restores_live_backend_and_allows_retry(
        validator, cache, factory, names, monkeypatch):
    import gabes
    from gabes.schemes import magneto_unified

    live = gabes.magneto_multilevel
    original_load = validator._load
    monkeypatch.setattr(validator, cache, None)

    def fail_scheme(name, path):
        if name == names[1]:
            raise RuntimeError("scheme fixture load failed")
        return original_load(name, path)

    monkeypatch.setattr(validator, "_load", fail_scheme)
    with pytest.raises(RuntimeError, match="scheme fixture load failed"):
        getattr(validator, factory)()
    assert getattr(validator, cache) is None
    assert gabes.magneto_multilevel is live
    assert magneto_unified.magneto_multilevel is live
    assert not any(name in sys.modules for name in ALIASES)
    monkeypatch.setattr(validator, "_load", original_load)
    backend, scheme = getattr(validator, factory)()
    assert scheme.magneto_multilevel is backend
    assert gabes.magneto_multilevel is live
    assert not any(name in sys.modules for name in ALIASES)


def test_optical_unit_test_reference_stays_private_and_usable():
    assembly_tests = importlib.import_module("test_magneto_optical_assembly")
    assert not any(name in sys.modules for name in ALIASES)
    assembly_tests._reference.cache_clear()
    module = assembly_tests._reference()
    assert module.build_multilevel_atom("Rb87", "D1").dimension == 34
    assert assembly_tests._reference() is module
    assert not any(name in sys.modules for name in ALIASES)
