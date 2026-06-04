"""
Scheme plugin contract.

A Scheme owns its physics *and* its UI surface, so the front-end stays generic:
it reads `param_schema()` and renders only those controls (no fixed slider wall).

Contract
  param_schema()           -> [ParamSpec]   knobs; `recompute=False` ones never
                                             trigger a resolve (navigate-only).
  presets()                -> [Preset]       one-click coherent parameter sets.
  info()                   -> str | None     reference / about markdown.
  extra_views()            -> [ExtraView]    optional heavy on-demand panels.
  compute(params)          -> dict           heavy solve; uses ONLY recompute knobs
                                             so the UI can cache it.
  observables(raw, params) -> dict           cheap; uses all knobs (incl. navigate)
                                             -> {"metrics":[...], "figure":fig,
                                                 "tables":[...]}.
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Callable, Optional


@dataclass(frozen=True)
class ParamSpec:
    """One UI control. Numeric slider unless `choices` is given (-> selectbox)."""
    name: str
    label: str
    group: str
    default: Any
    vmin: Optional[float] = None
    vmax: Optional[float] = None
    step: Optional[float] = None
    unit: str = ""
    recompute: bool = True
    help: str = ""
    choices: Optional[tuple] = None
    advanced: bool = False        # render inside an "Advanced" expander
    endpoints: Optional[tuple] = None   # (left, right) caption under a slider


@dataclass(frozen=True)
class Preset:
    """A named set of parameter values applied to the sidebar in one click."""
    name: str
    values: dict
    help: str = ""
    icon: str = "⚡"


@dataclass
class ExtraView:
    """An optional heavy panel rendered behind a button (e.g. FWM full scan)."""
    key: str
    description: str
    compute: Callable[[dict], dict]      # heavy -> picklable dict (UI caches it)
    render: Callable[[dict], Any]        # dict -> matplotlib Figure


class Scheme(ABC):
    name: str = ""
    cluster: str = ""
    title: str = ""
    caption: str = ""
    cache_version: str = "1"
    defaults_version: str = "1"

    @abstractmethod
    def param_schema(self) -> list:
        ...

    def presets(self) -> list:
        return []

    def recommended_defaults(self, params: dict):
        """Labelled default slider presets for the current selection, as an
        ordered {button_label: {param_name: value}} dict — the UI renders one
        button per entry. None if the scheme offers no such defaults."""
        return None

    def info(self):
        return None

    def extra_views(self) -> list:
        return []

    @abstractmethod
    def compute(self, params: dict) -> dict:
        ...

    @abstractmethod
    def observables(self, raw: dict, params: dict) -> dict:
        ...

    # ---- helpers shared by the UI ----
    def recompute_keys(self):
        return tuple(s.name for s in self.param_schema() if s.recompute)

    def defaults(self):
        return {s.name: s.default for s in self.param_schema()}
