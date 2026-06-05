"""
Scheme registry. The front-end iterates this; adding a scheme = one import +
one entry here (no UI change). FWM today; AT/EIT/CPT/SAS/Hanle/… land next.
"""
from .base import ExtraView, ParamSpec, Preset, Scheme
from .absorption import ODScheme, LambdaScheme
from .sas import SASScheme
from .magneto import MagnetoScheme
from .fwm import FWMScheme

# Each dropdown entry is one engine cluster; same-engine schemes are merged behind
# a single entry that carries an in-panel `view` selector + per-regime presets:
#   SASScheme()     → OD / SAS         (pump power = 0 recovers OD)
#   LambdaScheme()  → EIT / AT / CPT   (the coupling-Ω_c regimes of one Λ system)
#   MagnetoScheme() → Hanle / EIA / NMOR  (transmission vs rotation readouts, one solve)
# (ODScheme stays in absorption.py as the single-2-level validation primitive.)
_SCHEMES = [
    SASScheme(),
    LambdaScheme(),
    MagnetoScheme(),
    FWMScheme(),
]

# Single-regime instances, resolvable by name for the physics tests and direct use,
# but kept out of the dropdown — the merged entries above present them.
_ALIASES = [
    LambdaScheme("eit"), LambdaScheme("at"), LambdaScheme("cpt"),
    MagnetoScheme("hanle"), MagnetoScheme("eia"), MagnetoScheme("nmor"),
]

REGISTRY = {s.name: s for s in (_SCHEMES + _ALIASES)}


def all_schemes():
    """Dropdown schemes in registration order (merged entries)."""
    return list(_SCHEMES)


def get(name):
    return REGISTRY[name]


__all__ = ["Scheme", "ParamSpec", "Preset", "ExtraView",
           "REGISTRY", "all_schemes", "get"]
