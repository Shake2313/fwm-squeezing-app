"""
Scheme registry. The front-end iterates this; adding a scheme = one import +
one entry here (no UI change). FWM today; AT/EIT/CPT/SAS/Hanle/… land next.
"""
from .base import ExtraView, ParamSpec, Preset, Scheme
from .absorption import ODScheme, LambdaScheme
from .sas import SASScheme
from .magneto import MagnetoScheme
from .fwm import FWMScheme

_SCHEMES = [
    ODScheme(),
    LambdaScheme("eit"),
    LambdaScheme("at"),
    LambdaScheme("cpt"),
    SASScheme(),
    MagnetoScheme("hanle"),
    MagnetoScheme("eia"),
    MagnetoScheme("nmor"),
    FWMScheme(),
]

REGISTRY = {s.name: s for s in _SCHEMES}


def all_schemes():
    """Schemes in registration order."""
    return list(_SCHEMES)


def get(name):
    return REGISTRY[name]


__all__ = ["Scheme", "ParamSpec", "Preset", "ExtraView",
           "REGISTRY", "all_schemes", "get"]
