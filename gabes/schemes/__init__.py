"""
Scheme registry. The front-end iterates this; adding a scheme = one import +
one entry here (no UI change). FWM today; AT/EIT/CPT/SAS/Hanle/… land next.
"""
from .base import ExtraView, ParamSpec, Preset, Scheme
from .absorption import ODScheme, LambdaScheme
from .sas import SASScheme
from .magneto import MagnetoScheme
from .fwm import FWMScheme

# SASScheme is the merged absorption scheme (OD at pump = 0, SAS with pump on),
# so the standalone ODScheme is no longer registered — its class stays in
# absorption.py as the single-2-level validation primitive the Λ tests use.
_SCHEMES = [
    SASScheme(),
    LambdaScheme("eit"),
    LambdaScheme("at"),
    LambdaScheme("cpt"),
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
