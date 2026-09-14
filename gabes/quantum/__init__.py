"""Opt-in quantum research; no experimental FWM squeezing prediction yet.

Finite temporal-mode channels use interleaved (x, p) quadratures, [x,p]=i,
and vacuum covariance I/2. Complex RF/Nambu blocks require an explicit
conversion and are deliberately not accepted as real quadrature channels.
"""

from .contracts import CONVENTION_ID

__all__ = ["CONVENTION_ID"]
