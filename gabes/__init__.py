"""
GABES — Generic Atomic Bloch Equation Solver.

Layered package:
  constants    physical constants + 85Rb D1 line data + basic field helpers
  core         physics-agnostic engine (super-operators, Liouvillian, Floquet solve)
  doppler      Maxwell velocity grid + Δ_eff batching + Doppler average
  atoms        AtomModel dataclass + level-scheme registry
  observables  χ → gain, squeezing, (later) absorption / OD
  schemes/     experiment plugins (FWM today; AT/EIT/CPT/SAS/Hanle/… next)

The Streamlit front-end (streamlit_app.py) and the fwm_obe.py CLI shim sit on top.
"""

__version__ = "0.1.0"
