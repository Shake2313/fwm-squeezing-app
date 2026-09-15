"""Select pre-hotfix evidence for the runtime absorption normalization.

The original shared-worktree capture used C_F² line weights. The isolated
parent-commit capture used p_F·C_F². This only selects comparison evidence;
it never changes the runtime response, normalizes outputs, or refits the gain.
"""
from pathlib import Path

from gabes.schemes import fwm


HERE = Path(__file__).resolve().parent


def snapshot_directory():
    """Return the matching frozen capture; reject unrecorded conventions."""
    factors = fwm._reference_population_factors()
    if factors == {2: 1.0, 3: 1.0}:
        return HERE / "before"
    if factors == {2: 5.0 / 12.0, 3: 7.0 / 12.0}:
        return HERE / "before_parent"
    raise RuntimeError(
        "No pre-hotfix snapshot for absorption population factors "
        f"{factors!r}; expected C_F² or p_F·C_F² for both ground manifolds")
