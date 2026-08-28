"""Reproducible cell-heating analysis workflow for Rydberg-EIT electrometry."""

from typing import Any


def run_analysis(*args: Any, **kwargs: Any):
    """Lazy public entry point; avoids preloading the CLI module under ``-m``."""
    from .workflow import run_analysis as _run_analysis
    return _run_analysis(*args, **kwargs)


__all__ = ["run_analysis"]
