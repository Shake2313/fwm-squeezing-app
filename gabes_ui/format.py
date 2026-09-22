"""Number-display rules for controls and readouts (docs/ui_redesign/decisions.md D4)."""
import re

_MAX_DECIMALS = 6
_DEFAULT_EXTRA = 2
_NUMERIC_TEXT = re.compile(r"^\s*[+\-−]?(?:\d|\.\d)")


def decimals(value):
    """Digits after the decimal point needed to show `value` exactly (capped)."""
    text = f"{abs(float(value)):.{_MAX_DECIMALS}f}".rstrip("0")
    return len(text.split(".", 1)[1]) if "." in text else 0


def slider_format(spec):
    """printf format for a numeric slider.

    An explicit `ParamSpec.format` wins. Integer sliders keep Streamlit's own
    format (None). Float sliders show as many decimals as the step needs,
    widened by up to `_DEFAULT_EXTRA` digits so a hand-set off-grid default
    (86.94 with step 0.1) is not displayed rounded, while a computed default
    (1326.2572434514343 with step 1) does not flood the label with digits.
    """
    explicit = getattr(spec, "format", None)
    if explicit:
        return explicit
    numbers = [v for v in (spec.default, spec.vmin, spec.vmax, spec.step) if v is not None]
    if not any(isinstance(v, float) for v in numbers):
        return None
    step_places = decimals(spec.step) if spec.step is not None else 0
    places = max(step_places, min(decimals(spec.default), step_places + _DEFAULT_EXTRA))
    return f"%.{places}f"


def looks_numeric(text):
    """True when a presentation string starts with a number (for mono display)."""
    return bool(_NUMERIC_TEXT.match(str(text)))
