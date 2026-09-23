"""Unit-aware parsing for typed numeric entry (docs/ui_redesign/decisions.md D3).

The ScrubField sends the raw text a person typed and lets Python turn it into a
number, so there is exactly one parser to reason about and to test — the
component's JavaScript never needs a copy of these tables.

Accepted entry: a number (``1.5``, ``-2e3``, ``1,234.5``, ``−7`` with a real
minus sign) optionally followed by a unit. A bare number means "in the knob's
own unit"; a unit converts, as long as it names the same physical quantity.
"""
import re

# Powers of ten, kept as exponents: a single 10**(a-b) is exact where a
# multiply-then-divide by two floats is not (0.5 mT would land on 500.00000000000006 µT).
_PREFIXES = {
    "T": 12, "G": 9, "M": 6, "k": 3,
    "d": -1, "c": -2, "m": -3, "µ": -6, "n": -9, "p": -12, "f": -15,
}

# A leading prefix is only stripped when what remains names a base quantity.
# That is what keeps dB from becoming deci-B and Torr from becoming tera-orr.
_BASES = {
    "Hz", "m", "s", "W", "T", "A", "V", "rad", "Torr", "g", "eV", "J", "N",
    "Pa", "bar", "K", "F", "S", "H", "mol", "cd", "lm", "sr", "Wb", "C",
}

_ALIASES = {
    "u": "µ", "μ": "µ",            # typed ASCII u / Greek mu -> micron sign
    "C": "°C", "degC": "°C", "celsius": "°C", "°c": "°C",
    "kelvin": "K",
    "deg": "°", "degree": "°", "degrees": "°",
    "gauss": "G", "tesla": "T",
    "ohm": "Ω",
}

# Same quantity, different base unit, as a power of ten: 1 G = 1e-4 T.
_BRIDGE = {("G", "T"): -4, ("T", "G"): 4}

_NUMBER = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?")


class UnitError(ValueError):
    """The text parsed as a number but its unit does not fit the knob."""


def _normalize_unit(unit):
    """Fold the spellings a keyboard can produce onto one canonical form."""
    text = (unit or "").strip().replace(" ", "").replace(" ", "")
    text = text.replace("μ", "µ").replace("µ", "µ")
    if text in _ALIASES:
        return _ALIASES[text]
    # A lone "u"/"μ" prefix in front of a base ("uW", "μs") normalizes too.
    if text[:1] == "u" and len(text) > 1:
        text = "µ" + text[1:]
    return text


def _head(base):
    """The quantity a compound unit leads with: 'W/cm²' -> 'W'."""
    return re.split(r"[/·*]", base, maxsplit=1)[0]


def _decompose(unit):
    """Split a unit into (power-of-ten onto its base, base). Unknown units are atomic."""
    text = _normalize_unit(unit)
    if not text:
        return 0, ""
    for prefix, exponent in _PREFIXES.items():
        if text.startswith(prefix) and len(text) > len(prefix):
            rest = text[len(prefix):]
            if _head(rest) in _BASES:
                return exponent, rest
    return 0, text


def convert(value, from_unit, to_unit):
    """``value`` given in ``from_unit``, expressed in ``to_unit``.

    Raises UnitError when the two units do not measure the same thing.
    """
    from_exp, from_base = _decompose(from_unit)
    to_exp, to_base = _decompose(to_unit)
    if from_base == to_base:
        return value * 10.0 ** (from_exp - to_exp)
    # Temperature is the one offset scale the app uses; an offset scale and a
    # prefix do not mix, so "mK" is refused rather than quietly mis-shifted.
    if {from_base, to_base} == {"°C", "K"} and from_exp == to_exp == 0:
        return value - 273.15 if from_base == "K" else value + 273.15
    bridge = _BRIDGE.get((from_base, to_base))
    if bridge is not None:
        return value * 10.0 ** (from_exp + bridge - to_exp)
    raise UnitError(f"expected {to_unit or 'a plain number'}")


def parse_entry(text, unit):
    """Turn typed text into a number in the knob's own ``unit``.

    Returns ``(value, None)`` on success and ``(None, message)`` on failure,
    where the message is short enough to sit under the control.
    """
    raw = (text or "").strip().replace("−", "-").replace("–", "-")
    raw = re.sub(r"(?<=\d),(?=\d\d\d\b)", "", raw)
    match = _NUMBER.match(raw)
    if not match:
        return None, "not a number"
    try:
        value = float(match.group(0))
    except ValueError:                                  # pragma: no cover
        return None, "not a number"
    typed_unit = raw[match.end():].strip()
    if not typed_unit:
        return value, None
    try:
        return convert(value, typed_unit, unit), None
    except UnitError as exc:
        return None, str(exc)


def clamp(value, vmin, vmax):
    """Pull a value inside the knob's range; report it when that happened."""
    if vmin is not None and value < vmin:
        return vmin, f"clamped to {vmin:g}"
    if vmax is not None and value > vmax:
        return vmax, f"clamped to {vmax:g}"
    return value, None


def decimals_for(fmt):
    """How many decimals a printf-style slider format shows."""
    match = re.search(r"\.(\d+)f", fmt or "")
    return int(match.group(1)) if match else 0


def display_text(value, fmt):
    """The string shown in the value cell.

    Keeps the format the slider would have used, unless a typed value carries
    more precision than that format can hold — an exact entry must not be
    silently rounded away in the one place that shows it.
    """
    if value is None:
        return ""
    if isinstance(value, int) or (isinstance(value, float) and value.is_integer()
                                  and not fmt):
        return str(int(value))
    if fmt:
        shown = fmt % value
        try:
            if float(shown) == float(value):
                return shown
        except ValueError:                              # pragma: no cover
            return shown
    for precision in range(6, 13):
        shown = f"{value:.{precision}g}"
        if float(shown) == float(value):
            return shown
    return f"{value:.12g}"                              # pragma: no cover
