"""Display rules for slider values and readout strings (gabes_ui.format)."""
import pytest

from gabes import schemes
from gabes.schemes.base import ParamSpec
from gabes_ui.format import decimals, looks_numeric, slider_format


def _spec(default, vmin, vmax, step, fmt=None):
    return ParamSpec("x", "x", "g", default, vmin, vmax, step, format=fmt)


@pytest.mark.parametrize("value, places", [
    (1.0, 0), (40, 0), (0.5, 1), (0.01, 2), (86.94, 2), (761.702, 3), (0.0001, 4),
])
def test_decimals_counts_only_needed_digits(value, places):
    assert decimals(value) == places


@pytest.mark.parametrize("spec, expected", [
    (_spec(40.0, 20.0, 200.0, 1.0), "%.0f"),      # "40", not "40.00"
    (_spec(0.5, 0.0, 2.0, 0.01), "%.2f"),
    (_spec(75.0, 0.5, 200.0, 0.5), "%.1f"),
    (_spec(86.94, 0.0, 100.0, 0.1), "%.2f"),      # off-grid default keeps its digits
    (_spec(761.702, 300.0, 2500.0, 0.001), "%.3f"),
    (_spec(1326.2572434514343, 1.0, 5000.0, 1.0), "%.2f"),  # computed default: step + 2
    (_spec(0.0, 0.0, 1.0, 0.0001, "%.4f"), "%.4f"),  # explicit format wins
    (_spec(1401, 401, 4001, 100), None),          # integer slider: Streamlit default
])
def test_slider_format(spec, expected):
    assert slider_format(spec) == expected


def test_every_registered_numeric_slider_formats_its_default():
    for scheme in schemes.REGISTRY.values():
        for sp in scheme.param_schema():
            if sp.choices is not None or sp.control != "auto" or sp.vmin is None:
                continue
            fmt = slider_format(sp)
            if fmt is not None:
                shown = fmt % sp.default
                places = int(fmt[2:-1])
                # Shown value is the default rounded to the displayed digits.
                assert abs(float(shown) - sp.default) <= 0.5 * 10 ** -places + 1e-12, (
                    scheme.name, sp.name, fmt)


@pytest.mark.parametrize("text, numeric", [
    ("18.87 MHz", True), ("-7.52 dB", True), ("−1304.16 to −1285.29 MHz", True),
    (".5", True), ("resolution-limited", False), ("crossover", False), ("", False),
])
def test_looks_numeric(text, numeric):
    assert looks_numeric(text) is numeric
