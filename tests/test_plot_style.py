"""Figure styling must not depend on whatever set the global matplotlib font.

Importing a lab script is enough to change it process-wide: the AutoOD GUI in
`references/` runs `plt.rcParams['font.family'] = 'Malgun Gothic'` at import,
and that face has no superscript digits — so a ⁸⁵Rb title silently lost its
isotope number in any process that had imported it (a full `pytest` run does).
`apply_gabes_plot_style` therefore pins its own family stack per figure.

    python tests/test_plot_style.py     # or: pytest tests/test_plot_style.py
"""
import sys
import warnings
from io import BytesIO
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.text import Text

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import schemes  # noqa: E402
from gabes.plot_style import FONT_STACK, apply_gabes_plot_style  # noqa: E402

# Every non-ASCII character GABES puts into a figure: isotope superscripts,
# Greek symbols, units, the minus sign and the en/em dashes.
GLYPH_SAMPLE = "⁸⁵Rb ⁸⁷Rb ¹³³Cs — Γ_eff, Ω_c, µW, 2π, 40 °C, −1.5 GHz"


def _glyph_warnings(fig):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fig.savefig(BytesIO(), format="png")
    return [str(w.message) for w in caught if "missing from font" in str(w.message)]


def _hostile_font(family="Malgun Gothic"):
    """Context manager mimicking a lab script that took over the global font."""
    return plt.rc_context({"font.family": family})


def test_style_pins_a_font_stack_that_carries_every_gabes_glyph():
    with _hostile_font():
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.plot([0, 1], [0, 1], label=GLYPH_SAMPLE)
        ax.set_title(GLYPH_SAMPLE)
        ax.set_xlabel(GLYPH_SAMPLE)
        ax.set_ylabel(GLYPH_SAMPLE)
        ax.legend()
        assert _glyph_warnings(fig), "fixture no longer reproduces the failure"

        apply_gabes_plot_style(fig)
        assert _glyph_warnings(fig) == []
        plt.close(fig)


def test_style_reaches_titles_labels_legends_and_tick_labels():
    with _hostile_font():
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.plot([0, 1], [0, 1], label="curve")
        ax.set_title("title")
        ax.legend()
        apply_gabes_plot_style(fig)

        families = {tuple(text.get_fontfamily()) for text in fig.findobj(Text)}
        assert families == {tuple(FONT_STACK)}
        fig.canvas.draw()                      # tick labels are built here
        tick_families = {
            tuple(label.get_fontfamily()) for label in ax.get_xticklabels()}
        assert tick_families == {tuple(FONT_STACK)}
        plt.close(fig)


def test_style_does_not_touch_global_rcparams():
    before = list(matplotlib.rcParams["font.family"])
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.plot([0, 1], [0, 1])
    apply_gabes_plot_style(fig)
    plt.close(fig)
    assert matplotlib.rcParams["font.family"] == before


def test_style_accepts_an_axes_directly():
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.plot([0, 1], [0, 1])

    assert apply_gabes_plot_style(ax) is ax
    assert not ax.spines["top"].get_visible()
    plt.close(fig)


def test_a_styled_scheme_figure_renders_isotope_labels_under_a_hostile_font():
    sas = schemes.get("sas")
    params = sas.defaults()
    params.update(species="⁸⁵Rb", line="D1", pump_power_mw=0.0, temp_c=45.0,
                  scan_points=401)
    raw = sas.compute(params)
    with _hostile_font(), warnings.catch_warnings():
        # The scheme lays the figures out before any styling runs, so building
        # them under the hostile font warns by construction; that is the point
        # of the fixture, and those warnings stay inside this block.
        warnings.simplefilter("ignore")
        view = sas.observables(raw, params)
        figures = [item["figure"] for item in view["figure_views"]]
        assert any("⁸⁵Rb" in fig.axes[0].get_title() for fig in figures)
        for fig in figures:
            assert _glyph_warnings(fig), "styling has not been applied yet"
            apply_gabes_plot_style(fig)
            assert _glyph_warnings(fig) == []
            plt.close(fig)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print("ok:", fn.__name__)
    print(f"\nPlot style OK ({len(fns)} tests): font stack pinned per figure, "
          "global rcParams untouched.")
