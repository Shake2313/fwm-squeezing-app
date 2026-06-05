"""Shared visual polish for GABES matplotlib figures."""

PALETTE = {
    "ink": "#0F172A",
    "muted": "#64748B",
    "border": "#DCE6EF",
    "grid": "#E6EDF5",
    "surface": "#FFFFFF",
    "cyan": "#0284C7",
    "blue": "#2563EB",
    "teal": "#0F766E",
    "warm": "#F97316",
    "rose": "#F43F5E",
    "violet": "#7C3AED",
}

LINE_CYCLE = [
    PALETTE["cyan"],
    PALETTE["rose"],
    PALETTE["teal"],
    PALETTE["warm"],
    PALETTE["violet"],
    PALETTE["blue"],
]

DEFAULT_LINE_COLORS = {
    "C0", "C1", "C2", "C3", "C4", "C5",
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
    "tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown",
}


def apply_gabes_plot_style(target):
    """Apply the quiet scientific-console style to a Figure or Axes."""
    fig = getattr(target, "figure", None)
    axes = [target] if fig is not None and not hasattr(target, "axes") else None

    if hasattr(target, "axes"):
        fig = target
        axes = list(target.axes)

    if fig is None or axes is None:
        return target

    fig.patch.set_facecolor(PALETTE["surface"])

    for ax in axes:
        ax.set_facecolor(PALETTE["surface"])
        ax.grid(True, color=PALETTE["grid"], linewidth=0.8, alpha=0.85)
        ax.set_axisbelow(True)

        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(PALETTE["border"])
            ax.spines[side].set_linewidth(0.9)

        ax.tick_params(colors=PALETTE["muted"], labelsize=9)
        ax.xaxis.label.set_color(PALETTE["muted"])
        ax.yaxis.label.set_color(PALETTE["muted"])
        ax.title.set_color(PALETTE["ink"])
        ax.title.set_fontweight("semibold")

        for idx, line in enumerate(ax.get_lines()):
            color = line.get_color()
            if color in DEFAULT_LINE_COLORS:
                line.set_color(LINE_CYCLE[idx % len(LINE_CYCLE)])
            line.set_linewidth(max(line.get_linewidth(), 1.8))

        legend = ax.get_legend()
        if legend is not None:
            frame = legend.get_frame()
            frame.set_facecolor(PALETTE["surface"])
            frame.set_edgecolor(PALETTE["border"])
            frame.set_linewidth(0.8)
            for text in legend.get_texts():
                text.set_color(PALETTE["ink"])

    fig.tight_layout()
    return target
