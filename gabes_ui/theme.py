"""Design tokens — the single source for UI colours and font stacks.

`.streamlit/config.toml` repeats the few values Streamlit's own widgets need
(primary, background, text, border); everything GABES draws itself reads
these tokens, either as CSS custom properties (`--g-<name>`) or directly when
it renders into an isolated iframe. Rationale: docs/ui_redesign/decisions.md D2.
"""
from pathlib import Path

ASSETS = Path(__file__).resolve().parents[1] / "assets"

LIGHT = {
    "ground": "#F5F7FA",
    "surface": "#FFFFFF",
    "rail": "#FBFCFD",
    "ink": "#0E1726",
    "ink-2": "#1F2D3D",
    "ink-3": "#3B4A5E",
    "muted": "#5B6B80",        # 5.3:1 on white (WCAG AA for body text)
    "line": "#E3E8EF",
    "line-soft": "#EEF2F6",
    "line-strong": "#CBD5E1",
    "fill-soft": "#EDF1F5",
    "accent": "#0284C7",       # fills, tracks, focus rings
    "accent-ink": "#0369A1",   # accent-coloured text (5.9:1 on white)
    "accent-soft": "#E6F2FA",
    "warn-ink": "#92400E",
    "warn-bg": "#FEF3C7",
    "warn-line": "#F6D48F",
    "ok-ink": "#047857",
    "danger-ink": "#B91C1C",
}

# Dark values exist only where GABES renders into iframes today; the full dark
# theme is phase P5 (theme.dark.* in config.toml plus these tokens).
DARK = {
    **LIGHT,
    "ground": "#0E1117",
    "surface": "#0E1117",
    "rail": "#11151C",
    "ink": "#F1F5F9",
    "ink-2": "#E2E8F0",
    "ink-3": "#CBD5E1",
    "muted": "#94A3B8",
    "line": "#334155",
    "line-soft": "#1E293B",
    "line-strong": "#475569",
    "fill-soft": "#1E293B",
    "accent-ink": "#7DD3FC",
    "accent-soft": "#0C2D48",
}

FONT_SANS = '"IBM Plex Sans", -apple-system, "Segoe UI", system-ui, sans-serif'
FONT_MONO = '"IBM Plex Mono", ui-monospace, Consolas, monospace'


def tokens(theme_base="light"):
    """Token dict for Streamlit's `theme.base` ("light" or "dark")."""
    return DARK if str(theme_base).lower() == "dark" else LIGHT


def css_variables(theme_base="light"):
    """`:root` block exposing every token as `--g-<name>` plus font stacks."""
    values = "".join(f"--g-{name}:{value};" for name, value in tokens(theme_base).items())
    return f":root{{{values}--g-font-sans:{FONT_SANS};--g-font-mono:{FONT_MONO};}}"


def app_css(theme_base="light"):
    """Token variables followed by the GABES stylesheet (assets/ui/gabes.css)."""
    sheet = (ASSETS / "ui" / "gabes.css").read_text(encoding="utf-8")
    return css_variables(theme_base) + "\n" + sheet
