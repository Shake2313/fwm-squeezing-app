"""Readout strip, status chips and plot captions (docs/ui_redesign/decisions.md D4, D8).

Pure string builders: Streamlit only receives the finished HTML, so the rules
here are unit-tested without a running app.
"""
import re
import math
from html import escape

from gabes.ui_metrics import partition_metrics, split_metric_value
from gabes_ui.format import looks_numeric

MAX_SECONDARY = 4


def scheme_label(title):
    """Compact display name; scheme identity and saved choices stay unchanged."""
    return re.sub(r"\s+\([^()]*\)$", "", str(title))


def present_metrics(metrics, *, scheme_name, raw=None):
    """Display-only qualifications; never rewrite model/export descriptors."""
    items = [dict(m) for m in metrics]
    if scheme_name != "fwm":
        return items
    gain_supported = bool((raw or {}).get("claim_gate", {}).get("quantitative_gain_supported"))
    for metric in items:
        if metric.get("label") == "Squeezing Degree [dB]" and metric.get("delta"):
            note = str(metric.pop("delta"))
            if "physical squeezing unavailable" in note:
                scope = ("Conditional NF" if "conditional passive completion" in note else
                         "NF unavailable" if "passive completion unavailable" in note else
                         "Legacy diagnostic")
                metric["status"] = dict(label="Noise interpretation",
                                        value=f"{scope} · physical squeezing unavailable", help=note)
            else:
                metric["delta"] = note
        if not gain_supported and metric.get("label") in {
                "Seed gain G_s", "Conjugate gain G_c", "Atomic probe gain G_s", "Atomic conjugate gain G_c"}:
            try:
                number = float(metric["value"])
                if math.isfinite(number):
                    metric["value"] = f"{number:.3g}"
                    metric["approximate"] = True
                    metric["help"] = "Qualitative model estimate; quantitative gain claim unsupported. " + str(metric.get("help") or "")
            except (ValueError, TypeError):
                pass
    return items

# Status values that read as a caveat get the warning tone. P4 lets a scheme
# state the tone explicitly; until then the wording decides.
_CAVEAT = re.compile(r"limited|required|unavailable|unresolved|invalid|fail|warn|not ", re.I)
_LEADING_SYMBOLS = re.compile(r"^[^\w(]+", re.U)


def partition_readout(metrics, hero_count=2, max_secondary=MAX_SECONDARY):
    """Split scheme metrics into (heroes, secondary, statuses, overflow).

    Status metrics (kind="status") never occupy a value cell; they become
    chips. Heroes follow `partition_metrics` over the remaining values, the
    next `max_secondary` stay visible and the rest go behind "all metrics".
    """
    items = list(metrics)
    statuses = [m for m in items if str(m.get("kind")).lower() == "status"]
    values = [m for m in items if str(m.get("kind")).lower() != "status"]
    heroes, rest = partition_metrics(values, hero_count=hero_count)
    return heroes, rest[:max_secondary], statuses, rest[max_secondary:]


def _value_html(metric):
    value = str(metric.get("value", ""))
    number, unit = split_metric_value(value, kind=metric.get("kind"))
    num_class = " g-num" if looks_numeric(value) else ""
    unit_html = f"<span class='g-unit'>{escape(unit)}</span>" if unit else ""
    approx = "<span class='g-estimate' aria-label='approximate'>≈ </span>" if metric.get("approximate") else ""
    return f"<div class='g-cell-value{num_class}'>{approx}{escape(number)}{unit_html}</div>"


def _cell_html(metric, role):
    label = str(metric.get("label", ""))
    _, unit = split_metric_value(metric.get("value", ""), kind=metric.get("kind"))
    if unit:
        label = re.sub(r"\s*\[" + re.escape(unit) + r"\]$", "", label)
    label = caption_html(label)
    title = escape(str(metric.get("help") or ""), quote=True)
    delta = metric.get("delta")
    sub = f"<div class='g-cell-sub'>{escape(str(delta))}</div>" if delta is not None else ""
    if metric.get("status"):
        sub += f"<div class='g-cell-status'>{_chip_html(metric['status'])}</div>"
    return (f"<div class='g-cell g-cell--{role}' title='{title}'>"
            f"<div class='g-cell-label'>{label}</div>{_value_html(metric)}{sub}</div>")


def status_tone(value):
    return "warn" if _CAVEAT.search(str(value)) else "neutral"


def _chip_html(metric):
    label = str(metric.get("label", ""))
    value = str(metric.get("value", ""))
    tip = f"{label} — {metric.get('help')}" if metric.get("help") else label
    return (f"<span class='g-chip g-chip--{status_tone(value)}' title='{escape(tip, quote=True)}' "
            f"aria-label='{escape(f'{label}: {value}', quote=True)}'>{escape(value)}</span>")


def strip_html(heroes, secondary, statuses):
    """One-line readout: hero cells, secondary cells, then status chips."""
    cells = [_cell_html(m, "hero") for m in heroes]
    cells += [_cell_html(m, "secondary") for m in secondary]
    if statuses:
        chips = "".join(_chip_html(m) for m in statuses)
        cells.append(f"<div class='g-cell g-cell--status'>{chips}</div>")
    return f"<section class='g-strip' aria-label='Key results'>{''.join(cells)}</section>"


def metrics_table_markdown(metrics):
    """Every metric with its help text, for the "all metrics" popover."""
    rows = ["| Metric | Value | Note |", "|---|---|---|"]
    for m in metrics:
        note = " ".join(str(x) for x in (m.get("delta"), m.get("help")) if x)
        cells = (str(m.get("label", "")), str(m.get("value", "")), note)
        rows.append("| " + " | ".join(c.replace("|", "\\|").replace("\n", " ") for c in cells) + " |")
    return "\n".join(rows)


def clean_choice_label(label):
    """Drop leading emoji/symbols a scheme put in a choice label ("🕳️ EIT dip" → "EIT dip")."""
    text = _LEADING_SYMBOLS.sub("", str(label)).strip()
    return text or str(label)


# Matplotlib titles stay ASCII (mathtext layout lock); the HTML caption can
# show the physics symbols they spell out.
_CAPTION_RULES = [
    (re.compile(r"\bOmega_(\w+)"), r"Ω<sub>\1</sub>"),
    (re.compile(r"\bOmega\b"), "Ω"),
    (re.compile(r"\bDelta\b"), "Δ"),
    (re.compile(r"\bdelta\b"), "δ"),
    (re.compile(r"\bGamma\b"), "Γ"),
    (re.compile(r"\bgamma\b"), "γ"),
    (re.compile(r"\b(?:eta|η)_(\w+)"), r"η<sub>\1</sub>"),
    (re.compile(r"\bG_(\w+)"), r"G<sub>\1</sub>"),
    (re.compile(r"\beta\b"), "η"),
    (re.compile(r"(\d)\s?uW\b"), r"\1 µW"),
    (re.compile(r"(\d)\s?uT\b"), r"\1 µT"),
    (re.compile(r"(\d)\s?deg\b"), r"\1°"),
    (re.compile(r"(\d) C\b"), r"\1 °C"),
    (re.compile(r"\^2\b"), "²"),
    (re.compile(r"-&gt;"), "→"),
    (re.compile(r"\s{2,}"), " "),
]


def caption_html(title):
    """Escaped, symbol-restored caption text for a figure title ("" if none)."""
    text = escape(" ".join(str(title or "").split()))
    for pattern, repl in _CAPTION_RULES:
        text = pattern.sub(repl, text)
    return text
