"""Group the existing diagnostic table without changing or omitting results."""
from collections import OrderedDict

import streamlit as st


def diagnostic_sections(markdown):
    groups = OrderedDict((name, []) for name in ("Claims", "Solver", "Geometry & coupling", "Passivity"))
    reasons = []
    for line in markdown.splitlines()[2:]:
        if not line.startswith("|") or "|" not in line[1:]:
            continue
        label, value = line[1:].rsplit("|", 1)[0].split("|", 1)
        label, value = label.strip(), value.strip()
        if label == "Claim-gate reasons":
            reasons = value.split("; ")
            continue
        if label in {"Model scope", "Validation level", "Quantitative gain claim", "Physical squeezing claim"}:
            group = "Claims"
        elif label.startswith(("Canonical", "Small-signal photon", "Bare-map", "Noise-trace", "Passive", "L1", "Legacy gain")):
            group = "Passivity"
        elif label.startswith(("2-D", "Probe / conjugate effective", "Coupling norm", "Combined residual", "Effective coupling", "Pump-depletion", "Small-signal peak", "Cross-coupling")) or "phase" in label.lower():
            group = "Geometry & coupling"
        else:
            group = "Solver"
        groups[group].append(line)
    return groups, reasons


def render_diagnostics(markdown):
    groups, reasons = diagnostic_sections(markdown)
    for title, rows in groups.items():
        if not rows and not (title == "Claims" and reasons):
            continue
        st.markdown(f"#### {title}")
        if rows:
            st.markdown("| Check | Value |\n|---|---|\n" + "\n".join(rows))
        if title == "Claims" and reasons:
            st.markdown("**Claim-gate reasons**")
            st.markdown("\n".join(f"- {reason}" for reason in reasons))
