"""
Figure: PPT-internal consistency diagnostic for the 2026-07-15 cell-heating deck.

IMPORTANT: this figure uses ONLY numbers transcribed from the PPT.
No GABES curve is re-plotted here (the model sweep lives in brief_summary.pdf),
so nothing in this figure is model-derived except two annotated optimum markers.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- PPT transcription -------------------------------------------------
# setpoint, actual temperature annotated in the PPT PSN note, sensitivity, err, PSN limit
rows = [
    ("room",  20.0, 12.5, 3.6, 11.2, False),   # actual T not annotated in deck
    ("40 C",  38.0,  6.68, 2.55, 6.72, True),
    ("45 C",  42.5, 17.6, 2.2, 17.6, True),
    ("50 C",  44.0, 23.8, 6.5, 25.0, True),
]
setp = [r[0] for r in rows]
Tact = np.array([r[1] for r in rows])
sens = np.array([r[2] for r in rows])
serr = np.array([r[3] for r in rows])
psn = np.array([r[4] for r in rows])
ratio = sens / psn

GABES_SLOPE_OPT = 42.0
GABES_SENS_OPT = 47.0

fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.5))

# ---- panel (a) ---------------------------------------------------------
ax = axes[0]
ax.axvspan(41.5, 43.0, color="#0F766E", alpha=0.10, zorder=0)
ax.axvline(GABES_SLOPE_OPT, color="#0F766E", ls=":", lw=1.3, zorder=1)
ax.axvline(GABES_SENS_OPT, color="#94a3b8", ls="--", lw=1.2, zorder=1)

o = np.argsort(Tact)
ax.errorbar(Tact[o], sens[o], yerr=serr[o], fmt="o-", color="#C2410C",
            ms=6, lw=1.6, capsize=3, zorder=3, label="measured sensitivity [PPT]")
ax.plot(Tact[o], psn[o], "s--", color="#17324D", ms=5, lw=1.4, mfc="white",
        zorder=3, label="PSN limit quoted in PPT")
ax.plot([Tact[0]], [sens[0]], "o", ms=13, mfc="none", mec="#C2410C",
        mew=1.4, zorder=4)
ax.annotate("room point:\nlikely earlier campaign",
            xy=(20.0, 12.5), xytext=(21.5, 25.5), fontsize=7.5, color="#C2410C",
            arrowprops=dict(arrowstyle="->", color="#C2410C", lw=0.9))
ax.text(GABES_SLOPE_OPT - 0.45, 33.2, "GABES slope opt. 42.0 C",
        fontsize=7.5, color="#0F766E", rotation=90, va="top", ha="center")
ax.text(GABES_SENS_OPT - 0.45, 33.2, "GABES cond. sens. opt. 47 C",
        fontsize=7.5, color="#64748b", rotation=90, va="top", ha="center")
ax.text(42.5, 20.6, "clearest EIT\n(actual 42.5 C)", fontsize=7.5, color="#0F766E",
        ha="center", va="bottom")

for lab, x, y in zip(setp, Tact, sens):
    if lab != "room":
        ax.annotate(f"setpoint {lab}", xy=(x, y), xytext=(0, -14),
                    textcoords="offset points", ha="center", fontsize=7,
                    color="#475569")
ax.set_xlabel("Actual / effective temperature from PPT PSN note [C]")
ax.set_ylabel("E-field sens. [nV cm$^{-1}$ Hz$^{-1/2}$]")
ax.set_title("(a) sensitivity and its own PSN floor track each other", fontsize=9)
ax.set_xlim(18, 49)
ax.set_ylim(0, 34)
ax.legend(fontsize=7.5, loc="upper left", framealpha=0.9)
ax.grid(alpha=0.25, lw=0.6)

# ---- panel (b) ---------------------------------------------------------
ax = axes[1]
xpos = np.arange(len(rows))
cols = ["#94a3b8", "#C2410C", "#17324D", "#17324D"]
ax.bar(xpos, ratio, color=cols, width=0.55, zorder=3)
ax.axhline(1.0, color="#0F766E", lw=1.4, ls="--", zorder=4)
ax.text(3.42, 1.02, "measured = own PSN floor", fontsize=7.5, color="#0F766E",
        ha="right", va="bottom")
for x, r in zip(xpos, ratio):
    ax.text(x, r + 0.03, f"{r:.2f}", ha="center", fontsize=8)
ax.set_xticks(xpos)
ax.set_xticklabels([f"{s}\n({t:.0f} C act.)" if s != "room" else "room\n(n/a)"
                    for s, t in zip(setp, Tact)], fontsize=8)
ax.set_ylabel("measured / PSN limit  [1]")
ax.set_ylim(0, 1.45)
ax.set_title("(b) no excess technical noise at any setpoint", fontsize=9)
ax.grid(axis="y", alpha=0.25, lw=0.6)

fig.suptitle("PPT-internal consistency check (2026-07-15 deck numbers only)",
             fontsize=10, y=1.005)
fig.tight_layout()
fig.savefig("ppt_consistency.pdf", bbox_inches="tight")
print("ratios:", np.round(ratio, 3))
print("PSN floor 38->44 C degradation:", round(psn[3] / psn[1], 2))
