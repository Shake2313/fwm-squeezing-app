"""
Build the v2 report figures from the ultra-fidelity frontier scan output
(analysis/squeezing_frontier.npz, produced by
.claude/skills/fwm-squeezing-frontier/scripts/scan_squeezing_frontier.py
--fidelity ultra). Run from the repo root:

    python docs/squeezing_report/make_v2_figs.py

Outputs into docs/squeezing_report/:
  * squeezing_frontier.png      -- copy of the scan's main 3-panel figure
  * squeezing_frontier_od.png   -- NEW: in-cell optical-depth map + xi(T) slices
                                   showing the gain-limited / loss-limited
                                   interior optimum that the ultra model adds.
"""
import shutil
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
SRC = REPO / "analysis" / "squeezing_frontier.npz"

d = np.load(SRC)
D, T = d["delta_ghz"], d["temp_c"]
Xi, OD = d["xi_dB"], d["in_cell_od"]
det_floor = float(d["det_floor_dB"])
xi_opt = float(d["xi_opt_dB"])
oiD, oiT = int(d["opt_iD"]), int(d["opt_iT"])

# main 3-panel figure: just copy it next to the report
shutil.copyfile(REPO / "analysis" / "squeezing_frontier.png",
                OUT / "squeezing_frontier.png")

fig, (axOD, axS) = plt.subplots(1, 2, figsize=(13, 4.8))
extent = [T[0], T[-1], D[0], D[-1]]

# (1) in-cell optical-depth map: transparency windows vs absorption notches
im = axOD.imshow(np.clip(OD, 0, 2.0), origin="lower", aspect="auto",
                 extent=extent, cmap="inferno")
fig.colorbar(im, ax=axOD, label="in-cell OD (clip 2.0)")
axOD.scatter([T[oiT]], [D[oiD]], color="cyan", marker="o", s=70, edgecolor="k",
             zorder=6, label=f"optimum ({xi_opt:.2f} dB)")
axOD.set_xlabel("T [C]"); axOD.set_ylabel("Delta [GHz]")
axOD.set_title("in-cell optical depth OD(Delta, T)")
axOD.legend(loc="upper left", fontsize=8)

# (2) xi(T) slices: low T gain-limited, high T loss-limited -> interior optimum
slices = [(-2.2, "transparency (optimum band)"),
          (-1.9, "frontier"),
          (-1.0, "loss-limited"),
          (+0.9, "Sim detuning"),
          (0.0, "one-photon resonance")]
for dv, lab in slices:
    i = int(np.argmin(np.abs(D - dv)))
    axS.plot(T, Xi[i, :], "-", lw=1.6, label=f"Delta={D[i]:+.1f} GHz ({lab})")
    j = int(np.nanargmin(Xi[i, :]))
    axS.scatter([T[j]], [Xi[i, j]], s=30, zorder=5, color=axS.lines[-1].get_color())
axS.axhline(det_floor, color="crimson", ls="--", lw=1.2,
            label=f"detection floor {det_floor:.2f} dB")
axS.set_xlabel("T [C]"); axS.set_ylabel("xi [dB]")
axS.set_title("xi vs T: gain-limited (low T) vs loss-limited (high T)")
axS.grid(alpha=0.3); axS.legend(fontsize=7, loc="lower left")

fig.tight_layout()
fig.savefig(OUT / "squeezing_frontier_od.png", dpi=130)
print("wrote", OUT / "squeezing_frontier.png")
print("wrote", OUT / "squeezing_frontier_od.png")
