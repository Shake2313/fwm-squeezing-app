"""
Build the v6 hardened-Ultra report figures. The v6 finite-loss frontier is the
SAME grid/solve as v4 (Delta in [-3.5,3.0] GHz x T in [60,150] C, coarse=81,
vstep=5) after the seeded sideband beat-sign correction, with the hardened Ultra
noise model intrinsic to the Ultra path (conjugate-arm linear OD + pump scatter).

    python docs/squeezing_report/make_v6_hardened_figs.py

Outputs into docs/squeezing_report/ (v6-suffixed):
  * squeezing_frontier_v6.png            -- copy of the v6 scan's 3-panel figure
  * squeezing_frontier_od_v6.png         -- in-cell OD + conjugate OD + xi(T) slices
  * squeezing_frontier_hardened_compare_v6.png -- v4 legacy vs v6 hardened xi(Delta,T)
"""
import shutil
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
V6 = REPO / "analysis" / "squeezing" / "squeezing_frontier_finite_loss_v6" / "squeezing_frontier.npz"
V4 = REPO / "analysis" / "squeezing" / "squeezing_frontier_finite_loss_v4" / "squeezing_frontier.npz"

d = np.load(V6)
D, T = d["delta_ghz"], d["temp_c"]
Xi, OD = d["xi_dB"], d["in_cell_od"]
OdConj = d["od_conj"] if "od_conj" in d.files else np.zeros_like(Xi)
det_floor = float(d["det_floor_dB"])
xi_opt = float(d["xi_opt_dB"])
oiD, oiT = int(d["opt_iD"]), int(d["opt_iT"])

# 1) copy the scan's own 3-panel figure
src_png = REPO / "analysis" / "squeezing" / "squeezing_frontier_finite_loss_v6" / "squeezing_frontier.png"
if src_png.exists():
    shutil.copyfile(src_png, OUT / "squeezing_frontier_v6.png")
    print("wrote", OUT / "squeezing_frontier_v6.png")

# 2) OD + conjugate OD map + xi(T) slices
fig, (axOD, axC, axS) = plt.subplots(1, 3, figsize=(17, 4.8))
extent = [T[0], T[-1], D[0], D[-1]]

im = axOD.imshow(np.clip(OD, 0, 2.0), origin="lower", aspect="auto",
                 extent=extent, cmap="inferno")
fig.colorbar(im, ax=axOD, label="probe in-cell OD (clip 2.0)")
axOD.scatter([T[oiT]], [D[oiD]], color="cyan", marker="o", s=70, edgecolor="k",
             zorder=6, label=f"optimum ({xi_opt:.2f} dB)")
axOD.set_xlabel("T [C]"); axOD.set_ylabel("Delta [GHz]")
axOD.set_title("probe in-cell OD(Delta, T)")
axOD.legend(loc="upper left", fontsize=8)

imc = axC.imshow(np.clip(OdConj, 0, 1.0), origin="lower", aspect="auto",
                 extent=extent, cmap="viridis")
fig.colorbar(imc, ax=axC, label="conjugate-arm linear OD (clip 1.0)")
axC.scatter([T[oiT]], [D[oiD]], color="crimson", marker="o", s=70, edgecolor="k",
            zorder=6)
axC.set_xlabel("T [C]"); axC.set_ylabel("Delta [GHz]")
axC.set_title("NEW hardened channel: conjugate OD(Delta, T)")

slices = sorted({round(float(D[oiD]), 1), 0.9, -2.2, -1.0, 0.0})
for dv in slices:
    i = int(np.argmin(np.abs(D - dv)))
    lab = "optimum" if i == oiD else ("Sim +0.9" if abs(dv - 0.9) < 1e-6 else
          ("v4 artifact" if abs(dv + 2.2) < 1e-6 else f"{D[i]:+.1f}"))
    axS.plot(T, Xi[i, :], "-", lw=1.6, label=f"Delta={D[i]:+.1f} ({lab})")
    j = int(np.nanargmin(Xi[i, :]))
    axS.scatter([T[j]], [Xi[i, j]], s=30, zorder=5, color=axS.lines[-1].get_color())
axS.axhline(det_floor, color="crimson", ls="--", lw=1.2,
            label=f"detection floor {det_floor:.2f} dB")
axS.set_xlabel("T [C]"); axS.set_ylabel("xi [dB]")
axS.set_title("hardened xi vs T")
axS.grid(alpha=0.3); axS.legend(fontsize=7, loc="lower left")
fig.tight_layout()
fig.savefig(OUT / "squeezing_frontier_od_v6.png", dpi=130)
print("wrote", OUT / "squeezing_frontier_od_v6.png")

# 3) v4 legacy vs v6 hardened xi(Delta,T) side by side
if V4.exists():
    d4 = np.load(V4)
    D4, T4, Xi4 = d4["delta_ghz"], d4["temp_c"], d4["xi_dB"]
    o4D, o4T = int(d4["opt_iD"]), int(d4["opt_iT"])
    xi4 = float(d4["xi_opt_dB"])
    fig2, (a4, a5) = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
    vmin = float(np.nanmin([np.nanmin(Xi4), np.nanmin(Xi)]))
    ext4 = [T4[0], T4[-1], D4[0], D4[-1]]
    for ax, XX, ex, oD, oT, DD, TT, xo, tag in [
            (a4, Xi4, ext4, o4D, o4T, D4, T4, xi4, "v4 legacy Ultra"),
            (a5, Xi, extent, oiD, oiT, D, T, xi_opt, "v6 hardened Ultra")]:
        im = ax.imshow(XX, origin="lower", aspect="auto", extent=ex,
                       cmap="viridis_r", vmin=vmin, vmax=0.0)
        ax.scatter([TT[oT]], [DD[oD]], color="cyan", marker="*", s=180,
                   edgecolor="k", zorder=6,
                   label=f"opt {xo:.2f} dB @ {DD[oD]:+.1f} GHz")
        ax.axhline(0.9, color="w", ls=":", lw=1.0)
        ax.set_xlabel("T [C]"); ax.set_title(tag)
        ax.legend(loc="upper left", fontsize=8)
    a4.set_ylabel("Delta [GHz]")
    fig2.colorbar(im, ax=[a4, a5], label="xi [dB]")
    fig2.savefig(OUT / "squeezing_frontier_hardened_compare_v6.png", dpi=130,
                 bbox_inches="tight")
    print("wrote", OUT / "squeezing_frontier_hardened_compare_v6.png")
