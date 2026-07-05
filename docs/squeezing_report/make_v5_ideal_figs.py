"""Build the v5 IDEAL source-limit (QE=1, loss=0) figure under the hardened Ultra
noise model. This is the missing v5 counterpart of v4's ideal-detection scan
(docs/squeezing_report/squeezing_frontier_ideal_qe1_v4.png): with the detection
floor removed, the only remaining limits are the hardened in-cell channels
(probe/conjugate arm linear absorption + pump scatter). The figure shows that the
hardened model caps the ideal optimum too -- it does NOT run away to the extreme
gain (Gs~636) that the gap-only v4-style objective still reaches on the same grid.

    python docs/squeezing_report/make_v5_ideal_figs.py

Output: docs/squeezing_report/squeezing_frontier_ideal_qe1_v5.png
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO))
from gabes import observables  # noqa: E402

d = np.load(REPO / "analysis" / "squeezing_frontier_ideal_v5" / "squeezing_frontier.npz")
D, T = d["delta_ghz"], d["temp_c"]
Xi, Gs, Gc, OD = d["xi_dB"], d["G_s"], d["G_c"], d["in_cell_od"]
Untrusted = d["edge"] | d["gap_bad"]
Xit = np.where(Untrusted, np.nan, Xi)

o = int(np.nanargmin(Xit))
oiD, oiT = np.unravel_index(o, Xi.shape)
xi_opt = float(Xi[oiD, oiT])

# gap-only (v4-style, NO excess noise) ideal remap on the SAME grid, for contrast.
S_ideal = observables.ideal_twin_beam_noise(Gs, Gc)
tau = np.exp(-np.maximum(OD, 0.0))
Xi_go = 10.0 * np.log10(np.maximum(tau * S_ideal + (1.0 - tau), 1e-30))
Xi_go_t = np.where(Untrusted, np.nan, Xi_go)
g = int(np.nanargmin(Xi_go_t))
giD, giT = np.unravel_index(g, Xi.shape)

fig, (axM, axF, axT) = plt.subplots(1, 3, figsize=(17, 4.8))
extent = [T[0], T[-1], D[0], D[-1]]

# (1) ideal hardened xi(Delta,T) heatmap + optimum
Xi_clip = np.clip(Xi, -18, 0)
im = axM.imshow(Xi_clip, origin="lower", aspect="auto", extent=extent, cmap="magma_r")
fig.colorbar(im, ax=axM, label="xi [dB] (clip -18..0)")
axM.scatter([T[oiT]], [D[oiD]], color="cyan", marker="o", s=95, edgecolor="k",
            zorder=6, label=f"ideal optimum ({xi_opt:.2f} dB)")
axM.axhline(0.9, color="w", ls=":", lw=1.0, label="Sim +0.9 GHz")
axM.set_xlabel("T [C]"); axM.set_ylabel("Delta [GHz]")
axM.set_title("Ideal (QE=1, loss=0) hardened xi(Delta,T)")
axM.legend(loc="upper right", fontsize=7)

# (2) best-achievable xi per T: hardened (capped) vs gap-only (runaway)
best_h = np.nanmin(Xit, axis=0)
best_go = np.nanmin(Xi_go_t, axis=0)
axF.plot(T, best_h, "o-", color="#1f77b4", ms=3, lw=1.5, label="hardened (v5)")
axF.plot(T, best_go, "s--", color="#d62728", ms=3, lw=1.2,
         label="gap-only (v4-style, no excess noise)")
axF.scatter([T[oiT]], [xi_opt], color="cyan", edgecolor="k", s=70, zorder=6)
axF.set_xlabel("T [C]"); axF.set_ylabel("best xi [dB] at that T")
axF.set_title("hardened caps the ideal optimum; gap-only runs away")
axF.grid(alpha=0.3); axF.legend(fontsize=8, loc="lower left")

# (3) xi vs T at the optimum Delta: hardened vs gap-only (same Gs, diverge via odc)
sl = oiD
axT.plot(T, Xi[sl, :], "o-", color="#2ca02c", lw=1.6, ms=3, label="hardened (v5)")
axT.plot(T, Xi_go[sl, :], "s--", color="#d62728", lw=1.2, ms=3, label="gap-only")
axT.scatter([T[oiT]], [xi_opt], color="gold", marker="o", s=80, edgecolor="k",
            zorder=6, label=f"ideal optimum {xi_opt:.2f} dB")
axT.set_title(f"xi vs T at Delta = {D[sl]:+.2f} GHz (ideal detection)")
axT.set_xlabel("T [C]"); axT.set_ylabel("xi [dB]")
axT.grid(alpha=0.3); axT.legend(fontsize=8, loc="upper right")

fig.tight_layout()
fig.savefig(OUT / "squeezing_frontier_ideal_qe1_v5.png", dpi=130)
print(f"ideal optimum (hardened)  : {xi_opt:.4f} dB @ Delta={D[oiD]:+.2f} GHz, "
      f"T={T[oiT]:.0f} C, Gs={Gs[oiD,oiT]:.4g}")
print(f"gap-only ideal on grid    : {Xi_go[giD,giT]:.4f} dB @ Delta={D[giD]:+.2f} GHz, "
      f"T={T[giT]:.0f} C, Gs={Gs[giD,giT]:.4g}  (hardened at same pt {Xi[giD,giT]:.4f} dB)")
print(f"wrote {OUT / 'squeezing_frontier_ideal_qe1_v5.png'}")
