"""
Merge the per-T-value chunk .npz files (analysis/squeezing_frontier_ideal_v4_chunks/
T_<T>/squeezing_frontier.npz) into one combined ideal (QE=1, loss=0) wide-dense
grid, matching v3's Delta in [-6,8] step 0.05, T in [50,220] step 2.5 range.
Each chunk covers the FULL Delta axis for ONE T value (chunked this way because
long-running single background processes were unreliable in this sandbox; see
SKILL.md / v4 report appendix). Recomputes the global optimum, gap-gate
guardrail, and efficient frontier over the merged 2D grid, and writes:
  analysis/squeezing_frontier_ideal_v4_merged.npz
  docs/squeezing_report/squeezing_frontier_ideal_qe1_v4.png
"""
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
CHUNK_DIR = REPO / "analysis" / "squeezing_frontier_ideal_v4_chunks"
OUT_NPZ = REPO / "analysis" / "squeezing_frontier_ideal_v4_merged.npz"
OUT_DIR = Path(__file__).resolve().parent

sys.path.insert(0, str(REPO))
from gabes.schemes import fwm  # noqa: E402

T_LIST = np.round(np.arange(50.0, 220.0 + 1e-9, 2.5), 4)
GAP_MIN, GAP_MAX = 0.5, 1.5

D_ref = None
nD = None
Xi = None
Gs = None
Gc = None
Gap = None
OD = None
Dlt = None
Edge = None
GapBad = None

for iT, T in enumerate(T_LIST):
    d = np.load(CHUNK_DIR / f"T_{T}" / "squeezing_frontier.npz")
    if D_ref is None:
        D_ref = d["delta_ghz"]
        nD = D_ref.size
        nT = T_LIST.size
        Xi = np.full((nD, nT), np.nan)
        Gs = np.full((nD, nT), np.nan)
        Gc = np.full((nD, nT), np.nan)
        Gap = np.full((nD, nT), np.nan)
        OD = np.zeros((nD, nT))
        Dlt = np.full((nD, nT), np.nan)
        Edge = np.zeros((nD, nT), dtype=bool)
        GapBad = np.zeros((nD, nT), dtype=bool)
    assert np.allclose(d["delta_ghz"], D_ref), f"Delta axis mismatch at T={T}"
    assert d["temp_c"].size == 1 and abs(float(d["temp_c"][0]) - T) < 1e-6
    Xi[:, iT] = d["xi_dB"][:, 0]
    Gs[:, iT] = d["G_s"][:, 0]
    Gc[:, iT] = d["G_c"][:, 0]
    Gap[:, iT] = d["gap"][:, 0]
    OD[:, iT] = d["in_cell_od"][:, 0]
    Dlt[:, iT] = d["delta_mhz"][:, 0]
    Edge[:, iT] = d["edge"][:, 0]
    GapBad[:, iT] = d["gap_bad"][:, 0]

print(f"merged grid: {nD} OPD x {T_LIST.size} T = {nD * T_LIST.size} points")
print(f"NaN/Inf check: Xi {int(np.sum(~np.isfinite(Xi)))}, "
      f"Gs {int(np.sum(~np.isfinite(Gs)))}, Gc {int(np.sum(~np.isfinite(Gc)))}")

eta = 1.0
det_floor = 10.0 * np.log10(max(1.0 - eta, 1e-300))
Untrusted = Edge | GapBad
Xi_trusted = np.where(Untrusted, np.nan, Xi)

n_gap_bad = int(np.sum(GapBad))
n_edge = int(np.sum(Edge))
n_below = int(np.sum(Xi[np.isfinite(Xi)] < det_floor - 1e-6))

oFlat = int(np.nanargmin(Xi_trusted))
oiD, oiT = np.unravel_index(oFlat, Xi.shape)
xi_opt = float(Xi[oiD, oiT])
opt_at_edge = bool(oiD in (0, nD - 1) or oiT in (0, T_LIST.size - 1))

print(f"\n================  MERGED BEST ACHIEVABLE (ideal QE=1, loss=0)  ================")
print(f"deepest xi on the full grid = {xi_opt:.4f} dB @ Delta={D_ref[oiD]:+.2f} GHz, "
      f"T={T_LIST[oiT]:.1f} C  (G_s={Gs[oiD,oiT]:.4g}, G_c={Gc[oiD,oiT]:.4g}, "
      f"gap={Gap[oiD,oiT]:+.4f}, in-cell OD={OD[oiD,oiT]:.4f}, "
      f"delta={Dlt[oiD,oiT]:.1f} MHz)")
print(f"optimum on grid boundary: {opt_at_edge}")
print(f"points below detection floor (must be 0, floor is -inf here so trivial): {n_below}")
print(f"gap_bad points (excluded, dtype=bool True count): {n_gap_bad} / {Xi.size}")
print(f"edge points (excluded): {n_edge} / {Xi.size}")

# Efficient frontier: lowest-T (per Delta) within eps of xi_opt, excluding untrusted.
eps = 0.1
reached_all = Xi <= xi_opt + eps
reached_ok = reached_all & (~Untrusted)
Tmin_iT = np.full(nD, -1, dtype=int)
for iD in range(nD):
    idx = np.where(reached_ok[iD, :])[0]
    if idx.size:
        Tmin_iT[iD] = int(idx[0])
cands = [(T_LIST[jT], Gs[iD, jT], iD, jT)
         for iD, jT in enumerate(Tmin_iT) if jT >= 0]
eff_iD = eff_jT = -1
if cands:
    _, _, eff_iD, eff_jT = min(cands, key=lambda c: (c[0], c[1]))
if eff_iD >= 0:
    print(f"\nefficient frontier (min-T, eps={eps}dB): Delta={D_ref[eff_iD]:+.2f} GHz, "
          f"T={T_LIST[eff_jT]:.1f} C, G_s={Gs[eff_iD,eff_jT]:.4g}, "
          f"xi={Xi[eff_iD,eff_jT]:.4f} dB (gap to opt +{Xi[eff_iD,eff_jT]-xi_opt:.4f})")
else:
    print("\nno trustworthy frontier point found")

np.savez(OUT_NPZ, delta_ghz=D_ref, temp_c=T_LIST, xi_dB=Xi, G_s=Gs, G_c=Gc,
         gap=Gap, in_cell_od=OD, delta_mhz=Dlt, edge=Edge, gap_bad=GapBad,
         det_floor_dB=det_floor, xi_opt_dB=xi_opt, opt_iD=oiD, opt_iT=oiT,
         eff_iD=eff_iD, eff_jT=eff_jT, eps=eps, n_gap_bad=n_gap_bad, n_edge=n_edge)
print(f"\nsaved: {OUT_NPZ}")

# ---- figure ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, (axM, axF, axT) = plt.subplots(1, 3, figsize=(17, 4.8))
extent = [T_LIST[0], T_LIST[-1], D_ref[0], D_ref[-1]]
Xi_clip = np.clip(Xi, -100, 0)
im = axM.imshow(Xi_clip, origin="lower", aspect="auto", extent=extent, cmap="magma_r")
fig.colorbar(im, ax=axM, label="xi [dB] (clip -100..0)")
axM.scatter([T_LIST[oiT]], [D_ref[oiD]], color="cyan", marker="o", s=90,
            edgecolor="k", zorder=6, label=f"global optimum ({xi_opt:.2f} dB)")
if eff_iD >= 0:
    axM.scatter([T_LIST[eff_jT]], [D_ref[eff_iD]], color="lime", marker="*", s=150,
                edgecolor="k", zorder=7, label="min-T frontier point")
axM.set_xlabel("T [C]"); axM.set_ylabel("Delta [GHz]")
axM.set_title("Ideal (QE=1, loss=0) xi(Delta,T) -- v4, corrected physics")
axM.legend(loc="upper right", fontsize=7)

best_per_T = np.nanmin(Xi_trusted, axis=0)
best_D_per_T = D_ref[np.nanargmin(np.where(np.isfinite(Xi_trusted), Xi_trusted, np.inf), axis=0)]
axF.plot(T_LIST, best_per_T, "o-", color="#1f77b4", ms=3, lw=1.2)
axF.set_xlabel("T [C]"); axF.set_ylabel("best xi [dB] at that T")
axF.set_title("Best achievable xi per temperature")
axF.grid(alpha=0.3)

sl = oiD
axT.plot(T_LIST, Xi[sl, :], "o-", color="#2ca02c", lw=1.4, ms=3)
axT.scatter([T_LIST[oiT]], [xi_opt], color="gold", marker="o", s=80, edgecolor="k",
            zorder=6, label=f"optimum {xi_opt:.2f} dB")
axT.set_title(f"xi vs T at Delta = {D_ref[sl]:+.2f} GHz")
axT.set_xlabel("T [C]"); axT.set_ylabel("xi [dB]")
axT.grid(alpha=0.3); axT.legend(fontsize=8)

fig.tight_layout()
fig.savefig(OUT_DIR / "squeezing_frontier_ideal_qe1_v4.png", dpi=130)
print(f"saved: {OUT_DIR / 'squeezing_frontier_ideal_qe1_v4.png'}")
