"""Build the figures for the consolidated English paper (squeezing_paper.tex).

Pure post-processing: reads archived scan data (npz) and re-renders it with
neutral publication labels. No solver calls, no engine imports.

    python docs/squeezing_report/make_paper_figs.py

Outputs into docs/squeezing_report/paper_figs/ :
  fig1_floor_theorem.png    -- analytic detection-floor curves (Stage I)
  fig2_stage1_plateau.png   -- copy of the archived Stage-I plateau map
  fig3_stage2_od.png        -- copy of the archived Stage-II in-cell OD figure
  fig4_formula_bug.png      -- old vs corrected noise formula on identical gains
  fig5_excess_noise.png     -- Stage-IV hardening: v4-era vs v5-era maps + conjugate OD
  fig6_final_frontier.png   -- final finite-loss frontier (Stage V)
  fig7_ideal_limit.png      -- final ideal-detection source limit + no-excess contrast
  fig8_geometry.png         -- TPD/angle detail and low-angle ridge
  fig9_tolerance.png        -- tolerance scans around the recommended point

Historical provenance: the original Stage-I/II scan archives were later
overwritten by small smoke-test grids, so figs 2-3 reuse the archived PNGs
produced at those stages (their labels are already neutral); everything else is
re-rendered from the surviving full-resolution npz archives.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ANA = REPO / "analysis"
OUT = HERE / "paper_figs"
OUT.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 9.5,
    "axes.labelsize": 9,
    "legend.fontsize": 7.5,
    "figure.dpi": 100,
})
DPI = 170


def s_ideal_corrected(gs, gc):
    """Corrected intensity-difference noise (photon-number difference is
    conserved in lossless parametric gain): S = (Gs-Gc)^2 / (Gs+Gc)."""
    return (gs - gc) ** 2 / (gs + gc)


def s_ideal_old(gs, gc):
    """Deprecated pre-correction formula (Stage III bug): suffers catastrophic
    cancellation for Gs ~ Gc >> 1 and has the wrong high-gain asymptote."""
    return (gs + gc - 2.0 * np.sqrt(np.maximum(gs * gc - 1.0, 0.0))) / (gs + gc)


def xi_db(s):
    return 10.0 * np.log10(np.maximum(s, 1e-300))


# ---------------------------------------------------------------- fig 1
def fig1_floor_theorem():
    gs = np.logspace(0.01, 4, 400)
    gc = gs - 1.0
    s_src = s_ideal_corrected(gs, gc)

    fig, (a, b) = plt.subplots(1, 2, figsize=(9.6, 3.6))
    etas = [0.8694, 0.92, 0.95, 0.99, 1.0]
    for eta in etas:
        s = eta * s_src + (1.0 - eta)
        lab = rf"$\eta={eta:g}$"
        (ln,) = a.plot(gs, xi_db(s), lw=1.6, label=lab)
        if eta < 1.0:
            a.axhline(xi_db(1.0 - eta), color=ln.get_color(), ls=":", lw=0.9)
    a.set_xscale("log")
    a.set_ylim(-45, 1)
    a.set_xlabel(r"probe gain $G_s$  (balanced twin beams, $G_c=G_s-1$)")
    a.set_ylabel(r"$\xi$ [dB]")
    a.set_title("(a) squeezing vs gain: collapse onto the detection floor")
    a.legend(loc="lower left")
    a.grid(alpha=0.3)

    eta_ax = np.linspace(0.80, 0.999, 300)
    b.plot(eta_ax, 10 * np.log10(1 - eta_ax), lw=1.8, color="#333333")
    marks = [(0.8694, "reference detection"), (0.92, "unit optical loss (QE only)")]
    for ev, lab in marks:
        b.scatter([ev], [10 * np.log10(1 - ev)], zorder=5, s=35,
                  label=f"{lab}: {10 * np.log10(1 - ev):.2f} dB")
    b.set_xlabel(r"total detection efficiency $\eta$")
    b.set_ylabel(r"floor $10\log_{10}(1-\eta)$ [dB]")
    b.set_title("(b) the floor moves only with efficiency")
    b.legend(loc="upper right")
    b.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUT / "fig1_floor_theorem.png", dpi=DPI)
    plt.close(fig)
    print("fig1 done")


# ---------------------------------------------------------------- figs 2-3 (archived)
def figs_2_3_archived():
    shutil.copyfile(HERE / "squeezing_map.png", OUT / "fig2_stage1_plateau.png")
    shutil.copyfile(HERE / "squeezing_frontier_od.png", OUT / "fig3_stage2_od.png")
    print("fig2, fig3 copied from stage archives")


# ---------------------------------------------------------------- fig 4
def fig4_formula_bug():
    """Bare source-noise formulas evaluated on IDENTICAL archived gain pairs
    (Gs, Gc) from the wide ideal-detection scan. This mirrors the diagnostic
    that exposed the bug: same data, two formulas, ~50 dB apart."""
    d = np.load(ANA / "squeezing_frontier_ideal_v4_merged.npz")
    D, T = d["delta_ghz"], d["temp_c"]
    gs, gc = d["G_s"], d["G_c"]

    xi_old = xi_db(s_ideal_old(gs, gc))
    xi_new = xi_db(s_ideal_corrected(gs, gc))

    fig, (a, b, c) = plt.subplots(1, 3, figsize=(13.2, 3.9))
    extent = [T[0], T[-1], D[0], D[-1]]

    im = a.imshow(np.clip(xi_old, -90, 0), origin="lower", aspect="auto",
                  extent=extent, cmap="magma_r")
    fig.colorbar(im, ax=a, label=r"$\xi_{\rm src}$ [dB]")
    hiD, hiT = np.unravel_index(int(np.argmax(gs)), gs.shape)
    a.scatter([T[hiT]], [D[hiD]], marker="o", s=60, facecolor="none",
              edgecolor="cyan", lw=1.6,
              label=f"highest gain: {xi_old[hiD, hiT]:.1f} dB (spurious)")
    a.set_xlabel(r"$T$ [$^\circ$C]"); a.set_ylabel(r"$\Delta$ [GHz]")
    a.set_title("(a) deprecated noise formula")
    a.legend(loc="upper left")

    im = b.imshow(np.clip(xi_new, -35, 0), origin="lower", aspect="auto",
                  extent=extent, cmap="magma_r")
    fig.colorbar(im, ax=b, label=r"$\xi_{\rm src}$ [dB]")
    b.scatter([T[hiT]], [D[hiD]], marker="o", s=60, facecolor="none",
              edgecolor="cyan", lw=1.6,
              label=f"same point: {xi_new[hiD, hiT]:.1f} dB")
    b.set_xlabel(r"$T$ [$^\circ$C]"); b.set_ylabel(r"$\Delta$ [GHz]")
    b.set_title("(b) corrected formula, identical gain data")
    b.legend(loc="upper left")

    iD = int(np.argmin(np.abs(D - (-2.10))))
    c.plot(T, xi_old[iD], "s--", ms=3, lw=1.2, color="#d62728",
           label="deprecated formula")
    c.plot(T, xi_new[iD], "o-", ms=3, lw=1.5, color="#1f77b4",
           label="corrected formula")
    k = int(np.argmin(np.abs(T - 142.5)))
    c.annotate(f"{xi_old[iD, k]:.1f} dB", (T[k], xi_old[iD, k]),
               textcoords="offset points", xytext=(8, -3), fontsize=8,
               color="#d62728")
    c.annotate(f"{xi_new[iD, k]:.1f} dB", (T[k], xi_new[iD, k]),
               textcoords="offset points", xytext=(-52, 6), fontsize=8,
               color="#1f77b4")
    c.set_xlabel(r"$T$ [$^\circ$C]"); c.set_ylabel(r"$\xi_{\rm src}$ [dB]")
    c.set_title(rf"(c) slice at $\Delta = {D[iD]:+.2f}$ GHz")
    c.grid(alpha=0.3); c.legend(loc="center left")

    fig.tight_layout()
    fig.savefig(OUT / "fig4_formula_bug.png", dpi=DPI)
    plt.close(fig)
    print(f"fig4 done  (max-gain point: old {xi_old[hiD, hiT]:.2f} dB, "
          f"new {xi_new[hiD, hiT]:.2f} dB; slice@142.5C old {xi_old[iD, k]:.2f}, "
          f"new {xi_new[iD, k]:.2f})")


# ---------------------------------------------------------------- fig 5
def fig5_excess_noise():
    """Stage-IV contrast at ideal detection (eta=1), both sides genuinely
    pre-beat-correction. Panel (a): the pre-hardening archive, analytically
    remapped from its stored eta=0.8694 readout to eta=1 (detection
    efficiency is pure post-processing, S_cell = (S - (1-eta))/eta).
    Panel (b): the archived Stage-IV ideal scan with the extended
    excess-noise model. The finite-loss Stage-IV archive was later
    overwritten with beat-corrected data, so the finite-detection numbers
    are quoted in the text from the stage record instead."""
    d4 = np.load(ANA / "squeezing_frontier_finite_loss_v4" / "squeezing_frontier.npz")
    d5 = np.load(ANA / "squeezing_frontier_ideal_v5" / "squeezing_frontier.npz")

    eta4 = float(d4["eta"])
    s4 = 10.0 ** (d4["xi_dB"] / 10.0)
    s4_cell = np.maximum((s4 - (1.0 - eta4)) / eta4, 1e-30)
    xi4_ideal = 10.0 * np.log10(s4_cell)

    fig, (a, b, c) = plt.subplots(1, 3, figsize=(13.2, 3.9))

    D4, T4 = d4["delta_ghz"], d4["temp_c"]
    o4 = int(np.nanargmin(xi4_ideal))
    o4D, o4T = np.unravel_index(o4, xi4_ideal.shape)
    im = a.imshow(np.clip(xi4_ideal, -37, 0), origin="lower", aspect="auto",
                  extent=[T4[0], T4[-1], D4[0], D4[-1]], cmap="viridis_r")
    fig.colorbar(im, ax=a, label=r"$\xi$ [dB]")
    a.scatter([T4[o4T]], [D4[o4D]], marker="*", s=150, color="cyan",
              edgecolor="k", zorder=6,
              label=(f"optimum {xi4_ideal[o4D, o4T]:.1f} dB, "
                     f"$G_s={d4['G_s'][o4D, o4T]:.0f}$"))
    a.axhline(0.9, color="w", ls=":", lw=1.0)
    a.set_xlabel(r"$T$ [$^\circ$C]"); a.set_ylabel(r"$\Delta$ [GHz]")
    a.set_title("(a) no excess-noise channels: runaway")
    a.legend(loc="upper left")

    D5, T5 = d5["delta_ghz"], d5["temp_c"]
    xi5 = d5["xi_dB"]
    o5D, o5T = int(d5["opt_iD"]), int(d5["opt_iT"])
    im = b.imshow(np.clip(xi5, -37, 0), origin="lower", aspect="auto",
                  extent=[T5[0], T5[-1], D5[0], D5[-1]], cmap="viridis_r")
    fig.colorbar(im, ax=b, label=r"$\xi$ [dB]")
    b.scatter([T5[o5T]], [D5[o5D]], marker="*", s=150, color="cyan",
              edgecolor="k", zorder=6,
              label=(f"optimum {float(d5['xi_opt_dB']):.1f} dB, "
                     f"$G_s={d5['G_s'][o5D, o5T]:.0f}$"))
    b.axhline(0.9, color="w", ls=":", lw=1.0)
    b.set_xlabel(r"$T$ [$^\circ$C]"); b.set_ylabel(r"$\Delta$ [GHz]")
    b.set_title("(b) extended excess-noise model: capped")
    b.legend(loc="upper left")

    im = c.imshow(np.clip(d5["od_conj"], 0, 1.0), origin="lower", aspect="auto",
                  extent=[T5[0], T5[-1], D5[0], D5[-1]], cmap="cividis")
    fig.colorbar(im, ax=c, label="optical depth (clipped at 1)")
    c.scatter([135], [-2.20], marker="o", s=60, facecolor="none",
              edgecolor="crimson", lw=1.6, label="former high-gain optimum")
    c.set_xlabel(r"$T$ [$^\circ$C]"); c.set_ylabel(r"$\Delta$ [GHz]")
    c.set_title("(c) conjugate-arm linear absorption")
    c.legend(loc="upper left")

    fig.tight_layout()
    fig.savefig(OUT / "fig5_excess_noise.png", dpi=DPI)
    plt.close(fig)
    print(f"fig5 done  (bare ideal opt {xi4_ideal[o4D, o4T]:.2f} dB @ "
          f"{D4[o4D]:+.2f}/{T4[o4T]:.0f}, Gs={d4['G_s'][o4D, o4T]:.0f}; "
          f"hardened ideal opt {float(d5['xi_opt_dB']):.2f} dB @ "
          f"{D5[o5D]:+.2f}/{T5[o5T]:.0f}, Gs={d5['G_s'][o5D, o5T]:.0f})")


# ---------------------------------------------------------------- fig 6
def fig6_final_frontier():
    d = np.load(ANA / "squeezing_frontier_finite_loss_v6" / "squeezing_frontier.npz")
    D, T, xi = d["delta_ghz"], d["temp_c"], d["xi_dB"]
    oiD, oiT = int(d["opt_iD"]), int(d["opt_iT"])
    floor = float(d["det_floor_dB"])

    fig, (a, b, c) = plt.subplots(1, 3, figsize=(13.2, 3.9))
    extent = [T[0], T[-1], D[0], D[-1]]

    im = a.imshow(xi, origin="lower", aspect="auto", extent=extent,
                  cmap="viridis_r", vmax=0.0)
    fig.colorbar(im, ax=a, label=r"$\xi$ [dB]")
    a.scatter([T[oiT]], [D[oiD]], marker="*", s=160, color="cyan", edgecolor="k",
              zorder=6, label=f"optimum {float(d['xi_opt_dB']):.2f} dB")
    a.axhline(0.9, color="w", ls=":", lw=1.0, label="reference expt. $+0.9$ GHz")
    a.set_xlabel(r"$T$ [$^\circ$C]"); a.set_ylabel(r"$\Delta$ [GHz]")
    a.set_title(r"(a) final $\xi(\Delta,T)$, realistic detection")
    a.legend(loc="upper left")

    im = b.imshow(np.clip(d["od_conj"], 0, 1.0), origin="lower", aspect="auto",
                  extent=extent, cmap="cividis")
    fig.colorbar(im, ax=b, label="optical depth (clipped at 1)")
    b.scatter([T[oiT]], [D[oiD]], marker="*", s=120, color="cyan", edgecolor="k",
              zorder=6)
    b.set_xlabel(r"$T$ [$^\circ$C]"); b.set_ylabel(r"$\Delta$ [GHz]")
    b.set_title("(b) conjugate-arm linear absorption")

    slices = [(float(D[oiD]), "optimum"), (1.4, "best blue lobe"),
              (0.9, "reference expt."), (-2.0, "Stage-IV optimum")]
    for dv, lab in slices:
        i = int(np.argmin(np.abs(D - dv)))
        c.plot(T, xi[i, :], lw=1.5, label=rf"$\Delta={D[i]:+.1f}$ GHz ({lab})")
    c.axhline(floor, color="crimson", ls="--", lw=1.1,
              label=f"detection floor {floor:.2f} dB")
    c.set_xlabel(r"$T$ [$^\circ$C]"); c.set_ylabel(r"$\xi$ [dB]")
    c.set_title(r"(c) $\xi(T)$ at representative detunings")
    c.grid(alpha=0.3); c.legend(loc="lower left")

    fig.tight_layout()
    fig.savefig(OUT / "fig6_final_frontier.png", dpi=DPI)
    plt.close(fig)
    print("fig6 done")


# ---------------------------------------------------------------- fig 7
def fig7_ideal_limit():
    d = np.load(ANA / "squeezing_frontier_ideal_v6" / "squeezing_frontier.npz")
    D, T = d["delta_ghz"], d["temp_c"]
    xi, gs, gc, od = d["xi_dB"], d["G_s"], d["G_c"], d["in_cell_od"]
    untrusted = d["edge"] | d["gap_bad"]
    xit = np.where(untrusted, np.nan, xi)
    o = int(np.nanargmin(xit))
    oiD, oiT = np.unravel_index(o, xi.shape)
    xi_opt = float(xi[oiD, oiT])

    # contrast: same grid evaluated WITHOUT the excess-noise channels
    tau = np.exp(-np.maximum(od, 0.0))
    xi_bare = xi_db(tau * s_ideal_corrected(gs, gc) + (1.0 - tau))
    xi_bare_t = np.where(untrusted, np.nan, xi_bare)

    fig, (a, b, c) = plt.subplots(1, 3, figsize=(13.2, 3.9))
    extent = [T[0], T[-1], D[0], D[-1]]

    im = a.imshow(np.clip(xi, -18, 0), origin="lower", aspect="auto",
                  extent=extent, cmap="magma_r")
    fig.colorbar(im, ax=a, label=r"$\xi$ [dB] (clipped at $-18$)")
    a.scatter([T[oiT]], [D[oiD]], marker="o", s=80, color="cyan", edgecolor="k",
              zorder=6, label=f"optimum {xi_opt:.2f} dB")
    a.axhline(0.9, color="w", ls=":", lw=1.0)
    a.set_xlabel(r"$T$ [$^\circ$C]"); a.set_ylabel(r"$\Delta$ [GHz]")
    a.set_title(r"(a) $\xi(\Delta,T)$, ideal detection ($\eta=1$)")
    a.legend(loc="upper left")

    b.plot(T, np.nanmin(xit, axis=0), "o-", ms=3, lw=1.5, color="#1f77b4",
           label="full excess-noise model")
    b.plot(T, np.nanmin(xi_bare_t, axis=0), "s--", ms=3, lw=1.2, color="#d62728",
           label="excess-noise channels removed")
    b.scatter([T[oiT]], [xi_opt], s=60, color="cyan", edgecolor="k", zorder=6)
    b.set_xlabel(r"$T$ [$^\circ$C]"); b.set_ylabel(r"deepest $\xi$ at each $T$ [dB]")
    b.set_title("(b) absorption caps the ideal optimum")
    b.grid(alpha=0.3); b.legend(loc="lower left")

    c.plot(T, xi[oiD, :], "o-", ms=3, lw=1.5, color="#2ca02c",
           label="full excess-noise model")
    c.plot(T, xi_bare[oiD, :], "s--", ms=3, lw=1.2, color="#d62728",
           label="excess-noise channels removed")
    c.scatter([T[oiT]], [xi_opt], s=70, color="gold", edgecolor="k", zorder=6)
    c.set_xlabel(r"$T$ [$^\circ$C]"); c.set_ylabel(r"$\xi$ [dB]")
    c.set_title(rf"(c) $\xi(T)$ at $\Delta={D[oiD]:+.2f}$ GHz")
    c.grid(alpha=0.3); c.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(OUT / "fig7_ideal_limit.png", dpi=DPI)
    plt.close(fig)
    print(f"fig7 done  (ideal opt {xi_opt:.2f} dB @ {D[oiD]:+.2f} GHz, {T[oiT]:.0f} C)")


# ---------------------------------------------------------------- fig 8
def fig8_geometry():
    dd = np.load(ANA / "squeezing_frontier_ideal_v6_tpd_angle_detail"
                 / "squeezing_tpd_angle_detail.npz")
    dr = np.load(ANA / "squeezing_frontier_ideal_v6_low_angle_ridge"
                 / "squeezing_low_angle_ridge.npz")

    ang_d, tpd_d = dd["angle_deg"], dd["tpd_mhz"]
    xi_d = np.where(dd["trusted"], dd["xi_dB"], np.nan)      # (D,T,ang,tpd)
    import warnings
    warnings.filterwarnings("ignore", "All-NaN")  # untrusted-only columns
    xi_at = np.nanmin(xi_d, axis=(0, 1))                     # (ang,tpd)

    fig, (a, b) = plt.subplots(1, 2, figsize=(10.2, 3.9))

    im = a.imshow(xi_at, origin="lower", aspect="auto",
                  extent=[tpd_d[0], tpd_d[-1], ang_d[0], ang_d[-1]],
                  cmap="magma_r")
    fig.colorbar(im, ax=a, label=r"$\xi$ [dB]")
    ia = int(np.argmin(np.abs(ang_d - 0.32)))
    jt = int(np.nanargmin(xi_at[ia]))
    a.scatter([tpd_d[jt]], [ang_d[ia]], color="w", edgecolor="k", s=55, zorder=6,
              label=f"default geometry: {xi_at[ia, jt]:.2f} dB")
    fa, ft = np.unravel_index(int(np.nanargmin(xi_at)), xi_at.shape)
    a.scatter([tpd_d[ft]], [ang_d[fa]], color="cyan", edgecolor="k", s=55,
              zorder=6, label=f"scan best: {xi_at[fa, ft]:.2f} dB")
    a.set_xlabel(r"two-photon detuning $\delta$ [MHz]")
    a.set_ylabel(r"pump-probe angle $\theta$ [$^\circ$]")
    a.set_title(r"(a) $\xi(\theta,\delta)$, minimized over $(\Delta,T)$")
    a.legend(loc="upper left")

    best_d = np.nanmin(xi_d, axis=(0, 1, 3))
    ang_r = dr["angle_deg"]
    xi_r = np.where(dr["trusted"], dr["xi_dB"], np.nan)
    best_r = np.nanmin(xi_r, axis=(0, 1, 3))
    b.plot(ang_d, best_d, "o-", ms=4, lw=1.5, color="#1f77b4",
           label="local detail scan")
    b.plot(ang_r, best_r, "s--", ms=4, lw=1.3, color="#9467bd",
           label="low-angle follow-up scan")
    b.axvline(0.32, color="k", ls=":", lw=1.0, label=r"default $0.32^\circ$")
    b.set_xlabel(r"pump-probe angle $\theta$ [$^\circ$]")
    b.set_ylabel(r"deepest accepted $\xi$ [dB]")
    b.set_title("(b) ideal-detection limit vs beam angle")
    b.grid(alpha=0.3); b.legend(loc="upper left")

    fig.tight_layout()
    fig.savefig(OUT / "fig8_geometry.png", dpi=DPI)
    plt.close(fig)
    print(f"fig8 done  (default-angle best {xi_at[ia, jt]:.2f} dB, "
          f"ridge best {np.nanmin(best_r):.2f} dB)")


# ---------------------------------------------------------------- fig 9
def fig9_tolerance():
    d = np.load(ANA / "squeezing_v6_tolerance" / "squeezing_v6_tolerance.npz")
    # center record: (Delta_GHz, T_C, delta_MHz, P_probe_W, loss_frac, xi_dB)
    center = d["center"]
    xi_c = float(center[-1])
    print("tolerance center xi:", xi_c)

    fig, axes = plt.subplots(2, 3, figsize=(12.6, 6.6))
    (a_opd, a_tpd, a_T), (a_p, a_loss, a_leg) = axes

    def band(ax):
        ax.axhspan(xi_c, xi_c + 0.5, color="tab:green", alpha=0.12, lw=0)
        ax.axhline(xi_c, color="k", lw=0.8, alpha=0.5)

    def fixed_best(ax, x, xf, xb, tf, xlab, title, center_x=None):
        band(ax)
        ax.plot(x, xf, "-", color="#1f77b4", lw=1.6, label=r"fixed $\delta$ (pure drift)")
        bad = ~tf.astype(bool)
        if bad.any():
            ax.scatter(np.asarray(x)[bad], np.asarray(xf)[bad], marker="x",
                       color="crimson", s=30, zorder=6, label="outside trust gate")
        ax.plot(x, xb, "--", color="#2ca02c", lw=1.4, label=r"re-optimized $\delta$")
        if center_x is not None:
            ax.axvline(center_x, color="k", ls=":", lw=0.9)
        ax.set_xlabel(xlab); ax.set_ylabel(r"$\xi$ [dB]")
        ax.set_title(title); ax.grid(alpha=0.3)

    fixed_best(a_opd, (d["opd_ghz"] - (-1.50)) * 1e3, d["xi_opd_fixed"],
               d["xi_opd_best"], d["trust_opd_fixed"],
               r"one-photon detuning drift $\Delta-\Delta_0$ [MHz]",
               "(a) one-photon detuning", 0.0)

    band(a_tpd)
    gate = (d["gap_tpd"] >= float(d["gap_min"])) & (d["gap_tpd"] <= float(d["gap_max"]))
    a_tpd.plot(d["tpd_mhz"], d["xi_tpd"], "-", color="#1f77b4", lw=1.2,
               label="fine scan")
    a_tpd.scatter(d["tpd_mhz"][~gate], d["xi_tpd"][~gate], marker="x",
                  color="crimson", s=14, zorder=6, label="outside trust gate")
    a_tpd.axvline(-280.0, color="k", ls=":", lw=0.9)
    a_tpd.set_xlim(-360, -200)
    a_tpd.set_ylim(xi_c - 0.6, xi_c + 3.0)
    a_tpd.set_xlabel(r"two-photon detuning $\delta$ [MHz]")
    a_tpd.set_ylabel(r"$\xi$ [dB]")
    a_tpd.set_title("(b) two-photon detuning (5 MHz grid)")
    a_tpd.grid(alpha=0.3)

    fixed_best(a_T, d["temp_c"], d["xi_temp_fixed"], d["xi_temp_best"],
               d["trust_temp_fixed"], r"cell temperature $T$ [$^\circ$C]",
               "(c) cell temperature", 110.0)

    band(a_p)
    a_p.plot(d["pprobe_factor"], d["xi_pprobe_fixed"], "o-", ms=3, lw=1.4,
             color="#1f77b4")
    a_p.set_xscale("log", base=2)
    a_p.set_ylim(xi_c - 0.6, xi_c + 3.0)
    a_p.axvline(1.0, color="k", ls=":", lw=0.9)
    a_p.set_xlabel(r"probe power / 8 $\mu$W")
    a_p.set_ylabel(r"$\xi$ [dB]")
    a_p.set_title("(d) probe (seed) power")
    a_p.grid(alpha=0.3)

    band(a_loss)
    a_loss.plot(100 * d["loss_frac"], d["xi_loss"], "-", lw=1.6, color="#1f77b4")
    a_loss.axvline(100 * float(d["loss0"]), color="k", ls=":", lw=0.9)
    a_loss.set_xlabel("post-cell optical loss [%]")
    a_loss.set_ylabel(r"$\xi$ [dB]")
    a_loss.set_title("(e) detection-path loss")
    a_loss.grid(alpha=0.3)

    a_leg.axis("off")
    handles, labels = a_opd.get_legend_handles_labels()
    a_leg.legend(handles, labels, loc="center left", fontsize=9,
                 title="line styles (panels a, c)")
    a_leg.text(0.02, 0.24,
               "shaded band: within +0.5 dB of the optimum\n"
               "dotted vertical line: recommended operating value",
               transform=a_leg.transAxes, fontsize=8.5, va="top")

    fig.tight_layout()
    fig.savefig(OUT / "fig9_tolerance.png", dpi=DPI)
    plt.close(fig)
    print("fig9 done")


if __name__ == "__main__":
    fig1_floor_theorem()
    figs_2_3_archived()
    fig4_formula_bug()
    fig5_excess_noise()
    fig6_final_frontier()
    fig7_ideal_limit()
    fig8_geometry()
    fig9_tolerance()
    print("all paper figures written to", OUT)
