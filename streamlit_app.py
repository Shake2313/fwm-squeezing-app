"""
Interactive front-end for the 85Rb D1 double-Λ four-wave-mixing model.

Run with:
    streamlit run streamlit_app.py

Backend physics lives in fwm_obe.py; this file only wires sliders to
compute_spectrum() and draws the gain / squeezing curves. Heavy solves are
memoised with st.cache_data, so dragging the two-photon-detuning (TPD) slider
is instant — it just navigates an already-computed curve — while changing the
one-photon detuning (OPD) or any cell parameter triggers one cached recompute.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")          # headless server backend (no GUI / Tk)
import matplotlib.pyplot as plt
import streamlit as st

import fwm_obe as fwm

st.set_page_config(page_title="85Rb FWM — Gain & Squeezing",
                   page_icon="🔬", layout="wide")

# Resolution presets: (coarse points across the window, velocity step m/s).
RESOLUTION = {
    "Fast  (~3 s)":     dict(coarse_points=121, velocity_step=5.0),
    "Balanced  (~6 s)": dict(coarse_points=181, velocity_step=4.0),
    "Fine  (~20 s)":    dict(coarse_points=301, velocity_step=2.0),
}
WINDOW_GHZ = 0.55          # half-width of the probe window around the (−) Raman line
TPD_LIMIT_MHZ = 500.0      # two-photon-detuning slider range (kept inside the window)


# ----------------------------------------------------------------------
# Cached compute layer
# ----------------------------------------------------------------------
@st.cache_data(show_spinner=False)
def focused_spectrum(D_GHz, T_K, P_pump_mW, P_probe_uW,
                     line_strength, loss_pct, coarse_points, velocity_step):
    """Probe scan in a ±WINDOW_GHZ window around the (−) Raman resonance."""
    center = fwm.branch_center_GHz(D_GHz, -1)
    return fwm.compute_spectrum(
        D_GHz,
        T=T_K,
        P_pump=P_pump_mW * 1e-3,
        P_probe=P_probe_uW * 1e-6,
        line_strength=line_strength,
        loss_frac=loss_pct / 100.0,
        coarse_points=coarse_points,
        fine_points=0,
        scan_min=center - WINDOW_GHZ,
        scan_max=center + WINDOW_GHZ,
        velocity_step=velocity_step,
        velocity_cutoff=3.0,
        branches=fwm.BRANCHES,
    )


@st.cache_data(show_spinner=False)
def full_spectrum(D_GHz, T_K, P_pump_mW, P_probe_uW, line_strength, loss_pct):
    """Original wide −8…12 GHz scan with fine windows at both Raman lines."""
    return fwm.compute_spectrum(
        D_GHz,
        T=T_K,
        P_pump=P_pump_mW * 1e-3,
        P_probe=P_probe_uW * 1e-6,
        line_strength=line_strength,
        loss_frac=loss_pct / 100.0,
        coarse_points=301,
        fine_points=401,
        velocity_step=2.0,
    )


def delta_axis_mhz(spectrum):
    """Probe detuning axis re-expressed as two-photon detuning δ on the (−) branch."""
    return (spectrum["probe_axis_GHz"] - spectrum["raman_center_minus_GHz"]) * 1e3


# ----------------------------------------------------------------------
# Sidebar — controls
# ----------------------------------------------------------------------
st.sidebar.title("Controls")

# Slider state lives in st.session_state so the preset button can overwrite it.
DEFAULTS = dict(opd=0.9, tpd=0.0, temp_c=121.0, pump_mw=600.0,
                probe_uw=10.0, loss_pct=0.0, ls=0.05)
for _k, _v in DEFAULTS.items():
    st.session_state.setdefault(_k, _v)

# 85Rb conditions optimised for squeezing in Sim, Kim & Moon,
# Sci. Rep. 15, 7727 (2025): T = 121 °C, Δ ≈ 0.9 GHz, δ = −8 MHz, pump 600 mW,
# probe-seed 8 µW (Fig. 3 squeezing run), optical loss ≈ 5.5 %. Gain ≈ 15,
# measured IDS −7.8 dB.  (Line-strength is a model calibration knob, not a
# paper value, so the preset leaves it untouched.)
OPTIMAL_85RB = dict(opd=0.9, tpd=-8.0, temp_c=121.0, pump_mw=600.0,
                    probe_uw=8.0, loss_pct=5.5)


def _apply_optimal():
    st.session_state.update(OPTIMAL_85RB)


st.sidebar.button(
    "⚡ Sim et al. 85Rb optimum", on_click=_apply_optimal,
    use_container_width=True,
    help="One click → squeezing-optimised 85Rb conditions from "
         "Sim, Kim & Moon, Sci. Rep. 15, 7727 (2025): "
         "Δ = 0.9 GHz, δ = −8 MHz, T = 121 °C, pump 600 mW, seed 8 µW, loss 5.5 %.",
)

st.sidebar.subheader("Detunings")
opd_ghz = st.sidebar.slider(
    "OPD — one-photon detuning Δ  [GHz]",
    min_value=-1.0, max_value=3.0, step=0.1, key="opd",
    help="ω_pump = ω(F=2→F'=3) + Δ.  Sets where the pump sits and recomputes the spectrum.",
)
tpd_mhz = st.sidebar.slider(
    "TPD — two-photon detuning δ  [MHz]",
    min_value=-TPD_LIMIT_MHZ, max_value=TPD_LIMIT_MHZ, step=1.0, key="tpd",
    help="ω_seed = ω_pump − ν_HF + δ.  Navigates the curve instantly (no recompute).",
)

st.sidebar.subheader("Cell & beams")
T_C = st.sidebar.slider("Temperature  [°C]", 60.0, 150.0, step=1.0, key="temp_c")
P_pump_mW = st.sidebar.slider("Pump power  [mW]", 50.0, 1200.0, step=10.0, key="pump_mw")
P_probe_uW = st.sidebar.slider("Seed / probe power  [µW]", 1.0, 200.0, step=1.0, key="probe_uw")

st.sidebar.subheader("Detection & scaling")
loss_pct = st.sidebar.slider("Loss after cell  [%]", 0.0, 50.0, step=0.5, key="loss_pct",
                             help="Folds into η = QE × (1 − loss).")
line_strength = st.sidebar.slider(
    "Line-strength factor", 0.01, 1.0, step=0.01, key="ls",
    help="Effective |d|² rescaling (Clebsch-Gordan lumping). Calibration knob — "
         "tune to reproduce the paper's gain ≈ 15; not a tabulated paper value.",
)
resolution = st.sidebar.selectbox("Resolution", list(RESOLUTION.keys()), index=1, key="res")

T_K = T_C + 273.15
res = RESOLUTION[resolution]


# ----------------------------------------------------------------------
# Compute (cached)
# ----------------------------------------------------------------------
with st.spinner(f"Solving Bloch equations (Δ = {opd_ghz:.1f} GHz)…"):
    spec = focused_spectrum(opd_ghz, T_K, P_pump_mW, P_probe_uW,
                            line_strength, loss_pct,
                            res["coarse_points"], res["velocity_step"])

op = fwm.operating_point(spec, tpd_mhz, branch=-1)
d_axis = delta_axis_mhz(spec)


# ----------------------------------------------------------------------
# Header + operating-point readout
# ----------------------------------------------------------------------
st.title("85Rb D1 double-Λ four-wave mixing")
st.caption("Seed/probe gain and intensity-difference squeezing vs two-photon detuning. "
           "OPD and cell parameters recompute (cached); TPD navigates instantly.")

c1, c2, c3 = st.columns(3)
c1.metric("Seed / probe gain  G_s", f"{op['G_s']:.2f}",
          help="Power gain of the seeded probe through the cell.")
c2.metric("Squeezing", f"{op['S_dB']:.2f} dB",
          delta="below shot noise" if op["S_dB"] < 0 else "above shot noise",
          delta_color="inverse")
c3.metric("Conjugate gain  G_c", f"{op['G_c']:.2f}",
          help="Generated conjugate power gain (drives the twin-beam squeezing).")


# ----------------------------------------------------------------------
# Plots — G_s and squeezing vs TPD δ, marker at operating point
# ----------------------------------------------------------------------
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
fig, (axG, axS) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)

axG.plot(d_axis, spec["G_s"], color="#1f77b4", lw=1.8)
axG.axvline(tpd_mhz, color="crimson", ls="--", lw=1.2)
axG.axhline(1.0, color="black", lw=0.6)
axG.scatter([tpd_mhz], [op["G_s"]], color="crimson", zorder=5)
axG.set_ylabel("Seed / probe gain  $G_s$")
axG.set_title(f"Δ = {opd_ghz:.1f} GHz,  T = {T_C:.0f} °C,  η = {spec['eta']:.3f}")
if np.nanmax(spec["G_s"]) > 50:
    axG.set_yscale("log")

axS.plot(d_axis, spec["S_dB"], color="#2ca02c", lw=1.8)
axS.axvline(tpd_mhz, color="crimson", ls="--", lw=1.2)
axS.axhline(0.0, color="black", lw=0.6)
axS.scatter([tpd_mhz], [op["S_dB"]], color="crimson", zorder=5)
axS.set_ylabel("Intensity-difference\nsqueezing  [dB]")
axS.set_xlabel("Two-photon detuning δ  [MHz]   (probe on the − Raman branch)")
axS.set_xlim(-TPD_LIMIT_MHZ, TPD_LIMIT_MHZ)

fig.tight_layout()
st.pyplot(fig)


# ----------------------------------------------------------------------
# Reference optimum + derived quantities + optional full scan
# ----------------------------------------------------------------------
with st.expander("Sim et al. (2025) optimum — reference values"):
    st.markdown(
        """
**Source:** G. Sim, H. Kim, H. S. Moon, *Sci. Rep.* **15**, 7727 (2025).
85Rb conditions optimised for intensity-difference squeezing — the
**⚡ Sim et al. 85Rb optimum** button loads the sliders below.

| Parameter | Paper value | Slider |
|---|---|---|
| One-photon detuning Δ | ≈ 0.9 GHz | ✅ |
| Two-photon detuning δ | −8 MHz | ✅ |
| Cell temperature | 121 °C | ✅ |
| Pump power | 600 mW | ✅ |
| Probe-seed power | 8 µW (squeezing run) | ✅ |
| Optical loss after cell | 5.5 % | ✅ (loss) |
| Cell length | 12.5 mm | fixed |
| Pump / seed waist w₀ | 530 / 330 µm (1/e² radius) | fixed geometry |
| Measured result | gain ≈ 15, IDS −7.8 dB | — |

ℹ️ **Line-strength factor** is a model calibration knob (effective |d|²
rescaling), not a tabulated paper value, so the preset leaves it as-is. Tune it
until the on-resonance gain matches the paper's ≈ 15.
        """
    )

with st.expander("Derived quantities"):
    st.markdown(
        f"""
| Quantity | Value |
|---|---|
| N(85Rb) | {spec['N_atoms']:.3e} /m³ |
| σ_v (1-D thermal) | {spec['sigma_v']:.1f} m/s |
| Velocity classes | {spec['n_velocity']} |
| Ω_pump / 2π | {spec['Op_A_2pi_GHz']:.3f} GHz |
| Ω_seed / 2π | {spec['Os_2pi_MHz']:.3f} MHz |
| (−) Raman line (probe axis) | {spec['raman_center_minus_GHz']:.3f} GHz |
| Detection η = QE·(1−loss) | {spec['eta']:.4f} |
| Operating probe detuning | {op['probe_GHz']:.4f} GHz |

Fixed: cell L = {fwm.L_CELL*1e3:.1f} mm · pump w₀ {fwm.W_PUMP*1e6:.0f} µm ·
seed w₀ {fwm.W_PROBE*1e6:.0f} µm · QE {fwm.QE_DETECTOR*100:.2f}% ·
responsivity {fwm.RESPONSIVITY_AW} A/W @ 795 nm · pump⊥probe at PBS.
        """
    )

with st.expander("Full −8…12 GHz probe scan (slow, both Raman branches)"):
    st.caption("The focused view above zooms on the (−) Raman line. "
               "This runs the original wide scan showing both branches.")
    if st.button("Run full scan"):
        with st.spinner("Running wide scan…"):
            full = full_spectrum(opd_ghz, T_K, P_pump_mW, P_probe_uW,
                                 line_strength, loss_pct)
        figF, (aG, aS) = plt.subplots(2, 1, figsize=(8.5, 6.4), sharex=True)
        aG.plot(full["probe_axis_GHz"], full["G_s"], color="#1f77b4", lw=1.4)
        aG.axhline(1.0, color="black", lw=0.6)
        aG.set_ylabel("Seed / probe gain  $G_s$")
        if np.nanmax(full["G_s"]) > 50:
            aG.set_yscale("log")
        aS.plot(full["probe_axis_GHz"], full["S_dB"], color="#2ca02c", lw=1.4)
        aS.axhline(0.0, color="black", lw=0.6)
        aS.set_ylabel("Squeezing [dB]")
        aS.set_xlabel(r"Probe detuning from $F=2\to F'=3$  [GHz]")
        for a in (aG, aS):
            a.axvline(opd_ghz, color="gray", ls=":", lw=0.8)
        figF.tight_layout()
        st.pyplot(figF)
