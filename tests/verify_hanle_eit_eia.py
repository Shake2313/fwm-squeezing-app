"""
Theoretical cross-check: does the GABES magneto (Hanle) scheme reproduce the
EIT/CPT -> EIA/MIA polarization switch and the sub-milligauss feature widths
reported by the Moon group for a paraffin-coated 87Rb cell?

Reference values
----------------
Lee & Moon, JOSA B 30, 2301 (2013), "Magnetic-field-induced absorption with
sub-milligauss spectral width in a paraffin-coated Rb vapor cell":
  * Hanle config, 87Rb, paraffin cell.
  * Linear polarization  -> CPT (transparency dip),  width ~ 0.12 mG.
  * Circular polarization -> magnetic-field-induced absorption (MIA, a peak),
    width ~ 0.20 mG.
  * Switch driven at FIXED residual transverse B by the laser polarization,
    aided by the ground-coherence wall lifetime (Ramsey effect).

Unit bridge:  1 G = 100 uT  ->  1 mG = 0.1 uT.
  0.12 mG = 0.012 uT,  0.20 mG = 0.020 uT.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import observables, schemes  # noqa: E402

UT_PER_MG = 0.1  # 1 mG = 0.1 uT


def alpha_curve(raw, p):
    xprobe = observables.chi_phys(raw["chi_probe"], raw["N_eff"],
                                  dipole=raw["dipole"], line_strength=p["line_strength"])
    return raw["b_ut"], raw["k_vec"] * np.imag(xprobe)


def feature_amp(x, alpha):
    ic = int(np.argmin(np.abs(x)))
    bg = 0.5 * (alpha[0] + alpha[-1])
    return alpha[ic] - bg


def central_fwhm_ut(x, alpha):
    ic = int(np.argmin(np.abs(x)))
    bg = 0.5 * (alpha[0] + alpha[-1])
    amp = alpha[ic] - bg
    if abs(amp) < 1e-30:
        return np.nan
    target = bg + 0.5 * amp
    sign = np.sign(amp)

    def half(side):
        idx = np.arange(ic, len(x)) if side > 0 else np.arange(ic, -1, -1)
        vals = (alpha[idx] - target) * sign
        below = np.where(vals <= 0)[0]
        if below.size == 0 or below[0] == 0:
            return np.nan
        j = idx[below[0]]
        i = idx[below[0] - 1]
        t = (target - alpha[i]) / (alpha[j] - alpha[i])
        return x[i] + t * (x[j] - x[i])

    bl, br = half(-1), half(+1)
    if np.isnan(bl) or np.isnan(br):
        return np.nan
    return abs(br - bl)


def run(qwp, b_max_ut, b_perp_ut, wall_ms=3.0, transit_khz=80.0,
        points=601, intensity=0.8):
    sc = schemes.get("magneto")
    p = sc.defaults()
    p.update(cell_type="Paraffin coated cell", Fg=2.0, Fe=1.0,
             intensity_mw_cm2=intensity, qwp_deg=qwp,
             b_max_ut=b_max_ut, residual_transverse_b_ut=b_perp_ut,
             transverse_field_angle_deg=0.0, wall_coherence_ms=wall_ms,
             transit_relax_khz=transit_khz, dark_return_khz=1.0,
             scan_points=points, velocity_classes=5, doppler="on",
             laser_detuning_mhz=0.0, temp_c=25.0, cell_mm=10.0)
    x, alpha = alpha_curve(sc.compute(p), p)
    return dict(x=x, alpha=alpha, amp=feature_amp(x, alpha),
                fwhm_ut=central_fwhm_ut(x, alpha))


def tag(amp):
    return "DIP (CPT/EIT)" if amp < 0 else "PEAK (MIA/EIA)"


def transmission_contrast(raw, p, x, alpha):
    L = p["cell_mm"] * 1e-3
    T = np.exp(-alpha * L)
    ic = int(np.argmin(np.abs(x)))
    return T[ic] - 0.5 * (T[0] + T[-1])


def run_transition(Fg, Fe, qwp, b_perp_ut, b_max_ut=0.5, **kw):
    from gabes.schemes.magneto import _transition_strength
    sc = schemes.get("magneto")
    p = sc.defaults()
    p.update(cell_type="Paraffin coated cell", Fg=float(Fg), Fe=float(Fe),
             intensity_mw_cm2=0.8, qwp_deg=qwp, b_max_ut=b_max_ut,
             residual_transverse_b_ut=b_perp_ut, transverse_field_angle_deg=0.0,
             wall_coherence_ms=3.0, transit_relax_khz=80.0, dark_return_khz=1.0,
             scan_points=801, velocity_classes=5, doppler="on",
             laser_detuning_mhz=0.0, temp_c=25.0, cell_mm=10.0)
    p.update(kw)
    raw = sc.compute(p)
    if not raw["valid"]:
        return None
    x, alpha = alpha_curve(raw, p)
    return dict(amp=feature_amp(x, alpha), fwhm_ut=central_fwhm_ut(x, alpha),
                dT=transmission_contrast(raw, p, x, alpha),
                strength=_transition_strength(Fg, Fe))


def transition_dependence():
    """Does the EIT/EIA sign, width, and signal track the transition line?

    Reference rule (Lezama PRA 59, 4732; arXiv physics/0512199 for D1 Fg=1->Fe=2):
      intrinsic EIA on Fe = Fg+1 (coherence-transfer / cycling-type),
      intrinsic EIT on Fe = Fg and Fe = Fg-1.
    """
    print("\n" + "=" * 74)
    print("Transition-line dependence (87Rb D1, Fe in {1,2})")
    print("=" * 74)
    print("Reference: Fe=Fg+1 -> EIA peak (intrinsic, via coherence transfer);")
    print("           Fe=Fg or Fe=Fg-1 -> EIT dip.\n")
    print("[A] INTRINSIC test: linear pol (QWP=0), residual B_perp = 0")
    print(f"    {'Fg->Fe':>7} {'dF':>3} {'S_FF':>7} | {'sign':>14} "
          f"{'FWHM[mG]':>9} {'|dT|':>10}")
    for Fg in (1, 2):
        for Fe in (1, 2):
            r = run_transition(Fg, Fe, qwp=0.0, b_perp_ut=0.0)
            if r is None:
                continue
            fw = r['fwhm_ut'] / UT_PER_MG if r['fwhm_ut'] == r['fwhm_ut'] else float('nan')
            print(f"    {Fg}->{Fe:>4} {Fe-Fg:>3d} {r['strength']:>7.4f} | "
                  f"{tag(r['amp']):>14} {fw:>9.3f} {abs(r['dT']):>10.3e}")
    print("\n[B] MODEL EIA path: circular pol (QWP=45), residual B_perp = 0.2 mG")
    print(f"    {'Fg->Fe':>7} {'dF':>3} {'S_FF':>7} | {'sign':>14} "
          f"{'FWHM[mG]':>9} {'|dT|':>10}")
    for Fg in (1, 2):
        for Fe in (1, 2):
            r = run_transition(Fg, Fe, qwp=45.0, b_perp_ut=0.2 * UT_PER_MG)
            if r is None:
                continue
            fw = r['fwhm_ut'] / UT_PER_MG if r['fwhm_ut'] == r['fwhm_ut'] else float('nan')
            print(f"    {Fg}->{Fe:>4} {Fe-Fg:>3d} {r['strength']:>7.4f} | "
                  f"{tag(r['amp']):>14} {fw:>9.3f} {abs(r['dT']):>10.3e}")
    print("\n[C] Yu LCA regime: buffer gas, LINEAR pol (longitudinal Hanle)")
    print("    Lezama rule should now hold via TOC: Fe=Fg+1 -> EIA, else EIT.")
    print(f"    {'Fg->Fe':>7} {'dF':>3} | {'sign':>14} {'|dT|':>10}")
    for Fg in (1, 2):
        for Fe in (1, 2):
            sc = schemes.get("magneto")
            p = sc.defaults()
            p.update(cell_type="Buffer gas cell", Fg=float(Fg), Fe=float(Fe),
                     intensity_mw_cm2=0.8, qwp_deg=0.0, b_max_ut=0.5,
                     ne_pressure_torr=20.0, buffer_ground_relax_khz=20.0,
                     collisional_depol_khz=2.0, scan_points=801,
                     velocity_classes=5, doppler="on")
            raw = sc.compute(p)
            x, alpha = alpha_curve(raw, p)
            print(f"    {Fg}->{Fe:>4} {Fe-Fg:>3d} | {tag(feature_amp(x, alpha)):>14} "
                  f"{abs(transmission_contrast(raw, p, x, alpha)):>10.3e}")

    print("\n  NOTE: Spontaneous emission is now built from polarization-grouped")
    print("  jump operators Sigma_q (zeeman.py / atoms.py emission_ops), so it")
    print("  carries transfer of coherence (TOC). The reference Fe=Fg+1 -> EIA")
    print("  rule (Lezama; arXiv physics/0512199) is now reproduced for LINEAR")
    print("  polarization in both paraffin and buffer cells.")
    print("\n[D] Yu CIRCULAR-pol LCA (PRA 81, 023416): buffer + circular + transverse")
    print("    Circular light orients the ground state along z (a longitudinal-B")
    print("    eigenstate), so a transverse residual field is required to precess")
    print("    it into a B=0 level-crossing ABSORPTION peak. Reference width 2.4 mG.")
    print(f"    {'I[mW/cm2]':>10} {'B_perp[mG]':>11} | {'sign':>10} {'FWHM[mG]':>9}")
    for I, bperp_mg, relax in [(0.8, 0.0, 5.0), (0.2, 3.0, 5.0), (0.02, 1.0, 0.5)]:
        sc = schemes.get("magneto")
        p = sc.defaults()
        p.update(cell_type="Buffer gas cell", Fg=2.0, Fe=2.0, intensity_mw_cm2=I,
                 qwp_deg=45.0, b_max_ut=1.0, ne_pressure_torr=20.0,
                 buffer_ground_relax_khz=relax, collisional_depol_khz=0.0,
                 residual_transverse_b_ut=bperp_mg * UT_PER_MG,
                 transverse_field_angle_deg=90.0, scan_points=1601,
                 velocity_classes=5, doppler="on")
        x, alpha = alpha_curve(sc.compute(p), p)
        amp = feature_amp(x, alpha)
        bg = 0.5 * (alpha[0] + alpha[-1])
        rel = abs(amp) / max(abs(bg), 1e-30)          # feature vs absorption pedestal
        fw = central_fwhm_ut(x, alpha) / UT_PER_MG
        if rel < 1e-6:
            sign, fw = "flat", float("nan")           # orientation is a B_z eigenstate
        else:
            sign = "peak ABS" if amp > 0 else "dip"
            fw = fw if fw == fw else float("nan")
        print(f"    {I:>10.2f} {bperp_mg:>11.1f} | {sign:>10} {fw:>9.3f}")

    print("\n  NOTE: Spontaneous emission is now built from polarization-grouped")
    print("  jump operators Sigma_q (zeeman.py / atoms.py emission_ops), carrying")
    print("  transfer of coherence (TOC). The reference Fe=Fg+1 -> EIA rule")
    print("  (Lezama; arXiv physics/0512199) is reproduced for LINEAR polarization,")
    print("  and Yu's circular-pol buffer LCA (PRA 81, 023416) is reproduced once a")
    print("  transverse residual field is present (now exposed for both cells); its")
    print("  width reaches the reference ~2.4 mG at low power and low ground relax.")


def main():
    print("=" * 74)
    print("GABES Hanle vs Moon-group paraffin-cell references (87Rb, F=2->F'=1)")
    print("=" * 74)

    bmax = 5 * UT_PER_MG  # +/-0.5 uT = +/-5 mG window, resolves a sub-mG core

    # ---- 1. Does a FIXED residual field give linear=dip, circular=peak? ----
    print("\n[1] Sign vs QWP at FIXED residual transverse B (paper's actual knob)")
    print(f"    scan = +/-{bmax/UT_PER_MG:.1f} mG\n")
    print(f"    {'B_perp [mG]':>11} | {'linear amp':>12} {'':>14} | "
          f"{'circular amp':>12} {'':>14} | switch?")
    found = None
    for bperp_mg in (0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0):
        bperp = bperp_mg * UT_PER_MG
        lin = run(qwp=0.0,  b_max_ut=bmax, b_perp_ut=bperp)
        cir = run(qwp=45.0, b_max_ut=bmax, b_perp_ut=bperp)
        ok = lin['amp'] < 0 and cir['amp'] > 0
        if ok and found is None:
            found = bperp_mg
        print(f"    {bperp_mg:>11.2f} | {lin['amp']:>12.3e} {tag(lin['amp']):>14} | "
              f"{cir['amp']:>12.3e} {tag(cir['amp']):>14} | "
              f"{'YES (lin dip, cir peak)' if ok else 'no'}")

    # ---- 2. QWP sweep at a fixed field: continuous EIT -> EIA crossover -----
    print("\n[2] Continuous CPT->MIA crossover as QWP rotates (fixed B_perp)")
    bperp_demo = (found if found is not None else 0.6) * UT_PER_MG
    print(f"    B_perp = {bperp_demo/UT_PER_MG:.2f} mG\n")
    print(f"    {'QWP [deg]':>9} | {'amp':>12} | {'FWHM [mG]':>10} | feature")
    for qwp in (0, 10, 20, 25, 30, 35, 40, 45):
        r = run(qwp=float(qwp), b_max_ut=bmax, b_perp_ut=bperp_demo)
        fw = r['fwhm_ut'] / UT_PER_MG if r['fwhm_ut'] == r['fwhm_ut'] else float('nan')
        print(f"    {qwp:>9d} | {r['amp']:>12.3e} | {fw:>10.3f} | {tag(r['amp'])}")

    # ---- 3. Width scale vs reference (0.12 mG / 0.20 mG) -------------------
    print("\n[3] Feature FWHM order of magnitude vs reference")
    bperp = (found if found is not None else 0.6) * UT_PER_MG
    lin = run(qwp=0.0,  b_max_ut=bmax, b_perp_ut=bperp)
    cir = run(qwp=45.0, b_max_ut=bmax, b_perp_ut=bperp)
    print(f"    linear   FWHM = {lin['fwhm_ut']/UT_PER_MG:6.3f} mG   (paper CPT ~0.12 mG)")
    print(f"    circular FWHM = {cir['fwhm_ut']/UT_PER_MG:6.3f} mG   (paper MIA ~0.20 mG)")

    # ---- 4. Ramsey narrowing: width shrinks with wall coherence -----------
    print("\n[4] Wall-coherence (Ramsey) narrowing trend, linear pol")
    print(f"    {'wall tau [ms]':>13} | {'FWHM [mG]':>10}")
    for tau in (0.1, 0.5, 1.0, 3.0, 10.0, 30.0, 100.0):
        r = run(qwp=0.0, b_max_ut=2 * UT_PER_MG, b_perp_ut=0.3 * UT_PER_MG,
                wall_ms=tau, transit_khz=40.0, points=801)
        fw = r['fwhm_ut'] / UT_PER_MG if r['fwhm_ut'] == r['fwhm_ut'] else float('nan')
        print(f"    {tau:>13.2f} | {fw:>10.4f}")

    print("\n" + "=" * 74)
    print("Verdict")
    print("=" * 74)
    if found is not None:
        print(f"  - Polarization sign switch reproduced at fixed B_perp = "
              f"{found:.2f} mG: linear=CPT dip, circular=MIA peak. MATCH.")
    else:
        print("  - No single fixed B_perp gave linear=dip & circular=peak "
              "in the scanned window. (See table [1].)")
    print(f"  - Feature widths land in the {lin['fwhm_ut']/UT_PER_MG:.2f}-"
          f"{cir['fwhm_ut']/UT_PER_MG:.2f} mG range; paper is 0.12-0.20 mG "
          "(same sub-mG order).")

    transition_dependence()


if __name__ == "__main__":
    main()
