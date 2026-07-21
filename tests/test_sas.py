"""
Physics checks for the merged absorption scheme (OD / SAS).

Data layer (gabes.species):
  6j/3j strengths reproduce the validated ⁸⁵Rb D1 CF2; Casimir HF energies match
  the known ground/excited splittings; correct hyperfine line counts.

Pump off (P = 0) → OD:
  reduces to linear Doppler-broadened absorption; ⁸⁵Rb D1 reproduces the AutoOD
  scale (integrated + peak) and the 49/25 F=3/F=2 manifold ratio.

Pump on → SAS:
  sharp Doppler-free features appear (Lamb dips + crossovers); crossovers are
  enhanced and grow as the transit rate falls — the hyperfine-pumping signature.

Generic Γ-unit hole-burning fallback: Voigt background, Lamb dip, crossover.

    python tests/test_sas.py     # or: pytest tests/test_sas.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import schemes, constants, observables, species, hyperfine  # noqa: E402
from gabes.schemes.sas import GENERIC  # noqa: E402
from gabes.schemes.absorption import ODScheme  # noqa: E402
from gabes.lineshape import narrowest_subdoppler, window_fwhm  # noqa: E402

G = constants.GAMMA
GMHZ = G / (2 * np.pi) / 1e6
SAS = schemes.get("sas")
_tz = getattr(np, "trapezoid", getattr(np, "trapz", None))
RB85_KEY, RB87_KEY, CS_KEY = (species.SPECIES_ORDER[1], species.SPECIES_ORDER[2],
                              species.SPECIES_ORDER[3])


def _params(**over):
    p = SAS.defaults()
    p.update(over)
    return p


def _spectrum(**over):
    raw = SAS.compute(_params(**over))
    return raw["scan"] / (2 * np.pi) / 1e9, raw["alpha_unit"], raw     # GHz, 1/m


def _metric_number(view, label):
    value = next(m["value"] for m in view["metrics"] if m["label"] == label)
    return float(value.split()[0])


# ---------------------------------------------------------------- data layer
def test_cf2_matches_validated_85rb_d1():
    for (Fg, Fe), ref in hyperfine.CF2.items():
        assert abs(species.cf2(Fg, Fe, 2.5, 0.5, 0.5) - ref) < 1e-9


def test_reduced_dipole_and_density_match_od():
    lam = constants.C_LIGHT / species.RB85.nu_D1
    d2 = species.reduced_dipole_sq(2 * np.pi * 5.75e6, lam, 0.5, 0.5)
    assert abs(d2 / hyperfine.DIPOLE_SQ - 1.0) < 1e-3
    assert abs(species.number_density(species.RB85, 363.15)
               / hyperfine.number_density(363.15) - 1.0) < 1e-6


def test_hf_energies_match_known_splittings():
    assert abs(species.hf_energy_mhz(species.RB85.A_S, 0, 2.5, 0.5, 3)
               - species.hf_energy_mhz(species.RB85.A_S, 0, 2.5, 0.5, 2) - 3035.732) < 0.01
    assert abs(species.hf_energy_mhz(species.RB87.A_S, 0, 1.5, 0.5, 2)
               - species.hf_energy_mhz(species.RB87.A_S, 0, 1.5, 0.5, 1) - 6834.683) < 0.01
    assert abs(species.hf_energy_mhz(species.CS133.A_S, 0, 3.5, 0.5, 4)
               - species.hf_energy_mhz(species.CS133.A_S, 0, 3.5, 0.5, 3) - 9192.632) < 0.01


def test_manifold_line_counts():
    assert len(species.build_manifold(species.RB85, "D2").omega) == 6
    assert len(species.build_manifold(species.RB87, "D2").omega) == 6
    assert len(species.build_manifold(species.CS133, "D2").omega) == 6
    assert len(species.build_manifold(species.RB85, "D1").omega) == 4


def test_decay_branching_sums_and_values():
    man = species.build_manifold(species.RB85, "D1")
    br = {}
    for (frm, to, rate) in man.atom.decay:
        if frm in man.atom.excited:
            br.setdefault(frm, {})[to] = rate / man.gamma
    for d in br.values():
        assert abs(sum(d.values()) - 1.0) < 1e-9
    e2 = man.atom.ground[-1] + 1                      # first excited (Fe=2)
    assert abs(br[e2][0] - 2.0 / 9.0) < 1e-6          # → Fg=2
    assert abs(br[e2][1] - 7.0 / 9.0) < 1e-6          # → Fg=3


# ---------------------------------------------------- pump off (OD) fidelity
def test_pump_off_reproduces_autood_85rb_d1():
    """Pump = 0, ⁸⁵Rb D1 reproduces the AutoOD-validated OD scheme to <1 %."""
    raw = SAS.compute(_params(species=RB85_KEY, line="D1", pump_power_mw=0.0,
                              temp_c=90.0, cell_mm=12.5))
    od = ODScheme()
    ro = od.compute({**od.defaults(), "model": "85Rb D1 hyperfine",
                     "temp_c": 90.0, "cell_mm": 12.5, "doppler": "on"})
    int_ratio = _tz(raw["alpha_unit"], raw["scan"]) / _tz(ro["alpha"], ro["scan"])
    peak_ratio = raw["alpha_unit"].max() / ro["alpha"].max()
    assert abs(int_ratio - 1.0) < 0.01
    assert abs(peak_ratio - 1.0) < 0.01


def test_pump_off_is_smooth_and_49_25():
    p = _params(species=RB85_KEY, line="D1", pump_power_mw=0.0, temp_c=45.0)
    raw = SAS.compute(p)
    x, a = raw["scan"] / (2 * np.pi) / 1e9, raw["alpha_unit"]
    f2 = a[x > 0.5].sum()
    f3 = a[x < -0.5].sum()
    assert 1.90 <= f3 / f2 <= 2.02                    # validated 49/25 ≈ 1.96
    assert np.abs(np.diff(a, 2)).max() / max(a.max(), 1e-9) < 0.02   # no features

    view = SAS.observables(raw, p, include_figures=False)
    assert view["hero_count"] == 1
    assert [metric["label"] for metric in view["metrics"]] == ["Peak OD"]
    assert view["metrics"][0]["tier"] == "hero"


# ------------------------------------------------------ pump on (SAS) physics
def test_pump_creates_subdoppler_features():
    _, a_off, _ = _spectrum(species=RB85_KEY, line="D2", pump_power_mw=0.0, temp_c=30.0)
    _, a_on, _ = _spectrum(species=RB85_KEY, line="D2", pump_power_mw=1.5, temp_c=30.0)

    def rough(a):
        return np.abs(np.diff(a, 2)).max() / max(a.max(), 1e-9)
    assert rough(a_on) > 5 * rough(a_off)             # sharp features only with pump
    assert a_on.sum() < a_off.sum()                   # pump reduces total absorption


def test_hyperfine_pumping_enhances_crossover():
    """Crossover transmission rises as the transit rate falls (pumping signature)."""
    co = 1.719                                        # ⁸⁵Rb D2 F=2 (2′×3′) crossover [GHz]
    Tco = []
    for gt in (2000.0, 100.0, 20.0):
        x, a, _ = _spectrum(species=RB85_KEY, line="D2", pump_power_mw=1.5,
                            temp_c=30.0, cell_mm=50.0, transit_khz=gt)
        i = int(np.argmin(np.abs(x - co)))
        Tco.append(observables.transmission(a, 0.05)[i])
    assert Tco[0] < Tco[1] < Tco[2]                   # smaller γ_t → stronger crossover
    x, a, raw = _spectrum(species=RB85_KEY, line="D2", pump_power_mw=1.5,
                          temp_c=30.0, cell_mm=50.0)
    Ttr = observables.transmission(a, 0.05)
    dips = [Ttr[int(np.argmin(np.abs(x - gx)))] for gx, _ in raw["markers"] if gx > 1.0]
    assert Ttr[int(np.argmin(np.abs(x - co)))] > max(dips)   # crossover > Lamb dips


def test_natural_rb_overlays_both_isotopes():
    raw = SAS.compute(_params(species="Rb (natural)", line="D2", temp_c=40.0))
    assert len(raw["markers"]) == 12                  # 6 (⁸⁵Rb) + 6 (⁸⁷Rb)
    labels = " ".join(lbl for _, lbl in raw["markers"])
    assert species.RB85.label in labels and species.RB87.label in labels


def test_recommended_defaults_od_and_sas():
    sets = SAS.recommended_defaults(dict(species=CS_KEY, line="D2"))
    assert set(sets) == {"OD default", "SAS default"}
    assert sets["OD default"]["pump_power_mw"] == 0.0          # OD = pump off
    assert sets["SAS default"]["pump_power_mw"] > 0.0          # SAS = pump on
    assert sets["OD default"]["temp_c"] == sets["SAS default"]["temp_c"]   # same cell/T
    d1 = SAS.recommended_defaults(dict(species=CS_KEY, line="D1"))
    assert d1["SAS default"] != sets["SAS default"]            # genuinely per-line


def test_observables_render_species():
    p = _params(species=CS_KEY, line="D2", temp_c=35.0)
    view = SAS.observables(SAS.compute(p), p)
    assert view["figure"] is not None
    figure_views = view["figure_views"]
    assert [item["label"] for item in figure_views] == [
        "Transmission", "Optical density"
    ]
    assert figure_views[0]["figure"] is view["figure"]
    assert all(len(item["figure"].axes) == 1 for item in figure_views)
    assert figure_views[0]["figure"].axes[0].get_ylabel() == "Transmission"
    assert figure_views[1]["figure"].axes[0].get_ylabel() == "Optical density"
    assert view["comparison"]["axis_index"] == 0
    assert view["comparison"]["x_unit"] == "GHz"
    assert view["comparison"]["raw_x_unit"] == "Arb. unit"
    assert view["hero_count"] == 1
    assert [m["label"] for m in view["metrics"][:3]] == [
        "Lock Slope", "Lock Detuning", "Peak OD"
    ]
    assert [m["label"] for m in view["metrics"] if m.get("tier") == "hero"] == [
        "Lock Slope"
    ]
    assert "Doppler FWHM" not in {m["label"] for m in view["metrics"]}
    assert "Buffer Gas Broadening" not in {m["label"] for m in view["metrics"]}


def test_default_lock_point_tracks_detected_subdoppler_feature():
    p = _params()
    raw = SAS.compute(p)
    x = raw["scan"] / (2 * np.pi) / 1e9
    alpha = raw["alpha_unit"] * p["line_strength"]
    T_trans = observables.transmission(alpha, p["cell_mm"] * 1e-3)
    feature_fwhm, feature_at = narrowest_subdoppler(x, T_trans)

    view = SAS.observables(raw, p, include_figures=False)
    lock_ghz = _metric_number(view, "Lock Detuning") / 1000.0

    assert np.isfinite(feature_fwhm)
    assert abs(lock_ghz - feature_at) <= feature_fwhm + 0.00011


def test_pump_on_unresolved_feature_keeps_lock_slope_hero():
    p = _params(pump_power_mw=0.01, temp_c=200.0, cell_mm=200.0)
    view = SAS.headless_observables(SAS.compute(p), p)
    heroes = [m for m in view["metrics"] if m.get("tier") == "hero"]
    assert [m["label"] for m in heroes] == ["Lock Slope"]
    assert any(m["label"] == "SAS status" for m in view["metrics"])


def test_mode_and_buffer_pressure_control_metric_hierarchy():
    for atom in (RB85_KEY, GENERIC):
        for pump_power in (0.0, 0.5):
            for pressure in (0.0, 1.0):
                p = _params(
                    species=atom,
                    pump_power_mw=pump_power,
                    ne_pressure_torr=pressure,
                    scan_points=401,
                )
                view = SAS.headless_observables(SAS.compute(p), p)
                labels = [metric["label"] for metric in view["metrics"]]
                heroes = [
                    metric["label"] for metric in view["metrics"]
                    if metric.get("tier") == "hero"
                ]

                assert view["hero_count"] == 1
                assert "Doppler FWHM" not in labels
                if pump_power > 0.0:
                    assert heroes == ["Lock Slope"]
                    assert labels[:3] == [
                        "Lock Slope", "Lock Detuning", "Peak OD"
                    ]
                else:
                    assert heroes == ["Peak OD"]
                    assert labels[0] == "Peak OD"
                    assert not any(label.startswith("Lock ") for label in labels)

                if pressure > 0.0:
                    assert labels[-1] == "Buffer Gas Broadening"
                    assert _metric_number(view, "Buffer Gas Broadening") > 0.0
                else:
                    assert "Buffer Gas Broadening" not in labels


# -------------------------------------------------------------- generic mode
def test_generic_lamb_dip_and_crossover():
    raw = SAS.compute(_params(species=GENERIC, transitions="single line", pump_power_mw=1.0))
    x = raw["scan"] / (2 * np.pi) / 1e6
    a = raw["alpha_unit"]
    ic, ish = int(np.argmin(np.abs(x))), int(np.argmin(np.abs(x - 15 * GMHZ)))
    assert a[ic] < 0.7 * a[ish]                       # sub-Doppler Lamb dip

    p2 = _params(species=GENERIC, transitions="two lines (crossover)",
                 splitting=60.0, pump_power_mw=1.0)
    raw2 = SAS.compute(p2)
    x2, a2 = raw2["scan"] / (2 * np.pi) / 1e6, raw2["alpha_unit"]
    assert a2[int(np.argmin(np.abs(x2)))] < a2[int(np.argmin(np.abs(x2 - 8 * GMHZ)))]

    # The broad Doppler-envelope flank is steeper for this case.  The lock
    # proxy must nevertheless stay on the detected Lamb-dip feature.
    T2 = observables.transmission(a2, p2["cell_mm"] * 1e-3)
    feature_i = int(np.argmin(np.abs(
        x2 - raw2["offsets"][0] / (2 * np.pi) / 1e6)))
    feature_fwhm = window_fwhm(x2, T2, feature_i)
    rendered = SAS.observables(raw2, p2)
    assert [item["label"] for item in rendered["figure_views"]] == [
        "Transmission", "Optical density"
    ]
    assert rendered["figure_views"][0]["figure"] is rendered["figure"]
    view2 = SAS.observables(raw2, p2, include_figures=False)
    assert view2["figure_views"] == []
    assert view2["comparison"]["x_unit"] == "MHz"
    assert view2["hero_count"] == 1
    assert [m["label"] for m in view2["metrics"][:3]] == [
        "Lock Slope", "Lock Detuning", "Peak OD"
    ]
    lock_mhz = _metric_number(view2, "Lock Detuning")
    assert abs(lock_mhz - x2[feature_i]) <= feature_fwhm + 0.11


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print("ok:", fn.__name__)
    print(f"\nAbsorption (OD/SAS) OK ({len(fns)} tests): data, pump-off AutoOD "
          "fidelity, sub-Doppler features, hyperfine pumping, generic.")
