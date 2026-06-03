"""
Physics checks for saturated absorption (SAS).

Realistic multilevel mode (gabes.species):
  data       : Wigner-6j strengths ∝ validated ⁸⁵Rb D1 CF2; Casimir HF energies
               reproduce the known ground/excited splittings; correct line counts.
  no pump    : reduces to the Doppler-broadened multi-line absorption — the ⁸⁵Rb
               F=3/F=2 manifold ratio matches the validated OD value 49/25.
  pump on    : sharp Doppler-free features appear (Lamb dips + crossovers).
  pumping    : crossovers are enhanced and grow as the transit rate falls — the
               hyperfine-optical-pumping signature a single-ground model can't give.

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

G = constants.GAMMA
GMHZ = G / (2 * np.pi) / 1e6
SAS = schemes.get("sas")
RB85_KEY, RB87_KEY, CS_KEY = species.SPECIES_ORDER[1], species.SPECIES_ORDER[2], species.SPECIES_ORDER[3]


def _params(**over):
    p = SAS.defaults()
    p.update(over)
    return p


def _spectrum(**over):
    raw = SAS.compute(_params(**over))
    return raw["scan"] / (2 * np.pi) / 1e9, raw["alpha_unit"], raw     # GHz, 1/m


# ---------------------------------------------------------------- data layer
def test_line_strengths_match_validated_cf2():
    """T = (2Fg+1)·line_strength reproduces the AutoOD-validated ⁸⁵Rb D1 CF2 (×9)."""
    I, Jg, Je = 2.5, 0.5, 0.5
    for (Fg, Fe), cf2 in hyperfine.CF2.items():
        T = (2 * Fg + 1) * species.line_strength(Fg, Fe, I, Jg, Je)
        assert abs(T - 9.0 * cf2) < 1e-9


def test_hf_energies_match_known_splittings():
    assert abs(species.hf_energy_mhz(species.RB85.A_S, 0, 2.5, 0.5, 3)
               - species.hf_energy_mhz(species.RB85.A_S, 0, 2.5, 0.5, 2) - 3035.732) < 0.01
    assert abs(species.hf_energy_mhz(species.RB87.A_S, 0, 1.5, 0.5, 2)
               - species.hf_energy_mhz(species.RB87.A_S, 0, 1.5, 0.5, 1) - 6834.683) < 0.01
    assert abs(species.hf_energy_mhz(species.CS133.A_S, 0, 3.5, 0.5, 4)
               - species.hf_energy_mhz(species.CS133.A_S, 0, 3.5, 0.5, 3) - 9192.632) < 0.01


def test_manifold_line_counts():
    # D2 each: 2 grounds × 3 allowed Fe (ΔF=0,±1) = 6 transitions; D1: 4.
    assert len(species.build_manifold(species.RB85, "D2").omega) == 6
    assert len(species.build_manifold(species.RB87, "D2").omega) == 6
    assert len(species.build_manifold(species.CS133, "D2").omega) == 6
    assert len(species.build_manifold(species.RB85, "D1").omega) == 4


def test_decay_branching_sums_and_values():
    """Σ_Fg branching = Γ per excited; ⁸⁵Rb D1 Fe=2 → Fg2:0.222, Fg3:0.778."""
    man = species.build_manifold(species.RB85, "D1")
    br = {}
    for (frm, to, rate) in man.atom.decay:
        if frm in man.atom.excited:
            br.setdefault(frm, {})[to] = rate / man.gamma
    for e, d in br.items():
        assert abs(sum(d.values()) - 1.0) < 1e-9
    e2 = man.atom.ground[-1] + 1                      # first excited index (Fe=2)
    assert abs(br[e2][0] - 2.0 / 9.0) < 1e-6          # → Fg=2
    assert abs(br[e2][1] - 7.0 / 9.0) < 1e-6          # → Fg=3


# ------------------------------------------------------------- species model
def test_no_pump_recovers_doppler_manifold():
    """No pump → smooth Doppler lines; ⁸⁵Rb D1 F=3/F=2 area ratio = 49/25."""
    x, a, raw = _spectrum(species=RB85_KEY, line="D1", temp_c=45.0, pump_rabi=0.02)
    f2 = a[x > 0.5].sum()
    f3 = a[x < -0.5].sum()
    assert 1.90 <= f3 / f2 <= 2.02                    # validated OD value 49/25 = 1.96
    rough = np.abs(np.diff(a, 2)).max() / max(a.max(), 1e-9)
    assert rough < 0.05                               # no sub-Doppler features


def test_pump_creates_subdoppler_features():
    _, a_off, _ = _spectrum(species=RB85_KEY, line="D2", temp_c=30.0, pump_rabi=0.02)
    _, a_on, _ = _spectrum(species=RB85_KEY, line="D2", temp_c=30.0, pump_rabi=3.0)

    def rough(a):
        return np.abs(np.diff(a, 2)).max() / max(a.max(), 1e-9)
    assert rough(a_on) > 5 * rough(a_off)             # sharp features only with pump
    assert a_on.sum() < a_off.sum()                   # pump reduces total absorption


def test_hyperfine_pumping_enhances_crossover():
    """Crossover transmission rises as the transit rate falls (pumping signature)."""
    co = 1.719                                        # ⁸⁵Rb D2 F=2 (2′×3′) crossover [GHz]
    T = []
    for gt in (2000.0, 100.0, 20.0):
        x, a, _ = _spectrum(species=RB85_KEY, line="D2", temp_c=30.0,
                            pump_rabi=3.0, transit_khz=gt)
        i = int(np.argmin(np.abs(x - co)))
        T.append(observables.transmission(a, 0.05)[i])
    assert T[0] < T[1] < T[2]                         # smaller γ_t → stronger crossover
    # and the crossover dominates the nearby individual Lamb dips
    x, a, raw = _spectrum(species=RB85_KEY, line="D2", temp_c=30.0, pump_rabi=3.0)
    Ttr = observables.transmission(a, 0.05)
    dips = [Ttr[int(np.argmin(np.abs(x - gx)))] for gx, _ in raw["markers"] if gx > 1.0]
    ico = int(np.argmin(np.abs(x - co)))
    assert Ttr[ico] > max(dips)                       # crossover taller than Lamb dips


def test_natural_rb_overlays_both_isotopes():
    raw = SAS.compute(_params(species="Rb (natural)", line="D2", temp_c=40.0))
    assert len(raw["markers"]) == 12                  # 6 (⁸⁵Rb) + 6 (⁸⁷Rb)
    labels = " ".join(lbl for _, lbl in raw["markers"])
    assert species.RB85.label in labels and species.RB87.label in labels


def test_observables_render_species():
    raw = SAS.compute(_params(species=CS_KEY, line="D2", temp_c=35.0))
    view = SAS.observables(raw, _params(species=CS_KEY, line="D2", temp_c=35.0))
    assert view["figure"] is not None
    assert any(m["label"] == "Peak OD" for m in view["metrics"])


# -------------------------------------------------------------- generic mode
def test_generic_no_pump_voigt():
    _, a, raw = _spectrum(species=GENERIC, transitions="single line", pump_rabi=0.01)
    x = raw["scan"] / (2 * np.pi) / 1e6
    above = x[a >= 0.5 * a.max()]
    env = above.max() - above.min()
    dopp = raw["dopp_fwhm"] / (2 * np.pi) / 1e6
    assert 0.95 <= env / dopp <= 1.05
    assert a[int(np.argmin(np.abs(x)))] >= a[int(np.argmin(np.abs(x - 15 * GMHZ)))]


def test_generic_lamb_dip_and_crossover():
    raw = SAS.compute(_params(species=GENERIC, transitions="single line", pump_rabi=2.0))
    x = raw["scan"] / (2 * np.pi) / 1e6
    a = raw["alpha_unit"]
    ic, ish = int(np.argmin(np.abs(x))), int(np.argmin(np.abs(x - 15 * GMHZ)))
    assert a[ic] < 0.7 * a[ish]                       # sub-Doppler Lamb dip

    raw2 = SAS.compute(_params(species=GENERIC, transitions="two lines (crossover)",
                               splitting=60.0, pump_rabi=2.0))
    x2, a2 = raw2["scan"] / (2 * np.pi) / 1e6, raw2["alpha_unit"]
    off = 30 * GMHZ
    assert a2[int(np.argmin(np.abs(x2)))] < a2[int(np.argmin(np.abs(x2 - 8 * GMHZ)))]


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print("ok:", fn.__name__)
    print(f"\nSAS OK ({len(fns)} tests): data layer, Doppler-free features, "
          "hyperfine pumping, natural-Rb overlay, generic fallback.")
