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

from gabes import core, schemes, constants, observables, species, hyperfine  # noqa: E402
from gabes.schemes.sas import (  # noqa: E402
    GENERIC,
    _basis_reset_pump_pops,
    _coated_ground_transfer,
    _lock_readout_metrics,
    _pump_detuning_axis,
    _pump_pops,
    _stationary_ground_populations,
    _velocity_correlated_average,
)
from gabes.schemes.absorption import ODScheme  # noqa: E402
from gabes.lineshape import narrowest_subdoppler, subdoppler_feature  # noqa: E402

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


def test_paraffin_checkbox_schema_is_advanced_and_recomputes():
    spec = {item.name: item for item in SAS.param_schema()}["paraffin_coated"]
    assert spec.default is False
    assert spec.advanced
    assert spec.control == "checkbox"
    assert spec.recompute
    assert GENERIC not in spec.visible_if["species"]


def test_basis_resets_span_the_legacy_thermal_transit_state():
    """Basis-conditioned reloads reproduce any incoherent transit reservoir."""
    gt = 2 * np.pi * 100e3
    bare = species.build_manifold(species.RB85, "D1", transit_rate=0.0)
    thermal = species.build_manifold(species.RB85, "D1", transit_rate=gt)
    pump = species.pump_hamiltonian(bare, 0.35 * bare.gamma)
    deff = bare.omega[0] + bare.gamma * np.array([-1.0, 0.0, 1.0])
    basis = _basis_reset_pump_pops(bare, pump, deff, gt)
    mixed = np.einsum("g,gdl->dl", bare.p_ground, basis)
    expected = _pump_pops(
        core.build_liouvillian(pump, thermal.atom), deff,
        thermal.atom.S_v, thermal.n_levels)
    np.testing.assert_allclose(mixed, expected, rtol=2e-11, atol=2e-12)

    transfer = _coated_ground_transfer(
        bare, basis, deff, deff[None, :], np.full(deff.size, 1.0 / deff.size))
    assert np.isfinite(transfer).all()
    assert np.all(transfer >= 0.0)
    np.testing.assert_allclose(transfer.sum(axis=1), 1.0,
                               rtol=0.0, atol=2e-15)
    cycle_rate = 2 * np.pi * 100e3
    stationary = _stationary_ground_populations(
        transfer, bare.p_ground, cycle_rate=cycle_rate)
    mapped = np.einsum("sgh,sh->sg", transfer, stationary)
    wall_rate = 1.0 / 25.1e-3
    residual = (cycle_rate * (mapped - stationary)
                + wall_rate * (bare.p_ground - stationary))
    assert np.max(np.abs(residual)) < 2e-10


def test_coated_pump_axis_far_wings_converge_without_edge_clamping():
    man = species.build_manifold(species.RB85, "D1", transit_rate=0.0)
    gamma = man.gamma
    gt = 2 * np.pi * 100e3
    pump = species.pump_hamiltonian(man, 0.5 * gamma)
    scan_detuning = np.array([
        man.omega.min() - 80 * gamma,
        man.omega.mean(),
        man.omega.max() + 80 * gamma,
    ])
    velocity_shift = np.linspace(-60 * gamma, 60 * gamma, 9)
    reached = scan_detuning[:, None] + velocity_shift[None, :]
    z = np.linspace(-2.0, 2.0, velocity_shift.size)
    weights = np.exp(-0.5 * z**2)
    weights /= weights.sum()

    def reservoir(axis):
        basis = _basis_reset_pump_pops(man, pump, axis, gt)
        transfer = _coated_ground_transfer(
            man, basis, axis, reached, weights)
        return _stationary_ground_populations(
            transfer, man.p_ground, cycle_rate=gt)

    legacy_axis = _pump_detuning_axis(man.omega, gamma)
    covered_axis = _pump_detuning_axis(
        man.omega, gamma, (reached.min(), reached.max()))
    padded_axis = _pump_detuning_axis(
        man.omega, gamma,
        (reached.min() - 30 * gamma, reached.max() + 30 * gamma))
    clamped = reservoir(legacy_axis)
    covered = reservoir(covered_axis)
    padded = reservoir(padded_axis)

    assert covered_axis[0] <= reached.min()
    assert covered_axis[-1] >= reached.max()
    np.testing.assert_allclose(covered, padded, rtol=0.0, atol=1e-6)
    assert np.max(np.abs(clamped - covered)) > 1e-3


def test_coated_reservoir_identity_map_returns_thermal_population():
    transfer = np.broadcast_to(np.eye(2), (4, 2, 2)).copy()
    thermal = np.array([5.0 / 12.0, 7.0 / 12.0])
    populations = _stationary_ground_populations(
        transfer, thermal, cycle_rate=2 * np.pi * 100e3)
    np.testing.assert_allclose(populations, np.broadcast_to(thermal, (4, 2)),
                               rtol=0.0, atol=2e-13)
    np.testing.assert_allclose(populations.sum(axis=1), 1.0,
                               rtol=0.0, atol=2e-15)

    zero_cycle = _stationary_ground_populations(
        np.broadcast_to(np.array([[[0.8, 0.1], [0.2, 0.9]]]), (4, 2, 2)),
        thermal, cycle_rate=0.0)
    np.testing.assert_allclose(zero_cycle, np.broadcast_to(thermal, (4, 2)),
                               rtol=0.0, atol=2e-15)


def test_velocity_average_keeps_pump_and_probe_in_the_same_class():
    prepared = np.array([[0.0, 1.0]])
    probe = np.array([[1.0, 0.0]])
    weights = np.array([0.5, 0.5])
    correlated = _velocity_correlated_average(prepared, probe, weights)
    factorized = (prepared @ weights) * (probe @ weights)
    np.testing.assert_array_equal(correlated, np.array([0.0]))
    assert factorized[0] == 0.25


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


def test_paraffin_coating_has_no_direct_pump_off_optical_loss():
    common = dict(species=RB85_KEY, line="D1", pump_power_mw=0.0,
                  temp_c=45.0, scan_points=401)
    plain = SAS.compute(_params(**common, paraffin_coated=False))
    coated = SAS.compute(_params(**common, paraffin_coated=True))
    np.testing.assert_array_equal(coated["scan"], plain["scan"])
    np.testing.assert_array_equal(coated["alpha_unit"], plain["alpha_unit"])
    assert coated["paraffin_t1_s"] > 0.0


def test_paraffin_population_memory_is_opt_in_and_changes_pumped_sas():
    common = dict(species=RB85_KEY, line="D1", pump_power_mw=0.5,
                  temp_c=45.0, scan_points=401)
    explicit_off = SAS.compute(_params(**common, paraffin_coated=False))
    missing = _params(**common)
    missing.pop("paraffin_coated")
    legacy = SAS.compute(missing)
    coated = SAS.compute(_params(**common, paraffin_coated=True))

    np.testing.assert_array_equal(legacy["alpha_unit"],
                                  explicit_off["alpha_unit"])
    assert np.isfinite(coated["alpha_unit"]).all()
    assert np.all(coated["alpha_unit"] >= 0.0)
    change = np.max(np.abs(coated["alpha_unit"] - explicit_off["alpha_unit"]))
    assert change > 1e-3 * explicit_off["alpha_unit"].max()


def test_paraffin_checkbox_is_inert_for_generic_toy():
    common = dict(species=GENERIC, transitions="single line", scan_points=401)
    plain = SAS.compute(_params(**common, paraffin_coated=False))
    coated = SAS.compute(_params(**common, paraffin_coated=True))
    np.testing.assert_array_equal(coated["scan"], plain["scan"])
    np.testing.assert_array_equal(coated["alpha_unit"], plain["alpha_unit"])


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
    assert [metric["label"] for metric in view["metrics"]] == [
        "Peak OD", "Gaussian Doppler FWHM"
    ]
    assert view["metrics"][0]["tier"] == "hero"
    broad = next(m for m in view["metrics"] if m["label"] == "Gaussian Doppler FWHM")
    assert "not a measured Voigt" in broad["help"]


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
    labels = [m["label"] for m in view["metrics"]]
    assert labels[:4] == [
        "SAS resolution", "Sub-Doppler FWHM", "Half-height edges", "Samples / FWHM"
    ]
    assert [m["label"] for m in view["metrics"] if m.get("tier") == "hero"] == [
        "SAS resolution"
    ]
    assert _metric_number(view, "Samples / FWHM") < 6.0
    assert _metric_number(view, "Scan-edge distance") > 0.0
    assert "Gaussian Doppler FWHM" in labels
    assert "Doppler FWHM" not in labels
    assert "Buffer Gas Broadening" not in {m["label"] for m in view["metrics"]}


def test_default_underresolved_feature_reports_diagnostics_instead_of_lock_hero():
    p = _params()
    raw = SAS.compute(p)
    x = raw["scan"] / (2 * np.pi) / 1e9
    alpha = raw["alpha_unit"] * p["line_strength"]
    T_trans = observables.transmission(alpha, p["cell_mm"] * 1e-3)
    feature = subdoppler_feature(x, T_trans)

    view = SAS.observables(raw, p, include_figures=False)
    labels = [m["label"] for m in view["metrics"]]
    heroes = [m["label"] for m in view["metrics"] if m.get("tier") == "hero"]

    assert feature.detected and feature.status == "resolution-limited"
    assert heroes == ["SAS resolution"]
    assert "Lock Slope" not in labels and "Lock Detuning" not in labels
    assert {"Sub-Doppler FWHM", "Half-height edges", "Samples / FWHM",
            "Scan-edge distance"}.issubset(labels)


def test_pump_on_unresolved_feature_makes_status_hero_and_hides_envelope_slope():
    p = _params(pump_power_mw=0.01, temp_c=200.0, cell_mm=200.0)
    view = SAS.headless_observables(SAS.compute(p), p)
    heroes = [m for m in view["metrics"] if m.get("tier") == "hero"]
    labels = [m["label"] for m in view["metrics"]]
    assert [(m["label"], m["value"]) for m in heroes] == [
        ("SAS resolution", "unresolved")
    ]
    assert "Lock Slope" not in labels and "Lock Detuning" not in labels


def test_legacy_subdoppler_wrapper_keeps_nan_nan_for_a_flat_trace():
    x = np.linspace(-1.0, 1.0, 11)
    width, center = narrowest_subdoppler(x, np.ones_like(x))
    assert np.isnan(width) and np.isnan(center)


def test_rb85_d1_fixed_window_2x_4x_readouts_converge_within_five_percent():
    """Reference-only refinement: freeze the baseline window before 2x/4x solves."""
    p = _params(species=RB85_KEY, line="D1")
    p.update(SAS.recommended_defaults(p)["SAS default"])
    baseline_points = int(p["scan_points"])

    def measure(points, search_window=None):
        q = {**p, "scan_points": points}
        raw = SAS.compute(q)
        x = raw["scan"] / (2 * np.pi) / 1e9
        T_trans = observables.transmission(
            raw["alpha_unit"] * q["line_strength"], q["cell_mm"] * 1e-3)
        feature = subdoppler_feature(x, T_trans, search_window=search_window)
        view = SAS.headless_observables(raw, q)
        return x, T_trans, feature, view

    x1, _, feature1, _ = measure(baseline_points)
    assert feature1.status == "resolution-limited"
    assert feature1.detected and feature1.samples_per_fwhm < 6.0
    # Both half-height edges must come from interpolation, not sample snapping.
    assert not np.any(np.isclose(x1, feature1.left_half_height, rtol=0.0, atol=1e-12))
    assert not np.any(np.isclose(x1, feature1.right_half_height, rtol=0.0, atol=1e-12))

    fixed_window = (
        feature1.center - feature1.fwhm,
        feature1.center + feature1.fwhm,
    )
    # (N-1) sets the sample intervals, so these are exactly 2x and 4x density.
    x2, T2, feature2, view2 = measure(2 * (baseline_points - 1) + 1, fixed_window)
    x4, T4, feature4, view4 = measure(4 * (baseline_points - 1) + 1, fixed_window)

    for feature in (feature2, feature4):
        assert feature.resolved
        assert fixed_window[0] <= feature.center <= fixed_window[1]
        assert feature.scan_edge_distance > feature.fwhm
    assert 1.9 <= feature4.samples_per_fwhm / feature2.samples_per_fwhm <= 2.1
    for view in (view2, view4):
        heroes = [m["label"] for m in view["metrics"] if m.get("tier") == "hero"]
        assert heroes == ["Lock Slope"]
        assert next(m["value"] for m in view["metrics"]
                    if m["label"] == "SAS resolution") == "resolved"

    lock2 = _lock_readout_metrics(
        x2, T2, 1000.0, search_window=fixed_window)
    lock4 = _lock_readout_metrics(
        x4, T4, 1000.0, search_window=fixed_window)
    slope2 = float(lock2[0]["value"].split()[0])
    slope4 = float(lock4[0]["value"].split()[0])

    assert abs(feature4.fwhm - feature2.fwhm) / feature4.fwhm <= 0.05
    assert abs(slope4 - slope2) / slope4 <= 0.05


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
                assert "Gaussian Doppler FWHM" in labels
                assert "Doppler FWHM" not in labels
                if pump_power > 0.0:
                    assert heroes == ["SAS resolution"]
                    assert labels[0] == "SAS resolution"
                    assert not any(label.startswith("Lock ") for label in labels)
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

    # This default grid does not resolve a local discriminator.  The broad
    # Doppler-envelope flank must not be substituted as a lock slope.
    rendered = SAS.observables(raw2, p2)
    assert [item["label"] for item in rendered["figure_views"]] == [
        "Transmission", "Optical density"
    ]
    assert rendered["figure_views"][0]["figure"] is rendered["figure"]
    view2 = SAS.observables(raw2, p2, include_figures=False)
    assert view2["figure_views"] == []
    assert view2["comparison"]["x_unit"] == "MHz"
    assert view2["hero_count"] == 1
    labels = [m["label"] for m in view2["metrics"]]
    heroes = [m["label"] for m in view2["metrics"] if m.get("tier") == "hero"]
    assert labels[0] == "SAS resolution"
    assert heroes == ["SAS resolution"]
    assert not any(label.startswith("Lock ") for label in labels)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print("ok:", fn.__name__)
    print(f"\nAbsorption (OD/SAS) OK ({len(fns)} tests): data, pump-off AutoOD "
          "fidelity, sub-Doppler features, hyperfine pumping, generic.")
