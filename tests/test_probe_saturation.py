"""Finite-power Natural-Rb D2 probe model and Gaussian propagation checks."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gabes import probe_saturation, schemes, species  # noqa: E402
from gabes.schemes.sas import GENERIC  # noqa: E402


SAS = schemes.get("sas")


def _params(**overrides):
    params = SAS.defaults()
    params.update(
        species="Rb (natural)",
        line="D2",
        pump_power_mw=0.0,
        temp_c=25.0,
        scan_points=401,
    )
    params.update(overrides)
    return params


def test_probe_controls_are_natural_rb_d2_only_and_generic_is_hidden():
    specs = {spec.name: spec for spec in SAS.param_schema()}

    assert specs["species"].choices == tuple(species.SPECIES_ORDER)
    assert GENERIC not in specs["species"].choices
    assert specs["transitions"].hidden
    assert specs["splitting"].hidden
    expected_gate = {"species": "Rb (natural)", "line": "D2"}
    assert specs["probe_power_mw"].visible_if == expected_gate
    assert specs["probe_waist_mm"].visible_if == expected_gate
    assert specs["probe_power_mw"].default == 0.0
    assert specs["probe_power_mw"].recompute
    assert specs["probe_waist_mm"].recompute


def test_autood_d2_saturation_scale_and_gaussian_waist_convention():
    manifold = species.build_manifold(species.RB85, "D2")
    saturation_intensity = (
        probe_saturation.closed_two_level_saturation_intensity(
            manifold.nu0, manifold.gamma))

    np.testing.assert_allclose(
        saturation_intensity / 10.0,
        1.6693356214912014,
        rtol=0.0,
        atol=1e-12,
    )
    narrow = probe_saturation.gaussian_peak_intensity(0.25, 1.0)
    equal_intensity = probe_saturation.gaussian_peak_intensity(1.0, 2.0)
    np.testing.assert_array_equal(narrow, equal_intensity)


def test_power_broadening_has_two_level_peak_and_width():
    detuning = np.linspace(-100.0, 100.0, 40_001)
    spacing = detuning[1] - detuning[0]
    hwhm = 1.0
    weak = (hwhm / np.pi) / (detuning**2 + hwhm**2)
    table = probe_saturation.power_broadened_profile_table(
        weak, spacing, hwhm, 8.0, np.array([0.0]))
    saturated = table[0]
    centre = detuning.size // 2
    right = centre + int(np.argmin(
        np.abs(saturated[centre:] - 0.5 * saturated[centre])))

    np.testing.assert_allclose(
        saturated[centre], weak[centre] / 9.0, rtol=2e-5)
    np.testing.assert_allclose(detuning[right], 3.0, rtol=0.0, atol=spacing)
    np.testing.assert_allclose(
        np.trapezoid(saturated, detuning), 1.0 / 3.0, rtol=0.025)


def test_gaussian_propagator_reduces_to_beer_lambert_without_saturation():
    x_grid = np.linspace(0.0, 12.0, 97)
    alpha_weak = np.array([10.0, 100.0])
    alpha_table = np.broadcast_to(alpha_weak, (x_grid.size, 2)).copy()

    transmission, effective_alpha = (
        probe_saturation.propagate_gaussian_spectrum(
            x_grid, alpha_table, alpha_weak, length_m=0.03))

    np.testing.assert_allclose(
        transmission, np.exp(-alpha_weak * 0.03), rtol=2e-13, atol=0.0)
    np.testing.assert_allclose(effective_alpha, alpha_weak, rtol=2e-13)


def test_gaussian_propagator_handles_exact_local_transparency_without_gain():
    x_grid = np.linspace(0.0, 12.0, 121)
    alpha_weak = np.array([0.15, 2.0])
    alpha_table = np.broadcast_to(alpha_weak, (x_grid.size, 2)).copy()
    alpha_table[0, 0] = 0.0
    alpha_table[17, 1] = 0.0

    transmission, effective_alpha = (
        probe_saturation.propagate_gaussian_spectrum(
            x_grid, alpha_table, alpha_weak, length_m=0.075))

    assert np.isfinite(transmission).all()
    assert np.isfinite(effective_alpha).all()
    assert np.all((0.0 <= transmission) & (transmission <= 1.0))
    max_local_alpha = np.maximum(alpha_table.max(axis=0), alpha_weak)
    assert np.all(effective_alpha <= max_local_alpha * (1.0 + 1e-13))


def test_gaussian_propagator_is_continuous_at_zero_cell_length():
    x_grid = np.linspace(0.0, 12.0, 121)
    alpha_weak = np.array([1.0])
    alpha_table = np.linspace(0.1, 1.0, x_grid.size)[:, None]

    zero_transmission, zero_alpha = (
        probe_saturation.propagate_gaussian_spectrum(
            x_grid, alpha_table, alpha_weak, length_m=0.0))
    tiny_transmission, tiny_alpha = (
        probe_saturation.propagate_gaussian_spectrum(
            x_grid, alpha_table, alpha_weak, length_m=1e-12))

    np.testing.assert_array_equal(zero_transmission, np.ones(1))
    np.testing.assert_array_equal(tiny_alpha, zero_alpha)
    np.testing.assert_allclose(
        tiny_transmission, np.exp(-zero_alpha * 1e-12), rtol=0.0, atol=0.0)


def test_zero_probe_power_is_exact_legacy_limit_even_when_new_keys_are_missing():
    explicit = _params(probe_power_mw=0.0, probe_waist_mm=4.5)
    legacy = dict(explicit)
    legacy.pop("probe_power_mw")
    legacy.pop("probe_waist_mm")

    explicit_raw = SAS.compute(explicit)
    legacy_raw = SAS.compute(legacy)

    for key in ("scan", "alpha_unit", "chi_real_unit"):
        np.testing.assert_array_equal(explicit_raw[key], legacy_raw[key])
    assert explicit_raw["probe_model"] == "weak_probe_legacy"
    assert explicit_raw["probe_alpha_table_unit"] is None


def test_equal_probe_power_per_waist_squared_gives_identical_local_tables():
    narrow = SAS.compute(_params(probe_power_mw=0.25, probe_waist_mm=1.0))
    wide = SAS.compute(_params(probe_power_mw=1.0, probe_waist_mm=2.0))

    np.testing.assert_array_equal(narrow["scan"], wide["scan"])
    np.testing.assert_array_equal(
        narrow["probe_intensity_log_grid"], wide["probe_intensity_log_grid"])
    np.testing.assert_array_equal(
        narrow["probe_alpha_table_unit"], wide["probe_alpha_table_unit"])


def test_finite_probe_matches_supplied_autood_natural_rb_d2_anchor():
    """Pinned independent result from references/AutoOD-NatRbD2."""
    low_od_params = _params(
        temp_c=25.0,
        cell_mm=75.0,
        probe_power_mw=0.5,
        probe_waist_mm=1.0,
        scan_points=2001,
    )
    low_od_raw = SAS.compute(low_od_params)
    low_transmission, low_effective_alpha = (
        probe_saturation.propagate_gaussian_spectrum(
            low_od_raw["probe_intensity_log_grid"],
            low_od_raw["probe_alpha_table_unit"],
            low_od_raw["alpha_unit"],
            low_od_params["cell_mm"] * 1e-3,
        ))
    low_natural_od = low_effective_alpha * low_od_params["cell_mm"] * 1e-3
    np.testing.assert_allclose(
        low_transmission.min(), 0.719396022587, rtol=3e-5, atol=0.0)
    np.testing.assert_allclose(
        low_natural_od.max(), 0.329343276529, rtol=2e-4, atol=0.0)

    high_od_params = _params(
        temp_c=90.0,
        cell_mm=12.5,
        probe_power_mw=0.5,
        probe_waist_mm=1.0,
        scan_points=2001,
    )
    high_od_raw = SAS.compute(high_od_params)
    _high_transmission, high_effective_alpha = (
        probe_saturation.propagate_gaussian_spectrum(
            high_od_raw["probe_intensity_log_grid"],
            high_od_raw["probe_alpha_table_unit"],
            high_od_raw["alpha_unit"],
            high_od_params["cell_mm"] * 1e-3,
        ))
    high_natural_od = high_effective_alpha * high_od_params["cell_mm"] * 1e-3

    # High-OD Tmin is ~1e-11 and therefore a brittle relative anchor; pin the
    # well-conditioned natural optical depth instead.
    np.testing.assert_allclose(
        high_natural_od.max(), 25.1689672059, rtol=2e-4, atol=0.0)


def test_finite_probe_reduces_effective_od_and_suppresses_undefined_dispersion():
    weak_params = _params(cell_mm=75.0, probe_power_mw=0.0)
    finite_params = _params(cell_mm=75.0, probe_power_mw=0.5)
    weak_raw = SAS.compute(weak_params)
    finite_raw = SAS.compute(finite_params)
    weak_view = SAS.headless_observables(weak_raw, weak_params)
    finite_view = SAS.headless_observables(finite_raw, finite_params)

    weak_peak = float(next(
        metric["value"] for metric in weak_view["metrics"]
        if metric["label"] == "Peak OD"))
    finite_peak = float(next(
        metric["value"] for metric in finite_view["metrics"]
        if metric["label"] == "Effective peak OD"))
    assert finite_peak < weak_peak
    finite_labels = {metric["label"] for metric in finite_view["metrics"]}
    assert "Peak probe saturation" in finite_labels
    assert "Peak phase shift" not in finite_labels
    assert [table["title"] for table in finite_view["tables"]] == [
        "Hyperfine lines", "Finite-probe model"]
    assert "Dispersion is unavailable" in finite_view["tables"][1]["markdown"]


def test_finite_pump_probe_approximation_stays_passive_and_saturates_monotonically():
    effective_peaks = []
    for probe_power in (0.05, 0.5):
        params = _params(
            pump_power_mw=0.5,
            probe_power_mw=probe_power,
            cell_mm=50.0,
        )
        raw = SAS.compute(params)
        assert np.isfinite(raw["probe_alpha_table_unit"]).all()
        assert np.all(raw["probe_alpha_table_unit"] >= 0.0)
        _transmission, effective_alpha = (
            probe_saturation.propagate_gaussian_spectrum(
                raw["probe_intensity_log_grid"],
                raw["probe_alpha_table_unit"],
                raw["alpha_unit"],
                params["cell_mm"] * 1e-3,
            ))
        effective_peaks.append(float(np.max(effective_alpha)))

    assert effective_peaks[1] < effective_peaks[0]


def test_passivity_guard_covers_strong_pump_probe_and_small_waist_boundaries():
    cases = (
        dict(
            pump_power_mw=0.5,
            probe_power_mw=0.5,
            probe_waist_mm=1.0,
            temp_c=40.0,
            cell_mm=75.0,
            transit_khz=100.0,
            paraffin_coated=True,
        ),
        dict(
            pump_power_mw=2.0,
            waist_mm=0.1,
            probe_power_mw=0.5,
            probe_waist_mm=1.0,
            temp_c=40.0,
            cell_mm=75.0,
            transit_khz=5.0,
            paraffin_coated=False,
        ),
        dict(
            pump_power_mw=2.0,
            waist_mm=0.1,
            probe_power_mw=2.0,
            probe_waist_mm=0.1,
            temp_c=20.0,
            cell_mm=200.0,
            transit_khz=5.0,
            paraffin_coated=False,
        ),
    )

    saw_geometry_warning = False
    for overrides in cases:
        params = _params(**overrides)
        raw = SAS.compute(params)
        table = raw["probe_alpha_table_unit"]
        assert np.isfinite(table).all()
        assert np.all(table >= 0.0)
        assert raw["probe_passivity_limited"]
        assert 0.0 <= raw["probe_residual_scale_min"] < 1.0

        _transmission, effective_alpha = (
            probe_saturation.propagate_gaussian_spectrum(
                raw["probe_intensity_log_grid"],
                table,
                raw["alpha_unit"],
                params["cell_mm"] * 1e-3,
            ))
        max_local_alpha = np.maximum(table.max(axis=0), raw["alpha_unit"])
        assert np.all(effective_alpha <= max_local_alpha * (1.0 + 1e-13))

        view = SAS.headless_observables(raw, params)
        labels = {metric["label"]: metric["value"]
                  for metric in view["metrics"]}
        assert labels["Finite-probe approximation"] == "Passivity-limited"
        saw_geometry_warning |= "Probe beam geometry" in labels

    assert saw_geometry_warning


def test_passivity_limited_probe_power_sweep_is_continuous_and_monotonic():
    base = dict(
        pump_power_mw=0.5,
        waist_mm=1.0,
        probe_waist_mm=1.0,
        temp_c=40.0,
        cell_mm=75.0,
        transit_khz=100.0,
        paraffin_coated=True,
    )
    effective_peaks = []
    for probe_power in (0.0, 0.005, 0.01, 0.02, 0.05, 0.1):
        params = _params(**base, probe_power_mw=probe_power)
        raw = SAS.compute(params)
        if probe_power == 0.0:
            effective_alpha = raw["alpha_unit"]
        else:
            _transmission, effective_alpha = (
                probe_saturation.propagate_gaussian_spectrum(
                    raw["probe_intensity_log_grid"],
                    raw["probe_alpha_table_unit"],
                    raw["alpha_unit"],
                    params["cell_mm"] * 1e-3,
                ))
        effective_peaks.append(float(np.max(effective_alpha)))

    assert np.all(np.diff(effective_peaks) <= 1e-12)


def test_passivity_limited_result_converges_with_intensity_radial_and_scan_grids(
        monkeypatch):
    params = _params(
        pump_power_mw=0.5,
        waist_mm=1.0,
        probe_power_mw=0.1,
        probe_waist_mm=1.0,
        temp_c=40.0,
        cell_mm=75.0,
        transit_khz=100.0,
        paraffin_coated=True,
        scan_points=401,
    )
    coarse = SAS.compute(params)
    assert coarse["probe_passivity_limited"]

    def propagated(raw, radial_order):
        return probe_saturation.propagate_gaussian_spectrum(
            raw["probe_intensity_log_grid"],
            raw["probe_alpha_table_unit"],
            raw["alpha_unit"],
            params["cell_mm"] * 1e-3,
            radial_order=radial_order,
        )[1]

    coarse_32 = propagated(coarse, 32)
    coarse_64 = propagated(coarse, 64)
    with monkeypatch.context() as patch:
        patch.setattr(probe_saturation, "SATURATION_LOG_STEP", 0.05)
        fine_intensity = SAS.compute(params)
    fine_32 = propagated(fine_intensity, 32)

    scale = float(np.max(coarse_64))
    assert np.max(np.abs(coarse_32 - coarse_64)) / scale < 1e-3
    assert np.max(np.abs(coarse_32 - fine_32)) / scale < 1e-3

    dense_scan = SAS.compute({**params, "scan_points": 801})
    dense_32 = propagated(dense_scan, 32)
    dense_on_coarse = np.interp(
        coarse["scan"], dense_scan["scan"], dense_32)
    assert abs(np.max(coarse_32) - np.max(dense_32)) / np.max(dense_32) < 1e-3
    assert np.max(np.abs(coarse_32 - dense_on_coarse)) / np.max(dense_32) < 0.01


def test_unsupported_line_ignores_stale_probe_session_value_exactly():
    inactive = _params(
        line="D1", probe_power_mw=1.0, probe_waist_mm=0.2)
    weak = {**inactive, "probe_power_mw": 0.0, "probe_waist_mm": 1.0}
    inactive_raw = SAS.compute(inactive)
    weak_raw = SAS.compute(weak)

    for key in ("scan", "alpha_unit", "chi_real_unit"):
        np.testing.assert_array_equal(inactive_raw[key], weak_raw[key])
    assert not inactive_raw["probe_saturation_supported"]
    assert inactive_raw["probe_model"] == "weak_probe_legacy"
