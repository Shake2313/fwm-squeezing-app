"""Jasperse seed-noise UI routing and preserved standalone legacy helpers.

The obsolete constant-N/gain-linear-a laws remain unit-tested only as legacy
helper algebra. Current readout uses classical seed RIN a and calculated shot b.
"""
from copy import deepcopy

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from gabes import observables
from gabes.schemes import fwm


@pytest.fixture
def seeded_spectrum():
    center = fwm.branch_center_GHz(0.9, -1)
    gain = np.array([1.0, 3.0, 8.0])
    matrices = np.empty((3, 2, 2), dtype=complex)
    matrices[:, 0, 0] = matrices[:, 1, 1] = np.sqrt(gain)
    matrices[:, 0, 1] = matrices[:, 1, 0] = np.sqrt(gain - 1.0)
    raw = {
        "D_GHz": 0.9,
        "probe_axis_GHz": center + np.array([-8.0, 0.0, 8.0])*1e-3,
        "raman_center_minus_GHz": center,
        "raman_center_plus_GHz": fwm.branch_center_GHz(0.9, 1),
        "T_canonical_small_signal": matrices,
        "input_seed_power_W": 8e-6,
        "G_s": np.array([0.8, 3.0, 8.0]),
        "G_c": np.array([0.0, 2.0, 7.0]),
        "gain_referred_noise_dB": np.array([-0.4, -4.7, -6.0]),
        "S_dB": np.array([-0.4, -4.7, -6.0]),
        "physical_squeezing_dB": None,
        "eta": 1.0,
        "phase_detail": fwm.PHASE_LEGACY,
        "model_fidelity": fwm.FIDELITY_FAST,
        "N_atoms": 1e18, "sigma_v": 200.0, "n_velocity": 3,
        "Op_A_2pi_GHz": 0.5, "Os_2pi_MHz": 0.1,
        "claim_gate": {"physical_squeezing_prediction": False},
        "squeezing_status": "microscopic atomic diffusion absent",
    }
    return raw


def _ideal_expected(raw, params):
    _, k_seed, _ = fwm.seeded_option_a_wavenumbers(
        raw["D_GHz"], raw["probe_axis_GHz"])
    b = 2*fwm.constants.HBAR*fwm.constants.C_LIGHT*k_seed/raw["input_seed_power_W"]
    a = params.get("seed_classical_rin_scaled", 0.0)*1e-14
    eta = params["after_cell_eff_pct"]*params["qe_pct"]/1e4
    gain = np.abs(raw["T_canonical_small_signal"][:, 0, 0])**2
    nf = 1-eta + eta*(1+a/b)/(2*gain-1)
    return b, nf


def test_current_seed_noise_schema_replaces_retired_models_and_efficiency_override():
    scheme = fwm.FWMScheme()
    specs = {sp.name: sp for sp in scheme.param_schema()}
    retired = {"excess_noise_mode", "excess_noise", "excess_noise_slope", "seed_fano",
               "detection_eff_pct", "gain_closure_enabled"}
    assert retired.isdisjoint(specs)
    a = specs["seed_classical_rin_scaled"]
    b = specs["seed_shot_noise_b"]
    assert a.default == 0.0 and a.recompute is False
    assert a.visible_if == {"mode": fwm.MODE_SEEDED}
    assert b.control == "derived_slider"
    names = [sp.name for sp in scheme.param_schema()]
    assert names.index("seed_shot_noise_b") == names.index("seed_classical_rin_scaled") + 1
    for name in ("seed_classical_rin_scaled", "after_cell_eff_pct", "qe_pct"):
        assert specs[name].recompute is False
        assert name not in scheme.recompute_keys()
        assert name not in scheme.extra_views()[0].param_keys
    for name in ("seed_shot_noise_b", "eta_probe_pct", "eta_conjugate_pct"):
        assert name not in scheme.defaults()
        assert name not in scheme.recompute_keys()
    reset = scheme.recommended_defaults(dict(scheme.defaults(),
        seed_classical_rin_scaled=5.0))[fwm.MODE_SEEDED]
    assert reset["seed_classical_rin_scaled"] == 0.0


def test_focused_and_full_heavy_solves_ignore_current_and_retired_noise_controls(monkeypatch):
    scheme = fwm.FWMScheme()
    params = scheme.defaults()
    changed = dict(params, seed_classical_rin_scaled=12.0,
                   after_cell_eff_pct=50.0, qe_pct=60.0,
                   excess_noise=4.0, excess_noise_slope=2.0, seed_fano=50.0)
    calls = []
    marker = {"computed": True}
    for target, compute in (("compute_spectrum", scheme.compute),
                            ("full_spectrum", scheme.extra_views()[0].compute)):
        monkeypatch.setattr(fwm, target, lambda *a, **kw:
                            calls.append((a, kw)) or marker)
        calls.clear()
        assert compute(params) is marker
        assert compute(changed) is marker
        assert calls[0] == calls[1]


def test_current_a_b_nf_follow_ideal_photon_statistics_without_mutating_cache(seeded_spectrum):
    import matplotlib.pyplot as plt
    scheme = fwm.FWMScheme()
    params = dict(scheme.defaults(), tpd=0.0, seed_classical_rin_scaled=2.0,
                  after_cell_eff_pct=80.0, qe_pct=90.0)
    original = deepcopy(seeded_spectrum)
    b, nf = _ideal_expected(seeded_spectrum, params)
    view = scheme.observables(seeded_spectrum, params)
    try:
        current = view["jasperse_readout"]
        np.testing.assert_allclose(current["NF"], nf, rtol=1e-13)
        np.testing.assert_allclose(current["NF_dB"], 10*np.log10(nf), atol=1e-13)
        np.testing.assert_allclose(current["seed_noise"]["shot_rin_per_Hz"], b)
        assert current["eta_after_cell_QE"] == pytest.approx(0.72)
        assert view["metrics"][0]["label"] == "Squeezing Degree [dB]"
        assert view["metrics"][0]["value"] == f"{10*np.log10(nf[1]):.2f} dB"
        axis = view["figure"].axes[1]
        np.testing.assert_allclose(axis.lines[0].get_ydata(), 10*np.log10(nf))
        assert axis.get_ylabel() == "Squeezing Degree [dB]"
        headless = scheme.headless_observables(seeded_spectrum, params)
        assert headless["metrics"] == view["metrics"]
        b_slider = view["sidebar_readouts"]["seed_shot_noise_b"]["value"]
        assert b_slider == pytest.approx(b[1]/1e-14)
    finally:
        plt.close(view["figure"])
    for name in ("G_s", "G_c", "S_dB", "gain_referred_noise_dB", "T_canonical_small_signal"):
        np.testing.assert_array_equal(seeded_spectrum[name], original[name])


def test_retired_constant_and_gain_linear_noise_never_change_current_nf(seeded_spectrum):
    scheme = fwm.FWMScheme()
    params = dict(scheme.defaults(), tpd=0.0, seed_classical_rin_scaled=1.5)
    baseline = scheme.headless_observables(seeded_spectrum, params)
    stale = dict(params, excess_noise=9.0, excess_noise_mode=fwm.NOISE_GAIN_PROPORTIONAL,
                 excess_noise_slope=7.0, seed_fano=30.0, detection_eff_pct=10.0)
    changed = scheme.headless_observables(seeded_spectrum, stale)
    np.testing.assert_array_equal(changed["jasperse_readout"]["NF"],
                                  baseline["jasperse_readout"]["NF"])
    assert changed["metrics"] == baseline["metrics"]


def test_off_grid_marker_uses_log_of_interpolated_nf_and_lies_on_curve(seeded_spectrum):
    import matplotlib.pyplot as plt
    params = dict(fwm.FWMScheme().defaults(), tpd=4.0, seed_classical_rin_scaled=1.0)
    _, nf = _ideal_expected(seeded_spectrum, params)
    expected_nf = 0.5*(nf[1]+nf[2])
    view = fwm.FWMScheme().observables(seeded_spectrum, params)
    try:
        point = view["passive_l1_operating_point"]
        assert point["NF"] == pytest.approx(expected_nf)
        assert point["NF_dB"] == pytest.approx(10*np.log10(expected_nf))
        axis = view["figure"].axes[1]
        assert axis.collections[0].get_offsets()[0, 1] == pytest.approx(point["NF_dB"])
        line = axis.lines[0]
        assert np.interp(4.0, line.get_xdata(), line.get_ydata()) == pytest.approx(point["NF_dB"])
    finally:
        plt.close(view["figure"])


def test_full_scan_current_seed_noise_uses_each_branch_input_power(seeded_spectrum):
    import matplotlib.pyplot as plt
    params = dict(fwm.FWMScheme().defaults(), seed_classical_rin_scaled=1.0,
                  after_cell_eff_pct=75.0, qe_pct=80.0)
    full = {"D_GHz": 0.9, "minus": seeded_spectrum,
            "plus": dict(seeded_spectrum, input_seed_power_W=4e-6)}
    figure = fwm.FWMScheme().extra_views()[0].render(full, params)
    try:
        for index, branch in enumerate(("minus", "plus")):
            _, nf = _ideal_expected(full[branch], params)
            np.testing.assert_allclose(figure.axes[1].lines[index].get_ydata(),
                                       10*np.log10(nf), atol=1e-13)
        assert figure.axes[1].get_ylabel() == "Squeezing Degree [dB]"
    finally:
        plt.close(figure)

@pytest.mark.parametrize("eta", [0.0, 0.4, 1.0])
def test_excess_noise_matches_source_formula_and_efficiency_limits(eta):
    gain = np.array([1.0, 3.0, 15.0])
    excess = 1.2
    expected = (1.0 - eta) + eta * (1.0 / (2.0 * gain - 1.0) + excess)

    actual = observables.gain_referred_noise_dB(
        gain, gain - 1.0, eta, excess_noise=excess)
    np.testing.assert_allclose(10.0 ** (actual / 10.0), expected, atol=1e-14)
    np.testing.assert_array_equal(
        actual,
        observables.intensity_difference_squeezing_dB(
            gain, gain - 1.0, eta, excess_noise=excess),
    )
    if eta == 0.0:
        np.testing.assert_array_equal(actual, np.zeros_like(gain))

@pytest.mark.parametrize("eta", [0.25, 0.9, 1.0])
def test_excess_noise_crosses_sql_at_the_expected_source_threshold(eta):
    below, threshold, above = [
        observables.gain_referred_noise_dB(3.0, 2.0, eta, excess_noise=value)
        for value in (0.799, 0.8, 0.801)
    ]
    assert below < 0.0
    assert threshold == pytest.approx(0.0, abs=1e-14)
    assert above > 0.0

def test_source_noise_adds_to_existing_linear_noise_without_clipping():
    baseline = np.array([-9.0, -0.3, 0.0, 4.0])
    expected = 10.0 ** (baseline / 10.0) + 0.6 * 2.0
    noisy = observables.add_source_excess_noise_dB(baseline, 0.6, 2.0)
    np.testing.assert_allclose(10.0 ** (noisy / 10.0), expected, atol=1e-14)
    assert np.all(noisy > 0.0)
    np.testing.assert_array_equal(baseline, [-9.0, -0.3, 0.0, 4.0])

def test_zero_excess_is_exact_and_zero_efficiency_removes_the_added_noise():
    baseline = np.array([-300.0, -12.3456789, -0.1, 0.0, 7.891234])
    np.testing.assert_array_equal(
        observables.add_source_excess_noise_dB(baseline, 0.87, 0.0), baseline)
    np.testing.assert_allclose(
        observables.add_source_excess_noise_dB(baseline, 0.0, 5.0), baseline,
        rtol=0.0, atol=1e-14)

    gain = np.array([1.0, 3.0, 15.0])
    old_result = observables.gain_referred_noise_dB(gain, gain - 1.0, 0.87)
    np.testing.assert_array_equal(
        observables.gain_referred_noise_dB(
            gain, gain - 1.0, 0.87, excess_noise=0.0), old_result)

@pytest.mark.parametrize("excess", [-0.01, np.nan, np.inf, -np.inf])
def test_invalid_source_excess_noise_is_rejected(excess):
    with pytest.raises(ValueError):
        observables.add_source_excess_noise_dB(-3.0, 0.87, excess)
    with pytest.raises(ValueError):
        observables.gain_referred_noise_dB(3.0, 2.0, 0.87, excess_noise=excess)
    with pytest.raises(ValueError):
        observables.intensity_difference_squeezing_dB(
            3.0, 2.0, 0.87, excess_noise=excess)

def test_source_noise_array_broadcasts_in_linear_units_and_keeps_zero_exact():
    baseline = np.array([[-8.0, -3.0, 0.0], [-7.0, -2.0, 1.0]])
    noise = np.array([0.0, 0.2, 0.7])
    eta = np.array([[0.5], [0.9]])
    actual = observables.add_source_excess_noise_dB(baseline, eta, noise)
    expected = 10.0 ** (baseline / 10.0) + eta * noise
    np.testing.assert_allclose(10.0 ** (actual / 10.0), expected, atol=1e-14)
    np.testing.assert_array_equal(
        observables.add_source_excess_noise_dB(baseline, eta, np.zeros(3)),
        baseline)
    with pytest.raises(ValueError):
        observables.add_source_excess_noise_dB(baseline, eta, np.ones(4))

@pytest.mark.parametrize("invalid", [-0.01, np.nan, np.inf, -np.inf])
def test_invalid_entry_in_source_noise_array_is_rejected(invalid):
    with pytest.raises(ValueError):
        observables.add_source_excess_noise_dB(
            np.array([-3.0, -4.0]), 0.9, np.array([0.0, invalid]))

def test_gain_proportional_noise_uses_only_amplification_above_unity():
    gain = np.array([0.0, 0.7, 1.0, 3.0, 10.0])
    np.testing.assert_allclose(
        observables.gain_proportional_excess_noise(gain, 0.03),
        [0.0, 0.0, 0.0, 0.06, 0.27],
    )
    np.testing.assert_array_equal(
        observables.gain_proportional_excess_noise(gain, 0.0),
        np.zeros_like(gain),
    )
    assert observables.gain_proportional_excess_noise(4.0, 0.03) == pytest.approx(0.09)
    np.testing.assert_array_equal(gain, [0.0, 0.7, 1.0, 3.0, 10.0])

@pytest.mark.parametrize("invalid", [-0.01, np.nan, np.inf, -np.inf])
def test_gain_proportional_noise_rejects_invalid_slope_and_gain(invalid):
    with pytest.raises(ValueError):
        observables.gain_proportional_excess_noise([1.0, 3.0], invalid)
    with pytest.raises(ValueError):
        observables.gain_proportional_excess_noise([1.0, invalid], 0.01)

def test_gain_proportional_noise_requires_a_single_slope():
    with pytest.raises((TypeError, ValueError)):
        observables.gain_proportional_excess_noise([1.0, 3.0], [0.01, 0.02])

@pytest.mark.parametrize("eta", [0.25, 0.9, 1.0])
def test_gain_dependent_noise_has_finite_optimum_and_recrosses_sql(eta):
    slope = 0.05
    optimum_gain = (1.0 + np.sqrt(2.0 / slope)) / 2.0
    sql_gain = 0.5 + 1.0 / slope
    gains = np.array([
        1.0, optimum_gain - 0.1, optimum_gain, optimum_gain + 0.1,
        sql_gain - 0.1, sql_gain, sql_gain + 0.1,
    ])
    excess = observables.gain_proportional_excess_noise(gains, slope)
    noise_db = observables.gain_referred_noise_dB(
        gains, gains - 1.0, eta, excess_noise=excess)
    expected = (1.0 - eta) + eta * (
        1.0 / (2.0 * gains - 1.0) + slope * (gains - 1.0))
    np.testing.assert_allclose(10.0 ** (noise_db / 10.0), expected, atol=1e-14)
    assert noise_db[0] == pytest.approx(0.0, abs=1e-14)
    assert noise_db[2] < min(noise_db[1], noise_db[3]) < 0.0
    assert noise_db[4] < 0.0
    assert noise_db[5] == pytest.approx(0.0, abs=1e-14)
    assert noise_db[6] > 0.0
