"""Source-referred excess noise and cached seeded-FWM readout checks."""

from copy import deepcopy

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from gabes import observables  # noqa: E402
from gabes.schemes import fwm  # noqa: E402


@pytest.fixture
def seeded_spectrum():
    """A cached spectrum with an existing, non-ideal noise contribution."""
    center = fwm.branch_center_GHz(0.9, -1)
    noise_db = np.array([-0.4, -4.7, -6.0])
    return {
        "D_GHz": 0.9,
        "probe_axis_GHz": center + np.array([-8.0, 0.0, 8.0]) * 1e-3,
        "raman_center_minus_GHz": center,
        "raman_center_plus_GHz": fwm.branch_center_GHz(0.9, 1),
        "G_s": np.array([0.8, 3.0, 8.0]),
        "G_c": np.array([0.0, 2.0, 7.0]),
        "gain_referred_noise_dB": noise_db,
        "S_dB": noise_db,
        "physical_squeezing_dB": None,
        "eta": 0.8,
        "phase_detail": fwm.PHASE_LEGACY,
        "model_fidelity": fwm.FIDELITY_FAST,
        "N_atoms": 1e18,
        "sigma_v": 200.0,
        "n_velocity": 3,
        "Op_A_2pi_GHz": 0.5,
        "Os_2pi_MHz": 0.1,
        "claim_gate": {"physical_squeezing_prediction": False},
        "squeezing_status": "unavailable: microscopic atomic diffusion absent",
    }


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


def test_excess_noise_control_is_seeded_only_and_does_not_invalidate_solve():
    scheme = fwm.FWMScheme()
    spec = {item.name: item for item in scheme.param_schema()}["excess_noise"]
    assert spec.label == "Excess Noise N"
    assert spec.group == "Detection & scaling"
    assert (spec.vmin, spec.vmax, spec.step) == (0.0, 5.0, 0.01)
    assert spec.default == 0.0
    assert not spec.hidden and not spec.advanced
    assert spec.visible_if == {"mode": fwm.MODE_SEEDED}
    assert spec.recompute is False
    assert "excess_noise" not in scheme.recompute_keys()
    assert scheme.defaults()["excess_noise"] == 0.0

    changed = dict(scheme.defaults(), excess_noise=3.0)
    assert scheme.recommended_defaults(changed)[fwm.MODE_SEEDED]["excess_noise"] == 0.0
    full_view = scheme.extra_views()[0]
    assert full_view.render_with_params is True
    assert "excess_noise" not in full_view.param_keys


def test_focused_and_full_heavy_solves_ignore_the_readout_noise(monkeypatch):
    scheme = fwm.FWMScheme()
    params = scheme.defaults()
    calls = []
    marker = {"computed": True}

    def record_solve(*args, **kwargs):
        calls.append((args, kwargs))
        return marker

    monkeypatch.setattr(fwm, "compute_spectrum", record_solve)
    assert scheme.compute(params) is marker
    assert scheme.compute(dict(params, excess_noise=4.0)) is marker
    assert calls[0] == calls[1]

    calls.clear()
    monkeypatch.setattr(fwm, "full_spectrum", record_solve)
    full_view = scheme.extra_views()[0]
    assert full_view.compute(params) is marker
    assert full_view.compute(dict(params, excess_noise=4.0)) is marker
    assert calls[0] == calls[1]


def test_focused_readout_noise_updates_curve_marker_and_metric_without_mutation(
        seeded_spectrum):
    import matplotlib.pyplot as plt

    scheme = fwm.FWMScheme()
    params = dict(scheme.defaults(), tpd=4.0, excess_noise=2.0)
    original = deepcopy(seeded_spectrum)
    expected = 10.0 * np.log10(
        10.0 ** (original["gain_referred_noise_dB"] / 10.0)
        + original["eta"] * params["excess_noise"])
    operating_x = original["raman_center_minus_GHz"] + params["tpd"] * 1e-3
    operating_noise = np.interp(operating_x, original["probe_axis_GHz"], expected)

    rendered = scheme.observables(seeded_spectrum, params)
    try:
        gain_axis, noise_axis = rendered["figure"].axes
        np.testing.assert_allclose(noise_axis.lines[0].get_ydata(), expected)
        np.testing.assert_array_equal(gain_axis.lines[0].get_ydata(), original["G_s"])
        assert float(noise_axis.collections[0].get_offsets()[0, 1]) == pytest.approx(
            operating_noise)
        assert rendered["metrics"][0]["value"] == f"{operating_noise:.2f} dB"
        assert "physical squeezing unavailable" in rendered["metrics"][0]["delta"]
        ymin, ymax = noise_axis.get_ylim()
        assert ymin <= 0.0 < np.min(expected) <= np.max(expected) < ymax

        repeated = scheme.headless_observables(seeded_spectrum, params)
        assert repeated["metrics"] == rendered["metrics"]
        zero = scheme.headless_observables(
            seeded_spectrum, dict(params, excess_noise=0.0))
        legacy_params = dict(params)
        legacy_params.pop("excess_noise")
        assert scheme.headless_observables(seeded_spectrum, legacy_params)["metrics"] == zero["metrics"]
        assert zero["metrics"][0]["value"] != rendered["metrics"][0]["value"]
    finally:
        plt.close(rendered["figure"])

    for key in ("G_s", "G_c", "gain_referred_noise_dB", "S_dB"):
        np.testing.assert_array_equal(seeded_spectrum[key], original[key])
    assert seeded_spectrum["physical_squeezing_dB"] is None
    assert seeded_spectrum["claim_gate"] == original["claim_gate"]


def test_full_scan_render_applies_current_noise_to_both_cached_branches(
        seeded_spectrum):
    import matplotlib.pyplot as plt

    full = {
        "D_GHz": 0.9,
        "minus": seeded_spectrum,
        "plus": dict(seeded_spectrum, eta=0.5),
    }
    original = deepcopy(full)
    view = fwm.FWMScheme().extra_views()[0]
    figures = []
    try:
        baseline = view.render(full)
        figures.append(baseline)
        noisy = view.render(full, {"excess_noise": 2.0})
        figures.append(noisy)
        reset = view.render(full, {"excess_noise": 0.0})
        figures.append(reset)
        for index, key in enumerate(("minus", "plus")):
            source = original[key]
            expected = 10.0 * np.log10(
                10.0 ** (source["gain_referred_noise_dB"] / 10.0)
                + source["eta"] * 2.0)
            np.testing.assert_allclose(noisy.axes[1].lines[index].get_ydata(), expected)
            np.testing.assert_array_equal(
                baseline.axes[1].lines[index].get_ydata(), source["gain_referred_noise_dB"])
            np.testing.assert_array_equal(
                reset.axes[1].lines[index].get_ydata(), source["gain_referred_noise_dB"])
            np.testing.assert_array_equal(full[key]["gain_referred_noise_dB"], source["gain_referred_noise_dB"])
            np.testing.assert_array_equal(noisy.axes[0].lines[index].get_ydata(), source["G_s"])
            assert noisy.axes[1].get_ylim()[1] > np.max(expected) > 0.0
    finally:
        for figure in figures:
            plt.close(figure)
