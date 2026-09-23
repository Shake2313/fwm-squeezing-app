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
    specs = {item.name: item for item in scheme.param_schema()}
    spec = specs["excess_noise"]
    assert spec.label == "Excess Noise N"
    assert spec.group == "Detection & scaling"
    assert (spec.vmin, spec.vmax, spec.step) == (0.0, 5.0, 0.01)
    assert spec.default == 0.0
    assert not spec.hidden and not spec.advanced
    assert spec.visible_if["mode"] == fwm.MODE_SEEDED
    constant_modes = spec.visible_if["excess_noise_mode"]
    if not isinstance(constant_modes, (tuple, list)):
        constant_modes = (constant_modes,)
    assert fwm.NOISE_CONSTANT in constant_modes
    assert fwm.NOISE_GAIN_PROPORTIONAL not in constant_modes
    assert spec.recompute is False
    assert "excess_noise" not in scheme.recompute_keys()
    assert scheme.defaults()["excess_noise"] == 0.0

    mode = specs["excess_noise_mode"]
    assert mode.control == "segmented"
    assert mode.choices == (fwm.NOISE_CONSTANT, fwm.NOISE_GAIN_PROPORTIONAL)
    assert mode.choice_labels == {
        fwm.NOISE_CONSTANT: "Constant Noise N",
        fwm.NOISE_GAIN_PROPORTIONAL: "Gain proportional Noise a",
    }
    assert mode.visible_if == {"mode": fwm.MODE_SEEDED}
    assert not mode.applies_defaults

    slope = specs["excess_noise_slope"]
    assert (slope.vmin, slope.vmax, slope.step) == (0.0, 1.0, 0.0001)
    assert slope.format == "%.4f"
    assert slope.visible_if == {
        "mode": fwm.MODE_SEEDED,
        "excess_noise_mode": fwm.NOISE_GAIN_PROPORTIONAL,
    }
    noise_defaults = {
        "excess_noise_mode": fwm.NOISE_CONSTANT,
        "excess_noise": 0.0,
        "excess_noise_slope": 0.0,
    }
    changed = dict(scheme.defaults(), excess_noise=3.0, excess_noise_slope=0.04,
                   excess_noise_mode=fwm.NOISE_GAIN_PROPORTIONAL)
    reset = scheme.recommended_defaults(changed)[fwm.MODE_SEEDED]
    full_view = scheme.extra_views()[0]
    assert full_view.render_with_params is True
    for key, default in noise_defaults.items():
        assert specs[key].group == "Detection & scaling"
        assert not specs[key].hidden and not specs[key].advanced
        assert specs[key].recompute is False
        assert key not in scheme.recompute_keys()
        assert key not in full_view.param_keys
        assert scheme.defaults()[key] == default
        assert reset[key] == default


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
    gain_mode = dict(params, excess_noise=4.0, excess_noise_slope=0.04,
                     excess_noise_mode=fwm.NOISE_GAIN_PROPORTIONAL)
    assert scheme.compute(gain_mode) is marker
    assert calls[0] == calls[1] == calls[2]

    calls.clear()
    monkeypatch.setattr(fwm, "full_spectrum", record_solve)
    full_view = scheme.extra_views()[0]
    assert full_view.compute(params) is marker
    assert full_view.compute(dict(params, excess_noise=4.0)) is marker
    assert full_view.compute(gain_mode) is marker
    assert calls[0] == calls[1] == calls[2]


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


def test_noise_modes_keep_their_values_and_ignore_the_inactive_control():
    gains = np.array([0.8, 3.0, 8.0])
    params = dict(fwm.FWMScheme().defaults(), excess_noise=2.0,
                  excess_noise_slope=0.04)
    original = deepcopy(params)
    assert fwm._seeded_excess_noise(gains, params) == 2.0
    assert fwm._seeded_excess_noise(
        gains, dict(params, excess_noise_slope=0.8)) == 2.0

    selected = dict(params, excess_noise_mode=fwm.NOISE_GAIN_PROPORTIONAL)
    expected = np.array([0.0, 0.08, 0.28])
    np.testing.assert_allclose(fwm._seeded_excess_noise(gains, selected), expected)
    np.testing.assert_array_equal(
        fwm._seeded_excess_noise(gains, selected),
        fwm._seeded_excess_noise(gains, dict(selected, excess_noise=5.0)),
    )
    assert fwm._seeded_excess_noise(
        gains, dict(selected, excess_noise_mode=fwm.NOISE_CONSTANT)) == 2.0
    assert fwm._seeded_excess_noise(gains, {"excess_noise": 2.0}) == 2.0
    assert fwm._seeded_excess_noise(
        gains, dict(params, excess_noise_mode=None)) == 2.0
    assert params == original


def test_gain_noise_off_grid_readout_uses_displayed_gain_and_keeps_marker_on_curve(
        seeded_spectrum):
    import matplotlib.pyplot as plt

    scheme = fwm.FWMScheme()
    params = dict(scheme.defaults(), tpd=4.0, excess_noise=4.0,
                  excess_noise_mode=fwm.NOISE_GAIN_PROPORTIONAL,
                  excess_noise_slope=0.04)
    original = deepcopy(seeded_spectrum)
    center = original["raman_center_minus_GHz"]
    probe_axis = original["probe_axis_GHz"]
    delta_axis = (probe_axis - center) * 1e3
    operating_x = center + params["tpd"] * 1e-3
    operating_gain = np.interp(operating_x, probe_axis, original["G_s"])
    operating_baseline = np.interp(
        operating_x, probe_axis, original["gain_referred_noise_dB"])
    operating_excess = params["excess_noise_slope"] * (operating_gain - 1.0)
    operating_noise = 10.0 * np.log10(
        10.0 ** (operating_baseline / 10.0) + original["eta"] * operating_excess)
    expected_samples = 10.0 * np.log10(
        10.0 ** (original["gain_referred_noise_dB"] / 10.0)
        + original["eta"] * params["excess_noise_slope"]
        * np.maximum(original["G_s"] - 1.0, 0.0))

    rendered = scheme.observables(seeded_spectrum, params)
    try:
        gain_axis, noise_axis = rendered["figure"].axes
        noise_x, noise_y = noise_axis.lines[0].get_data()
        np.testing.assert_array_equal(gain_axis.lines[0].get_xdata(), delta_axis)
        np.testing.assert_array_equal(gain_axis.lines[0].get_ydata(), original["G_s"])
        assert len(noise_x) == len(delta_axis) + 1
        np.testing.assert_allclose(np.interp(delta_axis, noise_x, noise_y), expected_samples)
        selected = np.isclose(noise_x, params["tpd"], rtol=0.0, atol=1e-10)
        assert np.count_nonzero(selected) == 1
        assert noise_y[selected][0] == pytest.approx(operating_noise)
        assert float(noise_axis.collections[0].get_offsets()[0, 1]) == pytest.approx(
            operating_noise)
        assert rendered["metrics"][0]["value"] == f"{operating_noise:.2f} dB"
        assert scheme.headless_observables(seeded_spectrum, params)["metrics"] == rendered["metrics"]
        markdown = "\n".join(table["markdown"] for table in rendered["tables"])
        assert (
            f"| Excess Noise N (linear SQL units) | {operating_excess:.6g} |"
            in markdown
        )
        assert (
            "| Detected excess noise ηN (linear SQL units) | "
            f"{original['eta'] * operating_excess:.6g} |" in markdown
        )
        repeated = scheme.headless_observables(
            seeded_spectrum, dict(params, excess_noise=0.0))
        assert repeated["metrics"] == rendered["metrics"]
        zero_slope = scheme.headless_observables(
            seeded_spectrum, dict(params, excess_noise_slope=0.0))
        zero_constant = scheme.headless_observables(
            seeded_spectrum, dict(params, excess_noise=0.0,
                                  excess_noise_mode=fwm.NOISE_CONSTANT))
        assert zero_slope["metrics"] == zero_constant["metrics"]
    finally:
        plt.close(rendered["figure"])

    for key in ("G_s", "G_c", "gain_referred_noise_dB", "S_dB"):
        np.testing.assert_array_equal(seeded_spectrum[key], original[key])


def test_full_gain_noise_render_uses_each_branch_gain_and_efficiency(seeded_spectrum):
    import matplotlib.pyplot as plt

    full = {
        "D_GHz": 0.9,
        "minus": seeded_spectrum,
        "plus": dict(seeded_spectrum, eta=0.5, G_s=np.array([1.0, 2.0, 4.0])),
    }
    original = deepcopy(full)
    params = {"excess_noise_mode": fwm.NOISE_GAIN_PROPORTIONAL,
              "excess_noise_slope": 0.1, "excess_noise": 4.0}
    figure = fwm.FWMScheme().extra_views()[0].render(full, params)
    try:
        for index, key in enumerate(("minus", "plus")):
            source = original[key]
            expected = 10.0 * np.log10(
                10.0 ** (source["gain_referred_noise_dB"] / 10.0)
                + source["eta"] * 0.1 * np.maximum(source["G_s"] - 1.0, 0.0))
            np.testing.assert_allclose(figure.axes[1].lines[index].get_ydata(), expected)
            np.testing.assert_array_equal(
                figure.axes[0].lines[index].get_ydata(), source["G_s"])
            np.testing.assert_array_equal(
                full[key]["gain_referred_noise_dB"], source["gain_referred_noise_dB"])
    finally:
        plt.close(figure)


def test_hidden_noise_sliders_retain_values_across_streamlit_reruns():
    """Exercise real widget cleanup without running the atomic solver.

    AppTest stubs Streamlit's component manager, so the ScrubField cannot
    mount and every numeric knob falls back to a slider — which is exactly the
    path this test is here to pin (gabes_ui/scrub.py).
    """
    import ast
    from pathlib import Path

    from streamlit.proto.WidgetStates_pb2 import WidgetStates
    from streamlit.testing.v1 import AppTest

    source = (Path(__file__).resolve().parents[1] / "streamlit_app.py").read_text(
        encoding="utf-8-sig")
    tree = ast.parse(source)
    initialization_names = {"specs", "defaults_version", "defaults_key"}
    initialization = "\n".join(
        ast.get_source_segment(source, node)
        for node in tree.body
        if (
            isinstance(node, ast.Assign)
            and any(isinstance(target, ast.Name) and target.id in initialization_names
                    for target in node.targets)
        ) or (
            isinstance(node, ast.If)
            and any(isinstance(child, ast.Name) and child.id == "defaults_key"
                    for child in ast.walk(node.test))
        )
    )
    # Run the application's actual initialization and parameter collection:
    # both widget cleanup and restored frontend slider values matter. The
    # trailing loop mirrors the tail of gabes_ui.controls.render_rail, which
    # collects the knobs that were not rendered this run.
    app = (
        "import streamlit as st\nfrom gabes.schemes import fwm\n"
        "from gabes_ui.controls import param_visible, render_param, skey\n"
        "\nscheme = fwm.FWMScheme()\n"
        + initialization
        + "\nparams = {}\n"
        "for sp in specs:\n"
        "    if sp.name in ('excess_noise_mode', 'excess_noise', 'excess_noise_slope') "
        "and param_visible(scheme.name, sp):\n"
        "        params[sp.name] = render_param(st, scheme.name, sp, scheme)\n"
        "for sp in specs:\n"
        "    if sp.name not in params:\n"
        "        params[sp.name] = st.session_state[skey(scheme.name, sp.name)]\n"
        "\nst.metric('Source excess noise', "
        "float(fwm._seeded_excess_noise(3.0, params)))\n"
    )
    app_test = AppTest.from_string(app).run()
    assert not app_test.exception
    from gabes_ui import scrub
    assert not scrub.available(), "the sliders below are the ScrubField fallback"
    assert app_test.slider, "numeric knobs must still render something drivable"

    def rerun(mode, **slider_values):
        # AppTest's ButtonGroup serializer assumes a multiselect and fails for
        # formatted single-selection segmented controls. Send their real
        # browser protobuf state instead, including the old visible slider.
        widget_state = WidgetStates()
        for slider in app_test.slider:
            if slider.key in slider_values:
                slider.set_value(slider_values[slider.key])
            widget_state.widgets.add().CopyFrom(slider._widget_state)
        selector = app_test.get("button_group")[0]
        selected = widget_state.widgets.add()
        selected.id = selector.id
        selected.int_array_value.data[:] = [
            (fwm.NOISE_CONSTANT, fwm.NOISE_GAIN_PROPORTIONAL).index(mode)]
        app_test._run(widget_state=widget_state)
        assert not app_test.exception

    def assert_restored_slider(key, expected):
        slider = app_test.slider(key=key)
        assert slider.value == expected
        # Session state alone can look correct while a recreated browser
        # slider shows its default zero and sends that zero back next time.
        assert slider.proto.set_value
        assert list(slider.proto.value) == [expected]

    rerun(fwm.NOISE_CONSTANT, fwm__excess_noise=2.0)
    assert float(app_test.metric[0].value) == 2.0

    rerun(fwm.NOISE_GAIN_PROPORTIONAL)
    assert app_test.session_state["fwm__excess_noise"] == 2.0
    rerun(fwm.NOISE_GAIN_PROPORTIONAL, fwm__excess_noise_slope=0.02)
    assert float(app_test.metric[0].value) == pytest.approx(0.04)

    rerun(fwm.NOISE_CONSTANT)
    assert_restored_slider("fwm__excess_noise", 2.0)
    assert app_test.session_state["fwm__excess_noise_slope"] == 0.02
    rerun(fwm.NOISE_CONSTANT, fwm__excess_noise=3.0)
    assert float(app_test.metric[0].value) == 3.0

    rerun(fwm.NOISE_GAIN_PROPORTIONAL)
    assert_restored_slider("fwm__excess_noise_slope", 0.02)
    assert app_test.session_state["fwm__excess_noise"] == 3.0
    assert float(app_test.metric[0].value) == pytest.approx(0.04)
