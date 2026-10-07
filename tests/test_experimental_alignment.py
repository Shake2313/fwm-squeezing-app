"""Affine visual alignment, scale reachability, and calibration safeguards."""
from dataclasses import FrozenInstanceError
import io
import warnings

import numpy as np
import pytest

from gabes.experimental_alignment import (
    AxisAlignment,
    estimate_axis_scale_correction,
    fit_trace_alignment,
    frame_alignment,
)
from gabes.experimental_csv import (
    CALIBRATION_ABSOLUTE_GAIN_OFFSET,
    load_experimental_csv,
)


def _spectrum(x):
    return (1.0 - 0.68 * np.exp(-((x + 3.1) / 0.33)**2)
            - 0.36 * np.exp(-((x + 1.35) / 0.5)**2)
            - 0.82 * np.exp(-((x - 1.65) / 0.26)**2)
            - 0.5 * np.exp(-((x - 3.05) / 0.4)**2))


def _trace(x, y, *, absolute=False):
    buffer = io.StringIO()
    np.savetxt(buffer, np.column_stack((x, 300.0 * y + 12.0)), delimiter=",")
    kwargs = dict(denoise=False)
    if absolute:
        kwargs.update(calibration_mode=CALIBRATION_ABSOLUTE_GAIN_OFFSET,
                      gain=300.0, offset=12.0)
    return load_experimental_csv(buffer.getvalue().encode(), **kwargs)


def test_frame_uses_endpoint_midpoint_for_irregular_samples():
    x = np.array([10.0, 10.001, 10.002, 12.0])
    alignment = frame_alignment(x, [-6.0, 6.0])
    mapped = alignment.transform_x(x)
    assert mapped[0] == pytest.approx(-5.4)
    assert mapped[-1] == pytest.approx(5.4)
    assert alignment.shift == -11.0
    with pytest.raises(FrozenInstanceError):
        alignment.scale = 1.0


@pytest.mark.parametrize("width", [1e-15, 0.02, 2e15])
def test_framing_has_no_fixed_scale_floor_or_ceiling(width):
    x = np.array([-width / 2.0, width / 2.0])
    result = frame_alignment(x, [-5.0, 5.0])
    np.testing.assert_allclose(result.transform_x(x), [-4.5, 4.5], rtol=1e-14)
    assert result.scale == pytest.approx(9.0 / width)


@pytest.mark.parametrize("raw_span", [1e-14, .2, 10.0, 1e14])
def test_decade_preconditioning_keeps_manual_multiplier_reachable(raw_span):
    result = estimate_axis_scale_correction([0.0, raw_span], [-5.0, 5.0])
    assert result.valid
    assert .1 <= result.target_manual_scale <= 3.0
    assert result.multiplier * result.target_manual_scale == pytest.approx(10.0 / raw_span)
    assert np.log10(result.multiplier) == pytest.approx(round(np.log10(result.multiplier)))


def test_decade_correction_is_span_only_and_stays_separate_from_manual_scale():
    correction = estimate_axis_scale_correction([100.0, 100.2], [-5.0, 5.0])
    assert correction.multiplier == 100.0
    alignment = AxisAlignment(correction.multiplier * 0.5, -100.1)
    first = alignment.transform_x(np.array([100.0, 100.2]))
    second = alignment.transform_x(np.array([100.0, 100.2]))
    np.testing.assert_array_equal(first, second)
    np.testing.assert_allclose(first, [-5.0, 5.0], atol=1e-12)


@pytest.mark.parametrize("x", [[1.0, 1.0], [0.0], [0.0, np.nan]])
def test_invalid_span_correction_fails_without_applying_a_multiplier(x):
    correction = estimate_axis_scale_correction(x, [-5.0, 5.0])
    assert not correction.valid
    assert correction.multiplier == 1.0
    assert correction.reason
    with pytest.raises(ValueError):
        frame_alignment(x, [-5.0, 5.0])


@pytest.mark.parametrize("reverse", [False, True])
def test_recovers_multipeak_seconds_axis_and_hundredfold_detector_signal(reverse):
    theory_x = np.linspace(-6.0, 6.0, 2401)
    raw_x = np.linspace(-.012, .012, 1201)
    sign = -1.0 if reverse else 1.0
    trace = _trace(raw_x, _spectrum(sign * 400.0 * raw_x + .7))
    initial = frame_alignment(raw_x, theory_x, reverse=reverse)
    fitted = fit_trace_alignment(trace, theory_x, _spectrum(theory_x), initial=initial)
    assert fitted.success and fitted.improved
    assert fitted.alignment.scale == pytest.approx(400.0, rel=.003)
    assert fitted.alignment.shift == pytest.approx(.7, abs=.01)
    assert fitted.alignment.reverse == reverse
    assert fitted.rmse < .003
    assert fitted.rmse < fitted.initial_rmse
    assert fitted.overlap_fraction > .99
    corrected = fitted.alignment.transform_y(trace.transmission)
    # The existing relative-extrema importer clips its robust floor/ceiling.
    # Affine alignment cannot reconstruct those clipped extrema.
    unclipped = (trace.transmission > 0.0) & (trace.transmission < 1.0)
    np.testing.assert_allclose(corrected[unclipped],
                              _spectrum(sign * 400.0 * raw_x + .7)[unclipped], atol=.003)


def test_absolute_transmission_is_preserved_even_when_signal_fit_requested():
    theory_x = np.linspace(-6.0, 6.0, 2401)
    raw_x = np.linspace(-.012, .012, 1201)
    trace = _trace(raw_x, _spectrum(400.0 * raw_x + .7), absolute=True)
    initial = frame_alignment(raw_x, theory_x)
    fitted = fit_trace_alignment(trace, theory_x, _spectrum(theory_x),
                                 initial=initial, fit_signal=True)
    assert fitted.success
    assert fitted.alignment.scale == pytest.approx(400.0, rel=.003)
    assert fitted.alignment.signal_scale == 1.0
    assert fitted.alignment.signal_offset == 0.0
    np.testing.assert_array_equal(fitted.alignment.transform_y(trace.transmission),
                                  trace.transmission)
    with pytest.raises(ValueError, match="cannot be inverted"):
        fit_trace_alignment(trace, theory_x, _spectrum(theory_x), invert=True)


def test_signal_fit_can_be_disabled_for_relative_data():
    theory_x = np.linspace(-6.0, 6.0, 1401)
    raw_x = np.linspace(-.012, .012, 601)
    trace = _trace(raw_x, _spectrum(400.0 * raw_x + .7))
    initial = AxisAlignment(400.0, .7, signal_scale=.82, signal_offset=.18)
    fitted = fit_trace_alignment(trace, theory_x, _spectrum(theory_x),
                                 initial=initial, fit_signal=False)
    assert fitted.alignment.signal_scale == .82
    assert fitted.alignment.signal_offset == .18


def test_fit_preserves_direction_and_does_not_mutate_arrays():
    theory_x = np.linspace(-6.0, 6.0, 1401)
    raw_x = np.linspace(-.012, .012, 601)
    measured = _spectrum(400.0 * raw_x + .7)
    trace = _trace(raw_x, measured)
    original_x = trace.detuning.copy()
    original_y = trace.transmission.copy()
    fitted = fit_trace_alignment(trace, theory_x, _spectrum(theory_x),
                                 initial=AxisAlignment(400.0, .7))
    assert fitted.success
    assert not fitted.alignment.reverse
    np.testing.assert_array_equal(trace.detuning, original_x)
    np.testing.assert_array_equal(trace.transmission, original_y)


def test_fit_recovers_when_initial_alignment_has_no_overlap():
    theory_x = np.linspace(-6.0, 6.0, 1401)
    raw_x = np.linspace(10.0 - .012, 10.0 + .012, 601)
    trace = _trace(raw_x, _spectrum(400.0 * (raw_x - 10.0) + .7))
    result = fit_trace_alignment(trace, theory_x, _spectrum(theory_x),
                                 initial=AxisAlignment(1.0, 1000.0))
    assert result.success and result.improved
    assert result.initial_rmse == np.inf
    assert result.alignment.scale == pytest.approx(400.0, rel=.003)
    assert result.alignment.shift == pytest.approx(-9.3, abs=.01)
    assert result.overlap_fraction > .99


def test_fit_handles_relative_detector_polarity_only_when_explicit():
    theory_x = np.linspace(-6.0, 6.0, 1401)
    raw_x = np.linspace(-.012, .012, 601)
    trace = _trace(raw_x, 1.0 - _spectrum(400.0 * raw_x + .7))
    result = fit_trace_alignment(trace, theory_x, _spectrum(theory_x), invert=True)
    assert result.success
    assert result.alignment.scale == pytest.approx(400.0, rel=.003)
    assert result.alignment.shift == pytest.approx(.7, abs=.01)


def test_isolated_signal_spikes_do_not_move_the_multipeak_fit():
    theory_x = np.linspace(-6.0, 6.0, 1401)
    raw_x = np.linspace(-.012, .012, 601)
    measured = _spectrum(400.0 * raw_x + .7)
    measured[[27, 105, 260, 375, 580]] += .6
    trace = _trace(raw_x, measured)
    result = fit_trace_alignment(trace, theory_x, _spectrum(theory_x))
    assert result.success
    assert result.alignment.scale == pytest.approx(400.0, rel=.005)
    assert result.alignment.shift == pytest.approx(.7, abs=.02)
    assert result.overlap_fraction >= .65


@pytest.mark.parametrize("seed", [1, 2, 3, 4, 5])
def test_unrelated_random_signal_cannot_claim_success_by_y_scale_collapse(seed):
    theory_x = np.linspace(-6.0, 6.0, 1401)
    raw_x = np.linspace(-.012, .012, 601)
    unrelated = np.random.default_rng(seed).normal(size=raw_x.size)
    result = fit_trace_alignment(_trace(raw_x, unrelated), theory_x, _spectrum(theory_x))
    assert not result.success
    assert not result.improved


@pytest.mark.parametrize("kind", ["ramp", "sine"])
def test_unrelated_smooth_signals_do_not_claim_a_spectral_match(kind):
    theory_x = np.linspace(-6.0, 6.0, 1401)
    raw_x = np.linspace(-.012, .012, 601)
    position = np.linspace(-1.0, 1.0, raw_x.size)
    unrelated = position if kind == "ramp" else np.sin(7.2 * position + .2)
    result = fit_trace_alignment(_trace(raw_x, unrelated), theory_x, _spectrum(theory_x))
    assert not result.success
    assert not result.improved


@pytest.mark.parametrize("signal_scale,signal_offset", [(1e308, 1e308), (1e307, 0.0)])
def test_overflowing_initial_signal_transform_is_rejected_without_runtime_warnings(
        signal_scale, signal_offset):
    theory_x = np.linspace(-6.0, 6.0, 1401)
    raw_x = np.linspace(-.012, .012, 601)
    trace = _trace(raw_x, _spectrum(400.0 * raw_x + .7))
    initial = AxisAlignment(400.0, .7, signal_scale=signal_scale, signal_offset=signal_offset)
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        result = fit_trace_alignment(trace, theory_x, _spectrum(theory_x), initial=initial)
    assert result.success
    assert result.initial_rmse > 1e100
    assert result.alignment.scale == pytest.approx(400.0, rel=.003)


def test_theory_requires_distinct_detuning_and_matching_arrays():
    theory_x = np.linspace(-6.0, 6.0, 101)
    trace = _trace(theory_x, _spectrum(theory_x))
    repeated = theory_x.copy()
    repeated[20] = repeated[19]
    with pytest.raises(ValueError, match="distinct"):
        fit_trace_alignment(trace, repeated, _spectrum(theory_x))
    with pytest.raises(ValueError, match="same-length"):
        fit_trace_alignment(trace, theory_x, _spectrum(theory_x)[:-1])


def test_flat_and_short_spectra_do_not_claim_alignment():
    theory_x = np.linspace(-6.0, 6.0, 101)
    flat_trace = _trace(np.linspace(0.0, 1.0, 101), np.ones(101), absolute=True)
    with pytest.raises(ValueError, match="flat imported"):
        fit_trace_alignment(flat_trace, theory_x, _spectrum(theory_x))
    short_trace = _trace(np.linspace(0.0, 1.0, 8), np.linspace(0.0, 1.0, 8))
    with pytest.raises(ValueError, match="at least 12"):
        fit_trace_alignment(short_trace, theory_x, _spectrum(theory_x))
    trace = _trace(theory_x, _spectrum(theory_x))
    with pytest.raises(ValueError, match="flat theoretical"):
        fit_trace_alignment(trace, theory_x, np.ones_like(theory_x))


@pytest.mark.parametrize("kwargs", [dict(scale=0.0), dict(scale=np.inf),
                                    dict(shift=np.nan), dict(signal_scale=-1.0)])
def test_invalid_affine_parameters_are_rejected(kwargs):
    with pytest.raises(ValueError):
        AxisAlignment(**kwargs)
