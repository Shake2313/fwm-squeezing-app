"""Independent analytic RF crossings and trustworthy missing-band semantics."""

from dataclasses import replace

import numpy as np
import pytest

from gabes.quantum.contracts import AnalysisFrequencyAxis, OpticalDetunings
from gabes.quantum.readout import IntensityDifferenceSpectrum
from gabes.quantum.spectrum_analysis import SpectralTrace, analyze_squeezing_bands


def _trace(frequency, ratio, **kwargs):
    return SpectralTrace(AnalysisFrequencyAxis.from_hz(frequency), ratio,
                         layer="source", provenance="independent analytic fixture", **kwargs)


def _valley(frequency):
    return .25 + ((np.asarray(frequency)-2.1e6)/1.3e6)**2


def test_known_sql_crossings_converge_then_local_brackets_bound_error():
    expected = 2.1e6 + np.array([-1., 1.])*1.3e6*np.sqrt(.75)
    errors = []
    for count in (17, 65, 257):
        frequency = np.linspace(.1e6, 4.e6, count)
        result = analyze_squeezing_bands(_trace(frequency, _valley(frequency)))
        assert result.status == "bands_found" and result.complete_domain
        assert len(result.bands) == 1
        band = result.bands[0]
        obtained = np.array([band.lower.frequency_hz, band.upper.frequency_hz])
        errors.append(np.max(abs(obtained-expected)))
        assert band.width_hz == band.observed_span_hz
        assert result.rf_refinement_status == "sampled_only"
    assert errors[1] < errors[0] and errors[2] < errors[1]

    frequency = np.linspace(.1e6, 4.e6, 17)
    refined = analyze_squeezing_bands(_trace(frequency, _valley(frequency)),
                                      crossing_evaluator=_valley,
                                      crossing_tolerance_hz=.05)
    assert refined.rf_refinement_status == "refined"
    assert refined.crossing_tolerance_hz == .05
    for edge, reference in zip((refined.bands[0].lower, refined.bands[0].upper), expected):
        assert edge.refinement_status == "refined"
        assert edge.bracket_width_hz <= .05
        assert edge.bracket_hz[0] <= reference <= edge.bracket_hz[1]
        assert abs(edge.frequency_hz-reference) <= .025
    # Refining crossings does not silently refine the minimum or whole RF grid.
    assert refined.minimum_ratio == float(min(_valley(frequency)))
    assert refined.bands[0].width_hz == pytest.approx(np.diff(expected)[0], abs=.05)


def test_target_db_uses_linear_ratio_threshold_and_interpolation():
    target = 10*np.log10(.5)
    frequency = np.linspace(.1e6, 4.e6, 20)
    result = analyze_squeezing_bands(_trace(frequency, _valley(frequency)),
                                     threshold_db=target,
                                     crossing_evaluator=_valley,
                                     crossing_tolerance_hz=.1)
    assert result.threshold_ratio == pytest.approx(.5)
    assert result.bands[0].lower.frequency_hz == pytest.approx(1.45e6, abs=.05)
    assert result.bands[0].upper.frequency_hz == pytest.approx(2.75e6, abs=.05)
    # R moves linearly from .25 to 1.25: SQL crossing lies 75% across this cell.
    linear = analyze_squeezing_bands(_trace([1., 2., 3.], [1.25, .25, 1.25]))
    assert linear.bands[0].lower.frequency_hz == pytest.approx(1.25)
    assert linear.bands[0].upper.frequency_hz == pytest.approx(2.75)


def test_disjoint_bands_never_bridge_a_technical_noise_peak_or_equality():
    peak = analyze_squeezing_bands(_trace([1., 2., 3., 4., 5.], [1.2, .5, 1.4, .4, 1.2]))
    assert len(peak.bands) == 2
    assert peak.bands[0].upper.frequency_hz < 3. < peak.bands[1].lower.frequency_hz
    assert all(band.width_hz is not None for band in peak.bands)

    plateau = analyze_squeezing_bands(_trace([1., 2., 3., 4., 5., 6.], [1.2, .5, 1., 1., .5, 1.2]))
    assert len(plateau.bands) == 2
    assert plateau.bands[0].upper.frequency_hz == pytest.approx(3.)
    assert plateau.bands[1].lower.frequency_hz == pytest.approx(4.)
    assert plateau.bands[0].upper.refinement_status == "exact_sample"

    isolated_equality = analyze_squeezing_bands(_trace([1., 2., 3.], [.5, 1., .5]))
    assert len(isolated_equality.bands) == 2
    assert isolated_equality.bands[0].upper.frequency_hz == isolated_equality.bands[1].lower.frequency_hz
    assert all(band.width_hz is None for band in isolated_equality.bands)


def test_flat_squeezing_and_scan_edge_crossings_remain_domain_limited():
    result = analyze_squeezing_bands(_trace([1e5, 1e6, 4e6], [.2, .2, .2]))
    band = result.bands[0]
    assert band.width_hz is None
    assert band.lower.kind == band.upper.kind == "domain"
    assert band.observed_span_hz == pytest.approx(3.9e6)
    assert result.rf_refinement_status == "no_crossings"
    extended = analyze_squeezing_bands(_trace([1e4, 1e5, 1e6, 4e6, 8e6], [.2]*5))
    assert extended.bands[0].width_hz is None
    assert extended.bands[0].observed_span_hz > band.observed_span_hz

    clipped = analyze_squeezing_bands(_trace([1., 2., 3.], [.2, .4, 1.4]))
    assert clipped.bands[0].lower.kind == "domain"
    assert clipped.bands[0].upper.kind == "crossing"
    assert clipped.bands[0].width_hz is None


@pytest.mark.parametrize("hole", [np.nan, np.inf, -np.inf])
def test_missing_samples_split_bands_and_do_not_invent_crossings(hole):
    result = analyze_squeezing_bands(_trace([1., 2., 3., 4., 5.], [1.2, .5, hole, .4, 1.2]))
    assert len(result.bands) == 2 and not result.complete_domain
    assert result.invalid_sample_count == 1
    assert result.bands[0].upper.kind == result.bands[1].lower.kind == "invalid"
    assert result.bands[0].upper.frequency_hz == pytest.approx(2.)
    assert result.bands[1].lower.frequency_hz == pytest.approx(4.)
    assert all(band.width_hz is None for band in result.bands)


def test_invalid_and_unconverged_masks_keep_separate_unavailability_counts():
    trace = _trace([1., 2., 3., 4., 5.], [1.2, .5, .3, .4, 1.2],
                   converged=[True, True, False, True, True])
    result = analyze_squeezing_bands(trace)
    assert result.invalid_sample_count == 0 and result.unconverged_sample_count == 1
    assert result.bands[0].upper.kind == result.bands[1].lower.kind == "unconverged"
    changed = analyze_squeezing_bands(replace(trace, valid=[True, True, False, True, True]))
    assert changed.invalid_sample_count == 1 and changed.unconverged_sample_count == 0
    assert changed.bands[0].upper.kind == "invalid"


def test_no_squeezing_is_distinct_from_partial_or_entirely_unavailable_data():
    result = analyze_squeezing_bands(_trace([1., 2., 3.], [1., 1.1, 1.]))
    assert result.status == "no_squeezing" and result.bands == ()
    assert result.minimum_ratio == 1.
    # Missing a stronger target does not mean there is no squeezing at SQL.
    target = analyze_squeezing_bands(_trace([1., 2., 3.], [.7]*3), threshold_db=-3.)
    assert target.status == "no_target_band" and target.minimum_ratio < 1.
    result = analyze_squeezing_bands(_trace([1., 2., 3.], [1., np.nan, 1.1]))
    assert result.status == "incomplete" and not result.complete_domain
    result = analyze_squeezing_bands(_trace([1., 2., 3.], [np.nan]*3))
    assert result.status == "unavailable" and result.minimum_ratio is None
    result = analyze_squeezing_bands(_trace([1., 2., 3.], [.4]*3,
                                            converged=[False]*3))
    assert result.status == "unavailable" and result.unconverged_sample_count == 3
    # A one-point RF observation cannot establish any intrinsic bandwidth.
    result = analyze_squeezing_bands(_trace([1e6], [0.]))
    assert result.bands[0].width_hz is None and result.bands[0].observed_span_hz == 0.


def test_common_electronic_lowpass_preserves_ratio_and_does_not_set_source_bandwidth():
    frequency = np.linspace(.1e6, 4e6, 67)
    axis = AnalysisFrequencyAxis.from_hz(frequency)
    # Shared |H(f)|^2 multiplies numerator and SQL at the same current plane.
    transfer_power = 1/(1+(frequency/.4e6)**2)
    sql = 2.7e-24*(1+frequency/4e6)
    source = SpectralTrace.from_psds(axis, _valley(frequency)*sql, sql,
                                     layer="source", provenance="analytic source PSD")
    detected = SpectralTrace.from_psds(axis, _valley(frequency)*sql*transfer_power,
                                       sql*transfer_power, layer="detected",
                                       provenance="shared electronic low-pass fixture")
    np.testing.assert_allclose(source.ratio, detected.ratio, rtol=3e-16)
    source_band = analyze_squeezing_bands(source).bands[0]
    detected_band = analyze_squeezing_bands(detected).bands[0]
    assert detected_band.width_hz == pytest.approx(source_band.width_hz, rel=1e-15)
    assert detected_band.lower.frequency_hz > .4e6
    assert source.layer != detected.layer and source.provenance != detected.provenance
    assert "no RBW/VBW convolution" in detected.measurement_filter_scope


def test_readout_adapter_preserves_absolute_units_and_explicit_electronics_choice():
    axis = AnalysisFrequencyAxis.from_hz([1e5, 1e6, 4e6])
    sql = np.array([1e-25, 2e-25, 4e-25])
    quantum_ratio = np.array([1.1, .5, 1.2])
    electronics = .6*sql
    spectrum = IntensityDifferenceSpectrum(axis, quantum_ratio*sql, sql, electronics,
                                            quantum_ratio, quantum_ratio+.6,
                                            [1e12, 2e12], 0., 0., "declared detector fixture")
    quantum = SpectralTrace.from_intensity_difference(spectrum, contribution="quantum")
    total = SpectralTrace.from_intensity_difference(spectrum)
    np.testing.assert_allclose(quantum.ratio, quantum_ratio)
    np.testing.assert_allclose(total.ratio, quantum_ratio+.6)
    assert quantum.layer == total.layer == "detected"
    assert "quantum intensity-difference PSD" in quantum.provenance
    assert "total intensity-difference PSD" in total.provenance
    assert analyze_squeezing_bands(quantum).status == "bands_found"
    assert analyze_squeezing_bands(total).status == "no_squeezing"
    # Positive total PSD cannot conceal an invalid underlying quantum PSD.
    corrupt = replace(spectrum, quantum_psd_A2_Hz=[-1e-26, 1e-25, 1e-25])
    with pytest.raises(ValueError, match="quantum PSD"):
        SpectralTrace.from_intensity_difference(corrupt)


def test_psd_contract_rejects_nonphysical_units_boundary_and_retains_missing_values():
    axis = AnalysisFrequencyAxis.from_hz([1., 2., 3.])
    kwargs = dict(layer="detected", provenance="unit-declared fixture")
    for sql in ([1e-25, 0., 1e-25], [1e-25, -1e-25, 1e-25]):
        with pytest.raises(ValueError, match="SQL"):
            SpectralTrace.from_psds(axis, [1e-25]*3, sql, **kwargs)
    with pytest.raises(ValueError, match="nonnegative"):
        SpectralTrace.from_psds(axis, [1e-25, -1e-28, 1e-25], [1e-25]*3, **kwargs)
    with pytest.raises(ValueError, match="real"):
        SpectralTrace.from_psds(axis, [1e-25+1e-30j]*3, [1e-25]*3, **kwargs)
    missing = SpectralTrace.from_psds(axis, [1e-25, np.nan, 1e-25],
                                      [2e-25, 2e-25, np.inf], **kwargs)
    assert analyze_squeezing_bands(missing).invalid_sample_count == 2


def test_trace_owns_immutable_arrays_and_rejects_wrong_axes_or_status_inputs():
    values = np.array([1.2, .5, 1.2])
    trace = _trace([1., 2., 3.], values)
    values[:] = 42
    assert trace.ratio[1] == .5
    with pytest.raises(ValueError):
        trace.ratio.setflags(write=True)
    with pytest.raises(ValueError):
        trace.valid.setflags(write=True)
    for kwargs in (dict(valid=[1, 0, 1]), dict(converged=[True]),
                   dict(measurement_filter_scope="")):
        with pytest.raises(ValueError):
            _trace([1., 2., 3.], [1.2, .5, 1.2], **kwargs)
    with pytest.raises(ValueError, match="positive"):
        _trace([0., 1., 2.], [1.2, .5, 1.2])
    with pytest.raises(TypeError, match="AnalysisFrequencyAxis"):
        SpectralTrace(OpticalDetunings(1e9, 1e6), values, "source", "wrong axis")
    with pytest.raises(ValueError, match="negative"):
        _trace([1., 2., 3.], [1.2, -.5, 1.2])


def test_refinement_failure_and_iteration_limit_do_not_claim_numerical_convergence():
    frequency = np.linspace(.1e6, 4e6, 17)
    trace = _trace(frequency, _valley(frequency))
    result = analyze_squeezing_bands(trace, crossing_evaluator=_valley,
                                      crossing_tolerance_hz=.01,
                                      max_refinement_iterations=1)
    assert result.rf_refinement_status == "unresolved"
    assert all(edge.refinement_status == "unresolved"
               for edge in (result.bands[0].lower, result.bands[0].upper))
    assert result.bands[0].width_hz is not None  # existence differs from precision
    with pytest.raises(ValueError, match="unavailable"):
        analyze_squeezing_bands(trace, crossing_evaluator=lambda f: np.nan,
                                crossing_tolerance_hz=.01)
    with pytest.raises(ValueError, match="tolerance"):
        analyze_squeezing_bands(trace, crossing_evaluator=_valley)
    with pytest.raises(ValueError, match="provider"):
        analyze_squeezing_bands(trace, crossing_tolerance_hz=.01)
    with pytest.raises(ValueError, match="threshold"):
        analyze_squeezing_bands(trace, threshold_db=1.)


def test_refinement_rejects_a_different_provider_even_at_exact_sample_crossings():
    frequency = np.linspace(.1e6, 4e6, 17)
    trace = _trace(frequency, _valley(frequency))
    with pytest.raises(ValueError, match="trace endpoints"):
        analyze_squeezing_bands(trace, crossing_evaluator=lambda f: 2.,
                                crossing_tolerance_hz=.1)
    exact = _trace([1., 2., 3.], [1., .5, 1.])
    with pytest.raises(ValueError, match="trace endpoints"):
        analyze_squeezing_bands(exact, crossing_evaluator=lambda f: 2.,
                                crossing_tolerance_hz=.1)
    # Agreement to absolute precision alone cannot authorize a changed side.
    with pytest.raises(ValueError, match="threshold bracket"):
        analyze_squeezing_bands(exact, crossing_evaluator=lambda f: 1.+1e-15,
                                crossing_tolerance_hz=.1)


def test_extreme_finite_ratios_keep_interpolated_edges_inside_brackets():
    axis = AnalysisFrequencyAxis.from_hz([1e5, 1e6, 1e7])
    trace = SpectralTrace(axis, [1e308, .5, 1e308], layer="source",
                          provenance="extreme dynamic range fixture")
    result = analyze_squeezing_bands(trace)
    band = result.bands[0]
    assert np.isfinite(band.width_hz) and band.width_hz >= 0
    for edge in (band.lower, band.upper):
        assert edge.bracket_hz[0] <= edge.frequency_hz <= edge.bracket_hz[1]
