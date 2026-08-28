"""CSV import/correction checks for the OD/SAS experimental overlay."""

from pathlib import Path

import numpy as np
import pytest

from gabes import experimental_csv as ecsv


ROOT = Path(__file__).resolve().parent.parent
REFERENCE_OD = ROOT / "references" / "AutoOD" / "ReferenceOD.csv"


def _csv_bytes(x, y, *, header=None):
    lines = [] if header is None else list(header)
    lines.extend(f"{a:.12g},{b:.12g}" for a, b in zip(x, y))
    return ("\n".join(lines) + "\n").encode("utf-8")


def _roughness(y):
    return float(np.std(np.diff(np.asarray(y), n=2)))


def test_parser_uses_only_a_b_and_skips_scope_metadata():
    payload = (
        "Oscilloscope export,,operator note\n"
        "Detuning,Channel 1,Channel 2\n"
        "3,30,999\n"
        "1,10,888\n"
        "2,20,777\n"
        "2,22,666\n"
        "4,40,555\n"
        "5,50,444\n"
        "NaN,12,333\n"
        "6,Inf,222\n"
        ",,\n"
    ).encode("utf-8-sig")

    trace = ecsv.load_experimental_csv(payload, denoise=False)

    np.testing.assert_allclose(trace.detuning, [1, 2, 3, 4, 5])
    np.testing.assert_allclose(trace.raw_signal, [10, 21, 30, 40, 50])
    diag = trace.import_diagnostics
    assert diag.valid_rows == 6
    assert diag.unique_points == 5
    assert diag.duplicate_rows_merged == 1
    assert diag.nonnumeric_rows == 1
    assert diag.nonfinite_rows == 2
    assert diag.blank_rows == 2
    assert diag.extra_column_rows == diag.total_rows
    assert np.all(np.isfinite(trace.transmission))
    assert np.all((0.0 <= trace.transmission) & (trace.transmission <= 1.0))


def test_cp949_metadata_and_first_numeric_row_are_preserved():
    text = "장비 메타데이터,채널\n-2,20\n-1,40\n0,60\n1,80\n2,100\n"
    trace = ecsv.load_experimental_csv(text.encode("cp949"), denoise=False)
    assert trace.import_diagnostics.encoding == "cp949"
    assert trace.detuning[0] == -2.0
    assert trace.raw_signal[0] == 20.0


def test_unsupported_bytes_in_ignored_columns_do_not_reject_valid_a_b():
    payload = b"scope metadata,ignored,\x80\n0,0,x\n1,1,x\n2,2,x\n3,3,x\n4,4,\x80\n"

    trace = ecsv.load_experimental_csv(payload, denoise=False)

    assert trace.import_diagnostics.encoding == "latin-1"
    np.testing.assert_allclose(trace.detuning, [0, 1, 2, 3, 4])
    np.testing.assert_allclose(trace.raw_signal, [0, 1, 2, 3, 4])


def test_utf16_scope_export_is_supported():
    text = "scope header,channel B,ignored\n0,0,x\n1,1,x\n2,2,x\n3,3,x\n4,4,x\n"

    trace = ecsv.load_experimental_csv(text.encode("utf-16"), denoise=False)

    assert trace.import_diagnostics.encoding == "utf-16"
    np.testing.assert_allclose(trace.detuning, [0, 1, 2, 3, 4])


def test_mv_floor_and_ceiling_map_to_zero_and_one():
    x = np.linspace(-5.0, 5.0, 201)
    signal_mv = np.full_like(x, 180.0)
    signal_mv[70:131] = 20.0
    trace = ecsv.load_experimental_csv(
        _csv_bytes(x, signal_mv, header=("scope header,text",)),
        denoise=False,
    )

    assert trace.correction_diagnostics.floor == pytest.approx(20.0)
    assert trace.correction_diagnostics.ceiling == pytest.approx(180.0)
    assert np.median(trace.transmission[80:120]) <= 0.01
    assert np.median(np.r_[trace.transmission[:50], trace.transmission[-50:]]) >= 0.99
    assert trace.transmission_label == ecsv.RELATIVE_TRANSMISSION_LABEL
    assert not trace.is_absolute_transmission


def test_default_extrema_calibration_is_explicitly_relative():
    x = np.linspace(-1.0, 1.0, 101)
    physical_transmission = np.linspace(0.70, 0.92, x.size)

    trace = ecsv.load_experimental_csv(
        _csv_bytes(x, physical_transmission), denoise=False
    )

    assert trace.correction_diagnostics.calibration_mode == (
        ecsv.CALIBRATION_RELATIVE_EXTREMA
    )
    assert trace.transmission_label == "Relative normalized transmission"
    assert trace.transmission[0] == 0.0
    assert trace.transmission[-1] == 1.0
    assert not np.allclose(trace.transmission, physical_transmission)


def test_gain_offset_absolute_calibration_reconstructs_transmission():
    x = np.linspace(-2.0, 2.0, 101)
    physical_transmission = 0.81 + 0.11 * np.cos(x)
    gain = 2.75
    offset = -0.13
    signal = offset + gain * physical_transmission

    trace = ecsv.load_experimental_csv(
        _csv_bytes(x, signal),
        denoise=False,
        calibration_mode=ecsv.CALIBRATION_ABSOLUTE_GAIN_OFFSET,
        gain=gain,
        offset=offset,
    )

    np.testing.assert_allclose(
        trace.transmission, physical_transmission, rtol=1e-9, atol=1e-12
    )
    assert trace.is_absolute_transmission
    assert trace.transmission_label == "Absolute transmission"
    assert trace.correction_diagnostics.calibration_source == "explicit gain/offset"


def test_dark_reference_absolute_calibration_reconstructs_transmission():
    x = np.linspace(-1.0, 1.0, 81)
    physical_transmission = 0.76 + 0.08 * np.cos(2.0 * x)
    dark, reference = 1.4, -3.2
    signal = dark + (reference - dark) * physical_transmission

    trace = ecsv.load_experimental_csv(
        _csv_bytes(x, signal),
        denoise=False,
        calibration_mode=ecsv.CALIBRATION_ABSOLUTE_DARK_REFERENCE,
        dark_signal=dark,
        reference_signal=reference,
    )

    np.testing.assert_allclose(
        trace.transmission, physical_transmission, rtol=1e-9, atol=1e-12
    )
    assert trace.correction_diagnostics.calibration_gain == pytest.approx(
        reference - dark
    )


@pytest.mark.parametrize(
    "kwargs,match",
    [
        (
            dict(
                calibration_mode=ecsv.CALIBRATION_ABSOLUTE_DARK_REFERENCE,
                dark_signal=0.0,
            ),
            "both dark_signal and reference_signal",
        ),
        (
            dict(
                calibration_mode=ecsv.CALIBRATION_ABSOLUTE_DARK_REFERENCE,
                dark_signal=1.0,
                reference_signal=1.0,
            ),
            "non-zero",
        ),
        (
            dict(
                calibration_mode=ecsv.CALIBRATION_ABSOLUTE_GAIN_OFFSET,
                gain=0.0,
                offset=2.0,
            ),
            "non-zero",
        ),
        (
            dict(
                calibration_mode=ecsv.CALIBRATION_ABSOLUTE_GAIN_OFFSET,
                gain=float("nan"),
                offset=0.0,
            ),
            "finite",
        ),
        (dict(gain=2.0), "absolute calibration mode"),
    ],
)
def test_absolute_calibration_rejects_missing_or_invalid_evidence(kwargs, match):
    payload = b"0,0\n1,1\n2,2\n3,3\n4,4\n"
    with pytest.raises(ecsv.ExperimentalCSVError, match=match):
        ecsv.load_experimental_csv(payload, denoise=False, **kwargs)


def test_absolute_calibration_preserves_out_of_range_samples():
    x = np.arange(5.0)
    physical_transmission = np.array([-0.01, 0.2, 0.5, 0.8, 1.02])
    signal = 4.0 + 3.0 * physical_transmission

    trace = ecsv.load_experimental_csv(
        _csv_bytes(x, signal),
        denoise=False,
        calibration_mode=ecsv.CALIBRATION_ABSOLUTE_GAIN_OFFSET,
        gain=3.0,
        offset=4.0,
    )

    np.testing.assert_allclose(trace.transmission, physical_transmission, atol=1e-12)
    assert trace.correction_diagnostics.out_of_range_samples == 2
    assert "preserved rather than clipped" in trace.correction_diagnostics.warning


def test_round_trip_sweep_preserves_branches_and_exact_direction_counts():
    forward_x = np.arange(5.0)
    reverse_x = np.arange(4.0, -1.0, -1.0)
    x = np.concatenate((forward_x, reverse_x))
    forward_y = 10.0 + forward_x
    reverse_y = 30.0 + reverse_x
    y = np.concatenate((forward_y, reverse_y))

    trace = ecsv.load_experimental_csv(_csv_bytes(x, y), denoise=False)
    diag = trace.import_diagnostics

    assert diag.sweep_reversal_count == 1
    assert diag.forward_sample_count == 5
    assert diag.reverse_sample_count == 5
    assert diag.sweep_branch_count == 2
    assert diag.directional_merge_warning is not None
    assert [branch.direction for branch in trace.sweep_branches] == [
        "forward", "reverse"
    ]
    np.testing.assert_array_equal(trace.sweep_branches[0].detuning, forward_x)
    np.testing.assert_array_equal(trace.sweep_branches[1].detuning, reverse_x)
    np.testing.assert_array_equal(trace.sweep_branches[0].raw_signal, forward_y)
    np.testing.assert_array_equal(trace.sweep_branches[1].raw_signal, reverse_y)
    assert sum(branch.sample_count for branch in trace.sweep_branches) == 10
    assert trace.import_diagnostics.unique_points == 5


def test_equal_x_steps_do_not_create_false_sweep_reversals():
    x = np.array([0.0, 1.0, 1.0, 2.0, 3.0, 4.0])
    y = np.arange(x.size, dtype=float)
    trace = ecsv.load_experimental_csv(_csv_bytes(x, y), denoise=False)

    assert trace.import_diagnostics.sweep_reversal_count == 0
    assert trace.import_diagnostics.forward_sample_count == x.size
    assert [branch.direction for branch in trace.sweep_branches] == ["forward"]


def test_sort_duplicates_and_large_gaps_are_deterministic():
    ascending = "0,0\n1,1\n2,2\n2,4\n3,3\n4,4\n100,5\n101,6\n102,7\n103,8\n104,9\n"
    descending = "\n".join(reversed(ascending.strip().splitlines())) + "\n"
    first = ecsv.load_experimental_csv(ascending.encode(), denoise=False)
    second = ecsv.load_experimental_csv(descending.encode(), denoise=False)

    np.testing.assert_allclose(first.detuning, second.detuning)
    np.testing.assert_allclose(first.raw_signal, second.raw_signal)
    assert first.raw_signal[2] == 3.0
    assert first.import_diagnostics.segment_count == 2


def test_denoising_improves_noisy_absorption_without_moving_its_centre():
    rng = np.random.default_rng(2313)
    base_x = np.linspace(-2.0, 2.0, 801)
    x = np.sort(base_x + rng.normal(0.0, 0.00035, base_x.size))
    true_signal = 180.0 - 145.0 * np.exp(-0.5 * ((x - 0.23) / 0.34) ** 2)
    measured = true_signal + rng.normal(0.0, 4.0, x.size)
    spike_indices = rng.choice(x.size, size=18, replace=False)
    measured[spike_indices] += rng.choice((-75.0, 75.0), size=spike_indices.size)
    payload = _csv_bytes(x, measured)

    raw = ecsv.load_experimental_csv(payload, denoise=False)
    corrected = ecsv.load_experimental_csv(payload, denoise=True)
    target = np.clip((true_signal - 35.0) / 145.0, 0.0, 1.0)
    raw_rmse = float(np.sqrt(np.mean((raw.transmission - target) ** 2)))
    corrected_rmse = float(np.sqrt(np.mean((corrected.transmission - target) ** 2)))

    assert corrected_rmse <= 0.75 * raw_rmse
    assert _roughness(corrected.transmission) < 0.45 * _roughness(raw.transmission)
    low_tail = corrected.denoised_signal <= np.quantile(
        corrected.denoised_signal, 0.10
    )
    corrected_centre = float(np.mean(corrected.detuning[low_tail]))
    assert abs(corrected_centre - 0.23) <= 2.0 * np.median(np.diff(x))
    assert corrected.correction_diagnostics.hampel_replacements > 0
    assert corrected.correction_diagnostics.smoothing_window >= 5


@pytest.mark.parametrize("width", range(1, 6))
def test_automatic_correction_preserves_one_to_five_sample_lamb_dips(width):
    x = np.linspace(-5.0, 5.0, 1001)
    signal_mv = np.full_like(x, 180.0)
    start = x.size // 2 - width // 2
    signal_mv[start:start + width] = 20.0

    trace = ecsv.load_experimental_csv(_csv_bytes(x, signal_mv), denoise=True)

    assert trace.correction_diagnostics.floor == pytest.approx(20.0)
    assert trace.correction_diagnostics.ceiling == pytest.approx(180.0)
    assert np.min(trace.transmission[start:start + width]) <= 0.01
    assert np.median(trace.transmission[:400]) >= 0.99


def test_sub_five_percent_gaussian_feature_sets_the_calibration_floor():
    x = np.linspace(-1.0, 1.0, 1001)
    signal_mv = 180.0 - 160.0 * np.exp(-0.5 * (x / 0.008) ** 2)

    trace = ecsv.load_experimental_csv(_csv_bytes(x, signal_mv), denoise=True)

    assert trace.correction_diagnostics.floor < 25.0
    assert trace.correction_diagnostics.ceiling > 175.0
    assert abs(trace.detuning[np.argmin(trace.transmission)]) <= 0.002
    assert np.min(trace.transmission) <= 0.01


@pytest.mark.parametrize("width", range(1, 6))
def test_noisy_narrow_dip_keeps_a_unit_transmission_baseline(width):
    rng = np.random.default_rng(8000 + width)
    x = np.linspace(-1.0, 1.0, 1001)
    signal_mv = 100.0 + rng.normal(0.0, 0.1, x.size)
    start = x.size // 2 - width // 2
    signal_mv[start:start + width] -= 80.0

    trace = ecsv.load_experimental_csv(_csv_bytes(x, signal_mv), denoise=True)

    baseline = np.r_[trace.transmission[:400], trace.transmission[-400:]]
    assert np.median(baseline) >= 0.98
    assert np.min(trace.transmission[start:start + width]) <= 0.02
    assert trace.correction_diagnostics.sparse_extrema
    assert trace.correction_diagnostics.warning is not None


def test_span_defining_singleton_survives_beside_a_broad_dip():
    x = np.linspace(-2.0, 2.0, 1001)
    signal_mv = 180.0 - 80.0 * np.exp(-0.5 * (x / 0.35) ** 2)
    singleton = 180
    signal_mv[singleton] = 20.0

    trace = ecsv.load_experimental_csv(_csv_bytes(x, signal_mv), denoise=True)

    assert trace.transmission[singleton] <= 0.02
    assert trace.correction_diagnostics.sparse_extrema


@pytest.mark.parametrize(
    "payload,match",
    [
        (b"heading,text\n0,1\n1,2\n2,3\n3,4\n", "at least 5"),
        (b"heading,text\nmetadata,only\n", "at least 5"),
        (b"0,7\n1,7\n2,7\n3,7\n4,7\n", "no usable"),
    ],
)
def test_invalid_or_low_information_csvs_fail_cleanly(payload, match):
    with pytest.raises(ecsv.ExperimentalCSVError, match=match):
        ecsv.load_experimental_csv(payload)


def test_file_and_row_limits(monkeypatch):
    with pytest.raises(ecsv.ExperimentalCSVError, match="larger"):
        ecsv.load_experimental_csv(b"x" * (ecsv.MAX_FILE_BYTES + 1))

    monkeypatch.setattr(ecsv, "MAX_CSV_ROWS", 4)
    with pytest.raises(ecsv.ExperimentalCSVError, match="row import limit"):
        ecsv.load_experimental_csv(b"0,0\n1,1\n2,2\n3,3\n4,4\n")


def test_manual_x_transform_and_reset_math():
    x = np.array([0.0, 2.0, 4.0])
    np.testing.assert_allclose(ecsv.transform_x(x), x)
    np.testing.assert_allclose(
        ecsv.transform_x(x, scale=2.0, shift=1.0), [-1.0, 3.0, 7.0]
    )
    np.testing.assert_allclose(
        ecsv.transform_x(x, scale=2.0, shift=1.0, reverse=True), [7.0, 3.0, -1.0]
    )
    with pytest.raises(ValueError, match="positive"):
        ecsv.transform_x(x, scale=0.0)


def test_manual_x_transform_uses_an_overflow_safe_midpoint():
    x = np.array([1.0e308, 1.1e308])
    np.testing.assert_allclose(ecsv.transform_x(x), x)
    with pytest.raises(ValueError, match="remain finite"):
        ecsv.transform_x(x, scale=1.0e9)


@pytest.mark.skipif(not REFERENCE_OD.exists(), reason="local AutoOD reference absent")
def test_reference_od_full_file_is_imported_and_smoothed():
    payload = REFERENCE_OD.read_bytes()
    raw = ecsv.load_experimental_csv(payload, denoise=False)
    corrected = ecsv.load_experimental_csv(payload, denoise=True)

    assert raw.import_diagnostics.valid_rows == 7678
    assert raw.import_diagnostics.ignored_rows == 0
    assert raw.import_diagnostics.unique_points == 2341
    assert raw.import_diagnostics.duplicate_rows_merged == 5337
    assert raw.detuning[0] == pytest.approx(-7.93)
    assert raw.detuning[-1] == pytest.approx(7.53)
    assert np.all(np.isfinite(corrected.transmission))
    assert np.all((0.0 <= corrected.transmission) & (corrected.transmission <= 1.0))
    assert _roughness(corrected.transmission) < 0.50 * _roughness(raw.transmission)

    # Flat-bottom minima can move by one quantized sample, so compare the centre
    # of each low-signal basin rather than a single unstable argmin sample.
    median_dx = float(np.median(np.diff(raw.detuning)))
    for low, high in ((-6.5, -4.0), (-3.5, -0.5), (-0.5, 2.5)):
        mask = (raw.detuning >= low) & (raw.detuning <= high)
        raw_y = raw.transmission[mask]
        corrected_y = corrected.transmission[mask]
        x_window = raw.detuning[mask]
        raw_centre = float(np.mean(x_window[raw_y <= np.quantile(raw_y, 0.20)]))
        corrected_centre = float(
            np.mean(x_window[corrected_y <= np.quantile(corrected_y, 0.20)])
        )
        assert abs(corrected_centre - raw_centre) <= 2.5 * median_dx
