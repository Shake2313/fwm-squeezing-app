import hashlib

import numpy as np
import pytest

from gabes.rydberg_experimental_csv import (
    RydbergCSVError,
    load_eit_csv,
    load_psd_csv,
    load_rf_sweep_csv,
)


def test_eit_loader_finds_header_merges_duplicates_and_records_provenance():
    payload = (
        "Instrument,scope-1\n"
        "Operator,test\n"
        "Probe detuning [kHz],Photodiode voltage [V],ignored\n"
        "1000,1.0,a\n"
        "0,2.0,b\n"
        "0,4.0,c\n"
        "-1000,1.5,d\n"
        "bad,row,e\n"
        "2000,nan,f\n"
    ).encode("utf-8")
    trace = load_eit_csv(payload, source_name="heating_40C.csv")
    assert np.allclose(trace.detuning_mhz, [-1.0, 0.0, 1.0])
    assert np.allclose(trace.signal, [1.5, 3.0, 1.0])
    assert trace.provenance.source_name == "heating_40C.csv"
    assert trace.provenance.sha256 == hashlib.sha256(payload).hexdigest()
    assert trace.provenance.duplicate_rows_merged == 1
    assert trace.provenance.nonfinite_rows == 1
    assert trace.provenance.ignored_rows >= 4  # metadata, header, and bad row
    assert not trace.detuning_mhz.flags.writeable


def test_no_header_unit_assumption_is_explicit():
    trace = load_eit_csv(b"-1 0.2\n0 0.8\n1 0.2\n")
    assert np.allclose(trace.detuning_mhz, [-1.0, 0.0, 1.0])
    assert any("No header" in item for item in trace.provenance.assumptions)
    assert any("assumed MHz" in item for item in trace.provenance.assumptions)


def test_rf_field_loader_standardizes_nanovolts_per_centimeter():
    trace = load_rf_sweep_csv(
        b"SAM field [nV/cm],IF amplitude [V]\n1,0.1\n2,0.2\n3,0.3\n"
    )
    assert trace.drive_quantity == "field"
    assert trace.drive_unit == "V/m"
    assert np.allclose(trace.field_v_m, np.array([1.0, 2.0, 3.0]) * 1e-7)


def test_rf_source_power_loader_keeps_dbm_semantics():
    trace = load_rf_sweep_csv(
        b"Source power [dBm],Beat response\n-80,1\n-70,2\n-60,3\n",
        drive_quantity="source_power",
    )
    assert trace.drive_unit == "dBm"
    assert np.allclose(trace.drive, [-80.0, -70.0, -60.0])
    with pytest.raises(AttributeError):
        _ = trace.field_v_m


def test_asd_loader_converts_frequency_and_noise_units():
    trace = load_psd_csv(
        b"Frequency [kHz],Noise ASD [nV/sqrt(Hz)]\n1,2\n2,3\n3,4\n"
    )
    assert trace.input_spectrum_kind == "asd"
    assert trace.base_unit == "V"
    assert np.allclose(trace.frequency_hz, [1000.0, 2000.0, 3000.0])
    assert np.allclose(trace.asd, np.array([2.0, 3.0, 4.0]) * 1e-9)
    assert np.allclose(trace.psd, trace.asd**2)


def test_psd_loader_converts_squared_units_without_double_squaring():
    trace = load_psd_csv(
        b"Frequency [Hz],Voltage PSD [uV^2/Hz]\n1,4\n2,9\n3,16\n",
        spectrum_kind="psd",
        spectral_unit="uV^2/Hz",
    )
    assert np.allclose(trace.psd, np.array([4.0, 9.0, 16.0]) * 1e-12)
    assert np.allclose(trace.asd, np.array([2.0, 3.0, 4.0]) * 1e-6)


def test_path_source_and_named_columns(tmp_path):
    path = tmp_path / "trace.csv"
    path.write_text(
        "time,detuning_mhz,detector\n0,-1,2\n1,0,3\n2,1,2\n",
        encoding="utf-8",
    )
    trace = load_eit_csv(
        path, x_column="detuning_mhz", signal_column="detector"
    )
    assert trace.provenance.source_name == str(path)
    assert np.allclose(trace.detuning_mhz, [-1.0, 0.0, 1.0])

    mixed_selector = load_eit_csv(
        b"custom_axis,custom_signal\n-1,2\n0,3\n1,2\n",
        x_column="custom_axis",
        signal_column=1,
    )
    assert np.allclose(mixed_selector.signal, [2.0, 3.0, 2.0])


def test_psd_rejects_negative_density():
    with pytest.raises(RydbergCSVError, match="non-negative"):
        load_psd_csv(
            b"Frequency [Hz],ASD [V/sqrt(Hz)]\n1,1\n2,-1\n3,1\n"
        )
