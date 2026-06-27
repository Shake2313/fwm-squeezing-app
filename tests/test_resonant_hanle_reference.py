import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from analysis import resonant_hanle_squeezing_reference as ref  # noqa: E402


def test_oe_probe_locked_source_places_probe_on_87rb_d1_f2_f2():
    cfg = {c.name: c for c in ref.make_source_configs()}["oe_probe_locked"]
    probe = ref.source_probe_offset_ghz(cfg.D_GHz, cfg.delta_code_mhz)
    target = ref.offset_from_fwm_origin_ghz(*ref.RB87_PROBE)

    assert np.isclose((probe - target) * 1e3, 0.0, atol=1e-9)
    assert np.isclose(-cfg.delta_code_mhz, 16.0)


def test_oe_double_locked_source_places_both_modes_on_87rb_d1():
    cfg = {c.name: c for c in ref.make_source_configs()}["oe_double_locked"]
    probe = ref.source_probe_offset_ghz(cfg.D_GHz, cfg.delta_code_mhz)
    conj = ref.source_conjugate_offset_ghz(cfg.D_GHz, cfg.delta_code_mhz)
    probe_target = ref.offset_from_fwm_origin_ghz(*ref.RB87_PROBE)
    conj_target = ref.offset_from_fwm_origin_ghz(*ref.RB87_CONJ)

    assert np.isclose((probe - probe_target) * 1e3, 0.0, atol=1e-9)
    assert np.isclose((conj - conj_target) * 1e3, 0.0, atol=1e-9)
    assert np.isclose(-cfg.delta_code_mhz, -25.631133, atol=1e-6)


def test_absolute_sensitivity_has_shot_noise_power_scaling():
    args = SimpleNamespace(
        responsivity_aw=1.0,
        detector_qe=1.0,
        probe_wavelength_nm=795.0,
        probe_path_eff=1.0,
        reference_path_eff=1.0,
    )
    low = ref.absolute_field_sensitivity(
        1.0, 0.0, 1e-6, 0.0, 1.0, 1.0, args)
    high = ref.absolute_field_sensitivity(
        1.0, 0.0, 4e-6, 0.0, 1.0, 1.0, args)

    assert np.isclose(
        high["sensitivity_pT_per_sqrtHz"],
        0.5 * low["sensitivity_pT_per_sqrtHz"])


def test_electronic_noise_penalizes_total_but_not_shot_limited_sensitivity():
    base_args = SimpleNamespace(
        responsivity_aw=1.0,
        detector_qe=1.0,
        probe_wavelength_nm=795.0,
        probe_path_eff=1.0,
        reference_path_eff=1.0,
        detector_dark_current_na=0.0,
        detector_current_noise_pa_sqrt_hz=0.0,
        balanced_electronic_noise_pa_sqrt_hz=0.0,
        measurement_bandwidth_hz=1.0,
    )
    noisy_args = SimpleNamespace(**vars(base_args))
    noisy_args.balanced_electronic_noise_pa_sqrt_hz = 5.0

    shot = ref.absolute_field_sensitivity(
        0.25, 1.0, 1e-6, 1e-6, 1.0, 1.0, base_args)
    noisy = ref.absolute_field_sensitivity(
        0.25, 1.0, 1e-6, 1e-6, 1.0, 1.0, noisy_args)

    assert noisy["sensitivity_pT_per_sqrtHz"] > shot["sensitivity_pT_per_sqrtHz"]
    assert np.isclose(
        noisy["shot_limited_pT_per_sqrtHz"],
        shot["shot_limited_pT_per_sqrtHz"])


def test_gaussian_intensity_and_equivalent_waist_are_inverse():
    power_w = 42e-6
    waist_m = 1.3e-3

    intensity = ref.gaussian_peak_intensity_mw_cm2(power_w, waist_m)
    got_waist_um = ref.gaussian_waist_um_for_peak_intensity(
        power_w, intensity)

    assert np.isclose(got_waist_um, waist_m * 1e6)


def test_parse_float_axis_accepts_comma_list_and_rejects_empty_axis():
    got = ref.parse_float_axis("0.2, 0.8, 3")

    assert np.allclose(got, [0.2, 0.8, 3.0])
    with pytest.raises(ValueError):
        ref.parse_float_axis(" , ")


def test_load_hanle_csv_converts_milligauss(tmp_path):
    csv = tmp_path / "hanle.csv"
    csv.write_text("coil_mG,voltage\n-10,1.0\n0,2.0\n10,1.0\n",
                   encoding="utf-8")

    got = ref.load_hanle_csv(csv, b_column="coil_mG",
                             signal_column="voltage", b_unit="mG")

    assert np.allclose(got["b_ut"], [-1.0, 0.0, 1.0])
    assert np.allclose(got["signal"], [1.0, 2.0, 1.0])
    assert got["n_points"] == 3


def test_hanle_calibration_fit_recovers_b_axis_and_affine_signal():
    model_b = np.linspace(-1.5, 1.5, 501)
    model_t = 0.992 + 0.004 * np.exp(-((model_b - 0.08) / 0.22) ** 2) \
        - 0.0015 * model_b
    data_b = np.linspace(-1.2, 1.2, 61)
    true_scale = 1.04
    true_offset = -0.07
    true_signal_offset = 0.3
    true_signal_scale = 2.4
    data_t = np.interp(true_scale * data_b + true_offset,
                       model_b, model_t)
    data_signal = true_signal_offset + true_signal_scale * data_t

    fit = ref.fit_hanle_calibration(
        data_b, data_signal, model_b, model_t,
        offset_span_ut=0.12, offset_points=49,
        scale_span=0.06, scale_points=61)

    assert np.isclose(fit["b_scale"], true_scale, atol=0.0021)
    assert np.isclose(fit["b_offset_ut"], true_offset, atol=0.0051)
    assert fit["normalized_rms"] < 1e-9
    assert fit["r2"] > 0.999999


if __name__ == "__main__":
    test_oe_probe_locked_source_places_probe_on_87rb_d1_f2_f2()
    test_oe_double_locked_source_places_both_modes_on_87rb_d1()
    test_absolute_sensitivity_has_shot_noise_power_scaling()
    test_electronic_noise_penalizes_total_but_not_shot_limited_sensitivity()
    test_gaussian_intensity_and_equivalent_waist_are_inverse()
    test_parse_float_axis_accepts_comma_list_and_rejects_empty_axis()
    print("resonant Hanle reference frequency-lock tests OK")
