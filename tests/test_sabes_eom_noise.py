"""Phenomenological EOM-sideband noise must be positive, explicit, and local."""
import sys
from dataclasses import replace
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from sabes import detection  # noqa: E402
from sabes.beamline import SetupSettings, build_source_chain  # noqa: E402
from sabes.noise import (ELEMENTARY_CHARGE_C,  # noqa: E402
                         build_eom_noise_budget)


def _budget(chain=None, rin_db=-100.0):
    return build_eom_noise_budget(
        chain=chain or build_source_chain(),
        gains=(15.0, 14.0),
        optics_transmission=0.945,
        responsivity_a_per_w=0.58,
        rin_db_per_hz=rin_db,
        analysis_bandwidth_hz=100e3,
    )


def test_every_post_etalon_mode_except_the_wanted_seed_enters_the_budget():
    chain = build_source_chain()
    budget = _budget(chain)
    expected = {line.label for line in chain.seed.lines
                if line.offset_hz != chain.seed_offset_hz and line.power_w > 0.0}
    assert {mode.label for mode in budget.modes} == expected
    assert "carrier" in expected
    assert all(mode.order != -1 for mode in budget.modes)


@pytest.mark.parametrize("baseline_db", (-20.0, -8.0, -0.01, 0.0))
def test_any_unwanted_mode_with_positive_rin_can_only_worsen_sub_sql_noise(
        baseline_db):
    budget = _budget()
    assert budget.unwanted_detector_power_w > 0.0
    assert budget.classical_rin_psd_a2_per_hz > 0.0
    assert budget.rin_loaded_db(baseline_db) > baseline_db
    assert budget.rin_loaded_db(baseline_db) > budget.shot_noise_only_db(
        baseline_db)


def test_zero_rin_keeps_only_the_unpaired_poisson_penalty():
    budget = _budget(rin_db=-float("inf"))
    assert budget.rin_per_hz == 0.0
    assert budget.classical_rin_psd_a2_per_hz == 0.0
    assert budget.rin_loaded_db(-8.0) == pytest.approx(
        budget.shot_noise_only_db(-8.0), abs=1e-12)
    assert budget.shot_noise_only_db(-8.0) > -8.0


def test_penalty_label_is_guarded_above_sql_even_though_absolute_psd_is_positive():
    budget = _budget()
    assert budget.rin_penalty_db(-8.0) > 0.0
    assert budget.classical_rin_psd_a2_per_hz > 0.0
    with pytest.raises(ValueError, match="sub-SQL"):
        budget.rin_penalty_db(+3.0)


def test_mode_fano_factor_is_exactly_the_i_squared_rin_law():
    budget = _budget()
    carrier = next(mode for mode in budget.modes if mode.label == "carrier")
    expected_excess = (carrier.photocurrent_a * budget.rin_per_hz
                       / (2.0 * ELEMENTARY_CHARGE_C))
    assert carrier.fano_factor - 1.0 == pytest.approx(expected_excess, rel=1e-12)
    assert carrier.classical_rin_psd_a2_per_hz == pytest.approx(
        carrier.photocurrent_a ** 2 * budget.rin_per_hz, rel=1e-12)


def test_fixed_rin_technical_excess_scales_as_unwanted_power_squared():
    chain = build_source_chain()
    scaled_lines = [
        line if line.offset_hz == chain.seed_offset_hz else line.scaled(2.0)
        for line in chain.seed.lines
    ]
    doubled = replace(chain, seed=chain.seed.with_lines(scaled_lines))
    assert _budget(doubled).classical_rin_psd_a2_per_hz == pytest.approx(
        4.0 * _budget(chain).classical_rin_psd_a2_per_hz, rel=1e-12)


def test_legacy_noise_api_names_remain_aliases():
    budget = _budget()
    assert budget.thermal_like_psd_a2_per_hz == (
        budget.classical_rin_psd_a2_per_hz)
    assert budget.thermal_like_excess_sql == budget.classical_rin_excess_sql
    assert budget.coherent_mixture_db(-8.0) == budget.shot_noise_only_db(-8.0)
    assert budget.loaded_db(-8.0) == budget.rin_loaded_db(-8.0)


def test_declared_rin_reports_power_and_field_amplitude_rms():
    budget = _budget(rin_db=-100.0)
    assert budget.fractional_intensity_rms == pytest.approx(
        (1e-10 * 100e3) ** 0.5, rel=1e-12)
    assert budget.fractional_field_amplitude_rms == pytest.approx(
        0.5 * budget.fractional_intensity_rms, rel=1e-12)


def test_detector_dc_headroom_and_sql_include_unwanted_probe_power():
    settings = SetupSettings()
    chain = build_source_chain(settings)
    geom = detection.geometry(chain, settings)
    readout = detection.readout(geom, chain, (15.0, 14.0), settings)
    probe = next(arm for arm in readout.arms if arm.name == "probe")
    conjugate = next(arm for arm in readout.arms if arm.name == "conjugate")

    assert probe.power_w == pytest.approx(
        chain.seed_power_w * 15.0 * geom.arms[0].optics_transmission
        + readout.eom_noise.unwanted_detector_power_w)
    assert conjugate.power_w == pytest.approx(
        chain.seed_power_w * 14.0 * geom.arms[1].optics_transmission)
    assert readout.total_power_w == pytest.approx(probe.power_w
                                                  + conjugate.power_w)
    assert readout.shot_noise_a_per_rthz ** 2 == pytest.approx(
        readout.eom_noise.sql_psd_a2_per_hz, rel=1e-12)
