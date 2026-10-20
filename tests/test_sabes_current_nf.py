"""SABES uses current seed-noise NF and preserves unavailable readouts."""
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from gabes.passive_readout import passive_l1_from_spectrum
from gabes.schemes import fwm
from sabes import bridge
import sabes_page


def _result():
    result = bridge.run(solve=False)
    center = fwm.branch_center_GHz(result.params["opd"], -1)
    amplifier = np.array([[2.0, np.sqrt(3.0)], [np.sqrt(3.0), 2.0]], complex)
    raw = {
        "D_GHz": result.params["opd"],
        "raman_center_minus_GHz": center,
        "probe_axis_GHz": center + np.array([-8.0, 0.0, 8.0])*1e-3,
        "T_canonical_small_signal": np.tile(amplifier, (3, 1, 1)),
        "G_s": np.full(3, 4.0), "G_c": np.full(3, 3.0),
        "input_seed_power_W": 8e-6,
        "eta": 1.0,
        "gain_referred_noise_dB": np.full(3, -20.0),
    }
    return replace(result, raw=raw)


def test_geometry_and_legacy_overrides_map_to_the_current_external_efficiency():
    result = bridge.run(solve=False)
    assert result.params["after_cell_eff_pct"] == pytest.approx(100-result.geometry.loss_pct)
    assert result.eta == pytest.approx(
        result.params["after_cell_eff_pct"]*result.params["qe_pct"]/1e4)
    old = bridge.run(solve=False, loss_pct=25.0, qe_pct=80.0)
    assert old.params["after_cell_eff_pct"] == 75.0
    assert old.eta == pytest.approx(0.6)
    explicit = bridge.run(solve=False, loss_pct=25.0, after_cell_eff_pct=90.0)
    assert explicit.params["after_cell_eff_pct"] == 90.0
    assert explicit.params["loss_pct"] == 10.0
    assert bridge.solve_key(old.params) == bridge.solve_key(result.params)


def test_page_noise_readout_is_cheap_current_and_matches_the_passive_api(monkeypatch):
    result = _result()
    params = dict(result.params, tpd=0.0, seed_classical_rin_scaled=2.0,
                  after_cell_eff_pct=75.0, qe_pct=80.0)
    # This page path must not invoke figures, headless tables, or the optional
    # pump-dressed absorption/QRT extraction.
    monkeypatch.setattr(fwm.FWMScheme, "headless_observables", lambda *a, **kw:
                        pytest.fail("full headless observables are not needed"))
    view = bridge.noise_observables(result.raw, params)
    current = replace(result, params=params, observables=view)
    seed = view["jasperse_readout"]["seed_noise"]
    expected = passive_l1_from_spectrum(dict(result.raw, eta=0.6),
                                       seed_fano=seed["seed_fano"])
    assert current.noise_factor == pytest.approx(expected["NF"][1])
    assert current.gain_referred_noise_db == pytest.approx(expected["NF_dB"][1])
    assert current.squeezing_db == current.gain_referred_noise_db
    assert current.gain_referred_noise_db != -20.0
    assert current.physical_squeezing_db is None


def test_invalid_current_nf_never_falls_back_to_finite_legacy_raw_noise():
    result = _result()
    unavailable = {"jasperse_readout": {}, "passive_l1_operating_point": {
        "NF": float("nan"), "NF_dB": float("nan")}}
    current = replace(result, observables=unavailable, readout=SimpleNamespace())
    assert current.gain_referred_noise_db is None
    assert current.noise_factor is None
    assert current.eom_shot_noise_only_db is None
    assert current.eom_rin_loaded_noise_db is None
    assert current.eom_rin_penalty_db is None
    assert sabes_page._noise_db_text(current.gain_referred_noise_db) == "unavailable"
    assert sabes_page._noise_db_text(float("nan")) == "unavailable"
    assert sabes_page._noise_db_text(-3.0) == "-3.00 dB"
    # Old fixtures with no current-view descriptor retain their explicit alias.
    assert result.gain_referred_noise_db == -20.0


def test_current_nonpassive_and_outside_points_keep_unavailable_status():
    result = _result()
    raw = dict(result.raw)
    raw["T_canonical_small_signal"] = 1.1*raw["T_canonical_small_signal"]
    current = replace(result, raw=raw,
        observables=bridge.noise_observables(raw, result.params))
    assert current.gain_referred_noise_db is None
    outside_params = dict(result.params, tpd=30.0)
    outside = replace(result, params=outside_params,
        observables=bridge.noise_observables(result.raw, outside_params))
    assert outside.gain_referred_noise_db is None
    assert outside.observables["passive_l1_operating_point"]["passive_l1_status"] == "L1_OUTSIDE_SAMPLED_SCAN"
