"""Jasperse equation limits and current UI absorption/noise separation."""
import json

import numpy as np
import pytest

from gabes.passive_readout import passive_l1_readout
from gabes.schemes import fwm
from gabes_ui.export import build_export_payload


@pytest.mark.parametrize("squeeze,fano", [(0.4, 1.0), (1.5, 1.0), (0.4, 3.0), (1.5, 3.0)])
def test_unequal_external_efficiencies_match_jasperse_eq5_and_seed_variance(squeeze, fano):
    cosine, sine = np.cosh(squeeze), np.sinh(squeeze)
    phase = np.exp(.37j)
    matrix = np.array([[cosine, sine*phase], [sine/phase, cosine]])
    ep, ec = .37, .82
    gain = cosine**2
    denominator = ep*gain+ec*(gain-1)
    coherent = 1+2*(gain-1)*(gain*(ep-ec)**2-ec**2)/denominator
    classical = (fano-1)*(ep*gain-ec*(gain-1))**2/denominator
    result = passive_l1_readout(matrix, seed_fano=fano, eta_probe=ep, eta_conjugate=ec)
    assert result["NF"] == pytest.approx(coherent+classical, rel=1e-12)
    assert result["NF_dB"] == pytest.approx(10*np.log10(result["NF"]), abs=1e-12)


@pytest.fixture(scope="module")
def physical_view():
    scheme = fwm.FWMScheme()
    params = dict(scheme.defaults(), resolution=fwm.FIDELITY_BALANCED,
                  after_cell_eff_pct=50.0, qe_pct=60.0)
    raw = scheme.compute(params)
    view = scheme.headless_observables(raw, params)
    return scheme, params, raw, view


def test_absorption_efficiencies_use_positive_loss_not_net_gain_or_repeated_nf_loss(physical_view):
    _, _, raw, view = physical_view
    assert not raw["gain_closure"]["applied"]
    assert raw["gain_closure"]["coupling_multiplier"] == 1
    assert view["jasperse_readout"]["eta_after_cell_QE"] == pytest.approx(.3)
    absorption = view["dressed_cell_absorption"]
    for arm in ("probe", "conjugate"):
        rates = absorption[arm]
        transmission = np.exp(-rates["alpha_abs_per_m"]*raw["cell_length_m"])
        assert rates["in_cell_absorption_transmission"] == pytest.approx(transmission)
        assert view["sidebar_readouts"][f"eta_{arm}_pct"]["value"] == pytest.approx(100*.3*transmission)
        assert rates["alpha_abs_per_m"] >= 0 and rates["alpha_em_per_m"] >= 0
    assert view["metrics"][0]["label"] == "Squeezing Degree [dB]"


def test_classical_seed_and_external_efficiencies_do_not_change_absorption_or_cached_mean(physical_view):
    scheme, params, raw, first = physical_view
    gains_before = raw["G_s"].copy()
    changed = scheme.headless_observables(
        raw, dict(params, seed_classical_rin_scaled=5.0, after_cell_eff_pct=80., qe_pct=70.))
    assert changed["dressed_cell_absorption"] == first["dressed_cell_absorption"]
    np.testing.assert_array_equal(raw["G_s"], gains_before)
    assert changed["jasperse_readout"]["eta_after_cell_QE"] == pytest.approx(.56)
    for arm in ("probe", "conjugate"):
        assert changed["sidebar_readouts"][f"eta_{arm}_pct"]["value"] == pytest.approx(
            first["sidebar_readouts"][f"eta_{arm}_pct"]["value"]*.56/.3)
    assert changed["passive_l1_operating_point"]["NF"] != first["passive_l1_operating_point"]["NF"]


def test_strict_json_preserves_absorption_emission_and_current_noise_units(physical_view):
    scheme, params, raw, view = physical_view
    payload = build_export_payload(
        scheme=scheme, raw=raw, params=params, view=view, plot_payloads=[],
        comparison_payload=None, cache_version=scheme.cache_version,
        readout_cache_version="jasperse-integration")
    decoded = json.loads(json.dumps(payload, allow_nan=False))
    absorption = decoded["view"]["dressed_cell_absorption"]
    assert absorption["emission_retained_separately"]
    assert absorption["already_in_atomic_transfer_do_not_reapply"]
    assert decoded["view"]["jasperse_readout"]["seed_noise"]["convention"].startswith("one-sided")
    assert decoded["view"]["sidebar_readouts"]["eta_conjugate_pct"]["value"] is not None


@pytest.mark.parametrize("zero_control", ["after_cell_eff_pct", "qe_pct"])
def test_zero_detected_power_has_no_normalized_noise_factor(physical_view, zero_control):
    scheme, params, raw, _ = physical_view
    view = scheme.headless_observables(raw, dict(params, **{zero_control: 0.}))
    readout = view["jasperse_readout"]
    assert np.all(np.isnan(readout["NF"]))
    assert np.all(readout["passive_l1_status"] == "L1_NO_DETECTED_BRIGHT_OUTPUT")
    assert not view["passive_l1_operating_point"]["passive_l1_available"]
    for arm in ("probe", "conjugate"):
        assert view["sidebar_readouts"][f"eta_{arm}_pct"]["value"] == 0.
