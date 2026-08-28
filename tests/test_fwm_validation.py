"""Validation-boundary tests for the seeded-FWM diagnostic path."""

import numpy as np
import pytest

from gabes import constants, hyperfine, observables
from gabes.schemes import fwm
from sabes import bridge


def _trace_row(n=fwm.N_LEVELS):
    row = np.zeros(n * n)
    row[np.arange(n) * (n + 1)] = 1.0
    return row


def test_thermal_transit_reset_is_trace_preserving_and_has_thermal_fixed_point():
    rate = 2.0 * np.pi * 100e3
    reset = fwm.thermal_transit_reset_superoperator(rate)
    target = np.diag([
        hyperfine.GROUND_POP[2],
        hyperfine.GROUND_POP[3],
        0.0,
        0.0,
    ]).reshape(-1)
    trace = _trace_row()

    assert np.max(np.abs(trace @ reset)) < 1e-9
    assert np.max(np.abs(reset @ target)) < 1e-9

    rho = np.diag([0.1, 0.2, 0.3, 0.4]).astype(complex)
    rho[0, 1] = 0.07j
    rho[1, 0] = -0.07j
    expected = rate * (target * np.trace(rho) - rho.reshape(-1))
    assert reset @ rho.reshape(-1) == pytest.approx(expected)


def test_collisional_atom_splits_transit_reset_from_collision_dephasing():
    temperature = 394.15
    density = hyperfine.number_density(temperature)
    transit = 2.0 * np.pi * 93e3
    atom = fwm.collisional_atom(temperature, density, transit_rate=transit)
    collision = constants.ground_coherence_dephasing(
        temperature, density, floor=0.0)

    assert atom.transit_reset_rate == pytest.approx(transit)
    assert atom.ground_collision_dephasing_rate == pytest.approx(collision)
    assert np.max(np.abs(_trace_row() @ atom.lindblad)) < 1e-7
    coherence_index = atom.rho_index(fwm.G1, fwm.G2)
    assert atom.lindblad[coherence_index, coherence_index].real == pytest.approx(
        -(transit + collision))


def test_gain_referred_name_is_primary_and_legacy_alias_is_exact():
    probe = np.array([1.0, 3.0])
    conjugate = np.array([0.0, 2.0])
    primary = observables.gain_referred_noise_dB(probe, conjugate, 0.87)
    legacy = observables.intensity_difference_squeezing_dB(
        probe, conjugate, 0.87)
    assert legacy == pytest.approx(primary)


def test_operating_point_accepts_primary_diagnostic_without_legacy_key():
    spectrum = {
        "probe_axis_GHz": np.array([0.9, 0.901]),
        "raman_center_minus_GHz": 0.9,
        "raman_center_plus_GHz": 0.901,
        "G_s": np.array([2.0, 3.0]),
        "G_c": np.array([1.0, 2.0]),
        "gain_referred_noise_dB": np.array([-1.0, -2.0]),
    }
    point = fwm.operating_point(spectrum, 0.0, branch=-1)
    assert point["gain_referred_noise_dB"] == pytest.approx(-1.0)
    assert point["physical_squeezing_dB"] is None


def test_seeded_claim_gate_blocks_quantitative_gain_and_physical_squeezing():
    gate = fwm.seeded_validation_claim_gate(
        canonical_mode_status="conditional: matched waist",
        commutator_defect_max=np.array([0.25]),
        eom_residual_carrier_power=2e-6,
    )
    assert gate["level"] == "MEAN_FIELD_DIAGNOSTIC"
    assert not gate["quantitative_gain_supported"]
    assert not gate["physical_squeezing_prediction"]
    assert not gate["physical_claims_allowed"]
    assert "QUANTITATIVE_GAIN_UNSUPPORTED" in gate["badges"]
    assert "PHYSICAL_SQUEEZING_UNAVAILABLE" in gate["badges"]
    assert any("EOM" in reason for reason in gate["reasons"])


def test_seeded_claim_gate_reports_unapplied_eom_status_even_at_zero_power():
    gate = fwm.seeded_validation_claim_gate(
        canonical_mode_status="explicit",
        commutator_defect_max=np.array([0.0]),
        eom_spectrum_status="unsupported",
        eom_spectrum_application="unapplied",
    )
    assert any("unsupported/unapplied" in reason for reason in gate["reasons"])


def test_sabes_result_exposes_physical_squeezing_as_explicitly_unavailable():
    raw = {
        "D_GHz": np.array([0.0]),
        "probe_axis_GHz": np.array([0.0]),
        "G_s": np.array([1.0]),
        "G_c": np.array([0.0]),
        "gain_referred_noise_dB": np.array([0.0]),
        "physical_squeezing_dB": None,
        "validation_level": "MEAN_FIELD_DIAGNOSTIC",
        "claim_gate": {"physical_claims_allowed": False},
    }
    result = bridge.SabesResult(
        settings=object(), chain=object(), geometry=object(),
        params={"opd": 0.0, "tpd": 0.0}, raw=raw)

    assert result.physical_squeezing_db is None
    assert result.gain_referred_noise_db == pytest.approx(0.0)
    assert result.validation_level == "MEAN_FIELD_DIAGNOSTIC"
    assert result.claim_gate == {"physical_claims_allowed": False}


def test_solver_and_parameter_provenance_do_not_claim_a_fit():
    solver = fwm.seeded_atomic_solver_provenance()
    pump_reference = fwm.pump_weak_response_reference_provenance()
    noncollinear = fwm.noncollinear_doppler_reference_provenance()
    provenance = fwm.seeded_parameter_provenance(
        line_strength=0.74,
        transit_rate=constants.GAMMA_GG,
    )
    assert solver["floquet_order"] == fwm.SEEDED_FLOQUET_ORDER == 3
    assert solver["floquet_modes"] == tuple(range(-3, 4))
    assert solver["state_equation"] == "full density-matrix Liouvillian null solve"
    assert not solver["pump_only_self_consistent_nullspace"]
    assert solver["pump_only_weak_response_reference_available"]
    assert solver["pump_only_reference_solver_id"] == pump_reference["solver_id"]
    assert solver["noncollinear_doppler_reference_available"]
    assert not solver["rate_sylvester_approximation"]
    assert pump_reference["supported_branches"] == (-1,)
    assert +1 in pump_reference["unsupported_branches"]
    assert pump_reference["reference_fields"] == "none"
    assert not pump_reference["production_default"]
    assert pump_reference["quantum_noise"] == "not implemented"
    assert noncollinear["supported_branches"] == (-1,)
    assert not noncollinear["production_default"]
    assert "never velocity shifted" in noncollinear["lab_optical_beat"]
    for name in ("ell_s", "kappa", "gamma_transit"):
        assert "not an uncertainty interval" in provenance[name]["range_kind"]
    assert "not refitted" in provenance["ell_s"]["calibration_status"]


def test_direct_and_scheme_seeded_residual_defaults_are_identical():
    assert fwm.FWMScheme().defaults()["line_strength"] == pytest.approx(
        fwm.SEEDED_REFERENCE_RESIDUAL)


def test_scheme_forwards_transit_and_eom_inputs_without_noise_inference(monkeypatch):
    captured = {}

    def fake_compute_spectrum(detuning, **kwargs):
        captured["detuning"] = detuning
        captured.update(kwargs)
        return {"captured": True}

    monkeypatch.setattr(fwm, "compute_spectrum", fake_compute_spectrum)
    scheme = fwm.FWMScheme()
    params = scheme.defaults()
    params.update(
        probe_uw=7.5,
        eom_residual_carrier_uw=2.0,
        eom_other_sidebands_uw=0.75,
        eom_seed_spectrum_provenance="test spectrum",
        eom_seed_spectrum_status="unsupported",
        eom_seed_spectrum_application="unapplied",
        transit_rate_khz=91.0,
        floquet_order=4,
    )
    assert scheme.compute(params) == {"captured": True}
    assert captured["P_probe"] == pytest.approx(7.5e-6)
    assert captured["eom_residual_carrier_power"] == pytest.approx(2.0e-6)
    assert captured["eom_other_sidebands_power"] == pytest.approx(0.75e-6)
    assert captured["eom_seed_spectrum_application"] == "unapplied"
    assert captured["transit_rate"] == pytest.approx(2.0 * np.pi * 91e3)
    assert captured["floquet_order"] == 4
    assert captured["enforce_floquet_convergence"] is True


def test_full_spectrum_view_forwards_the_same_transit_reset_rate(monkeypatch):
    captured = {}

    def fake_full_spectrum(*args, **kwargs):
        captured.update(kwargs)
        return {"captured": True}

    monkeypatch.setattr(fwm, "full_spectrum", fake_full_spectrum)
    params = fwm.FWMScheme().defaults()
    params["transit_rate_khz"] = 91.0
    params["floquet_order"] = 4
    params["resolution"] = fwm.FIDELITY_BALANCED
    params["seeded_angle_deg"] = 0.45
    params.update(
        eom_residual_carrier_uw=2.0,
        eom_other_sidebands_uw=0.75,
        eom_seed_spectrum_provenance="full-view test spectrum",
        eom_seed_spectrum_status="unsupported",
        eom_seed_spectrum_application="unapplied",
    )
    result = fwm.FWMScheme().extra_views()[0].compute(params)
    assert result == {"captured": True}
    assert captured["transit_rate"] == pytest.approx(2.0 * np.pi * 91e3)
    assert captured["floquet_order"] == 4
    assert captured["phase_detail"] == fwm.PHASE_FINE
    assert captured["model_fidelity"] == fwm.FIDELITY_BALANCED
    assert captured["pump_probe_angle_deg"] == pytest.approx(0.45)
    assert captured["velocity_step"] == fwm.FWM_FIDELITY[
        fwm.FIDELITY_BALANCED]["velocity_step"]
    assert captured["eom_residual_carrier_power"] == pytest.approx(2.0e-6)
    assert captured["eom_other_sidebands_power"] == pytest.approx(0.75e-6)
    assert captured["eom_seed_spectrum_provenance"] == "full-view test spectrum"
    assert captured["eom_seed_spectrum_status"] == "unsupported"
    assert captured["eom_seed_spectrum_application"] == "unapplied"


def test_sabes_bridge_preserves_the_full_eom_power_ledger():
    result = bridge.run(solve=False)
    params = result.params
    assert params["probe_uw"] == pytest.approx(params["seed_wanted_sideband_uw"])
    assert params["eom_residual_carrier_uw"] == pytest.approx(
        result.chain.eom_residual_carrier_power_w * 1e6, abs=5.1e-4)
    assert params["eom_other_sidebands_uw"] == pytest.approx(
        result.chain.eom_other_sidebands_power_w * 1e6, abs=5.1e-4)
    assert params["eom_seed_spectrum_status"] == "unsupported"
    assert params["eom_seed_spectrum_application"] == "unapplied"


@pytest.mark.parametrize("override_key", ["probe_uw", "seed_wanted_sideband_uw"])
def test_sabes_wanted_sideband_aliases_share_one_solver_and_cache_value(override_key):
    result = bridge.run(solve=False, **{override_key: 12.0})
    assert result.params["probe_uw"] == pytest.approx(12.0)
    assert result.params["seed_wanted_sideband_uw"] == pytest.approx(12.0)
    assert ("probe_uw", 12.0) in bridge.solve_key(result.params)


def test_sabes_rejects_conflicting_wanted_sideband_aliases():
    with pytest.raises(ValueError, match="same wanted EOM sideband power"):
        bridge.run(
            solve=False, probe_uw=11.0, seed_wanted_sideband_uw=12.0)


@pytest.mark.parametrize(
    "override",
    [
        {"eom_seed_spectrum_status": "supported"},
        {"eom_seed_spectrum_application": "applied"},
    ],
)
def test_sabes_cannot_promote_unapplied_eom_metadata(override):
    with pytest.raises(ValueError, match="EOM"):
        bridge.run(solve=False, **override)


def test_direct_fwm_cannot_promote_unapplied_eom_metadata():
    with pytest.raises(ValueError, match="EOM spectrum status"):
        fwm.compute_spectrum(
            0.9, eom_seed_spectrum_status="supported",
            eom_seed_spectrum_application="applied")


def test_sabes_eom_ratios_are_recomputed_after_power_quantization():
    result = bridge.run(
        solve=False,
        probe_uw=8.0004,
        eom_residual_carrier_uw=0.0046,
        eom_other_sidebands_uw=0.0026,
    )
    params = result.params
    assert params["probe_uw"] == pytest.approx(8.0)
    assert params["eom_residual_carrier_uw"] == pytest.approx(0.005)
    assert params["eom_other_sidebands_uw"] == pytest.approx(0.003)
    assert params["eom_residual_carrier_to_wanted_ratio"] == pytest.approx(
        0.005 / 8.0)
    assert params["eom_other_sidebands_to_wanted_ratio"] == pytest.approx(
        0.003 / 8.0)
