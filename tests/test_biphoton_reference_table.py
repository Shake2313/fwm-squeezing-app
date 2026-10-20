"""Biphoton literature rows keep source, filter and input planes separate."""
from types import SimpleNamespace

from gabes.schemes import fwm


def test_reference_rate_uses_source_plane_and_filter_setting_is_not_a_pass():
    scheme = fwm.FWMScheme()
    raw = {
        "topology": SimpleNamespace(name=fwm.TOPOLOGY_RB87_TELECOM),
        "predictive": True,
        "velocity_converged": True,
        "pair_rate_cps": 380.0,
        "od_estimate": 112.0,
        "phase_match_weight": 1.0,
        "energy_mismatch_hz": 0.0,
    }
    params = {"pump_biphoton_uw": 10.0, "filter_bandwidth_mhz": 300.0}
    stats = {
        "pair_rate_cps": 124.6,
        "g2_peak": 44.0,
        "raw_g2_peak": 44.0,
        "detected_fwhm_ns": 0.56,
    }
    table = scheme._reference_validation_table(raw, params, stats)
    assert "| Source pair rate / pump | 38000 cps/mW |" in table
    assert "| Filter-passed pair rate / pump | 12460 cps/mW |" in table
    bandwidth_row = next(line for line in table.splitlines()
                         if line.startswith("| Bandwidth setting |"))
    assert "| setting only |" in bandwidth_row
    assert "PASS" not in bandwidth_row

    # An arbitrary filter change must not alter the unfiltered anchor row.
    stats["pair_rate_cps"] = 38.0
    params["filter_bandwidth_mhz"] = 30.0
    changed = scheme._reference_validation_table(raw, params, stats)
    assert "| Source pair rate / pump | 38000 cps/mW |" in changed
    assert "| Filter-passed pair rate / pump | 3800 cps/mW |" in changed
