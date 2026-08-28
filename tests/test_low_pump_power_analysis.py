import pytest

from analysis.squeezing.low_pump_power import scan_low_pump_power as low_pump


def test_transit_scaled_atom_accepts_current_fwm_builder_api():
    captured = {}
    sentinel = object()

    def atom_builder(T, density=None, *, transit_rate=None):
        captured.update(T=T, density=density, transit_rate=transit_rate)
        return sentinel

    wrapped = low_pump.transit_scaled_atom(300e-6, atom_builder)
    result = wrapped(
        low_pump.T_TRANSIT_REF,
        density=1.25,
        transit_rate=0.0,
    )

    expected = low_pump.constants.GAMMA_GG * low_pump.W_TRANSIT_REF / 300e-6
    assert result is sentinel
    assert captured["T"] == pytest.approx(low_pump.T_TRANSIT_REF)
    assert captured["density"] == pytest.approx(1.25)
    assert captured["transit_rate"] == pytest.approx(expected)
