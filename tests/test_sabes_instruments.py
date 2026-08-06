"""
Virtual-instrument checks (Stage E of the optical-table plan).

    python tests/test_sabes_instruments.py   # or: pytest tests/...

Two things matter most here. First, the instruments must be honest about which
of their output is the model's physics and which is a re-expression of a power
spectral density the model already knows -- a plot that looks like an
oscilloscope invites being read as evidence. Second, the package must stay free
of layout and UI imports, because the point of putting it in its own folder was
that it can be lifted out.
"""
import ast
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from sabes import instruments as I, layout  # noqa: E402
from sabes.beamline import SetupSettings, build_source_chain  # noqa: E402
from sabes.instruments import base  # noqa: E402

PACKAGE = ROOT / "sabes" / "instruments"


@pytest.fixture(scope="module")
def chain():
    return build_source_chain()


@pytest.fixture(scope="module")
def seed_offset():
    return -SetupSettings().eom_frequency_hz


def _beam(chain, part, x, y):
    reading = layout.beam_at(part, x, y, chain)
    assert reading is not None, f"no beam at ({x}, {y})"
    return reading.beam


# ------------------------------------------------------------ independence

def test_the_package_imports_nothing_from_the_layout_or_the_ui():
    """It lives in its own folder so it can be lifted out; keep it that way."""
    forbidden = {"streamlit", "matplotlib"}
    for path in PACKAGE.glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = {a.name.split(".")[0] for a in node.names}
            elif isinstance(node, ast.ImportFrom):
                names = {(node.module or "").split(".")[0]}
                assert "layout" not in (node.module or ""), path.name
            else:
                continue
            assert not (names & forbidden), f"{path.name} imports {names}"


# ------------------------------------------------------------ base types

def test_a_trace_rejects_series_that_do_not_match_its_axis():
    with pytest.raises(ValueError, match="points against"):
        base.Trace(x=np.arange(4), series={"y": np.arange(3)})


def test_a_reading_rejects_an_unknown_provenance():
    with pytest.raises(ValueError, match="unknown provenance"):
        base.Reading("x", provenance="probably-fine")


def test_quantum_efficiency_matches_the_detector_the_paper_used():
    """0.58 A/W at 795 nm is where the 90.45 % QE comes from."""
    assert base.quantum_efficiency(0.58, 795.0) == pytest.approx(0.9045,
                                                                 abs=1e-3)


def test_shot_noise_follows_the_square_root_of_power():
    quiet = base.shot_noise_a_per_rthz(0.58, 1e-4)
    loud = base.shot_noise_a_per_rthz(0.58, 4e-4)
    assert loud / quiet == pytest.approx(2.0, rel=1e-9)


# ------------------------------------------------------------ optical

def test_the_wavemeter_shows_the_etalon_chain_doing_its_job(chain, seed_offset):
    """The clearest single demonstration of why the seed path exists.

    Dropped either side of the filters, the carrier falls by nearly three orders
    of magnitude relative to the wanted sideband.
    """
    meter = I.Wavemeter(wanted_offset_hz=seed_offset)
    after_eom = meter.measure(_beam(chain, layout.PART1, 350, 190))
    after_chain = meter.measure(_beam(chain, layout.PART1, 366, 292))

    before = after_eom.quantity("Carrier : wanted").value
    after = after_chain.quantity("Carrier : wanted").value
    assert before > 30.0                      # per cent, straight out of the EOM
    assert after < 0.5                        # the paper's three-etalon figure
    assert before / after > 100.0
    assert after_chain.provenance == I.COMPUTED


def test_the_wavemeter_warns_when_the_carrier_would_kill_the_squeezing(chain,
                                                                      seed_offset):
    """Sim et al. still saw squeezing at 25 %; past that it is gone."""
    raw = I.Wavemeter(wanted_offset_hz=seed_offset).measure(
        _beam(chain, layout.PART1, 350, 190))
    assert any("past the point" in w for w in raw.warnings)


def test_a_power_meter_cannot_tell_the_wanted_line_from_the_rest(chain,
                                                                seed_offset):
    """Which is exactly why purity needs its own instrument."""
    beam = _beam(chain, layout.PART1, 350, 190)
    reading = I.PowerMeter(wanted_offset_hz=seed_offset).measure(beam)
    total = reading.quantity("Total power").value / 1e3
    assert total == pytest.approx(beam.total_power_w, rel=1e-9)
    assert reading.quantity("In the wanted line").value < 90.0
    assert any("overstates" in w for w in reading.warnings)


def test_the_profiler_reports_the_size_here_not_the_waist(chain):
    """Confusing the two is how a beam gets clipped on something that looked
    comfortable in a design note."""
    reading = I.BeamProfiler().measure(_beam(chain, layout.PART2, 745, 214))
    here = reading.quantity("1/e² radius here").value
    waist = reading.quantity("Waist w₀").value
    assert here > waist
    assert reading.quantity("Rayleigh range").value > 0
    assert reading.trace is not None and len(reading.trace.x) > 20


def test_the_profiler_profile_is_a_normalised_gaussian(chain):
    trace = I.BeamProfiler().measure(_beam(chain, layout.PART2, 300, 256)).trace
    peak = trace.series["intensity"].max()
    assert peak == pytest.approx(1.0, abs=1e-9)
    assert trace.series["intensity"][0] < 0.01


# ------------------------------------------------------------ electrical

def test_the_photodiode_turns_power_into_current_and_its_noise(chain):
    detector = I.Photodiode(responsivity_a_per_w=0.58, nep_w_per_rthz=3.3e-12)
    signal = detector.convert(chain.seed.scaled(14.0))
    assert signal.dc_current_a == pytest.approx(
        0.58 * signal.optical_power_w, rel=1e-9)
    assert signal.electronic_noise_a_per_rthz == pytest.approx(0.58 * 3.3e-12)
    assert signal.clearance_db == pytest.approx(
        20 * np.log10(signal.shot_noise_a_per_rthz
                      / signal.electronic_noise_a_per_rthz), abs=1e-9)


def test_the_photodiode_warns_about_saturation(chain):
    detector = I.Photodiode(saturation_w=1e-6)
    reading = detector.measure(chain.seed.scaled(14.0))
    assert any("saturation" in w for w in reading.warnings)


def test_the_scope_scan_is_computed_and_the_time_series_is_not(chain):
    """The distinction the whole package is built around."""
    signal = I.Photodiode().convert(chain.seed.scaled(14.0))
    scope = I.Oscilloscope()
    scan = scope.scan(signal, np.linspace(-1, 1, 32),
                      np.full(32, 1e-4))
    series = scope.timeseries(signal, duration_s=1e-5)
    assert scan.provenance == I.COMPUTED
    assert series.provenance == I.SYNTHESISED
    assert "not simulated" in series.note or "says nothing" in series.note


def test_the_scope_scan_converts_power_to_volts_through_the_head(chain):
    signal = I.Photodiode(responsivity_a_per_w=0.5,
                          transimpedance_v_per_a=1e5).convert(chain.seed)
    scan = I.Oscilloscope().scan(signal, [0.0, 1.0], [1e-5, 2e-5])
    assert scan.trace.series["signal"][1] == pytest.approx(2e-5 * 0.5 * 1e5)


def test_the_time_series_is_reproducible_and_matches_its_own_density(chain):
    signal = I.Photodiode().convert(chain.seed.scaled(14.0))
    scope = I.Oscilloscope()
    first = scope.timeseries(signal, duration_s=2e-4, seed=7)
    again = scope.timeseries(signal, duration_s=2e-4, seed=7)
    assert np.array_equal(first.trace.series["signal"],
                          again.trace.series["signal"])
    # The sampled spread must agree with the density it was drawn from.
    expected_v = (signal.total_noise_a_per_rthz()
                  * np.sqrt(signal.bandwidth_hz)
                  * signal.transimpedance_v_per_a)
    assert np.std(first.trace.series["signal"]) == pytest.approx(expected_v,
                                                                 rel=0.1)


# ------------------------------------------------------------ analyser

def test_the_analyser_stacks_the_three_levels_in_the_right_order(chain):
    signal = I.Photodiode().convert(chain.seed.scaled(14.0))
    reading = I.SpectrumAnalyzer().analyze(signal, -7.0, total_power_w=2e-3)
    snl = reading.quantity("Shot-noise level").value
    squeezed = reading.quantity("Squeezed level").value
    floor = reading.quantity("Amplifier floor").value
    assert squeezed == pytest.approx(snl - 7.0, abs=1e-9)
    assert floor < squeezed < snl
    assert reading.provenance == I.SYNTHESISED


def test_the_analyser_clearance_follows_the_power_it_was_given(chain):
    """Quoting a shot-noise level from one power and a clearance from another
    would be quietly inconsistent, and it was, before this was fixed."""
    signal = I.Photodiode().convert(chain.seed)
    quiet = I.SpectrumAnalyzer().analyze(signal, -7.0, total_power_w=1e-4)
    loud = I.SpectrumAnalyzer().analyze(signal, -7.0, total_power_w=1e-3)
    assert loud.quantity("Clearance").value > quiet.quantity("Clearance").value
    # Shot noise goes as sqrt(P), so a decade of power is 10 dB of clearance in
    # noise power -- and 20 dB for the two decades one might casually assume.
    assert (loud.quantity("Clearance").value
            - quiet.quantity("Clearance").value) == pytest.approx(10.0, abs=0.1)


def test_the_amplifier_floor_eats_into_the_observable_squeezing(chain):
    """What the measurement can show, not what the atoms produced."""
    signal = I.Photodiode().convert(chain.seed)
    reading = I.SpectrumAnalyzer().analyze(signal, -8.0, total_power_w=2e-4)
    observable = reading.quantity("Observable squeezing").value
    assert -8.0 < observable < 0.0
    assert any("above the amplifier floor" in w for w in reading.warnings)


def test_plenty_of_power_recovers_the_full_squeezing(chain):
    signal = I.Photodiode().convert(chain.seed)
    reading = I.SpectrumAnalyzer().analyze(signal, -8.0, total_power_w=5e-3)
    assert reading.quantity("Observable squeezing").value == pytest.approx(
        -8.0, abs=0.3)
    assert not reading.warnings


def test_the_squeezing_is_flat_unless_a_measured_corner_is_supplied(chain):
    """A real trace degrades towards DC. Inventing that shape would dress a
    guess up as a measurement, so the corner defaults to zero."""
    signal = I.Photodiode().convert(chain.seed)
    flat = I.SpectrumAnalyzer().analyze(signal, -8.0, total_power_w=5e-3)
    squeezed = flat.trace.series["squeezed"] - flat.trace.series["shot-noise level"]
    assert np.allclose(squeezed, -8.0, atol=1e-9)

    shaped = I.SpectrumAnalyzer(technical_corner_hz=1e6).analyze(
        signal, -8.0, total_power_w=5e-3)
    spoiled = (shaped.trace.series["squeezed"]
               - shaped.trace.series["shot-noise level"])
    assert spoiled[0] > spoiled[-1]          # worse at low frequency


def test_the_analyser_trace_carries_every_level_it_reports(chain):
    signal = I.Photodiode().convert(chain.seed)
    trace = I.SpectrumAnalyzer().analyze(signal, -7.0,
                                         total_power_w=2e-3).trace
    assert set(trace.series) == {"shot-noise level", "squeezed", "observed",
                                 "amplifier floor"}
    assert trace.x_unit == "MHz" and trace.y_unit == "dBm"


# ------------------------------------------------------------ honesty rule

def test_only_the_two_reconstructed_readings_are_marked_synthesised(chain,
                                                                    seed_offset):
    beam = _beam(chain, layout.PART1, 350, 190)
    signal = I.Photodiode().convert(beam)
    computed = [
        I.PowerMeter(wanted_offset_hz=seed_offset).measure(beam),
        I.BeamProfiler().measure(beam),
        I.Wavemeter(wanted_offset_hz=seed_offset).measure(beam),
        I.Photodiode().measure(beam),
        I.Oscilloscope().scan(signal, [0, 1], [1e-4, 1e-4]),
    ]
    assert all(r.provenance == I.COMPUTED for r in computed)
    synthesised = [
        I.Oscilloscope().timeseries(signal),
        I.SpectrumAnalyzer().analyze(signal, -7.0, total_power_w=1e-3),
    ]
    assert all(r.synthesised and r.note for r in synthesised)


def test_the_page_labels_synthesised_readings_on_screen():
    source = (ROOT / "sabes_page.py").read_text(encoding="utf-8")
    assert "reading.synthesised" in source
    assert "SYNTHESISED" in source


def test_the_page_routes_every_instrument_and_gates_the_electrical_ones():
    """Oscilloscope and analyser must go through a photodiode, as on the bench."""
    import sabes_page

    source = (ROOT / "sabes_page.py").read_text(encoding="utf-8")
    router = source[source.index("def _read_instrument"):
                    source.index("def _scope_reading")]
    for name in sabes_page.INSTRUMENTS:
        assert f'"{name}"' in router, name
    # The electrical branch is reached only after a detector is built.
    assert router.index("_detector_from_calibration") < router.index("Oscilloscope")
    assert "detector.convert(beam)" in router


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
