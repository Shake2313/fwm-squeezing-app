"""Smooth raw QRT controls: ordering, SI units, phase rotation and sparse lift."""

from dataclasses import replace

import numpy as np
import pytest

from analysis.grand_challenge.reference.segmented_qrt import segmented_qrt
from analysis.grand_challenge.reference.smooth_qrt import smooth_qrt
from analysis.grand_challenge.reference.transport_qrt import characteristic_qrt
from gabes import core
from gabes.quantum.reservoirs import CollapseChannel, ExplicitReservoirs
from gabes.quantum.transport import BallisticPath, poisson_beam_spectrum


@pytest.fixture(autouse=True)
def single_thread():
    with core.blas_single_thread():
        yield


def _case():
    lower = np.array([[0., 1.], [0., 0.]])
    return dict(
        h0=np.array([[0., .12+.07j], [.12-.07j, .6]]),
        h1=np.array([[.1, .37-.18j], [.37+.18j, -.2]]),
        envelope=lambda t: np.exp(-((t-.57)/.31)**2)+.15*np.sin(4*t),
        reservoirs=ExplicitReservoirs(2, (
            CollapseChannel("emission", np.sqrt(.8)*lower, "decay", "declared toy"),
            CollapseChannel("dephasing", np.sqrt(.13)*np.diag([1., -1.]),
                            "dephasing", "declared toy"))),
        boundary_state=np.array([[.63, .11+.07j], [.11-.07j, .37]]),
        duration_s=1.3,
        readouts=np.array([lower+(.2+.13j)*np.eye(2)+.3j*lower.T,
                           .4*lower-.7j*lower.T+np.diag([.2, -.15j])]),
        frequencies_rad_s=np.array([[.6, -.9], [-1.3, .25], [0., 0.]]))


def _direct(case, *, rtol=5e-13, atol=5e-15):
    """Old centered-source QRT in seconds, with unequal phases in the readout."""
    path = BallisticPath([0, 0, 0], [0, 0, 1], case["duration_s"], "smooth QRT test")

    def hamiltonian(age, position):
        return case["h0"]+case["envelope"](age)*case["h1"]

    outputs = []
    for frequencies in case["frequencies_rad_s"]:
        common = float(np.mean(frequencies))

        def readout(age, position):
            return np.exp(1j*(frequencies-common)*age)[:, None, None]*case["readouts"]

        outputs.append(characteristic_qrt(
            path, hamiltonian, case["reservoirs"], case["boundary_state"], readout,
            np.array([common]), rtol=rtol, atol=atol, intervals=64))
    result = {key: np.concatenate([out[key] for out in outputs])
              for key in ("greater", "lesser", "mean_pulse")}
    result["exit_state"] = outputs[0]["exit_state"]
    return result


def _relative_errors(actual, expected):
    """Each frequency is checked separately, so DC cannot hide carrier error."""
    errors = {}
    for key in ("greater", "lesser", "mean_pulse", "exit_state"):
        values = [(actual[key], expected[key])] if key == "exit_state" else zip(actual[key], expected[key])
        errors[key] = max(np.linalg.norm(value-reference)/max(np.linalg.norm(reference), np.finfo(float).tiny)
                          for value, reference in values)
    return errors


@pytest.mark.parametrize("amplitude", [0., .73, -1.2])
def test_constant_envelope_matches_segmented_exponential(amplitude):
    case = dict(_case(), envelope=lambda t: amplitude)
    actual = smooth_qrt(**case)
    expected = segmented_qrt(
        (case["h0"]+amplitude*case["h1"])[None], case["reservoirs"], case["boundary_state"],
        [case["duration_s"]], case["readouts"][None], case["frequencies_rad_s"])
    for key in ("greater", "lesser", "raw_greater", "raw_lesser", "mean_pulse", "exit_state"):
        np.testing.assert_allclose(actual[key], expected[key], rtol=2e-11, atol=3e-14)


def test_smooth_noncommuting_drive_matches_centered_physical_time_qrt():
    case = _case()
    assert np.linalg.norm(case["h0"]@case["h1"]-case["h1"]@case["h0"]) > .1
    actual, expected = smooth_qrt(**case), _direct(case)
    errors = _relative_errors(actual, expected)
    assert max(errors.values()) < 2e-10, errors
    assert np.linalg.norm(actual["mean_pulse"]) > .1
    assert np.max(abs(actual["greater"].imag)) > .01
    assert np.linalg.norm(actual["greater"]-actual["lesser"]) > .01
    for key in ("greater", "lesser"):
        np.testing.assert_allclose(actual[key], actual[key].conj().swapaxes(-1, -2), atol=2e-15)
        assert np.linalg.eigvalsh(actual[key]).min() > -2e-12


def test_seconds_rescaling_preserves_envelope_clock_and_pulse_units():
    case = _case()
    expected = smooth_qrt(**case)
    scale = 7e-9
    sampled = []

    def scaled_envelope(t):
        sampled.append(t)
        return case["envelope"](t/scale)

    actual = smooth_qrt(**dict(
        case, h0=case["h0"]/scale, h1=case["h1"]/scale, envelope=scaled_envelope,
        duration_s=case["duration_s"]*scale, frequencies_rad_s=case["frequencies_rad_s"]/scale,
        reservoirs=ExplicitReservoirs(2, tuple(
            replace(channel, operator=channel.operator/np.sqrt(scale))
            for channel in case["reservoirs"].channels))))
    assert min(sampled) == 0.
    assert max(sampled) == pytest.approx(case["duration_s"]*scale, rel=1e-14, abs=0.)
    for key, power in (("greater", 2), ("lesser", 2), ("raw_greater", 2),
                       ("raw_lesser", 2), ("mean_pulse", 1), ("exit_state", 0)):
        np.testing.assert_allclose(actual[key]/scale**power, expected[key], rtol=2e-12, atol=5e-15)


def test_frequency_batch_shares_one_density_and_matches_individual_rows():
    case = _case()
    actual = smooth_qrt(**case)
    for f, omega in enumerate(case["frequencies_rad_s"]):
        single = smooth_qrt(**dict(case, frequencies_rad_s=omega[None]))
        for key in ("greater", "lesser", "raw_greater", "raw_lesser", "mean_pulse"):
            np.testing.assert_allclose(actual[key][f], single[key][0], rtol=2e-11, atol=3e-14)
        np.testing.assert_allclose(actual["exit_state"], single["exit_state"], atol=3e-14)
    # Three rows cost 4 + 3*(16+8+2) = 82, rather than 3*30 = 90.
    assert actual["qrt_block_dimension"] == 82
    assert max(actual["affine_generator_nnz"]) < 82**2/4
    assert actual["stored_time_points"] == 1


def test_identity_offset_changes_raw_moments_and_mean_but_not_connected_noise():
    case = _case()
    expected = smooth_qrt(**case)
    offset = np.array([.4+.2j, -.7j])
    actual = smooth_qrt(**dict(case, readouts=case["readouts"]+offset[:, None, None]*np.eye(2)))
    w, duration = case["frequencies_rad_s"], case["duration_s"]
    mean_shift = offset*duration*np.exp(.5j*w*duration)*np.sinc(w*duration/(2*np.pi))
    np.testing.assert_allclose(actual["mean_pulse"]-expected["mean_pulse"], mean_shift, atol=3e-14)
    for key in ("greater", "lesser"):
        np.testing.assert_allclose(actual[key], expected[key], rtol=2e-11, atol=4e-14)
    assert np.linalg.norm(actual["raw_greater"]-expected["raw_greater"]) > .1


def test_identity_readout_retains_poisson_arrival_number_noise():
    case = dict(_case(), readouts=np.array([np.eye(2), (1+.3j)*np.eye(2)]))
    actual = smooth_qrt(**case)
    duration, w = case["duration_s"], case["frequencies_rad_s"]
    mean = np.array([1., 1+.3j])*duration*np.exp(.5j*w*duration)*np.sinc(w*duration/(2*np.pi))
    raw = mean[:, :, None]*mean[:, None, :].conj()
    poisson = poisson_beam_spectrum(actual, 3., source="declared test Poisson arrivals")
    np.testing.assert_allclose(actual["mean_pulse"], mean, atol=4e-14)
    for key in ("greater", "lesser"):
        np.testing.assert_allclose(actual[key], 0., atol=5e-14)
        np.testing.assert_allclose(actual["raw_"+key], raw, rtol=2e-11, atol=5e-14)
        np.testing.assert_allclose(poisson[key], 3*raw, rtol=2e-11, atol=2e-13)


def test_damped_lowering_has_exact_finite_lorentzian_window():
    lower = np.array([[0., 1.], [0., 0.]])
    reservoirs = ExplicitReservoirs(2, (
        CollapseChannel("emission", lower, "decay", "unit decay test"),))
    omega, duration = np.array([.2, .7, 1.1]), 3.
    actual = smooth_qrt(np.diag([0., .7]), np.zeros((2, 2)), lambda t: np.sin(t),
                        reservoirs, np.diag([1., 0.]), duration, lower[None], omega[:, None])
    s = -.5+1j*(omega-.7)
    expected = 2*((np.expm1(duration*s)-duration*s)/s**2).real
    np.testing.assert_allclose(actual["greater"][:, 0, 0], expected, rtol=2e-11)
    np.testing.assert_allclose(actual["lesser"], 0., atol=2e-15)
    np.testing.assert_allclose(actual["mean_pulse"], 0., atol=2e-15)


def _short_rb_case():
    from gabes.fwm_quantum.field import reduced_readout_operators
    from gabes.fwm_quantum.model import reduced_pump_system
    from gabes.schemes.fwm import OMEGA_HF

    h0, reservoirs = reduced_pump_system(0., 2*np.pi*.9e9)
    driven, _ = reduced_pump_system(2*np.pi*60e6, 2*np.pi*.9e9)
    duration = 40e-9
    return dict(h0=h0, h1=driven-h0, envelope=lambda t: np.exp(-((t/duration-.45)/.3)**2),
                reservoirs=reservoirs, boundary_state=np.diag([5/12, 7/12, 0., 0.]),
                duration_s=duration, readouts=reduced_readout_operators(),
                frequencies_rad_s=np.array([[-OMEGA_HF+2*np.pi*.4e6,
                                             OMEGA_HF-2*np.pi*.7e6]]))


def test_short_smooth_physical_rb85_d1_against_centered_qrt():
    case = _short_rb_case()
    actual = smooth_qrt(**case, rtol=2e-11, atol=2e-14)
    expected = _direct(case, rtol=2e-11, atol=2e-14)
    errors = _relative_errors(actual, expected)
    assert max(errors.values()) < 2e-8, errors
    tighter = smooth_qrt(**case, rtol=2e-12, atol=2e-15)
    refinement = _relative_errors(actual, tighter)
    assert max(refinement.values()) < 2e-8, refinement
    assert actual["qrt_block_dimension"] == 90
    assert np.linalg.norm(actual["mean_pulse"]) > 1e-14


@pytest.mark.parametrize("changes, message", [
    ({"h0": np.eye(3)}, "h0"),
    ({"h1": np.eye(3)}, "h1"),
    ({"h0": np.full((2, 2), np.nan)}, "finite"),
    ({"h1": [[0., 1.], [0., 0.]]}, "Hermitian"),
    ({"readouts": np.zeros((0, 2, 2))}, "readouts"),
    ({"readouts": np.full((2, 2, 2), np.nan)}, "finite"),
    ({"frequencies_rad_s": [0., 1.]}, "frequencies_rad_s"),
    ({"frequencies_rad_s": [[0., 1j]]}, "real"),
    ({"frequencies_rad_s": [[0., np.inf]]}, "finite"),
    ({"frequencies_rad_s": np.zeros((0, 2))}, "frequencies_rad_s"),
    ({"boundary_state": np.diag([1.1, -.1])}, "physical"),
    ({"boundary_state": np.diag([.7, .7])}, "physical"),
    ({"boundary_state": np.array([[.5, .2j], [.2j, .5]])}, "physical"),
    ({"boundary_state": np.eye(3)/3}, "boundary_state"),
    ({"duration_s": 0.}, "duration_s"),
    ({"duration_s": -1.}, "duration_s"),
    ({"duration_s": np.inf}, "finite"),
    ({"duration_s": 1j}, "real"),
    ({"duration_s": [1.]}, "scalar"),
    ({"max_step_s": 0.}, "max_step_s"),
    ({"max_step_s": np.nan}, "finite"),
    ({"rtol": 0.}, "rtol"),
    ({"rtol": 1.}, "rtol"),
    ({"atol": -1.}, "atol"),
    ({"envelope": lambda t: np.nan}, "finite"),
    ({"envelope": lambda t: 1j}, "real"),
    ({"envelope": lambda t: [1.]}, "scalar"),
])
def test_invalid_inputs_are_rejected(changes, message):
    with pytest.raises(ValueError, match=message):
        smooth_qrt(**dict(_case(), **changes))


@pytest.mark.parametrize("changes, message", [
    ({"envelope": 1.}, "callable"), ({"reservoirs": []}, "ExplicitReservoirs"),
])
def test_invalid_callback_and_reservoir_types_are_rejected(changes, message):
    with pytest.raises(TypeError, match=message):
        smooth_qrt(**dict(_case(), **changes))


def test_callback_is_validated_during_propagation():
    with pytest.raises(ValueError, match="envelope.*finite"):
        smooth_qrt(**dict(_case(), envelope=lambda t: 1. if t < .4 else np.nan))
