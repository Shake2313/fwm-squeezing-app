"""Independent segmented QRT: operator ordering, phases and boundary memory."""

from dataclasses import replace

import numpy as np
import pytest

from analysis.grand_challenge.reference.segmented_qrt import segmented_qrt
from analysis.grand_challenge.reference.transport_qrt import characteristic_qrt
from gabes import core
from gabes.quantum.reservoirs import CollapseChannel, ExplicitReservoirs
from gabes.quantum.transport import BallisticPath, poisson_beam_spectrum


@pytest.fixture(autouse=True)
def single_thread():
    with core.blas_single_thread():
        yield


def _two_level_case():
    lower = np.array([[0., 1.], [0., 0.]])
    h = np.array([[[0., .37+.18j], [.37-.18j, .6]],
                  [[.1, -.21+.31j], [-.21-.31j, -.4]]])
    reservoirs = ExplicitReservoirs(2, (
        CollapseChannel("emission", np.sqrt(.8)*lower, "decay", "declared toy"),
        CollapseChannel("dephasing", np.sqrt(.13)*np.diag([1., -1.]),
                        "dephasing", "declared toy")))
    rho = np.array([[.63, .11+.07j], [.11-.07j, .37]])
    first = lower + (.2+.13j)*np.eye(2) + .3j*lower.T
    second = .4*lower - .7j*lower.T + np.diag([.2, -.15j])
    readouts = np.array([[first, second],
                         [.7j*second + .1*np.eye(2), (.8-.2j)*first]])
    return dict(hamiltonians=h, reservoirs=reservoirs, boundary_state=rho,
                durations_s=np.array([.37, .63]), readouts=readouts,
                frequencies_rad_s=np.array([[.6, -.9], [-1.3, .25], [0., 0.]]))


def _direct_characteristic(case):
    """Original solve_ivp QRT, with port phases embedded in its callback.

    The original API takes one common omega per row. Giving port j the
    readout exp(i*(w_j-common)*t)*O_j makes its physical Y_j exactly ours.
    It integrates centered QRT sources directly in physical time, independently
    of the new raw-source rotating blocks and end-of-pulse mean subtraction.
    """
    edges = np.cumsum(case["durations_s"])
    duration = float(edges[-1])
    path = BallisticPath([0, 0, 0], [0, 0, 1], duration, "segmented QRT test")

    def segment(age):
        return min(int(np.searchsorted(edges, age, side="right")), len(edges)-1)

    def hamiltonian(age, position):
        return case["hamiltonians"][segment(age)]

    results = []
    for frequencies in case["frequencies_rad_s"]:
        common = float(np.mean(frequencies))

        def readout(age, position):
            return (np.exp(1j*(frequencies-common)*age)[:, None, None]
                    * case["readouts"][segment(age)])

        results.append(characteristic_qrt(
            path, hamiltonian, case["reservoirs"], case["boundary_state"],
            readout, np.array([common]), rtol=5e-13, atol=5e-15, intervals=64))
    output = {key: np.concatenate([result[key] for result in results])
              for key in ("greater", "lesser", "mean_pulse")}
    output["exit_state"] = results[0]["exit_state"]
    return output


def _relative_errors(actual, expected):
    return {key: np.linalg.norm(actual[key]-expected[key])/np.linalg.norm(expected[key])
            for key in ("greater", "lesser", "mean_pulse", "exit_state")}


def _phase_integral(omega, start, duration):
    return (duration*np.exp(1j*omega*(start+duration/2))
            * np.sinc(omega*duration/(2*np.pi)))


def _closed_system_operator_integrals(case):
    """Exact Heisenberg operator integral from each Hamiltonian's eigenvectors.

    There are no QRT source, triangle, covariance or augmented-generator
    equations in this analytic control. Global U is an n by n unitary.
    """
    n = case["boundary_state"].shape[0]
    frequencies = case["frequencies_rad_s"]
    pulses = np.zeros(frequencies.shape+(n, n), complex)
    unitary = np.eye(n, dtype=complex)
    start = 0.
    for h, operators, duration in zip(case["hamiltonians"], case["readouts"], case["durations_s"]):
        energies, vectors = np.linalg.eigh(h)
        for f, omega in enumerate(frequencies):
            for j, operator in enumerate(operators):
                local = vectors.conj().T @ operator @ vectors
                beat = omega[j] + energies[:, None]-energies[None, :]
                integral = vectors @ (local*_phase_integral(beat, 0., duration)) @ vectors.conj().T
                pulses[f, j] += np.exp(1j*omega[j]*start)*unitary.conj().T @ integral @ unitary
        step = (vectors*np.exp(-1j*energies*duration)) @ vectors.conj().T
        unitary = step @ unitary
        start += duration
    rho = case["boundary_state"]
    mean = np.einsum("fjab,ba->fj", pulses, rho)
    raw_greater = np.empty(frequencies.shape+(frequencies.shape[1],), complex)
    raw_lesser = np.empty_like(raw_greater)
    for f, operators in enumerate(pulses):
        for j, left in enumerate(operators):
            for k, right in enumerate(operators):
                raw_greater[f, j, k] = np.trace(left @ right.conj().T @ rho)
                raw_lesser[f, j, k] = np.trace(right.conj().T @ left @ rho)
    outer = mean[:, :, None]*mean[:, None, :].conj()
    return dict(greater=raw_greater-outer, lesser=raw_lesser-outer,
                raw_greater=raw_greater, raw_lesser=raw_lesser, mean_pulse=mean,
                exit_state=unitary @ rho @ unitary.conj().T)


def test_changed_hamiltonians_complex_readouts_and_unequal_signed_frequencies():
    case = _two_level_case()
    actual = segmented_qrt(**case)
    expected = _direct_characteristic(case)
    errors = _relative_errors(actual, expected)
    assert max(errors.values()) < 2e-9, errors
    assert np.linalg.norm(actual["mean_pulse"]) > .1
    assert np.max(abs(actual["greater"].imag)) > .01
    assert np.linalg.norm(actual["greater"]-actual["lesser"]) > .01
    for key in ("greater", "lesser"):
        np.testing.assert_allclose(actual[key], actual[key].conj().swapaxes(-1, -2), atol=2e-15)
        assert np.linalg.eigvalsh(actual[key]).min() > -2e-13
    assert actual["qrt_block_dimension"] == 30
    assert actual["qrt_exponentials"] == 6
    assert actual["density_exponentials"] == 2


@pytest.mark.parametrize("fast", [False, True])
def test_exact_closed_system_operator_integral_including_fast_ghz_phases(fast):
    case = _two_level_case()
    case["reservoirs"] = ExplicitReservoirs(2, ())
    if fast:
        # Thousands of carrier cycles; the analytic control never time-steps.
        case["durations_s"] *= 2e-6
        case["hamiltonians"] *= 1e6
        case["frequencies_rad_s"] = 2*np.pi*np.array([
            [3.0357e9, -3.0353e9], [-2.7181e9, 1.419e9], [0., 0.]])
    actual = segmented_qrt(**case)
    expected = _closed_system_operator_integrals(case)
    # Compare each frequency separately so the DC row cannot hide carrier errors.
    for key in ("greater", "lesser", "raw_greater", "raw_lesser", "mean_pulse"):
        for value, reference in zip(actual[key], expected[key]):
            error = np.linalg.norm(value-reference)/np.linalg.norm(reference)
            assert error < (3e-7 if fast else 2e-12), (key, error)
    np.testing.assert_allclose(actual["exit_state"], expected["exit_state"], rtol=2e-13, atol=2e-14)


def test_segment_subdivision_preserves_all_prior_source_and_phase_memory():
    case = _two_level_case()
    expected = segmented_qrt(**case)
    split = dict(case)
    split["hamiltonians"] = np.repeat(case["hamiltonians"], 2, axis=0)
    split["readouts"] = np.repeat(case["readouts"], 2, axis=0)
    split["durations_s"] = (case["durations_s"][:, None]*[.21, .79]).ravel()
    actual = segmented_qrt(**split)
    for key in ("greater", "lesser", "raw_greater", "raw_lesser", "mean_pulse", "exit_state"):
        np.testing.assert_allclose(actual[key], expected[key], rtol=3e-13, atol=2e-15)


def test_disjoint_segment_readouts_retain_nonzero_cross_correlations():
    lower = np.array([[0., 1.], [0., 0.]])
    readouts = np.zeros((2, 2, 2, 2), complex)
    readouts[0, 0] = lower
    readouts[1, 1] = (1+.4j)*lower
    case = dict(hamiltonians=np.zeros((2, 2, 2)), reservoirs=ExplicitReservoirs(2, ()),
                boundary_state=np.diag([.7, .3]), durations_s=np.array([.4, .6]),
                readouts=readouts, frequencies_rad_s=np.array([[.8, -.3]]))
    actual = segmented_qrt(**case)
    expected = _closed_system_operator_integrals(case)
    assert abs(actual["greater"][0, 0, 1]) > .1
    for key in ("greater", "lesser"):
        np.testing.assert_allclose(actual[key], expected[key], rtol=3e-14, atol=2e-15)


def test_piecewise_identity_shifts_change_only_the_global_mean_and_raw_moments():
    case = _two_level_case()
    expected = segmented_qrt(**case)
    offsets = np.array([[.4+.2j, -.7j], [.2-.5j, -.1+.3j]])
    shifted = dict(case, readouts=case["readouts"]+offsets[:, :, None, None]*np.eye(2))
    actual = segmented_qrt(**shifted)
    mean_shift = np.zeros_like(actual["mean_pulse"])
    start = 0.
    for duration, offset in zip(case["durations_s"], offsets):
        mean_shift += offset*_phase_integral(case["frequencies_rad_s"], start, duration)
        start += duration
    np.testing.assert_allclose(actual["mean_pulse"]-expected["mean_pulse"], mean_shift, atol=4e-15)
    outer = actual["mean_pulse"][:, :, None]*actual["mean_pulse"][:, None, :].conj()
    for key in ("greater", "lesser"):
        np.testing.assert_allclose(actual[key], expected[key], rtol=3e-13, atol=4e-15)
        np.testing.assert_allclose(actual["raw_"+key], actual[key]+outer, atol=2e-15)


def test_identity_has_zero_internal_noise_and_exact_poisson_raw_moments():
    case = _two_level_case()
    case["readouts"] = np.broadcast_to(np.eye(2), (2, 2, 2, 2)).copy()
    actual = segmented_qrt(**case)
    duration = case["durations_s"].sum()
    mean = _phase_integral(case["frequencies_rad_s"], 0., duration)
    raw = mean[:, :, None]*mean[:, None, :].conj()
    np.testing.assert_allclose(actual["mean_pulse"], mean, rtol=2e-14, atol=2e-15)
    poisson = poisson_beam_spectrum(actual, 3., source="independent Poisson identity pulses")
    for key in ("greater", "lesser"):
        np.testing.assert_allclose(actual[key], 0., atol=5e-15)
        np.testing.assert_allclose(actual["raw_"+key], raw, rtol=2e-14, atol=5e-15)
        np.testing.assert_allclose(poisson[key], 3*raw, rtol=2e-14, atol=2e-14)
    assert poisson["mean_occupancy"] == pytest.approx(3*duration)


def test_damped_lowering_readout_has_the_exact_finite_lorentzian_window():
    lower = np.array([[0., 1.], [0., 0.]])
    reservoirs = ExplicitReservoirs(2, (
        CollapseChannel("emission", lower, "decay", "unit decay test"),))
    omega = np.array([.2, .7, 1.1])
    duration = 3.
    actual = segmented_qrt(np.diag([0., .7])[None], reservoirs, np.diag([1., 0.]),
                           [duration], lower[None, None], omega[:, None])
    s = -.5+1j*(omega-.7)
    expected = 2*((np.expm1(duration*s)-duration*s)/s**2).real
    np.testing.assert_allclose(actual["greater"][:, 0, 0], expected, rtol=3e-14)
    np.testing.assert_allclose(actual["lesser"], 0., atol=2e-15)
    np.testing.assert_allclose(actual["mean_pulse"], 0., atol=2e-15)


def test_seconds_rescaling_preserves_physical_pulse_units():
    case = _two_level_case()
    expected = segmented_qrt(**case)
    scale = 7e-9
    rescaled = dict(case, hamiltonians=case["hamiltonians"]/scale,
                    durations_s=case["durations_s"]*scale,
                    frequencies_rad_s=case["frequencies_rad_s"]/scale,
                    reservoirs=ExplicitReservoirs(2, tuple(
                        replace(channel, operator=channel.operator/np.sqrt(scale))
                        for channel in case["reservoirs"].channels)))
    actual = segmented_qrt(**rescaled)
    for key, power in (("greater", 2), ("lesser", 2), ("mean_pulse", 1), ("exit_state", 0)):
        np.testing.assert_allclose(actual[key]/scale**power, expected[key], rtol=3e-13, atol=3e-15)


def test_short_physical_rb85_d1_protocol_against_full_density_solve_ivp():
    from gabes.fwm_quantum.field import reduced_readout_operators
    from gabes.fwm_quantum.model import reduced_pump_system
    from gabes.schemes.fwm import OMEGA_HF

    h0, reservoirs = reduced_pump_system(2*np.pi*60e6, 2*np.pi*.9e9)
    h1, _ = reduced_pump_system(2*np.pi*43e6, 2*np.pi*.91e9)
    operators = reduced_readout_operators()
    case = dict(hamiltonians=np.array([h0, h1]), reservoirs=reservoirs,
                boundary_state=np.diag([5/12, 7/12, 0., 0.]),
                durations_s=np.array([.15e-9, .21e-9]),
                readouts=np.array([operators, operators*np.array([.8+.2j, .7-.3j])[:, None, None]]),
                frequencies_rad_s=np.array([[-OMEGA_HF+2*np.pi*.4e6,
                                             OMEGA_HF-2*np.pi*.7e6]]))
    actual = segmented_qrt(**case)
    expected = _direct_characteristic(case)
    errors = _relative_errors(actual, expected)
    assert max(errors.values()) < 3e-9, errors
    assert actual["qrt_block_dimension"] == 90
    assert np.linalg.norm(actual["mean_pulse"]) > 1e-14


@pytest.mark.parametrize("durations", [[0., .5], [.5, 0.], [0., 0.], [-.1, 1.],
                                        [np.nan, 1.], [np.inf, 1.], [1.], [[.4, .6]],
                                        [.4+1j, .6]])
def test_zero_negative_or_invalid_segment_durations_are_rejected(durations):
    with pytest.raises(ValueError, match="durations_s"):
        segmented_qrt(**dict(_two_level_case(), durations_s=durations))


@pytest.mark.parametrize("changes, message", [
    ({"hamiltonians": np.zeros((0, 2, 2))}, "hamiltonians"),
    ({"hamiltonians": np.array([[[0, 1], [0, 0]]]*2)}, "Hermitian"),
    ({"readouts": np.zeros((2, 0, 2, 2))}, "readouts"),
    ({"readouts": np.full((2, 2, 2, 2), np.nan)}, "finite"),
    ({"frequencies_rad_s": [0., 1.]}, "frequencies_rad_s"),
    ({"frequencies_rad_s": [[0., 1j]]}, "real"),
    ({"frequencies_rad_s": [[0., np.inf]]}, "finite"),
    ({"frequencies_rad_s": np.zeros((0, 2))}, "frequencies_rad_s"),
    ({"boundary_state": np.diag([1.1, -.1])}, "physical"),
    ({"boundary_state": np.diag([.7, .7])}, "physical"),
    ({"boundary_state": np.array([[.5, .2j], [.2j, .5]])}, "physical"),
    ({"boundary_state": np.eye(3)/3}, "boundary_state"),
])
def test_invalid_dimensions_frequencies_and_boundary_states_are_rejected(changes, message):
    with pytest.raises(ValueError, match=message):
        segmented_qrt(**dict(_two_level_case(), **changes))
