"""Optical sideband, phase, vacuum and detector conventions of the PSR readout."""
import numpy as np
import pytest

from gabes.psr_homodyne import psr_homodyne_readout, psr_output_homodyne_readout
from gabes.psr_noise import psr_gaussian_readout
from gabes.psr_spectrum import PSR_NAMBU_MODES, PSR_RF_FRAME
from gabes.quantum.contracts import GeneratorFrequencyAxis
from gabes.quantum.traveling import NambuTransfer


def paired_transfer(positive, greater=None, lesser=None, frequencies_hz=(1e6,)):
    """Provide an analytic positive-f transfer and its adjoint at negative f."""
    frequency = np.array(frequencies_hz, dtype=float)
    shape = (frequency.size, 4, 4)
    positive = np.broadcast_to(positive, shape)
    greater = np.zeros(shape) if greater is None else np.broadcast_to(greater, shape)
    lesser = np.zeros(shape) if lesser is None else np.broadcast_to(lesser, shape)
    swap = [2, 3, 0, 1]

    def reflected(block):
        return block[::-1].conj()[:, swap][:, :, swap]

    return NambuTransfer(
        GeneratorFrequencyAxis(2*np.pi*np.r_[-frequency[::-1], frequency], PSR_RF_FRAME),
        [1, 1, -1, -1], PSR_NAMBU_MODES,
        np.concatenate((reflected(positive), positive)),
        np.concatenate((reflected(lesser), greater)),
        np.concatenate((reflected(greater), lesser)),
    )


def test_passive_frequency_asymmetric_sidebands_return_vacuum_at_every_phase():
    eta = np.array([.7, .4, .7, .15])
    t = np.diag(np.sqrt(eta)*np.exp(1j*np.array([.1, -.7, -.1, -.53])))
    positive_order = np.diag((1-eta)*[1, 1, 0, 0])
    negative_order = np.diag((1-eta)*[0, 0, 1, 1])
    transfer = paired_transfer(t, positive_order, negative_order, frequencies_hz=(100e3, 2e6))
    result = psr_homodyne_readout(transfer, detection_efficiency=.78)
    np.testing.assert_allclose(result["frequencies_hz"], [100e3, 2e6])
    np.testing.assert_allclose(result["source_sideband_covariance"],
                               np.tile(np.eye(4)/2, (2, 1, 1)), atol=1e-14)
    np.testing.assert_allclose(result["noise_a_db"], 0., atol=1e-13)
    np.testing.assert_allclose(result["noise_b_db"], 0., atol=1e-13)
    assert not result["angle_defined"].any()
    assert np.isnan(result["min_angle_deg"]).all()


def test_flat_shear_matches_canonical_phase_scan_with_one_external_loss():
    kappa = 1.7
    # b'=(1−iκ/2)b + iκ b†/2 gives x'=x+κp and p'=p.
    u, v = 1-0.5j*kappa, 0.5j*kappa
    t = np.eye(4, dtype=complex)
    t[np.ix_([1, 3], [1, 3])] = [[u, v], [v.conjugate(), u.conjugate()]]
    transfer = paired_transfer(t)
    phase = [0., 30., 78., 90., 210.]
    for eta in (1., .56, 0.):
        result = psr_homodyne_readout(transfer, detection_efficiency=eta, phase_degrees=phase)
        expected = psr_gaussian_readout(kappa, 0., detection_efficiency=eta, phase_degrees=phase)
        for key in ("source_covariance", "covariance", "noise_a", "noise_b", "minimum_db", "maximum_db"):
            np.testing.assert_allclose(result[key][0], expected[key], atol=1e-12)
        if eta > 0:
            assert result["angle_defined"][0]
            assert result["min_angle_deg"][0] == pytest.approx(expected["min_angle_deg"])


def test_positive_but_unresolved_squeezing_is_rejected_until_loss_restores_precision():
    r = 9.5
    u, v = np.cosh(r), np.exp(.23j)*np.sinh(r)
    t = np.eye(4, dtype=complex)
    t[np.ix_([1, 3], [1, 3])] = [[u, v], [v.conjugate(), u]]
    transfer = paired_transfer(t)
    with pytest.raises(ValueError, match="unresolved|unrepresentable|nonpositive"):
        psr_homodyne_readout(transfer)
    eta = .7
    result = psr_homodyne_readout(transfer, detection_efficiency=eta)
    assert result["minimum_variance"][0] == pytest.approx(1-eta + eta*np.exp(-2*r), abs=4e-8)


def test_phase_and_noise_follow_the_complex_product_of_two_distinct_sideband_transmissions():
    r = .61
    # A lossless pair source followed by different complex transmissions.
    # <a+ a−> = h+ h− sinh(r)cosh(r). The LO reference changes the angle twice.
    upper = np.sqrt(.7)*np.exp(-.26j)
    lower = np.sqrt(.23)*np.exp(.91j)
    t = np.eye(4, dtype=complex)
    c, s = np.cosh(r), np.sinh(r)
    t[np.ix_([1, 3], [1, 3])] = [[upper*c, upper*s], [lower.conj()*s, lower.conj()*c]]
    greater = np.diag([0., 1-abs(upper)**2, 0., 0.])
    lesser = np.diag([0., 0., 0., 1-abs(lower)**2])
    transfer = paired_transfer(t, greater, lesser)
    phi_lo, efficiency = .31, .82
    result = psr_homodyne_readout(transfer, detection_efficiency=efficiency, carrier_phase_rad=phi_lo)

    assert result["source_sideband_covariance"][0, 0, 0] == pytest.approx(.5 + abs(upper)**2*s*s)
    assert result["source_sideband_covariance"][0, 2, 2] == pytest.approx(.5 + abs(lower)**2*s*s)
    base = 1 + (abs(upper)**2+abs(lower)**2)*s*s
    pair = upper*lower*s*c
    theta = np.deg2rad(result["phase_degrees"]) + phi_lo
    expected_a = base + 2*(pair*np.exp(-2j*theta)).real
    expected_b = base - 2*(pair*np.exp(-2j*theta)).real
    np.testing.assert_allclose(result["noise_a"][0], 1 + efficiency*(expected_a-1), atol=1e-13)
    np.testing.assert_allclose(result["noise_b"][0], 1 + efficiency*(expected_b-1), atol=1e-13)
    angle = (np.degrees(.5*np.angle(pair) - phi_lo) + 90.) % 180
    assert result["min_angle_deg"][0] == pytest.approx(angle)
    assert result["minimum_variance"][0] == pytest.approx(1+efficiency*(base-2*abs(pair)-1))


def test_orphan_frequency_and_inconsistent_reflected_noise_are_rejected():
    base = paired_transfer(np.eye(4))
    single = NambuTransfer(GeneratorFrequencyAxis([2*np.pi*1e6], PSR_RF_FRAME),
                           base.signs, base.mode_labels, base.transfer[1:], base.noise_greater[1:],
                           base.noise_lesser[1:])
    with pytest.raises(ValueError, match="reflected RF partner"):
        psr_homodyne_readout(single)

    # Each signed transfer is physical by itself. Only one has added thermal
    # noise, which is inconsistent with the reflected ordering of the same modes.
    noise = base.noise_greater.copy()
    noise[0] = np.eye(4)*.25
    bad = NambuTransfer(base.frequency_axis, base.signs, base.mode_labels, base.transfer, noise, noise)
    assert bad.audit()["passed"]
    with pytest.raises(ValueError, match="conjugation"):
        psr_homodyne_readout(bad)


def test_dc_is_skipped_and_a_different_generator_frame_is_rejected():
    t = np.broadcast_to(np.eye(4), (3, 4, 4))
    noise = np.zeros_like(t)
    args = ([1, 1, -1, -1], PSR_NAMBU_MODES, t, noise, noise)
    omega = 2*np.pi*np.array([-5e5, 0., 5e5])
    transfer = NambuTransfer(GeneratorFrequencyAxis(omega, PSR_RF_FRAME), *args)
    result = psr_homodyne_readout(transfer)
    np.testing.assert_allclose(result["frequencies_hz"], [5e5])
    with pytest.raises(ValueError, match="rotating frame"):
        psr_homodyne_readout(NambuTransfer(GeneratorFrequencyAxis(omega, "optical_beat"), *args))


def test_output_x_reference_preserves_existing_readout_exactly():
    r, angle = .41, .63
    t = np.eye(4, dtype=complex)
    t[np.ix_([1, 3], [1, 3])] = [
        [np.cosh(r), np.exp(1j*angle)*np.sinh(r)],
        [np.exp(-1j*angle)*np.sinh(r), np.cosh(r)]]
    transfer = paired_transfer(t)
    before = psr_homodyne_readout(transfer, detection_efficiency=.81)
    after = psr_output_homodyne_readout(transfer, [1., 0.], detection_efficiency=.81)
    for key, value in before.items():
        if isinstance(value, np.ndarray):
            np.testing.assert_array_equal(value, after[key])
        elif key != "convention":
            assert value == after[key]


def test_elliptic_output_projection_preserves_whole_phase_scan_under_common_phase():
    # Physical output modes: carrier u and independent orthogonal mode v.
    # An elliptic unitary mixer and arbitrary common phase follow a lossy pair
    # source. This rotates the unequal noise orderings as well as the signal.
    r = .48
    source = np.eye(4, dtype=complex)
    source[np.ix_([1, 3], [1, 3])] = [[np.cosh(r), 1j*np.sinh(r)],
                                      [-1j*np.sinh(r), np.cosh(r)]]
    eta = np.array([.8, .65, .8, .4])
    transmission = np.sqrt(eta)*np.exp(1j*np.array([.1, .27, -.1, -.51]))
    base = paired_transfer(transmission[:, None]*source,
                           np.diag((1-eta)*[1, 1, 0, 0]),
                           np.diag((1-eta)*[0, 0, 1, 1]))
    theta, elliptic_phase = .6, .7
    u = np.array([np.cos(theta), np.exp(1j*elliptic_phase)*np.sin(theta)])
    v = np.array([-np.exp(-1j*elliptic_phase)*np.sin(theta), np.cos(theta)])
    unitary = np.column_stack((u, v))
    expected = psr_homodyne_readout(base, detection_efficiency=.73, carrier_phase_rad=.19)
    for phase in (0., .81, -2.7):
        physical = np.exp(1j*phase)*unitary
        transform = np.zeros((4, 4), complex)
        transform[:2, :2], transform[2:, 2:] = physical, physical.conj()
        changed = NambuTransfer(base.frequency_axis, base.signs, base.mode_labels,
            transform @ base.transfer,
            transform @ base.noise_greater @ transform.conj().T,
            transform @ base.noise_lesser @ transform.conj().T)
        actual = psr_output_homodyne_readout(changed, np.exp(1j*phase)*u,
                                            detection_efficiency=.73, carrier_phase_rad=.19)
        for key in ("source_covariance", "covariance", "noise_a", "noise_b", "minimum_db", "maximum_db", "min_angle_deg"):
            np.testing.assert_allclose(actual[key], expected[key], atol=2e-13)
        assert actual["polarization_phase_chart"] == "x"
        assert np.vdot(actual["canonical_output_jones_xy"], actual["orthogonal_jones_xy"]) == pytest.approx(0., abs=1e-15)


def test_y_polarization_chart_and_invalid_output_jones():
    transfer = paired_transfer(np.eye(4))
    result = psr_output_homodyne_readout(transfer, [0., np.exp(.81j)])
    assert result["polarization_phase_chart"] == "y"
    np.testing.assert_allclose(result["noise_a_db"], 0., atol=1e-13)
    assert not result["angle_defined"].any()
    for jones in ([0., 0.], [1., 1.], [np.nan, 0.], [1.]):
        with pytest.raises(ValueError, match="Jones"):
            psr_output_homodyne_readout(transfer, jones)


def test_driven_atom_output_readout_is_invariant_to_common_optical_phase():
    from gabes.psr_full_obe import solve_psr_pump
    from gabes.psr_spectrum import psr_local_generator
    from gabes.quantum.traveling import constant_segment

    jones = np.array([.8, .6*np.exp(.3j)])
    previous = None
    for phase in (0., .7, -1.5):
        rotated = np.exp(1j*phase)*jones
        state = solve_psr_pump(400., intensity_mw_cm2=12.5, density_m3=1e16,
                               ground_relaxation_hz=1e5, jones_xy=rotated)
        local = psr_local_generator(state, [-4e5, 4e5], method="optimized")
        transfer = constant_segment(local, .001)
        actual = psr_output_homodyne_readout(transfer, rotated, detection_efficiency=.83)
        if previous is not None:
            for key in ("covariance", "noise_a", "noise_b", "min_angle_deg"):
                np.testing.assert_allclose(actual[key], previous[key], atol=2e-11)
        previous = actual
