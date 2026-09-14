"""Independent thermal flux, first-exit chords and marked-Poisson identities."""

from dataclasses import FrozenInstanceError, replace

import numpy as np
import pytest

from gabes.constants import KB, MASS_85RB
from gabes.quantum.contracts import GeneratorFrequencyAxis
from gabes.quantum.inflow import (
    EntryPhaseQuadrature, average_marked_poisson, entry_phase_factors,
    maxwell_box_inflow, uniform_entry_phases,
)
from gabes.quantum.transport import BallisticPath


LO = np.array([-.0002, -.00015, -.00625])
HI = -LO
TEMPERATURE = 394.15
DENSITY = 1e18


def inflow(*, power=3, seed=0, lower=LO, upper=HI, **kwargs):
    parameters = dict(temperature_K=TEMPERATURE, mass_kg=MASS_85RB,
                      density_m3=DENSITY, source="independent conditional test inputs")
    parameters.update(kwargs)
    return maxwell_box_inflow(lower, upper, points_per_face_power=power, seed=seed, **parameters)


@pytest.fixture(scope="module")
def fine():
    return inflow(power=14)


def test_six_faces_have_analytic_positive_maxwell_flux_without_occupancy_rescaling():
    q = inflow(power=5)
    sigma = np.sqrt(KB*TEMPERATURE/MASS_85RB)
    lengths = HI-LO
    total_area = 2*(lengths[0]*lengths[1]+lengths[0]*lengths[2]+lengths[1]*lengths[2])
    assert len(q.rate_s_inverse) == 6*2**5
    assert set(q.face_index) == set(range(6))
    for face in range(6):
        selected = q.face_index == face
        normal = face//2
        area = np.prod(np.delete(lengths, normal))
        expected_rate = DENSITY*area*sigma/np.sqrt(2*np.pi)
        np.testing.assert_allclose(q.rate_s_inverse[selected], expected_rate/2**5, rtol=3e-15)
        assert q.rate_s_inverse[selected].sum() == pytest.approx(expected_rate, rel=3e-15)
        assert np.all((1 if face % 2 == 0 else -1)*q.velocity_m_s[selected, normal] > 0)
    assert q.total_arrival_rate_s_inverse == pytest.approx(DENSITY*total_area*sigma/np.sqrt(2*np.pi), rel=3e-15)
    assert q.equilibrium_atom_number == pytest.approx(DENSITY*np.prod(lengths), rel=3e-15)
    assert q.mean_occupancy == pytest.approx(q.rate_s_inverse@q.residence_time_s, rel=3e-15)
    # A finite quadrature reports its defect rather than forcing the target nV.
    assert q.mean_occupancy != q.equilibrium_atom_number


def test_incoming_normal_is_flux_weighted_rayleigh_and_tangents_are_maxwell(fine):
    sigma = np.sqrt(KB*TEMPERATURE/MASS_85RB)
    for face in range(6):
        rows = fine.face_index == face
        velocity = fine.velocity_m_s[rows]/sigma
        normal = face//2
        tangent = [i for i in range(3) if i != normal]
        inward = (1 if face % 2 == 0 else -1)*velocity[:, normal]
        assert inward.mean() == pytest.approx(np.sqrt(np.pi/2), abs=.001)
        assert np.mean(inward**2) == pytest.approx(2., abs=.004)
        np.testing.assert_allclose(velocity[:, tangent].mean(axis=0), 0., atol=.001, rtol=0)
        np.testing.assert_allclose(np.mean(velocity[:, tangent]**2, axis=0), 1., atol=.003, rtol=0)
        assert np.mean(np.sum(velocity**2, axis=1)) == pytest.approx(4., abs=.006)


def test_chords_end_at_the_first_box_exit_and_include_oblique_exits():
    q = inflow(power=5, seed=9)
    tangential_exit_count = 0
    for index in range(len(q.rate_s_inverse)):
        entry, velocity = q.entry_position_m[index], q.velocity_m_s[index]
        face = int(q.face_index[index])
        normal = face//2
        assert entry[normal] == (LO[normal] if face % 2 == 0 else HI[normal])
        candidates = []
        for component, speed in enumerate(velocity):
            if speed > 0:
                candidates.append((HI[component]-entry[component])/speed)
            elif speed < 0:
                candidates.append((LO[component]-entry[component])/speed)
            else:
                candidates.append(np.inf)
        expected = min(candidates)
        assert q.residence_time_s[index] == pytest.approx(expected, rel=2e-14)
        path = q.path(index)
        assert isinstance(path, BallisticPath)
        np.testing.assert_array_equal(path.entry_position_m, entry)
        np.testing.assert_array_equal(path.velocity_m_s, velocity)
        locations = path.position(np.linspace(0., expected, 9))
        assert np.all(locations >= LO-1e-17) and np.all(locations <= HI+1e-17)
        exit_component = int(np.argmin(candidates))
        exit_bound = HI[exit_component] if velocity[exit_component] > 0 else LO[exit_component]
        assert locations[-1, exit_component] == pytest.approx(exit_bound, abs=1e-17)
        tangential_exit_count += exit_component != normal
    assert tangential_exit_count > len(q.rate_s_inverse)//4


def test_occupation_recovers_nV_uniform_position_and_bulk_maxwell_moments(fine):
    sigma = np.sqrt(KB*TEMPERATURE/MASS_85RB)
    length = HI-LO
    result = fine.occupation_moments()
    assert result["occupancy_over_nV"] == pytest.approx(1., abs=.004)
    np.testing.assert_allclose(result["velocity_mean_m_s"]/sigma, 0., atol=.003, rtol=0)
    np.testing.assert_allclose(result["velocity_second_m2_s2"]/sigma**2, np.eye(3), atol=.004, rtol=0)
    assert result["speed_mean_m_s"]/(np.sqrt(8/np.pi)*sigma) == pytest.approx(1., abs=.004)
    assert result["speed_fourth_m4_s4"]/(15*sigma**4) == pytest.approx(1., abs=.006)
    np.testing.assert_allclose(result["position_centered_mean_m"]/length, 0., atol=.001, rtol=0)
    normalized_position = result["position_centered_second_m2"]/np.outer(length, length)
    np.testing.assert_allclose(normalized_position, np.eye(3)/12, atol=.0004, rtol=0)
    # Normalization remains analytic nV, not the quadrature's computed count.
    raw_weight = fine.rate_s_inverse*fine.residence_time_s/fine.equilibrium_atom_number
    assert result["occupancy_over_nV"] == pytest.approx(raw_weight.sum(), rel=3e-15)


def test_independent_seed_refinement_reduces_combined_occupation_moment_error():
    sigma = np.sqrt(KB*TEMPERATURE/MASS_85RB)
    length = HI-LO
    def error(q):
        values = q.occupation_moments()
        return max(abs(values["occupancy_over_nV"]-1),
            np.max(abs(values["velocity_second_m2_s2"]/sigma**2-np.eye(3))),
            np.max(abs(12*values["position_centered_second_m2"]/np.outer(length, length)-np.eye(3))))
    coarse, refined = [], []
    for seed in (0, 1):
        coarse.append(error(inflow(power=6, seed=seed)))
        refined.append(error(inflow(power=14, seed=seed)))
    assert max(refined) < .004
    assert np.mean(refined) < np.mean(coarse)/5


@pytest.mark.parametrize("change,rate_factor,time_factor,number_factor", [
    ({"density_m3": 3*DENSITY}, 3., 1., 3.),
    ({"temperature_K": 4*TEMPERATURE}, 2., .5, 1.),
    ({"mass_kg": 4*MASS_85RB}, .5, 2., 1.),
])
def test_density_temperature_and_mass_follow_analytic_flux_and_chord_scaling(change, rate_factor, time_factor, number_factor):
    base, changed = inflow(), inflow(**change)
    np.testing.assert_array_equal(base.entry_position_m, changed.entry_position_m)
    np.testing.assert_array_equal(base.face_index, changed.face_index)
    np.testing.assert_allclose(changed.rate_s_inverse, base.rate_s_inverse*rate_factor, rtol=3e-15)
    np.testing.assert_allclose(changed.residence_time_s, base.residence_time_s*time_factor, rtol=3e-15)
    assert changed.mean_occupancy == pytest.approx(base.mean_occupancy*number_factor, rel=3e-15)
    assert changed.equilibrium_atom_number == pytest.approx(base.equilibrium_atom_number*number_factor, rel=3e-15)


def test_uniform_geometry_scaling_obeys_area_volume_and_transit_units():
    base, larger = inflow(), inflow(lower=3*LO, upper=3*HI)
    np.testing.assert_allclose(larger.entry_position_m, 3*base.entry_position_m, rtol=3e-14, atol=1e-18)
    np.testing.assert_array_equal(larger.velocity_m_s, base.velocity_m_s)
    np.testing.assert_allclose(larger.rate_s_inverse, 9*base.rate_s_inverse, rtol=3e-15)
    # Near a face, subtracting entry from exit magnifies position roundoff.
    np.testing.assert_allclose(larger.residence_time_s, 3*base.residence_time_s, rtol=3e-13)
    assert larger.mean_occupancy == pytest.approx(27*base.mean_occupancy, rel=3e-15)
    translated = inflow(lower=LO+[.002, .003, .04], upper=HI+[.002, .003, .04])
    np.testing.assert_allclose(translated.residence_time_s, base.residence_time_s, rtol=5e-12)
    np.testing.assert_allclose(translated.occupation_moments()["position_centered_second_m2"],
                               base.occupation_moments()["position_centered_second_m2"], rtol=1e-12, atol=1e-20)


def test_seed_replays_every_boundary_mark_and_factory_arrays_are_owned_immutable():
    lower, upper = LO.copy(), HI.copy()
    q, replay = inflow(lower=lower, upper=upper, seed=17), inflow(seed=17)
    lower[:] = -99.
    upper[:] = 99.
    for field in ("lower_corner_m", "upper_corner_m", "entry_position_m", "velocity_m_s",
                  "residence_time_s", "rate_s_inverse", "face_index"):
        np.testing.assert_array_equal(getattr(q, field), getattr(replay, field))
        with pytest.raises(ValueError):
            getattr(q, field).setflags(write=True)
    assert not np.array_equal(q.velocity_m_s, inflow(seed=18).velocity_m_s)
    with pytest.raises(FrozenInstanceError):
        q.seed = 99


@pytest.mark.parametrize("change", [
    {"temperature_K": 0}, {"temperature_K": np.nan}, {"mass_kg": -1},
    {"density_m3": 0}, {"density_m3": np.inf}, {"source": ""},
])
def test_missing_nonphysical_or_unsourced_thermal_inputs_are_rejected(change):
    with pytest.raises(ValueError):
        inflow(**change)


@pytest.mark.parametrize("kwargs", [
    {"power": -1}, {"power": 21}, {"power": 1.5}, {"power": True},
    {"seed": -1}, {"seed": 2**32}, {"seed": 1.5}, {"seed": True},
    {"lower": [0., 0.]}, {"upper": LO}, {"upper": [np.inf, 1., 1.]},
])
def test_geometry_and_quadrature_controls_are_validated(kwargs):
    with pytest.raises(ValueError):
        inflow(**kwargs)


def test_direct_dataclass_replacement_cannot_bypass_the_physical_boundary_contract():
    q = inflow(power=1)
    for changed in (
        {"rate_s_inverse": -q.rate_s_inverse},
        {"rate_s_inverse": 2*q.rate_s_inverse},
        {"residence_time_s": 10*q.residence_time_s},
        {"face_index": np.full(len(q.rate_s_inverse), 99)},
        {"face_index": q.face_index.astype(float)+.5},
        {"velocity_m_s": -q.velocity_m_s},
        {"entry_position_m": q.entry_position_m+.1},
        {"points_per_face_power": 2},
    ):
        with pytest.raises(ValueError):
            replace(q, **changed)
    values = np.array(q.velocity_m_s, copy=True)
    copied = replace(q, velocity_m_s=values)
    values[:] = 0.
    np.testing.assert_array_equal(copied.velocity_m_s, q.velocity_m_s)
    with pytest.raises(ValueError):
        copied.velocity_m_s.setflags(write=True)
    for index in (-1, len(q.rate_s_inverse), .5, True):
        with pytest.raises(ValueError):
            q.path(index)


@pytest.mark.parametrize("offset", [0., .37, 1e20])
def test_uniform_common_phase_has_zero_nonzero_resolved_fourier_moments(offset):
    phases = uniform_entry_phases(8, offset_rad=offset)
    assert len(np.unique(phases.phase_rad)) == 8
    np.testing.assert_allclose(phases.probability, np.full(8, 1/8), atol=0, rtol=0)
    for harmonic in (-3, -1, 1, 2, 3):
        assert abs(phases.probability@np.exp(1j*harmonic*phases.phase_rad)) < 2e-14
    for values in (phases.phase_rad, phases.probability):
        with pytest.raises(ValueError):
            values.setflags(write=True)


def test_common_entry_phase_retains_fixed_spatial_offsets_and_harmonic_relations():
    path = inflow(power=0).path(2)
    vectors = np.array([[1000., 200., -30.], [-400., 50., 10.], [70., -120., 40.]])
    offsets = np.array([.2, -.3, .7])
    static = vectors@path.entry_position_m+offsets
    for phase in uniform_entry_phases(8).phase_rad:
        actual = entry_phase_factors(path, phase, vectors, [1, 1, -1], offsets_rad=offsets)
        np.testing.assert_allclose(actual, np.exp(1j*(static+np.array([1, 1, -1])*phase)), atol=2e-15)
        assert actual[0]*actual[1].conjugate() == pytest.approx(np.exp(1j*(static[0]-static[1])), abs=2e-15)
        assert actual[0]*actual[2] == pytest.approx(np.exp(1j*(static[0]+static[2])), abs=2e-15)
    # Age evolution is not randomized or injected into this entry-only helper.
    different_speed = replace(path, velocity_m_s=2*path.velocity_m_s, residence_time_s=path.residence_time_s/2)
    np.testing.assert_array_equal(entry_phase_factors(path, .4, vectors, [1, 1, -1]),
                                 entry_phase_factors(different_speed, .4, vectors, [1, 1, -1]))


def test_phase_contract_rejects_independent_phases_and_invalid_probabilities():
    path = inflow(power=0).path(0)
    for kwargs in (
        {"common_phase_rad": [.1, .2]}, {"common_phase_rad": 1j},
        {"temporal_harmonics": [1.5, -1]}, {"wavevectors_rad_m": np.zeros((2, 2))},
        {"offsets_rad": [0.]},
    ):
        arguments = dict(common_phase_rad=.1, wavevectors_rad_m=np.zeros((2, 3)), temporal_harmonics=[1, -1])
        arguments.update(kwargs)
        with pytest.raises(ValueError):
            entry_phase_factors(path, **arguments)
    for weights in ([.4, .4], [-1., 2.], [0., 1.], [np.nan, 1.]):
        with pytest.raises(ValueError):
            EntryPhaseQuadrature([0., 1.], weights, "invalid phase fixture")
    for count in (0, True, 1.5):
        with pytest.raises(ValueError):
            uniform_entry_phases(count)


AXIS = GeneratorFrequencyAxis([0., 1e5], "declared atomic wavepacket Fourier axis")


def packet(path, mean, *, greater=None, lesser=None, axis=AXIS):
    mean = np.asarray(mean, complex)
    shape = (len(axis.omega_rad_s), mean.shape[-1], mean.shape[-1])
    return {"frequency_axis": axis, "mean_pulse": mean,
        "greater": np.zeros(shape, complex) if greater is None else greater,
        "lesser": np.zeros(shape, complex) if lesser is None else lesser,
        "residence_time_s": path.residence_time_s}


def test_marked_poisson_averages_per_mark_raw_moments_before_any_mean_subtraction():
    q, phases = inflow(power=0), uniform_entry_phases(8)
    g = np.array([[2., .2j], [-.2j, 1.]])
    l = np.array([[1., -.3j], [.3j, .7]])
    amplitudes = np.array([1., 2*np.exp(.4j)])
    frequency_scale = np.array([1., .7])
    def factory(path, phase):
        mean = path.residence_time_s*frequency_scale[:, None]*amplitudes*np.exp(1j*phase)
        factor = path.residence_time_s**2*frequency_scale[:, None, None]**2
        return packet(path, mean, greater=factor*g*(1+.2*np.cos(phase)),
                      lesser=factor*l*(1+.1*np.sin(phase)))
    result = average_marked_poisson(q, phases, factory, source="analytic two-readout per-mark identity")
    rate_second_lifetime = q.rate_s_inverse@(q.residence_time_s**2)
    multiplier = rate_second_lifetime*frequency_scale[:, None, None]**2
    number = multiplier*np.outer(amplitudes, amplitudes.conj())
    np.testing.assert_allclose(result["poisson_number"], number, rtol=4e-15, atol=0)
    np.testing.assert_allclose(result["greater"], multiplier*g+number, rtol=4e-15, atol=0)
    np.testing.assert_allclose(result["lesser"], multiplier*l+number, rtol=4e-15, atol=0)
    np.testing.assert_allclose(result["greater"]-result["lesser"], multiplier*(g-l), rtol=5e-14, atol=0)
    assert np.linalg.norm(result["arrival_weighted_conditional_mean_pulse"]) < 1e-14*np.mean(q.residence_time_s)
    assert np.linalg.norm(result["poisson_number"]) > 0
    assert result["mean_occupancy"] == q.mean_occupancy
    assert result["arrival_rate_s_inverse"] == q.total_arrival_rate_s_inverse
    assert result["path_count"] == 6 and result["phase_count"] == 8
    assert "zero-cyclic" in result["scope"] and "no optical SQL" in result["scope"]


def test_uniform_entry_phase_keeps_same_atom_spatial_cross_terms_and_count_noise():
    q, phases = inflow(power=1), uniform_entry_phases(8)
    vectors = np.array([[0., 0., 0.], [1000., -400., 25.]])
    offsets = np.array([0., .4])
    def factory(path, phase):
        amplitude = np.array([1., 2.])*entry_phase_factors(path, phase, vectors, [1, 1], offsets_rad=offsets)
        return packet(path, np.tile(path.residence_time_s*amplitude, (2, 1)))
    result = average_marked_poisson(q, phases, factory, source="same-entry-phase spatial coherence fixture")
    expected = np.zeros((2, 2), complex)
    for rate, tau, position in zip(q.rate_s_inverse, q.residence_time_s, q.entry_position_m):
        amplitude = np.array([1., 2.])*np.exp(1j*(vectors@position+offsets))
        expected += rate*tau**2*np.outer(amplitude, amplitude.conj())
    np.testing.assert_allclose(result["greater"], np.broadcast_to(expected, (2, 2, 2)), rtol=4e-15)
    assert abs(expected[0, 1]) > .1*expected[0, 0]
    np.testing.assert_array_equal(result["internal_greater"], np.zeros((2, 2, 2)))
    np.testing.assert_array_equal(result["internal_lesser"], np.zeros((2, 2, 2)))
    # Missing the mean term leaves zero PSD for deterministic single-atom
    # pulses, despite their strictly positive Poisson counting spectrum.
    assert np.linalg.norm(result["greater"]) > 0
    assert np.linalg.norm(result["arrival_weighted_conditional_mean_pulse"]) < 1e-14*np.mean(q.residence_time_s)


@pytest.mark.parametrize("failure", ["missing_mean", "duration", "shape", "negative", "nonhermitian", "grid", "frame"])
def test_invalid_marked_packets_are_rejected_without_dropping_or_resampling(failure):
    q, phases = inflow(power=0), uniform_entry_phases(2)
    seen = []
    def factory(path, phase):
        seen.append((path, phase))
        output = packet(path, np.ones((2, 2))*path.residence_time_s)
        if len(seen) == 2:
            if failure == "missing_mean":
                del output["mean_pulse"]
            elif failure == "duration":
                output["residence_time_s"] *= 2
            elif failure == "shape":
                output["greater"] = np.zeros((2, 1, 1))
            elif failure == "negative":
                output["greater"] = -np.broadcast_to(np.eye(2)*1e-30, (2, 2, 2))
            elif failure == "nonhermitian":
                output["greater"] = np.array([[[1., .2j], [0., 1.]]]*2)
            elif failure == "grid":
                output["frequency_axis"] = GeneratorFrequencyAxis([0., 2e5], AXIS.frame)
            elif failure == "frame":
                output["frequency_axis"] = GeneratorFrequencyAxis(AXIS.omega_rad_s, "different rotating frame")
        return output
    with pytest.raises((ValueError, KeyError)):
        average_marked_poisson(q, phases, factory, source="invalid packet control")
    assert len(seen) == 2


def test_callback_exception_remains_visible_and_cannot_be_replaced_with_a_fit():
    calls = []
    def broken(path, phase):
        calls.append((path, phase))
        raise RuntimeError("unphysical characteristic fixture")
    with pytest.raises(RuntimeError, match="unphysical characteristic"):
        average_marked_poisson(inflow(power=0), uniform_entry_phases(2), broken, source="failed callback fixture")
    assert len(calls) == 1
