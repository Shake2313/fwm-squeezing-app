"""Sign, unit, geometry and source-plane checks shared by the PSR models."""
import json
from pathlib import Path

import numpy as np
import pytest

from gabes import constants, psr_conventions as contract


PAPERS = Path(__file__).resolve().parents[1] / "references" / "psr_experiments" / "papers"


def paper(paper_id):
    return json.loads((PAPERS / f"{paper_id}.json").read_text(encoding="utf-8"))


def ellipse(field):
    """Independent azimuth, ellipticity and power of any nonzero transverse field."""
    x, y = field
    power = abs(x)**2 + abs(y)**2
    return np.array([.5*np.arctan2(2*(x.conjugate()*y).real, abs(x)**2 - abs(y)**2),
                     .5*np.arcsin(2*(x.conjugate()*y).imag/power), power])


def test_qwp_stencil_and_independent_jones_perturbations():
    offset = .002
    np.testing.assert_allclose(ellipse(contract.qwp_jones(np.degrees(offset))),
                               [offset, -offset, 1.], atol=1e-13)
    for psi, epsilon in ((0., offset), (-offset, 0.), (.3, -.45)):
        field = contract.jones_from_angles(psi, epsilon)
        np.testing.assert_allclose(ellipse(field), [psi, epsilon, 1.], atol=1e-13)
        pol = contract.jones_to_circular(field)
        np.testing.assert_allclose(contract.circular_to_jones(pol[1], pol[-1]), field)
    with pytest.raises(ValueError, match="unit norm"):
        contract.jones_to_circular([2., 0.])


def test_stokes_rates_from_independent_cartesian_amplitude_perturbation():
    field = contract.jones_from_angles(.41, .23)
    derivative = np.array([-.11 + .32j, .2 - .46j])
    epsilon, dpsi, depsilon, alpha = contract.stokes_rates(field, derivative)
    step = 1e-5
    finite = (ellipse(field+step*derivative)-ellipse(field-step*derivative))/(2*step)
    assert epsilon == pytest.approx(.23)
    np.testing.assert_allclose([dpsi, depsilon, -alpha], finite, rtol=1e-8, atol=1e-9)
    with pytest.raises(FloatingPointError):
        contract.stokes_rates(field, [np.nan, 0.])


@pytest.mark.parametrize("sign", [-1, 1])
def test_circular_input_preserves_power_and_flags_angular_pole(sign):
    field = contract.jones_from_angles(.27, sign*np.pi/4)
    epsilon, dpsi, depsilon, alpha = contract.stokes_rates(field, -1.5*field)
    assert epsilon == pytest.approx(sign*np.pi/4)
    assert np.isnan([dpsi, depsilon]).all()
    assert alpha == pytest.approx(3.)
    pol = contract.jones_to_circular(field)
    assert abs(pol[sign]) == pytest.approx(1.)
    assert abs(pol[-sign]) < 1e-15


def test_velocity_shift_and_peak_intensity_unit_conversion():
    optical_hz = constants.C_LIGHT / 800e-9
    effective = contract.atom_frame_detuning_mhz(125., velocity_m_s=80., optical_frequency_hz=optical_hz)
    assert effective == pytest.approx(25.)
    peak = contract.gaussian_peak_intensity_mw_cm2(np.pi/2, intensity_radius_1e2_mm=1.)
    assert peak == pytest.approx(100.)
    assert contract.gaussian_peak_intensity_mw_cm2(1., intensity_radius_1e2_mm=None) is None


def test_literature_case_keeps_raw_conditions_planes_and_unique_experiment():
    original = paper("ries_2003_pra")
    case = contract.benchmark_case(original, "fig4_phase_scan")
    conditions = {q['quantity']: q for q in case['conditions']}
    assert conditions['analysis_frequency']['frequency_hz'] == 5e6
    assert conditions['analysis_frequency']['omega_rad_s'] == pytest.approx(2*np.pi*5e6)
    assert conditions['resolution_bandwidth']['bandwidth_hz'] == 1e6
    assert conditions['video_bandwidth']['bandwidth_hz'] == 100.
    assert 'maximum absorption' in conditions['optical_detuning']['basis']
    assert 'unspecified' in conditions['pump_beam_size']['basis']
    results = {q['quantity']: q for q in case['results']}
    assert results['squeezing_measured']['noise_db_sql'] == -.85
    assert results['squeezing_measured']['plane'] == 'detector'
    assert results['squeezing_loss_corrected']['noise_db_sql'] == -1.23
    assert results['squeezing_loss_corrected']['plane'] == 'cell_output_inferred'
    assert case['line'] == 'D2'
    assert '85Rb_residue' in case['isotope_composition']
    assert not any('noise_db_sql' in q for q in original['experiments'][0]['results'])
    with pytest.raises(ValueError, match="expected one experiment"):
        contract.benchmark_case(original, "absent_case")


def test_joint_scan_and_missing_beam_width_are_preserved():
    scan = contract.benchmark_case(paper("mikhailov_2009_jmo"), "fig3_density_power")
    assert scan['kind'] == 'sweep'
    assert [(q['quantity'], q['frequency_hz']) for q in scan['conditions']] == [('analysis_frequency', 1.2e6)]
    vortex = contract.benchmark_case(paper("zhang_2013_ol"), "vortex_optimum")
    beam = next(q for q in vortex['conditions'] if q['quantity'] == 'pump_beam_size')
    assert beam['value'] is None and beam['value_type'] == 'not_reported'


def test_reported_correction_remains_attached_to_its_author_claim_and_caveat():
    lima = paper('lima_2024_pra')
    case = contract.benchmark_case(lima, 'best_fig5')
    assert case['publication'] == lima['publication'] and case['paper_title'] == lima['title']
    assert case['paper_ai_notes'] == lima['ai_notes']
    assert any('81%' in note for note in case['paper_ai_notes'])
    assert case['notes'] == lima['experiments'][0]['notes']
    corrected = next(q for q in case['results'] if q['quantity'] == 'squeezing_detection_corrected')
    assert corrected['noise_db_sql'] == -6.4
    assert corrected['plane'] == 'inferred_cell_output_as_reported'


def test_digital_filter_bandwidth_preserves_its_meaning():
    case = contract.benchmark_case(paper('valente_2015_pra'), 'gaussianity_test')
    width = next(q for q in case['conditions'] if q['quantity'] == 'analysis_bandwidth')
    assert width['bandwidth_hz'] == 400_000.
    assert width['basis'] == 'digital band-pass width'
    assert 'omega_rad_s' not in width


def test_source_and_sensor_plane_noise_depths_keep_their_signs():
    li = contract.benchmark_case(paper('li_2022_josab'), 'amor_operating_point')
    sensor = next(q for q in li['results'] if q['quantity'] == 'squeezing_before_sensor')
    assert (sensor['noise_db_sql'], sensor['plane']) == (-1.9, 'magnetometer_input')
    horrom = contract.benchmark_case(paper('horrom_2013_jmo'), 'narrow_eit_fig3')
    before_filter = [q for q in horrom['results'] if q['plane'] == 'eit_input']
    assert {q['quantity']: q['noise_db_sql'] for q in before_filter} == {
        'squeezing_before_filter': -2., 'anti_squeezing_before_filter': 8.,
    }


@pytest.mark.parametrize("azimuth,ellipticity", [(-.6, -.23), (0., 0.), (.7, .31)])
def test_jones_polarization_angles_roundtrip_and_phase_invariance(azimuth, ellipticity):
    from gabes.psr_conventions import jones_polarization_angles, jones_from_angles

    for phase in (0., .8, -2.4):
        jones = np.exp(1j*phase)*jones_from_angles(azimuth, ellipticity)
        psi, epsilon = jones_polarization_angles(jones)
        assert psi == pytest.approx(azimuth, abs=1e-14)
        assert epsilon == pytest.approx(ellipticity, abs=1e-14)


def test_circular_polarization_has_undefined_azimuth():
    from gabes.psr_conventions import jones_polarization_angles, jones_from_angles

    for handedness in (-1, 1):
        psi, epsilon = jones_polarization_angles([1/np.sqrt(2), handedness*1j/np.sqrt(2)])
        assert np.isnan(psi)
        assert epsilon == pytest.approx(handedness*np.pi/4)
