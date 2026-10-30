"""Classical atomic PSR slopes, polarization conventions and applicability."""
import numpy as np
import pytest

from gabes import magneto_multilevel, psr


def _propagated_stokes(qwp_deg, chi_plus, chi_minus, k, distance):
    pol = magneto_multilevel._polarization(qwp_deg)
    plus = pol[1] * np.exp(.5j*k*chi_plus*distance)
    minus = pol[-1] * np.exp(.5j*k*chi_minus*distance)
    x, y = (plus+minus)/np.sqrt(2), 1j*(plus-minus)/np.sqrt(2)
    power = abs(x)**2 + abs(y)**2
    angle = .5*np.arctan2(2*np.real(x.conjugate()*y), abs(x)**2-abs(y)**2)
    ellipse = .5*np.arcsin(2*np.imag(x.conjugate()*y)/power)
    return np.array([angle, ellipse, power])


def test_local_stokes_rates_match_independent_jones_propagation():
    k, plus, minus, qwp = 8e6, 2e-5+3e-6j, -1e-5+1e-6j, 2.0
    epsilon, azimuth, ellipse, alpha = psr._stokes_rates(qwp, plus, minus, k)
    step = 1e-7
    left = _propagated_stokes(qwp, plus, minus, k, -step)
    right = _propagated_stokes(qwp, plus, minus, k, step)
    derivative = (right-left)/(2*step)
    zero = _propagated_stokes(qwp, plus, minus, k, 0.)
    assert epsilon == pytest.approx(np.deg2rad(-qwp), abs=1e-15)
    np.testing.assert_allclose([azimuth, ellipse, alpha],
                               [derivative[0], derivative[1], -derivative[2]/zero[2]], rtol=2e-8)


def test_common_susceptibility_attenuates_without_rotating_or_changing_ellipse():
    _, azimuth, ellipse, alpha = psr._stokes_rates(3., 1e-5+2e-6j, 1e-5+2e-6j, 8e6)
    assert abs(azimuth) < 1e-12 and abs(ellipse) < 1e-12
    assert alpha == pytest.approx(16.)


@pytest.fixture(scope="module")
def baseline():
    return psr.compute_psr_response([-1000., -950., -900.])


def test_zero_field_psr_is_nonzero_odd_in_ellipticity_and_within_declared_screens(baseline):
    assert np.all(baseline["numerical_valid"])
    assert np.all(baseline["approximation_valid"])
    assert np.all(baseline["passive_noise_valid"])
    assert np.all(np.abs(baseline["shear"]) > .1)
    rates = baseline["stencil"]["azimuth_rate_per_m"]
    np.testing.assert_allclose(rates[:, 1], -rates[:, 2], atol=1e-11, rtol=1e-8)
    np.testing.assert_allclose(rates[:, 3], -rates[:, 4], atol=1e-11, rtol=1e-8)
    np.testing.assert_allclose(rates[:, 0], 0., atol=1e-11)
    np.testing.assert_allclose(baseline["transmission"], np.exp(-baseline["optical_depth"]))
    k = 2*np.pi*baseline["frequency_hz"]/psr.constants.C_LIGHT
    np.testing.assert_allclose(baseline["alpha_per_m"], k*baseline["stencil"]["chi_probe"][:, 0].imag)
    assert "Langevin diffusion" in baseline["model_scope"]
    assert not any(key in baseline for key in ("squeezing_db", "covariance", "noise_variance"))


def test_zero_drive_has_no_psr_even_though_linear_absorption_remains():
    result = psr.compute_psr_response([-200., 200.], intensity_mw_cm2=0.,
                                      density_m3=1e16, cell_length_m=.01)
    np.testing.assert_allclose(result["shear"], 0., atol=1e-10)
    np.testing.assert_allclose(result["ellipticity_drift_length"], 0., atol=1e-10)
    assert np.all(result["optical_depth"] > 0)
    assert np.all(result["derivative_converged"])


def test_density_and_length_enter_local_and_integrated_quantities_once():
    common = dict(doppler=False, density_m3=1e16, cell_length_m=.01)
    original = psr.compute_psr_response([-200.], **common)
    density = psr.compute_psr_response([-200.], **dict(common, density_m3=2e16))
    length = psr.compute_psr_response([-200.], **dict(common, cell_length_m=.02))
    for key in ("g_per_m", "alpha_per_m", "ellipticity_drift_per_m", "shear", "optical_depth"):
        np.testing.assert_allclose(density[key], 2*original[key], rtol=1e-7, atol=1e-10)
    for key in ("g_per_m", "alpha_per_m", "ellipticity_drift_per_m"):
        np.testing.assert_allclose(length[key], original[key], rtol=1e-10, atol=1e-10)
    for key in ("shear", "optical_depth", "ellipticity_drift_length"):
        np.testing.assert_allclose(length[key], 2*original[key], rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(density["transmission"], original["transmission"]**2)
    vacuum = psr.compute_psr_response([-200.], **dict(common, density_m3=0.))
    np.testing.assert_array_equal(vacuum["shear"], [0.])
    np.testing.assert_array_equal(vacuum["optical_depth"], [0.])
    np.testing.assert_array_equal(vacuum["transmission"], [1.])
    assert vacuum["conditional_model_applicable"][0]


def test_ellipticity_half_step_refinement_and_qwp_angle_are_distinct(baseline):
    finer = psr.compute_psr_response(baseline["detunings_mhz"], ellipticity_step_deg=.025)
    np.testing.assert_allclose(finer["g_per_m"], baseline["g_per_m"], rtol=3e-6)
    assert np.all(finer["derivative_error"] < baseline["derivative_error"])
    np.testing.assert_allclose(baseline["stencil"]["ellipticity_rad"],
                               -np.deg2rad(baseline["stencil"]["qwp_deg"]), atol=1e-15)


@pytest.mark.parametrize("isotope,line,levels", [
    ("Rb87", "D1", (8, 8)), ("Rb87", "D2", (8, 16)),
    ("Rb85", "D1", (12, 12)), ("Rb85", "D2", (12, 24)),
])
def test_shared_atom_data_supports_each_rubidium_d_line(isotope, line, levels):
    result = psr.compute_psr_response([100.], isotope=isotope, line=line,
                                     doppler=False, density_m3=1e16, cell_length_m=.01)
    assert tuple(result["level_counts"].values()) == levels
    assert result["conditions"]["isotope"] == isotope and result["conditions"]["line"] == line
    assert result["numerical_valid"][0]
    assert result["alpha_per_m"][0] >= 0
    assert result["g_per_m"].shape == (1,)


def test_long_cell_is_flagged_without_clipping_estimated_psr(baseline):
    extended = psr.compute_psr_response(baseline["detunings_mhz"], cell_length_m=.75)
    np.testing.assert_allclose(extended["shear"], 10*baseline["shear"], rtol=1e-10)
    assert np.all(extended["numerical_valid"]) and not np.any(extended["approximation_valid"])
    assert all(any("uniform-pump" in item for item in reasons) for reasons in extended["approximation_reasons"])
    assert all(status == "outside declared approximation screens" for status in extended["approximation_status"])


def test_numerical_failure_is_separate_from_applicability(monkeypatch):
    def fail(*args, **kwargs):
        raise np.linalg.LinAlgError("singular synthetic solver")
    monkeypatch.setattr(magneto_multilevel, "multilevel_response", fail)
    result = psr.compute_psr_response([0.])
    assert not result["numerical_valid"][0]
    assert result["approximation_status"] == ["numerical failure"]
    assert result["approximation_reasons"] == [[]]
    assert "atomic solve failed" in result["numerical_reasons"][0][0]
    assert not result["passive_noise_valid"][0]


def test_negative_absorption_is_retained_but_cannot_enter_passive_noise_model(monkeypatch):
    original = magneto_multilevel.multilevel_response
    def active_fixture(*args, **kwargs):
        result = original(*args, **kwargs)
        for key in ("chi_plus", "chi_minus", "chi_probe"):
            result[key] = result[key].real - 1e-8j
        return result
    monkeypatch.setattr(magneto_multilevel, "multilevel_response", active_fixture)
    result = psr.compute_psr_response([-900.], doppler=False)
    assert result["numerical_valid"][0]
    assert result["alpha_per_m"][0] < 0 and result["transmission"][0] > 1
    assert not result["passive_noise_valid"][0] and not result["conditional_model_applicable"][0]
    assert "negative absorption" in result["passive_noise_reasons"][0][0]


@pytest.mark.parametrize("settings", [dict(intensity_mw_cm2=-1.), dict(cell_length_m=-1.),
    dict(ground_relaxation_hz=0.), dict(ellipticity_step_deg=0.), dict(ellipticity_step_deg=10.),
    dict(velocity_order=True), dict(density_m3=-1.), dict(isotope="Cs133")])
def test_invalid_physical_or_numerical_settings_are_rejected(settings):
    with pytest.raises(ValueError):
        psr.compute_psr_response([0.], **settings)
