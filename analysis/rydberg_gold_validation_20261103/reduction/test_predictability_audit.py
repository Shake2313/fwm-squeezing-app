"""Focused checks of reduction decisions, leakage boundaries and domain guards."""
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

SPEC = importlib.util.spec_from_file_location("predictability_audit", Path(__file__).with_name("predictability_audit.py"))
audit = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = audit
SPEC.loader.exec_module(audit)


def product_surface():
    # Nonlinear T dependence remains exactly separable at observed nodes.
    return {(t,p): np.exp(.002*(t-20)**2-.1*t)*(p/30)**-.3
            for t in [15.,20.,25.,30.,35.] for p in [20.,30.,40.]}


def test_nonlinear_but_separable_surface_reconstructs_non_cross_nodes():
    values=product_surface()
    model,decision=audit.select_from_pilot(values)
    error,rows=audit.evaluate(model,values,model.fitted_coordinates)
    assert decision['accepted']
    assert len(rows)==8
    assert error['max_absolute_relative_error']<1e-12
    # This tests a nonlinear factorization, not a claim that temperature is linear.
    assert not np.isclose(np.log(values[15.,30.])+np.log(values[25.,30.]),2*np.log(values[20.,30.]))
    diagnostics=audit.axis_diagnostics(values)
    assert diagnostics[0]['max_absolute_log_curvature']>0
    assert diagnostics[0]['decreasing_segment_count']>0


def test_cross_interaction_requires_direct_refinement():
    values={key:value*np.exp(.0006*(key[0]-20)*(key[1]-30))
            for key,value in product_surface().items()}
    _,decision=audit.select_from_pilot(values,tolerance=.01)
    assert not decision['accepted']
    assert decision['selected_action']=='branch_aware_refinement'
    assert decision['pilot_check']['max_absolute_relative_error']>.01


def test_surrogate_fit_does_not_use_non_cross_truth():
    values=product_surface()
    other=dict(values)
    other[35.,40.]*=100
    a=audit.CrossSlice.fit(values)
    b=audit.CrossSlice.fit(other)
    assert a==b
    # Validation must catch corruption even though it cannot leak into the fit.
    error,_=audit.evaluate(b,other,b.fitted_coordinates)
    assert error['max_absolute_relative_error']>.9


def test_no_extrapolation_and_invalid_physical_scalars_rejected():
    values=product_surface()
    model=audit.CrossSlice.fit(values)
    with pytest.raises(ValueError,match='extrapolation'):
        model.predict(36.,30.)
    with pytest.raises(ValueError,match='extrapolation'):
        model.predict(20.,19.)
    broken=dict(values);broken[20.,30.]=0
    with pytest.raises(ValueError,match='positive'):
        audit.CrossSlice.fit(broken)
    del broken[20.,30.]
    with pytest.raises(ValueError,match='incomplete'):
        audit.CrossSlice.fit(broken)


def test_shape_and_duplicates_fail_instead_of_silent_average():
    rows=[dict(temperature_c=t,coupling_power_mw=p,sensitivity_nv_cm_sqrt_hz=value)
          for (t,p),value in product_surface().items()]
    with pytest.raises(ValueError,match='Duplicate'):
        audit.surface(rows+[rows[0]])
    with pytest.raises(ValueError,match='Incomplete'):
        audit.surface(rows[:-1])


def test_discontinuous_controls_are_not_averaged_into_scalar_predictor():
    values=product_surface()
    rows=[];branches=[]
    for (t,p),value in values.items():
        label='central' if t<25 else 'sideband'
        rows.append(dict(temperature_c=t,coupling_power_mw=p,winning_branch=label))
        branches.append(dict(temperature_c=t,coupling_power_mw=p,
            central_sensitivity_nv_cm_sqrt_hz=value*(1.0001 if label=='sideband' else 1),
            sideband_sensitivity_nv_cm_sqrt_hz=value*(1.0001 if label=='central' else 1)))
    transitions=audit.branch_transitions(rows,branches)
    assert len(transitions)==3
    assert all(row['at_least_one_endpoint_near_tie'] for row in transitions)
    model=audit.CrossSlice.fit(values)
    assert isinstance(model.predict(25.,30.),float)
    assert not hasattr(model,'probe_power_uw')


def test_seed_guards_depend_on_axes_not_argmin_or_response():
    plan=audit.guard_plan(list(np.arange(15.,36.)),list(np.arange(20.,40.1,2.5)))
    assert plan['fit_count']==29
    assert plan['guard_count']==8
    assert plan['initial_count']==37
    assert not set(map(tuple,plan['fit_coordinates'])) & set(map(tuple,plan['guard_coordinates']))
    assert [15.,20.] in plan['guard_coordinates']
    assert [35.,40.] in plan['guard_coordinates']


def test_model_fingerprint_change_retires_old_reduction(tmp_path):
    path=tmp_path/'model.py';path.write_text('old')
    fingerprint={'model.py':audit.digest(path)}
    assert not audit.check_fingerprint(fingerprint,tmp_path)
    path.write_text('new physics')
    assert audit.check_fingerprint(fingerprint,tmp_path)==['model.py']


def test_parameter_change_retires_reduction_without_a_source_edit():
    params=dict(temp_c=20.,probe_power_uw=6.,coupling_power_mw=30.,lo_rabi_mhz=3.7,
        doppler='on',effective_opacity_participation=.001,coupling_rabi_mhz=3.,temperature_model='Linked')
    registry=dict(model_tag='old',fixed_model_settings_sha256=audit.fixed_settings_fingerprint(params),
        expected_critical_source_sha256={'engine':'same'},
        domain=dict(temperature_c=[15.,35.],probe_power_uw=[2.,10.],coupling_power_mw=[20.,40.],
                    lo_rabi_mhz=[2.,6.],opacity_scenarios=[.0005,.001,.0015]))
    assert not audit.reuse_rejections(registry,dict(params,temp_c=25.),{'engine':'same'},'old')
    assert 'fixed_physics_or_noise_settings_changed' in audit.reuse_rejections(
        registry,dict(params,coupling_rabi_mhz=3.5),{'engine':'same'},'old')
    assert 'untested_opacity_scenario' in audit.reuse_rejections(
        registry,dict(params,effective_opacity_participation=.0012),{'engine':'same'},'old')
    assert 'model_tag_changed' in audit.reuse_rejections(registry,params,{'engine':'same'},'new')
