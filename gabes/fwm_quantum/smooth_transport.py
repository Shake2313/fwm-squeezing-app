"""Continuous Gaussian pump for the existing moving reduced-Rb characteristic."""

import math
import numpy as np

from .model import reduced_pump_system
from .normalization import reduced_dipoles
from .transport import reduced_ballistic_problem
from ..quantum.contracts import readonly_array
from ..quantum.smooth_transport import smooth_wavepacket


def smooth_rb_problem(inputs,geometry,path,analysis_axis,*,entry_phase_rad=0.,
                       boundary_state=None,pump_center_xy_m=(0.,0.),convention='uniform-zeeman-rms'):
    """Reuse the fixed physical ports/validation; do not use its frozen H."""
    fixed=reduced_ballistic_problem(inputs,geometry,path,analysis_axis,segments=1,
        entry_phase_rad=entry_phase_rad,boundary_state=boundary_state,
        pump_center_xy_m=pump_center_xy_m,convention=convention)
    dipoles=reduced_dipoles(convention)
    dv=fixed['metadata']['one_photon_doppler_rad_s']
    h0,reservoirs=reduced_pump_system(0.,dv,transition_scales=dipoles.transition_scales)
    peak=dipoles.pump_rabi_rad_s(inputs.pump_power_W,inputs.pump_waist_m)
    hpeak,_=reduced_pump_system(peak,dv,transition_scales=dipoles.transition_scales)
    center=readonly_array(pump_center_xy_m,real=True)
    x,y=map(float,path.entry_position_m[:2]-center)
    vx,vy=map(float,path.velocity_m_s[:2])
    inverse_waist_sq=1/inputs.pump_waist_m**2
    def envelope(age):
        if not 0<=age<=path.residence_time_s:
            raise ValueError('age outside declared characteristic')
        # Exactly r_perp(a)=r_in,perp+v_perp*a; avoid rebuilding Cartesian
        # arrays at every ODE evaluation. No envelope or velocity averaging.
        return math.exp(-((x+vx*age)**2+(y+vy*age)**2)*inverse_waist_sq)
    metadata=dict(fixed['metadata'])
    metadata['pump_envelope']='continuous constant-waist transverse Gaussian'
    metadata['pump_peak_rabi_rad_s']=peak
    metadata['pump_amplitudes_rad_s']=[peak*envelope(t) for t in np.linspace(0.,path.residence_time_s,17)]
    metadata['pump_sample_ages_s']=np.linspace(0.,path.residence_time_s,17).tolist()
    return {'h0':h0,'h1':hpeak-h0,'envelope':envelope,'reservoirs':reservoirs,
        'boundary_state':fixed['boundary_state'],'duration_s':path.residence_time_s,
        'readouts':fixed['readouts'][0],'frequencies_rad_s':fixed['frequencies_rad_s'],
        'drives':fixed['drives'][0],'drive_frequencies_rad_s':fixed['drive_frequencies_rad_s'],
        'metadata':metadata}


def smooth_rb_wavepacket(inputs,geometry,path,analysis_axis,*,rtol=1e-10,atol=1e-13,
                         max_step_s=None,**kwargs):
    problem=smooth_rb_problem(inputs,geometry,path,analysis_axis,**kwargs)
    result=smooth_wavepacket(**{k:v for k,v in problem.items() if k!='metadata'},
        rtol=rtol,atol=atol,max_step_s=max_step_s)
    return {**result,'metadata':problem['metadata'],'analysis_axis':analysis_axis}
