"""One-time, manually checked foundation records; edit papers/*.json thereafter."""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
SOURCES = {s['paper_id']: s for s in json.loads((ROOT/'research_foundations_sources.json').read_text())}
META = json.loads((ROOT/'research_foundations_metadata.json').read_text(encoding='utf-8'))

def q(name, value, unit, page, basis='', plane='not_specified', kind='reported', **extra):
    return dict(quantity=name, value=value, unit=unit, value_type=kind, basis=basis,
                reference=None, plane=plane, source_id='main', locator=page,
                evidence='Reported for this configuration.', notes='', **extra)

def finding(fid, claim, page, category, implication):
    return dict(finding_id=fid, claim=claim, source_id='main', locator=page,
                category=category, model_implication='⟦AI⟧ '+implication)

def exp(eid, label, conditions, results=None, kind='operating_point', **extra):
    return dict(experiment_id=eid, label=label, kind=kind, conditions=conditions,
                results=results or [], notes=[], **extra)

def paper(pid, title, authors, year, journal, volume, pages, doi, arxiv, study, topics,
          species='87Rb', line='D1', inclusion='primary', config='single_pass_psr',
          measurement='orthogonal_vacuum_quadrature', seed='orthogonal_vacuum',
          reason='Direct PSR mechanism, noise or generation evidence.', source_pid=None, version='published'):
    s=dict(SOURCES[source_pid or pid]); s.pop('paper_id'); s.update(source_id='main',access='full_text',version=version,accessed_on='2026-10-30')
    return dict(schema_version=1,paper_id=pid,title=title,authors=authors,year=year,
        publication=dict(journal=journal,volume=volume,issue='',pages=pages,doi=doi,arxiv_id=arxiv,status='published'),
        legacy_column=None, study_type=study, topics=topics,
        scope=dict(species=species,line=line,configuration=config,measurement=measurement,
                   seed_mode=seed,inclusion=inclusion,reason=reason,isotope_composition='not_reported'),
        sources=[s],experiments=[],findings=[],ai_notes=[])

records=[]
p=paper('matsko_2002_pra','Vacuum squeezing in atomic media via self-rotation',
 ['A. B. Matsko','I. Novikova','G. R. Welch','D. Budker','D. F. Kimball','S. M. Rochester'],2002,
 'Physical Review A','66','043815','10.1103/PhysRevA.66.043815','quant-ph/0112072','theory',
 ['foundational_theory','classical_shear','absorption','hyperfine_structure','Doppler','buffer_gas','wall_coating'],line='D1/D2',version='arxiv_v2_2002-07-17')
for line,density,fig in [('D1',1e12,4),('D2',2e11,5)]:
 loc=f'PDF p. 8, Fig. {fig} and Sec. VI'
 p['experiments'].append(exp('model_'+line.lower(),f'{line} model example, Fig. {fig}',[
 q('pump_power',10,'mW',loc,plane='model_input'),q('pump_beam_size',.03,'cm',loc,'reported diameter; no radius convention inferred'),
 q('cell_length',10,'cm',loc,plane='model_medium'),q('atomic_density',density,'cm^-3',loc,plane='model_medium')],kind='model_case',line=line))
p['findings']=[finding('shear','Small-ellipticity rotation produces a phase-dependent vacuum variance through gL.',
 'PDF pp. 2-4, Eqs. (1), (12)-(17), (25)','theory','Retain the ellipticity and quadrature convention with the shear coefficient.'),
 finding('validity','The optimized Rb estimate uses large gL, small absorption and separately treated ground hyperfine manifolds.',
 'PDF pp. 7-8, Sec. VI','approximation','Do not treat the optimized prediction as measured data or as a finite-frequency atomic-noise calculation.')]
p['ai_notes']=['⟦AI⟧ Published metadata; evidence locators refer to the archived arXiv v2. Later noise-inclusive papers are required for comparison.','⟦AI⟧ Source PDF p. 4, Eq. (29) defines s = sqrt(Vvac/Vmin), while the subsequent text uses 10 log10(s). This author-specific dB convention is retained and is not imported as standard measured noise-power dB.'];records.append(p)

p=paper('ries_2003_pra','Experimental vacuum squeezing in rubidium vapor via self-rotation',
 ['J. Ries','B. Brezger','A. I. Lvovsky'],2003,'Physical Review A','68','025801','10.1103/PhysRevA.68.025801','quant-ph/0303109',
 'experiment',['early_demonstration','D2','homodyne','excess_noise','spatial_mode'],line='D2',version='arxiv_v1_2003-03-18')
p['scope']['isotope_composition']='nearly_pure_87Rb_with_85Rb_residue'
loc='PDF p. 3, Fig. 4 and Results and discussion'
p['experiments']=[exp('fig4_phase_scan','Phase scan at 5 MHz',[
 q('temperature',70,'degC',loc,plane='vapor_cell'),q('cell_length',75,'mm','PDF p. 2, Experimental apparatus'),
 q('pump_beam_size',600,'um','PDF p. 2, Experimental apparatus','focal diameter; intensity convention unspecified'),
 q('analysis_frequency',5,'MHz',loc),q('optical_detuning',.35,'GHz',loc,'zero is maximum absorption; no hyperfine-reference conversion'),
 q('resolution_bandwidth',1,'MHz',loc),q('video_bandwidth',100,'Hz',loc),
 q('detector_qe',91,'percent','PDF p. 3, detection paragraph','nominal diode QE')],[
 q('squeezing_measured',.85,'dB',loc,'positive reduction below SQL',plane='detector',uncertainty={'value':.05,'unit':'dB','kind':'reported','confidence':'not_reported'}),
 q('squeezing_loss_corrected',1.23,'dB',loc,'positive reduction; optical loss, diode inefficiency and electronic noise corrected',plane='cell_output_inferred',uncertainty={'value':.07,'unit':'dB','kind':'reported','confidence':'not_reported'})])]
p['findings']=[finding('nonminimum','Noise exceeds a minimum-uncertainty description; absorption alone cannot explain the excess.',
 'PDF p. 3, Results and discussion','noise','Include atomic excess noise separately from detector loss.'),
 finding('frequency','Best noise suppression occurs around 5-10 MHz; lower-frequency fluorescence degrades it.',
 'PDF pp. 3-4, Figs. 5-6','detection','Add analysis frequency before comparing a phase-only model to this measurement.')]
p['ai_notes']=['⟦AI⟧ The stated 35 mW is laser output, not an independently established cell-input operating-point power.'];records.append(p)

p=paper('josse_2003_job','Polarization squeezing in a four-level system',
 ['V. Josse','A. Dantan','A. Bramati','M. Pinard','E. Giacobino'],2003,'Journal of Optics B: Quantum and Semiclassical Optics','5','S513-S522',
 '10.1088/1464-4266/5/4/358','quant-ph/0302142','theory',['cavity','cross_kerr','atomic_noise','stokes'],
 species='four_level_X_model',line='model_transition',inclusion='context_only',config='cavity_X_system',
 reason='Cavity mechanism comparison; not a free-space warm-Rb PSR operating point.',source_pid='josse_2003_pra',version='arxiv_v1_2003-02-19')
p['sources'][0]['title_in_source']='Polarization squeezing in a 4-level system'
p['findings']=[finding('mechanism','Optical-pumping self-rotation adds low-frequency noise; crossed Kerr coupling enables higher-frequency vacuum squeezing.',
 'PDF pp. 12-16, Secs. 5.1-5.3','noise','Do not infer quantum-noise suppression from static rotation alone.'),
 finding('cavity','The calculation includes cavity feedback and polarization-switching thresholds.',
 'PDF pp. 5-7, Secs. 3.3-3.4; pp. 17-18, Sec. 6','scope','Use as a mechanism reference, not a numerical benchmark for single-pass vapor.')];records.append(p)

p=paper('lezama_2008_pra','Numerical investigation of the quantum fluctuations of optical fields transmitted through an atomic medium',
 ['A. Lezama','P. Valente','H. Failache','M. Martinelli','P. Nussenzveig'],2008,'Physical Review A','77','013806',
 '10.1103/PhysRevA.77.013806','0710.0609','theory',['heisenberg_langevin','diffusion_matrix','Zeeman','Doppler','noise_spectrum','entanglement'],
 species='angular_momentum_transition_models',line='model_transition',version='arxiv_v1_2007-10-02')
p['findings']=[finding('langevin','Two-polarization propagation includes Zeeman-resolved atomic Langevin diffusion and field spectra.',
 'PDF pp. 3-5, Sec. II, Eqs. (1)-(20); Appendix','noise_model','Reuse atomic operators while deriving a consistent frequency-dependent diffusion matrix.'),
 finding('doppler','Velocity classes enter both response and noise propagation.',
 'PDF pp. 5-6, Sec. II.A, Eqs. (21)-(29)','Doppler','Average the response and noise kernels before propagation; do not average squeezing in dB.'),
 finding('pump','Mean-field depletion is neglected despite propagation of quantum fluctuations.',
 'PDF p. 4, Sec. II, paragraph after Eqs. (5)-(6)','approximation','Add pump evolution separately when extending to appreciable absorption.')];records.append(p)

p=paper('mikhailov_2009_jmo','Vacuum squeezing via polarization self-rotation and excess noise in hot Rb vapors',
 ['Eugeniy E. Mikhailov','Arturo Lezama','Thomas W. Noel','Irina Novikova'],2009,'Journal of Modern Optics','56','1985-1992',
 '10.1080/09500340903159503','0903.3156','mixed',['atomic_excess_noise','excited_hyperfine','Doppler','frequency','density_power_map'],version='arxiv_v1_2009-03-18')
p['publication']['issue']='18-19';p['scope']['isotope_composition']='isotopically_enriched_87Rb'
p['experiments']=[exp('fig1_phase_scan','Figure 1(c) phase scan',[
 q('cell_length',75,'mm','PDF p. 2, Experimental Setup'),q('temperature',59,'degC','PDF p. 2, Fig. 1(c)'),
 q('atomic_density',2.3e11,'cm^-3','PDF p. 2, Fig. 1(c)'),q('analysis_frequency',1.36,'MHz','PDF p. 2, Fig. 1(c)'),
 q('resolution_bandwidth',100,'kHz','PDF p. 2, Fig. 1(c)'),q('video_bandwidth',30,'Hz','PDF p. 2, Fig. 1(c)'),
 q('pump_beam_size',.175,'mm','PDF p. 2, Experimental Setup','FWHM diameter',uncertainty={'value':.015,'unit':'mm','kind':'reported','confidence':'not_reported'})]),
 exp('fig3_density_power','Joint density/power sweep at F=2 to F-prime=2',[
 q('analysis_frequency',1.2,'MHz','PDF p. 3, Fig. 3')],kind='sweep')]
p['findings']=[finding('island','Squeezing occupies a small region in density and power; excess noise dominates outside it.',
 'PDF p. 3, Fig. 3 and discussion','noise','Use null-squeezing regions as constraints, alongside extrema.'),
 finding('hyperfine','Including both excited hyperfine levels changes the predicted squeezing qualitatively.',
 'PDF p. 4, Fig. 4 and Numerical simulations','atomic_structure','Keep off-resonant excited levels in response and Langevin calculations.'),
 finding('velocity','Far-detuned velocity groups substantially affect low-frequency noise.',
  'PDF pp. 4-5, Fig. 5 and Doppler discussion','Doppler','Check velocity convergence for noise independently of classical rotation.')]
p['ai_notes']=['⟦AI⟧ Figure 3 is retained as a joint sweep; no artificial Cartesian samples or digitized extrema.'];records.append(p)

p=paper('rochester_2001_pra','Self-rotation of resonant elliptically polarized light in collision-free rubidium vapor',
 ['S. M. Rochester','D. S. Hsiung','D. Budker','R. Y. Chiao','D. F. Kimball','V. V. Yashchuk'],2001,'Physical Review A','63','043814',
 '10.1103/PhysRevA.63.043814',None,'mixed',['classical_self_rotation','hyperfine_structure','optical_pumping','propagation','Doppler'],
 species='natural_Rb',line='D1/D2',measurement='classical_rotation_and_absorption',seed='elliptically_polarized_classical_field',
 reason='Classical atomic response foundation; no optical squeezing measurement.')
p['scope']['isotope_composition']='natural_abundance'
p['experiments']=[exp('classical_apparatus','Collision-free classical rotation apparatus',[
 q('cell_length',7.5,'cm','PDF p. 7, Sec. IV'),q('cell_diameter',2.5,'cm','PDF p. 7, Sec. IV'),
 q('residual_magnetic_field',.1,'mG','PDF p. 7, Sec. IV','approximate bound in every direction',kind='upper_bound')],kind='apparatus')]
p['findings']=[finding('hyperfine','Doppler broadening suppresses isolated-transition effects, allowing multiple-transition mechanisms to dominate.',
 'PDF pp. 4-5, Secs. II.D-II.E','atomic_structure','Retain excited hyperfine manifolds when validating the reused classical response.'),
 finding('slices','Significant absorption is handled by integrating intensity, rotation and ellipticity through thin slices.',
 'PDF pp. 6-7, Sec. III','propagation','An incremental extension can reuse the atomic solve along a longitudinal intensity trajectory.')];records.append(p)

p=paper('novikova_2015_pra','Excess optical quantum noise in atomic sensors',
 ['Irina Novikova','Eugeniy E. Mikhailov','Yanhong Xiao'],2015,'Physical Review A','91','051804(R)',
 '10.1103/PhysRevA.91.051804','1410.3810','mixed',['magnetometry','atomic_excess_noise','readout_quadrature','wall_coating'],
 inclusion='context_only',config='coated_cell_nmor_with_psr',measurement='polarization_rotation_noise_and_sensitivity',
 reason='Noise/readout tradeoff in a sensor; not an independent optimized PSR source benchmark.',version='arxiv_v1_2014-10-14')
p['scope']['isotope_composition']='enriched_87Rb'
p['experiments']=[exp('coated_cell','Paraffin-coated NMOR apparatus',[
 q('cell_length',75,'mm','PDF p. 2, apparatus paragraph'),q('cell_diameter',25,'mm','PDF p. 2, apparatus paragraph')],kind='apparatus')]
p['findings']=[finding('noise_vs_signal','Off-resonant coupling alters polarization noise even when its mean rotation contribution is negligible.',
 'PDF pp. 2-3, Figs. 2-4','noise','Fit or validate mean signal and covariance as distinct observables.'),
 finding('sensitivity','Rotating detection toward the squeezed quadrature reduces the magnetometric response as well.',
 'PDF pp. 3-4, Figs. 4-5','readout','Optimize signal-to-noise, rather than squeezing alone, for a future sensor mode.')];records.append(p)

p=paper('lanning_2017_pra','Gaussian-beam-propagation theory for nonlinear optics involving an analytical treatment of orbital-angular-momentum transfer',
 ['R. Nicholas Lanning','Zhihao Xiao','Mi Zhang','Irina Novikova','Eugeniy E. Mikhailov','Jonathan P. Dowling'],2017,'Physical Review A','96','013830',
 '10.1103/PhysRevA.96.013830','1702.01095','theory',['spatial_modes','paraxial_propagation','Born_approximation','OAM'],
 species='general_nonlinear_medium',line='not_applicable',inclusion='context_only',config='paraxial_nonlinear_propagation',
 measurement='classical_spatial_mode_structure',seed='varies',reason='Reusable propagation method; not a complete PSR quantum-noise theory.',version='arxiv_v2_2017-04-25')
p['sources'][0]['title_in_source']='Gaussian Beam-Propagation Theory for Nonlinear Optics - Featuring an Exact Treatment of Orbital Angular Momentum Transfer'
p['findings']=[finding('born','The Green-function mode expansion uses a first-order Born approximation and a weakly modified pump.',
 'PDF p. 2, Sec. II.A, Eqs. (2)-(3); p. 7, conclusion','spatial_modes','A candidate spatial propagation layer must retain its perturbative limit.'),
 finding('classical','The developed propagation treatment is classical and uses an idealized lossless, dispersionless medium.',
 'PDF p. 2, Sec. II.A','approximation','Do not use this alone as a Langevin covariance or strong-absorption model.')];records.append(p)

p=paper('wang_2026_lsa','Large-scale array of squeezed light and synchronization using atomic vapor',
 ['Lin Wang','Xichang Zhang','Konstantin Manannikov','Nir Davidson','Ying Hu','Dongdong Hao','Yanhong Xiao'],2026,'Light: Science & Applications','15','365',
 '10.1038/s41377-026-02428-w','2605.28316','experiment',['spatial_array','coherence_diffusion','wall_coating','synchronization','anti_squeezing'],
 inclusion='context_only',config='coated_cell_collective_coherence_array',measurement='polarization_quadrature_and_array_response',
 reason='Related coherence-enhanced polarization squeezing; requires atomic transport between beams beyond independent-velocity PSR.')
p['scope']['isotope_composition']='enriched_87Rb'
loc='PDF p. 9, Materials and methods'
p['experiments']=[exp('array30_fig2','30-beam array, Fig. 2',[
 q('channel_count',30,'1',loc),q('pump_power',2,'mW',loc,'per optical channel',plane='cell_input'),
 q('pump_beam_size',.505,'mm',loc,'1/e2 waist radius at cell center'),q('temperature',63.7,'degC',loc),
 q('cell_length',7.5,'cm',loc),q('analysis_frequency',180,'kHz','PDF p. 4, Fig. 2 caption'),
 q('resolution_bandwidth',10,'kHz','PDF p. 4, Fig. 2 caption')],[
 q('squeezing_measured',2.03,'dB','PDF pp. 3-4, Fig. 2; p. 1 abstract','positive reduction below SQL',plane='detector',uncertainty={'value':.02,'unit':'dB','kind':'reported','confidence':'not_reported'})])]
p['findings']=[finding('shared_coherence','Moving atoms share ground-state coherence between channels in a paraffin-coated cell.',
 'PDF pp. 2-3, principle; p. 9, methods','transport','Model inter-beam transport and wall memory before applying this array result.'),
 finding('scope','Squeezing-level synchronization is measured; direct dynamic phase-locking evidence is not established.',
  'PDF pp. 5-6, Synchronization section','scope','Keep phase locking and shared noise extrema separate in validation claims.')];records.append(p)

for p in records:
    (ROOT/'papers'/f"{p['paper_id']}.json").write_text(json.dumps(p,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(f'Wrote {len(records)} foundation and context papers.')
