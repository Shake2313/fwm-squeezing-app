"""Initial manually reviewed curation. Daily edits belong in papers/*.json.

Every value below was checked against the indicated archived PDF page. Source
discovery metadata is kept separately from the canonical scientific evidence.
"""
from __future__ import annotations

import html
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCES = json.loads((ROOT / "research_sources.json").read_text(encoding="utf-8"))
AS_OF = "2026-11-03"
POLICY = {
    "gold_reference_paper_id": "ju_2026_preprint",
    "status": "active_temporary",
    "origin": "direct_user_request",
    "adopted_on": AS_OF,
    "scope": "Main experimental geometry, transitions, operating conditions and readout reference.",
    "request": "Ju et al. (2026)을 당분간 Gold reference(실험 구도, 조건 등 메인으로 삼는 논문)으로 두겠습니다.",
}


def q(name, value, unit, page, basis="", *, kind="reported", uncertainty=None,
      reference=None, plane="not_specified", note=""):
    item = dict(quantity=name, value=value, unit=unit, value_type=kind,
                basis=basis, reference=reference, plane=plane, source_id="main",
                locator="PDF " + page, evidence=basis or name.replace("_", " "), notes=note)
    if uncertainty is not None:
        item["uncertainty"] = {"value": uncertainty, "unit": unit,
                               "basis": "Uncertainty as printed; confidence level only if stated in basis."}
    return item


def interval(name, low, high, unit, page, basis="", **kwargs):
    return dict(q(name, None, unit, page, basis, kind="range", **kwargs), lower=low, upper=high)


def missing(name, unit, page, note):
    return q(name, None, unit, page, note, kind="not_reported")


def exp(eid, label, conditions, results, *, kind="operating_point", notes=None, **overrides):
    return dict(experiment_id=eid, label=label, kind=kind, conditions=conditions,
                results=results, notes=notes or [], **overrides)


def finding(fid, claim, page, category):
    return dict(finding_id=fid, claim=claim, locator="PDF " + page, source_id="main",
                category=category, model_implication="")


def add(key, species, line, configuration, topics, experiments, findings, *,
        pid=None, study="experiment", inclusion="primary", reason="", group="", notes=None):
    source = SOURCES[key]
    cr, ax = source.get("crossref", {}), source.get("arxiv_metadata", {})
    date = cr.get("published", {}).get("date-parts", [[int(ax.get("citation_date", ["2000"])[0][:4])]])[0]
    authors = [" ".join(filter(None, (a.get("given"), a.get("family")))) for a in cr.get("author", [])]
    authors = authors or ax.get("citation_author", [])
    title = cr.get("title", ax.get("citation_title"))[0]
    url = source["source_url"]
    version = "arXiv preprint, accessed " + AS_OF if "arxiv.org" in url else "archived publisher/author manuscript"
    if key == "sedlacek_2012_natphys":
        version = "arXiv:1205.4461; preprint title Quantum Assisted Electrometry using Bright Atomic Resonances"
    if key == "ju_2026_preprint":
        version = "arXiv:2606.04354v1"
    if key == "schlossberger_2026_jscheme":
        version = "arXiv:2603.08678v1; publication metadata updated separately"
    if key == "yang_2024_pra":
        version = "published Physical Review A 109, 032609"
    paper = dict(schema_version=1, paper_id=pid or key, title=html.unescape(title),
        authors=authors, year=date[0], study_type=study, topics=topics,
        publication=dict(journal=html.unescape(next(iter(cr.get("container-title", [])), "arXiv")),
            volume=cr.get("volume", ""), issue=cr.get("issue", ""),
            pages=cr.get("page", cr.get("article-number", "")), doi=source.get("doi"),
            arxiv_id=source.get("arxiv_id"), status="published" if cr else "preprint",
            published_date="-".join(str(v).zfill(2) if i else str(v) for i, v in enumerate(date))),
        scope=dict(species=species, line=line, configuration=configuration,
            measurement="optical_EIT_readout_of_electric_field", seed_mode="classical_probe",
            inclusion=inclusion, reason=reason or "Rydberg-EIT 전기장 검출의 실험 조건·판독·잡음 근거."),
        apparatus_group=group,
        sources=[dict(source_id="main", url=url, path=source["path"], sha256=source["sha256"],
            pages=source["pages"], access="full_text", version=version, accessed_on=AS_OF,
            access_method=source.get("access_method", "public_pdf_download"),
            access_note="Numerical fields and findings manually checked at their PDF locators; bibliography may describe the later published edition.")],
        experiments=experiments, findings=findings, ai_notes=notes or [])
    if key == "ju_2026_preprint":
        paper["scope"]["reference_role"] = "gold_reference"
        paper["reference_designation"] = POLICY
    (ROOT / "papers" / (paper["paper_id"] + ".json")).write_text(
        json.dumps(paper, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


add("mohapatra_2007_prl", "85Rb addressed; isotope purity not specified", "5S-5P-nS/nD", "two_photon_EIT",
    ["foundation", "optical_depth", "dc_screening", "transient_response"], [
    exp("apparatus", "초기 열증기 Rydberg-EIT 장치", [
        q("cell_length", 75, "mm", "p. 1", "Physical vapor-cell length"),
        q("probe_power", 1, "uW", "p. 1"),
        q("probe_beam_radius", 0.4, "mm", "p. 1", "1/e^2 intensity radius"),
        q("coupling_beam_radius", 0.8, "mm", "p. 1", "1/e^2 intensity radius"),
        q("coupling_power", 200, "mW", "p. 1", "Available power up to this value", kind="upper_bound"),
        q("temperature", "room_temperature", "1", "p. 1", "No numerical cell temperature assigned"),
    ], [interval("eit_linewidth", 2, 4, "MHz", "p. 2", "Varies with laser power and transition")], kind="apparatus"),
    exp("charge_response", "벽 전하·DC 차폐 관측", [], [
        q("applied_dc_field_without_observed_response", 100, "V/cm", "p. 3", "External DC field screened in this cell", kind="upper_bound"),
        q("rf_frequency", 90, "MHz", "p. 4, Fig. 6", "AC Stark response test; not an adjacent-Rydberg resonant LO"),
    ], kind="scan", notes=["가한 DC 장의 상한과 최소 검출장·감도를 구분."]),
], [finding("screening", "벽 전하가 정적 외부 전기장을 차폐하며, 전기장 반전 뒤 과도 EIT 응답을 관찰.", "pp. 3-4, Figs. 5-6", "screening"),
    finding("doppler", "속도 평균한 ladder 감수율로 probe transmission을 계산. EIT 위치는 probe·coupling detuning에 의존.", "pp. 2-3, Eq. (1)", "optical_depth")],
    inclusion="context_only", reason="전기장 감응 Rydberg-EIT의 광학 검출·벽 전하 기초 실험.", group="Durham_Adams")

add("sedlacek_2012_natphys", "87Rb", "5S(F2)-5P(F3)-53D5/2-54P3/2", "bright_resonance_amplitude_readout",
    ["foundation", "AT_splitting", "weak_field", "Doppler"], [
    exp("weak_field", "밝은 공명에 의한 약한 MW 검출", [
        q("rf_frequency", 14.233, "GHz", "p. 5", "Observed microwave transition frequency"),
        q("cell_length", 10, "cm", "p. 8, Methods", "Physical cell length"),
        q("effective_interaction_length", 7.5, "cm", "p. 8, Methods", "Focused coupling limits the interaction length"),
        q("temperature", "room_temperature", "1", "p. 8, Methods"),
        q("probe_rabi_over_2pi", 6, "MHz", "p. 8, Methods", "Average Rabi frequency"),
        q("coupling_rabi_over_2pi", 2, "MHz", "p. 8, Methods", "Average Rabi frequency"),
    ], [q("field_sensitivity", 30, "uV/cm/sqrt(Hz)", "p. 1, Abstract", "Optical amplitude readout; laser spectral bandwidth limited", kind="approximate"),
        q("minimum_detectable_field", 8, "uV/cm", "p. 1, Abstract", "Observed field amplitude; not Hz-normalized", kind="approximate")]),
], [finding("bright", "마이크로파가 EIT 창 안의 밝은 공명을 바꾸며, 약한 장에서는 광학 진폭으로 전기장을 검출.", "pp. 5-6, Figs. 2-3", "readout"),
    finding("axis", "Probe 주파수 스캔에서 AT 분리에는 광학 파장비·Doppler 조건이 작용.", "p. 5", "frequency_axis")], group="Oklahoma_Shaffer",
    notes=["⟦AI⟧ 출판 논문 서지와 제목이 다른 공개 preprint를 원문 근거로 사용. 판본 출처 보존."])

add("sedlacek_2013_prl", "87Rb", "5S(F2)-5P(F3)-53D5/2-54P3/2", "polarization_vector_readout",
    ["polarization", "Zeeman", "AT_splitting"], [
    exp("vector", "편광 회전 기반 벡터 MW 측정", [
        q("cell_dimensions", "10 x 10 x 30", "mm", "p. 3", "Cuboidal cell dimensions; propagation length not inferred"),
        q("temperature", 45, "degC", "p. 4"),
        q("rf_frequency", 14.233, "GHz", "p. 4"),
        q("probe_power", 15, "uW", "p. 4"), q("coupling_power", 11, "mW", "p. 4"),
        q("probe_beam_size", 200, "um", "p. 4", "Spot size; radius/diameter definition not specified in this passage"),
        q("coupling_beam_size", 65, "um", "p. 4", "Spot size; radius/diameter definition not specified in this passage"),
    ], [q("polarization_angle_resolution", 0.5, "degree", "p. 4", "Resolution in both zeta and phi", kind="approximate"),
        q("at_splitting", 39.36, "MHz", "p. 4", "Probe-axis peak splitting", uncertainty=0.06)])
], [finding("manifold", "편광에 따른 3준위·4준위 경로 비율로 전기장 방향을 추정. 진폭 예측에는 52상태·광펌핑·Doppler 계산 사용.", "pp. 2-4, Figs. 2-4", "multilevel")], group="Oklahoma_Shaffer")

add("holloway_2014_ieee", "85Rb", "5S-5P-nD-(n+1)P", "AT_field_calibration",
    ["SI_traceability", "AT_splitting", "field_mapping", "frequency_axis"], [
    exp("17ghz", "17.04 GHz AT 보정 운전점", [
        q("cell_length", 75, "mm", "p. 6"), q("cell_diameter", 25, "mm", "p. 6"),
        q("rf_frequency", 17.04, "GHz", "p. 7, Fig. 12"),
        q("probe_power", 175, "nW", "p. 7", "Input power for this frequency", plane="cell_input"),
        q("coupling_power", 30, "mW", "p. 7"),
    ], [], notes=["p. 6의 22 mW 예시와 p. 7의 30 mW 측정점을 합치지 않음."]),
    exp("68ghz", "68.64 GHz 별도 운전점", [
        q("rf_frequency", 68.64, "GHz", "p. 8, Fig. 15"),
        q("probe_power", 250, "nW", "p. 8", "Input power for this frequency", plane="cell_input"),
        q("coupling_power", 24, "mW", "p. 8"),
    ], []),
], [finding("traceability", "AT 주파수 분리와 전이 쌍극자를 연결해 원자 위치의 RF 장을 보정. 여러 Rydberg 전이로 주파수 범위를 확장.", "pp. 4-8", "calibration"),
    finding("local_field", "셀 및 주변 물질의 RF 장 변화는 원자 위치의 장과 외부 입사장을 구분해 검토.", "pp. 9-10", "field_plane")], group="NIST_Raithel")

add("simons_2016_apl", "85Rb", "5S-5P-23D/24D-22F/23F", "detuned_AT_extrapolation",
    ["AT_splitting", "rf_detuning", "SI_traceability"], [
    exp("detuning_scan", "182·208 GHz RF detuning 측정", [
        q("probe_power", 30, "uW", "p. 4"), q("coupling_power", 53, "mW", "p. 4"),
        q("probe_beam_fwhm", 75, "um", "p. 4", "FWHM, not 1/e^2 radius"),
        q("coupling_beam_fwhm", 140, "um", "p. 4", "FWHM"),
    ], [q("detection_improvement_factor", 2, "1", "p. 5", "At least about twofold lower measurable field through detuned fits", kind="approximate")],
        kind="scan", notes=["208.30 GHz와 182.149 GHz는 별도 전이 측정. Hz 정규화 감도 값으로 해석하지 않음."]),
    exp("208ghz", "208.30 GHz 전이", [q("rf_frequency", 208.30, "GHz", "p. 4", "23D5/2-22F7/2")], []),
    exp("182ghz", "182.149 GHz 전이", [q("rf_frequency", 182.149, "GHz", "p. 4", "24D5/2-23F7/2")], []),
], [finding("generalized_rabi", "RF detuning에 따른 일반화 Rabi 분리를 맞춰 공명 분리를 외삽하면 직접 분리되지 않는 약한 장도 추정 가능.", "pp. 4-5, Eq. (2), Table I", "calibration")], group="NIST_Raithel")

add("kumar_2017_scirep", "133Cs", "6S-6P-52D5/2-53P3/2", "optical_homodyne",
    ["homodyne", "photon_shot_noise", "transit", "collisions", "power_broadening"], [
    exp("sensitivity", "Mach-Zehnder 광학 homodyne 검출", [
        q("cell_length", 4, "cm", "p. 2"),
        q("temperature", "room_temperature", "1", "p. 6", "Room-temperature sensitivity data"),
        q("beam_diameter_for_noise_estimate", 1, "mm", "p. 6", "Interaction volume used in atomic-noise estimate; distinct from beam-size scans"),
    ], [q("field_sensitivity", 5, "uV/cm/sqrt(Hz)", "p. 6", "Measured homodyne readout sensitivity", kind="approximate"),
        q("sql_field_sensitivity", 9, "nV/cm/sqrt(Hz)", "p. 6", "Calculated atom projection-noise estimate; not measured optical sensitivity", kind="approximate")]),
], [finding("tradeoff", "Probe 세기·밀도·빔 크기 변화에서 power·collision·transit broadening을 분리 비교.", "pp. 3-6, Figs. 3-5", "broadening"),
    finding("readout_noise", "MZI로 probe 기술잡음을 줄인 뒤 광자 샷노이즈가 판독 감도를 제한.", "pp. 6-8", "photon_shot_noise")], group="Oklahoma_Shaffer")

add("kumar_2017_oe", "133Cs", "6S-6P-52D5/2-53P3/2", "frequency_modulation_spectroscopy",
    ["frequency_modulation", "residual_amplitude_modulation", "photon_shot_noise"], [
    exp("fm", "RAM 보정 FM 분광 검출", [q("rf_frequency", 5.047, "GHz", "p. 3, Fig. 1", "Approximate Rydberg RF transition", kind="approximate")], [
        q("field_sensitivity", 3, "uV/cm/sqrt(Hz)", "p. 1, Abstract; p. 11", "FM readout; abstract supplies Hz^-1/2 convention", kind="approximate"),
        q("minimum_detectable_field", 1.8, "uV/cm", "p. 11", "With matched filtering; field amplitude", kind="approximate"),
    ], notes=["본문 p. 11의 Hz 지수 인쇄와 초록 정의를 함께 기록. matched-filter 최소장과 감도 구분."]),
], [finding("ram", "FM spectroscopy에 능동 RAM 보정을 더해 광학 기술잡음을 낮춤.", "pp. 4-8", "readout"),
    finding("filter", "Matched filtering은 약한 신호 가시성을 높이나 최저장 SNR은 광자 샷노이즈로 제한.", "pp. 10-12", "noise")], group="Oklahoma_Shaffer")

add("simons_2019_apl", "133Cs", "6S-6P-34D5/2-35P3/2", "rf_mixer_phase_readout",
    ["heterodyne", "phase", "RF_LO", "IF"], [
    exp("mixer", "RF mixer 위상 검출 장치", [
        q("rf_frequency", 19.626, "GHz", "p. 2"),
        q("probe_power", 49, "uW", "p. 2"), q("coupling_power", 60.6, "mW", "p. 2"),
        q("probe_beam_fwhm", 425, "um", "p. 2", "FWHM"),
        q("coupling_beam_fwhm", 620, "um", "p. 2", "FWHM"),
    ], []),
], [finding("phase", "LO와 SIG의 간섭을 원자가 혼합해 IF probe 변조를 만들며, 이 위상으로 RF 위상·전파거리를 읽음.", "pp. 2-3, Figs. 2-4", "heterodyne")], group="NIST_Boulder_mixer")

add("gordon_2019_aipadv", "133Cs", "6S-6P-34D5/2-35P3/2", "rf_mixer_weak_signal",
    ["heterodyne", "weak_field", "frequency_resolution", "IF"], [
    exp("weak_mixer", "90 kHz IF 약한 RF 신호 검출", [
        q("lo_frequency", 19.626000, "GHz", "p. 2, Abstract"),
        q("signal_frequency", 19.626090, "GHz", "p. 2, Abstract"),
        q("intermediate_frequency", 90, "kHz", "p. 2, Abstract"),
    ], [q("field_sensitivity", 79, "uV/m/sqrt(Hz)", "p. 2, Abstract", "Published-edition sensitivity", kind="approximate"),
        q("minimum_detectable_field", 46, "uV/m", "p. 2, Abstract", "Measured local field amplitude", uncertainty=2, kind="approximate"),
        q("frequency_discrimination", 0.1, "Hz", "p. 2, Abstract", "Equal-strength neighboring signals resolved")]),
], [finding("resolution", "RF carrier와 IF의 역할을 분리. 좁은 lock-in 대역으로 서로 가까운 신호를 판별.", "pp. 2-5", "frequency_resolution")],
    group="NIST_Boulder_mixer", notes=["⟦AI⟧ 출판본 PDF 초록의 79 µV/m/√Hz를 채택. 다른 preprint·후속 표의 값으로 덮어쓰지 않음."])

add("jing_2020_natphys", "133Cs", "6S-6P-47D5/2-48P3/2", "atomic_superheterodyne",
    ["heterodyne", "weak_field", "phase", "frequency_resolution", "RF_LO"], [
    exp("superhet", "Atomic superhet 감도 운전점", [
        q("intermediate_frequency", 150, "kHz", "p. 2, Fig. 2", "LO-SIG offset used for the comparison"),
        q("resolution_bandwidth", 1, "Hz", "p. 2, Fig. 2", "Spectrum analyzer RBW"),
        q("temperature", "room_temperature", "1", "p. 3"),
    ], [q("field_sensitivity", 55, "nV/cm/sqrt(Hz)", "p. 1", "Measured superheterodyne sensitivity"),
        q("minimum_detectable_field", 780, "pV/cm", "p. 1", "Long-integration minimum field; not 1-Hz sensitivity")]),
], [finding("lo_optimum", "LO dressing으로 약한 SIG에 선형 응답. 최적 LO는 EIT 폭과 관련되며 무제한 증가하지 않음.", "pp. 2-3, Fig. 2", "RF_LO"),
    finding("phase_frequency", "동일 장치가 RF 위상과 주파수도 검출. 감도·최소 검출장·주파수 정밀도는 서로 다른 지표.", "pp. 3-4, Figs. 3-4", "metric_definition")], group="Shanxi_Jing_Zhang")

add("prajapati_2021_apl", "85Rb", "5S-5P-nD; D1_ground_repump", "repumped_atomic_mixer",
    ["repumping", "heterodyne", "photon_shot_noise", "population"], [
    exp("repump_off", "동일 장치 repump off", [
        q("probe_power", 550, "uW", "p. 4"), q("coupling_power", 70, "mW", "p. 4"),
        q("repump_power", 0, "mW", "p. 4", "Repump absent"),
        q("lockin_time_constant", 1, "s", "p. 4"),
    ], [q("field_sensitivity", 5, "uV/m/sqrt(Hz)", "p. 4", "Authors equate their 1-s acquisition field threshold to sensitivity", kind="approximate")]),
    exp("repump_on", "동일 장치 repump on", [
        q("probe_power", 550, "uW", "p. 4"), q("coupling_power", 70, "mW", "p. 4"),
        q("repump_power", 10, "mW", "p. 3", "D1 repump before entering cell", plane="cell_input"),
    ], [q("field_sensitivity", 3, "uV/m/sqrt(Hz)", "p. 4", "Authors' 1-s acquisition convention", kind="approximate")]),
    exp("noise_budget", "광자 샷노이즈 판독 조건", [
        q("probe_power", 560, "uW", "p. 5", "Noise measurement input", plane="cell_input"),
        q("detected_probe_power", 220, "uW", "p. 5", "Power on the detector", plane="photodetector"),
    ], [q("shot_noise_dominated_frequency", 300, "kHz", "p. 5", "Above this approximate frequency the traces reach photon shot noise", kind="approximate")], kind="noise_characterization"),
], [finding("population", "D1 repump가 EIT 신호를 키우며 동일 정도의 power broadening 없이 판독 응답을 향상.", "pp. 2-4", "repumping"),
    finding("noise", "낮은 판독 주파수에서는 1/f laser noise, 높은 주파수에서는 광자 샷노이즈를 확인.", "p. 5, Fig. 6", "noise")], group="NIST_Boulder_repump")

add("meyer_2021_pra", "generic_alkali; 85Rb example", "three_level_ladder_EIT", "EIT_quantum_sensitivity_theory",
    ["theory", "optical_depth", "SQL", "Doppler", "propagation"], [], [
    finding("od_definition", "식 (18)의 OD는 공명·small-signal 기준의 원자수·쌍극자·길이 정규화. 열증기 Doppler 평균은 별도로 취급.", "pp. 4-6, Eqs. (18)-(20)", "optical_depth"),
    finding("finite_od", "큰 OD에서 probe 감쇠와 길이 방향 Rabi 변화가 감도 최적점에 영향. Optically thin 결과의 적용 범위 구분.", "pp. 4-6, Secs. III.C-III.D", "propagation"),
], study="theory", inclusion="context_only", reason="실험 감도 비교에 필요한 OD·Doppler·광학 판독 SQL 이론.", group="ARL_Meyer")

add("cai_2022_photonics", "85Rb", "5S-5P-70S1/2-70P3/2", "detuned_atomic_superheterodyne",
    ["heterodyne", "rf_detuning", "sensitivity_definition", "noise_spectrum"], [
    exp("zero_detuning_signal", "공명 LO 직접 응답 평가", [
        q("rf_frequency", 10.68, "GHz", "p. 5", "70S-70P transition", kind="approximate"),
        q("rf_detuning", 0, "MHz", "p. 12"),
    ], [q("field_sensitivity", 49.84, "nV/cm/sqrt(Hz)", "p. 12", "Direct field modulation calibration", uncertainty=0.72)]),
    exp("minus6_signal", "−6 MHz LO 직접 응답 평가", [q("rf_detuning", -6, "MHz", "p. 12", "RF detuning, not optical probe detuning")], [
        q("field_sensitivity", 35.75, "nV/cm/sqrt(Hz)", "p. 12", "Direct field modulation calibration", uncertainty=0.51)]),
    exp("minus6_noise", "−6 MHz 잡음 스펙트럼 평가", [
        q("rf_detuning", -6, "MHz", "p. 12"), q("acquisition_time", 2, "s", "p. 12"),
    ], [q("field_sensitivity", 12.50, "nV/cm/sqrt(Hz)", "p. 12", "Noise-spectrum fit and RF-optical gain calibration", uncertainty=0.04),
        q("minimum_detectable_field", 176.78, "pV/cm", "p. 12", "5000-s minimum discernible field", uncertainty=0.57)],
        kind="noise_characterization", notes=["35.75와 12.50은 평가 방법이 다름. 동일 결과로 합치지 않음."]),
], [finding("detuning", "비공명 RF dressing에서 응답 기울기 증가를 관측. 직접 변조와 잡음 스펙트럼 감도 수치는 평가 방식별로 다름.", "pp. 10-12, Figs. 9-10", "sensitivity_definition")], group="WIPM_Liu")

add("su_2022_oe", "87Rb", "5S(F2)-5P(F3)-33D3/2,5/2", "thermal_EIT_optimization",
    ["optical_depth", "temperature", "polarization", "Zeeman", "power_broadening"], [
    exp("temperature_power_scan", "열증기 OD·편광·probe 세기 최적화", [
        interval("temperature", 27, 65, "degC", "p. 2", "Experimental cell-temperature range"),
        q("coupling_power", 27, "mW", "p. 2", "Largest available coupling power", kind="upper_bound"),
        q("beam_diameter", 0.81, "mm", "p. 2", "Full width at e^-2 maximum for both optical beams", kind="approximate"),
    ], [q("eit_peak_height", 13, "percent", "p. 1, Abstract", "Optimized transmission peak height; not an absolute field sensitivity")], kind="scan"),
], [finding("od_optimization", "Probe 세기와 열증기 광학밀도를 함께 최적화하면 EIT contrast가 증가. 자기장·편광 조건의 영향도 측정.", "pp. 1-3, 8-10", "optical_depth")],
    inclusion="context_only", reason="Gold 실험의 온도·OD·EIT contrast 비교용 광학 최적화 실험.", group="NCCU_Chen")

add("berweger_2022_pra", "133Cs", "6S-6P-Rydberg_ladder", "optical_frequency_comb_self_heterodyne",
    ["frequency_comb", "AT_splitting", "optical_heterodyne", "transient_response"], [
    exp("comb", "Probe comb 동시 스펙트럼 판독", [
        q("probe_wavelength", 852, "nm", "p. 1"), q("coupling_wavelength", 509, "nm", "p. 1"),
    ], [q("eit_linewidth", 5, "MHz", "p. 1, Abstract", "Resolved EIT linewidth below 5 MHz", kind="upper_bound"),
        q("field_sensitivity", 2.3, "uV/cm/sqrt(Hz)", "p. 1, Abstract", "Self-heterodyne optical comb readout"),
        q("minimum_detectable_field", 66, "uV/cm", "p. 1, Abstract", "Measured RF field amplitude")], kind="apparatus"),
], [finding("comb_parallel", "여러 comb tooth의 투과를 self-heterodyne으로 동시에 측정해 레이저 스캔을 줄이고 RF pulse에 의한 AT 변화를 판독.", "pp. 1, 5-10", "readout")],
    pid="dixon_2023_prapplied", group="Quantum_Valley_Shaffer")

add("wang_2023_oe", "133Cs", "6S-6P-47D5/2-48P3/2", "flat_top_atomic_superheterodyne_noise",
    ["noise_spectrum", "transit", "beam_profile", "atom_number", "SQL"], [
    exp("flat_top_noise", "Flat-top 빔 직경별 잡음 비교", [
        q("beam_diameter_quantum_noise_regime", 2, "mm", "p. 1, Abstract", "Quantum-noise-dominated regime for diameter <=2 mm", kind="upper_bound"),
        q("readout_frequency_quantum_noise_regime", 70, "kHz", "p. 1, Abstract", "Above this frequency in the stated diameter regime", kind="lower_bound"),
        q("lo_rabi_over_2pi", 10.68, "MHz", "p. 2", "Local microwave Rabi frequency", kind="approximate"),
    ], [], kind="noise_characterization"),
], [finding("participation", "빔에 참여하는 모든 원자는 잡음에 기여하지만 신호 기여 원자는 일부. 양자잡음 지배 관측과 이상적 SQL 달성은 별도.", "pp. 1, 5-7", "atom_number"),
    finding("transit", "Flat-top 직경과 판독 주파수를 바꾸어 transit 및 기술잡음의 영향을 구분.", "pp. 3-7", "transit")], group="Shanxi_Jing_Zhang")

add("zhang_2024_quantum_scaling", "133Cs", "6S-6P-47D5/2-48P3/2", "interaction_length_scaling",
    ["atom_number", "SQL", "transit", "heterodyne"], [], [
    finding("scaling", "상호작용 길이를 바꾸는 실험에서 신호 진폭의 N 비례·잡음 진폭의 sqrt(N) 비례 및 감도의 N^-1/2 경향을 조사.", "pp. 1, 5-9", "atom_number"),
    finding("nonideal", "유한 빔·비균일 세기·transit noise가 이상적 원자수 스케일링에 주는 영향을 논의.", "pp. 6-10", "transit")],
    pid="zhang_2023_epjqt", study="mixed", group="Shanxi_Jing_Zhang",
    notes=["⟦AI⟧ 길이 스캔의 그림 점은 수치 digitization 없이 finding·원문 위치로 연결. 독립 스칼라 감도 표본 생성 없음."])

add("yang_2024_pra", "87Rb model", "four_level_Rydberg_ladder", "phase_diffusion_noise_theory",
    ["theory", "laser_frequency_noise", "Doppler", "temperature", "noise_spectrum"], [
    exp("theory_parameters", "잡음 이론 수치 계산 기준", [
        q("temperature", 300, "K", "p. 4", "Default numerical calculation temperature"),
        q("transit_relaxation_over_2pi", 500, "kHz", "p. 3, Eq. (6)", "Model transit relaxation"),
        q("excited_decay_over_2pi", 6.07, "MHz", "p. 3, Eq. (6)", "Model intermediate-state spontaneous rate"),
        q("rydberg_decay_over_2pi", 10, "kHz", "p. 3, Eq. (6)", "Model Rydberg spontaneous rate"),
    ], [], kind="theory_prediction"),
], [finding("phase_diffusion", "Laser phase diffusion으로 원자 coherence 잡음 PSD와 감도 한계를 계산. Laser linewidth·Doppler·온도 의존을 포함.", "pp. 2-8", "laser_frequency_noise")],
    study="theory", inclusion="context_only", reason="실제 laser frequency noise·Doppler 감도 한계 모델 비교.", group="SIOM_Yang_Cheng")

add("tu_2024_sciadv", "87Rb laser_cooled", "5S-5P-Rydberg_MW_ladder", "cold_atom_heterodyne",
    ["cold_atoms", "SQL", "noise_budget", "duty_cycle"], [
    exp("cold_sensor", "냉각 원자 SQL 접근 실험", [
        q("temperature", 200, "uK", "p. 3", "Atomic-cloud temperature", kind="approximate"),
        q("atom_number", 5.2e5, "1", "p. 1", "Laser-cooled ensemble", kind="approximate"),
        q("repetition_rate", 100, "Hz", "pp. 1, 3"),
        q("loading_time", 7.3, "ms", "p. 3"), q("detection_window", 2.7, "ms", "p. 3"),
    ], [q("field_sensitivity", 10.0, "nV/cm/sqrt(Hz)", "p. 1, Abstract", "Duty-cycled cold-atom measurement"),
        q("sensitivity_to_sql_ratio", 2.6, "1", "p. 1", "Measured sensitivity above theoretical SQL"),
        q("minimum_detectable_field", 540, "pV/cm", "p. 1", "Long-integration minimum field")]),
], [finding("cold_vs_vapor", "냉각·가동률·유효 원자수를 포함한 잡음 예산으로 SQL과 비교. 열증기 연속 판독과 장치 조건이 다름.", "pp. 1-5", "SQL")],
    inclusion="context_only", reason="열증기 Gold와 구분한 냉각 원자 SQL·가동률 비교 실험.", group="SCNU_Tu_Yan")

add("prajapati_2024_bandwidth", "133Cs", "6S-6P-Rydberg_ladder", "probe_comb_bandwidth",
    ["bandwidth", "frequency_comb", "heterodyne", "IF"], [
    exp("comb_bandwidth", "Probe comb 판독 대역폭", [
        q("lo_frequency", 4.0, "GHz", "p. 2, Fig. 1"),
        q("comb_tooth_spacing", 1, "MHz", "p. 2, Fig. 1"),
        q("comb_optical_span", 200, "MHz", "p. 2, Fig. 1", "Optical probe-comb span; distinct from RF instantaneous bandwidth"),
        q("acquisition_time", 100, "ms", "p. 3"),
    ], [q("instantaneous_bandwidth", 12, "MHz", "p. 1; p. 3", "RF modulation sensitivity bandwidth", uncertainty=1),
        q("field_sensitivity", 9, "mV/m/sqrt(Hz)", "p. 3", "Full-acquisition-time-normalized field sensitivity", kind="approximate")]),
    exp("single_frequency_comparison", "단일 주파수 probe 대역폭 대조", [q("probe_power", 19, "uW", "p. 3")], [
        q("instantaneous_bandwidth", 400, "kHz", "p. 3", "Approximate 3-dB point", kind="approximate")]),
], [finding("bandwidth", "광 comb로 판독 대역을 넓히면서 감도·광학 span·RF modulation bandwidth를 각각 평가.", "pp. 2-3, Fig. 2", "bandwidth")],
    pid="artusio_glimpse_2024_preprint", group="NIST_frequency_comb")

add("multichannel_2025_natcom", "133Cs", "6S-6P-7S-Rydberg_three_photon", "coupling_comb_multichannel",
    ["frequency_comb", "multichannel", "three_photon", "tuning_range"], [
    exp("multichannel", "Coupling comb 다중 Rydberg 채널", [
        q("probe_wavelength", 895, "nm", "p. 2"), q("dressing_wavelength", 636, "nm", "p. 2"),
        interval("rf_tuning_range", 1, 40, "GHz", "p. 1, Abstract", "Coverage through several Rydberg channels; not simultaneous IF bandwidth"),
    ], [q("rydberg_channel_count", 7, "1", "p. 1, Abstract", "Up to seven states", kind="upper_bound")], kind="apparatus"),
], [finding("multiplex", "Coupling comb로 여러 Rydberg 상태를 빠르게 준비·전환하여 다중 RF 주파수를 검출.", "pp. 1-4", "multichannel")],
    pid="prajapati_2025_natcom", group="NIST_frequency_comb")

add("ju_2026_preprint", "85Rb addressed; isotope purity not specified", "5S1/2(F3)-5P3/2(F4)-40D5/2-39F7/2", "counterpropagating_EIT_atomic_superheterodyne",
    ["gold_reference", "heterodyne", "photon_shot_noise", "temperature", "magnetic_compensation", "laser_locking"], [
    exp("gold_operating_point", "Gold 실험 메인 운전점", [
        q("cell_length", 50, "mm", "p. 9, Fig. 5 discussion", "Vapor-cell length"),
        q("temperature", 20, "degC", "p. 9, Fig. 5 discussion", "Room-temperature operation", kind="approximate"),
        q("probe_wavelength", 780, "nm", "pp. 2-3", "5S(F3)-5P(F4)", kind="approximate"),
        q("probe_transition_text", "5S1/2(F=3) -> 5P3/2(F'=4)", "1", "p. 2, Sec. 2", "Experimental schematic in the main text; used as the nominal reference"),
        q("probe_transition_figure", "5S1/2(F=2) -> 5P3/2(F'=3)", "1", "p. 3, Fig. 1(a)", "Figure labels conflict with the main-text hyperfine labels", note="⟦AI⟧ Source inconsistency preserved. Nominal reference follows Sec. 2 text."),
        q("coupling_wavelength", 480, "nm", "pp. 2-3", "5P-40D", kind="approximate"),
        q("probe_power", 6, "uW", "p. 3", "Reported probe power"),
        q("coupling_power", 30, "mW", "p. 3", "Reported coupling power"),
        q("beam_diameter", 0.15, "mm", "p. 3", "Reported beam diameter; 1/e^2/FWHM convention unspecified"),
        missing("beam_diameter_definition", "1", "p. 3", "Paper says diameter; no explicit Gaussian radius convention assigned"),
        q("optical_geometry", "counter_propagating", "1", "p. 3", "Probe and coupling counter-propagate"),
        q("optical_polarizations", "linear_orthogonal", "1", "p. 3"),
        q("detector_qe", 95, "percent", "p. 2", "Balanced photodetector quantum efficiency", kind="approximate"),
        q("lo_rabi_over_2pi", 3.7, "MHz", "p. 2", "Microwave local oscillator"),
        q("lo_frequency", 37.07104, "GHz", "p. 5"),
        q("intermediate_frequency", 40, "kHz", "pp. 2, 5", "LO-SIG frequency offset"),
        missing("probe_detuning", "MHz", "pp. 2-5", "Absolute +40 MHz optical detuning is not specified for this reported operating point"),
        missing("operating_optical_depth", "1", "pp. 3, 7, 9", "No absolute measured operating-point OD inserted from AutoOD"),
    ], [q("field_sensitivity", 12.5, "nV/cm/sqrt(Hz)", "pp. 1, 7-8", "Measured superheterodyne sensitivity; RF-field RMS conversion defined in Eq. (1)", uncertainty=0.8),
        q("psn_field_sensitivity", 11.2, "nV/cm/sqrt(Hz)", "pp. 7, 10", "Calculated probe-photon-shot-noise floor; field convention of Eqs. (1)-(2)"),
        q("eit_linewidth", 1.6, "MHz", "pp. 1, 4-5", "Narrow compensated spectrum; its probe power is not specified. Fig. 2(b) gives approximately 3.0 MHz at 6 uW; do not bind 1.6 MHz to the 6 uW operating point.")],
        notes=["사용자 2026-11-03 채택 Gold reference. 실험 메인 구도·조건의 기준.",
               "⟦AI⟧ 논문은 85Rb 전이를 사용. 순수 동위원소 셀이라는 추가 가정은 보류."]),
    exp("temperature_model", "온도 의존 PSN 계산", [
        q("cell_length", 50, "mm", "p. 9"), q("beam_diameter", 0.15, "mm", "p. 9", "Same stated diameter convention"),
    ], [q("predicted_optimum_temperature", 33, "degC", "p. 9, Fig. 5(a)", "Calculated PSN optimum; distinct from measured ~20 C operation", kind="approximate")], kind="theory_prediction"),
], [finding("readout", "감도 식 (2)는 동작점의 투과 광자수와 응답 기울기를 함께 사용. Probe 파워 증가의 광자수 이득과 EIT power broadening이 경쟁.", "pp. 7-9, Eq. (2), Figs. 4-5", "sensitivity_definition"),
    finding("coherence", "자기장 보상·두 레이저의 ULE cavity 고정·balanced 검출로 EIT 폭과 판독 잡음을 낮춤.", "pp. 2-5", "coherence_control"),
    finding("hyperfine_label_discrepancy", "본문 2절은 5S F=3 → 5P F′=4, 그림 1(a)는 F=2 → F′=3으로 표기되어 서로 다름.", "pp. 2-3, Sec. 2 and Fig. 1(a)", "source_consistency"),
    finding("linewidth_power_scope", "⟦AI⟧ Fig. 2(b)의 보상 후 폭은 0.5 µW에서 약 1.6 MHz, 6 µW에서 약 3.0 MHz. Fig. 2(a)의 1.6 MHz 스펙트럼에는 프로브 파워가 명시되지 않아 6 µW 앵커로 사용할 수 없음.", "p. 5, Fig. 2(a-b)", "source_consistency"),
    finding("psn_curve_shape_scope", "⟦AI⟧ Fig. 5(b)의 빨간 SNL 곡선은 파워 증가에 따라 단조 감소하지만 측정 곡선은 6 µW 부근에서 최소. 본문의 최적점 설명과 그림의 PSN 곡선을 구분해야 함.", "pp. 8-10, Figs. 4(d), 5(b)", "source_consistency"),
    finding("optical_scan_identity", "⟦AI⟧ Fig. 2·3 주파수 축은 스캔한 레이저와 고정 레이저의 절대 디튜닝을 명시하지 않음. 이전 학위논문의 probe scan은 다른 장치의 근거.", "pp. 3-6, Figs. 2-3", "source_consistency")], group="PNU_Moon_KRISS",
    notes=["⟦AI⟧ Gold 메인 전이는 본문 2절 표기 채택. 그림 1(a)의 다른 초미세 표기는 별도 근거로 보존.",
           "⟦AI⟧ 논문이 별도 Supplementary Materials를 언급. 현재 보관·검수한 PDF는 main article 15쪽이며 별도 보충문서는 확보하지 않음."])

add("continuous_2026_comphys", "133Cs", "6S-6P-Zeeman_tuned_Rydberg", "magnetically_tunable_superheterodyne",
    ["tuning_range", "Zeeman", "heterodyne"], [
    exp("tuning_campaign", "Zeeman 연속 주파수 조정 캠페인", [
        interval("rf_tuning_range", 1, 40, "GHz", "pp. 1-2", "Carrier tuning coverage, not instantaneous signal bandwidth"),
    ], [q("field_sensitivity", 65, "nV/cm/sqrt(Hz)", "pp. 1-2", "Reported sensitivity upper envelope over the tuning campaign", kind="upper_bound")],
        kind="scan", notes=["전체 구간 상한 보고. 모든 주파수의 개별 측정점이나 동일 운전점으로 확장하지 않음."]),
], [finding("zeeman", "자기장으로 Rydberg 전이를 이동시켜 이산 공명 사이의 carrier 주파수 조정을 구현.", "pp. 2-6", "tuning_range")], pid="yao_2026_comphys", group="NUDT_Yao_Fu")

add("schlossberger_2026_jscheme", "133Cs", "J_scheme_6P1/2-6S1/2-7P-Rydberg", "three_photon_J_scheme",
    ["three_photon", "repumping", "heterodyne", "linewidth"], [
    exp("narrow_line", "낮은 광파워의 좁은 EIT", [
        q("probe_power", 1.3, "uW", "p. 2, Table I"), q("dressing_power", 17, "uW", "p. 2, Table I"),
        q("coupling_power", 460, "uW", "p. 2, Table I"),
        q("probe_beam_fwhm", 1.6, "mm", "p. 2, Table I", "FWHM"),
    ], [q("eit_linewidth", 1.3, "MHz", "p. 1, Abstract", "Low-laser-power linewidth")],
        notes=["좁은 선폭용 Table I 조건. 아래 감도 운전점에 같은 파워를 자동 상속하지 않음."]),
    exp("linked_jscheme", "연결된 J-scheme RF 판독", [q("rf_frequency", 4.7, "GHz", "p. 1, Abstract")], [
        q("field_sensitivity", 27, "uV/m/sqrt(Hz)", "p. 1, Abstract", "Conventional heterodyne J-scheme")]),
    exp("population_repump", "분리된 hyperfine population 판독", [q("rf_frequency", 4.7, "GHz", "p. 1, Abstract")], [
        q("field_sensitivity", 39, "uV/m/sqrt(Hz)", "p. 1, Abstract", "Population-repump readout")]),
], [finding("architecture", "직접 diode 파장으로 가능한 J 모양 세 광자 결합과 별도 hyperfine 인구수 판독을 비교.", "pp. 1-4", "three_photon")], group="NIST_Boulder_Jscheme")

add("mems_2026_micronano", "133Cs", "6S(F4)-6P(F5)-48D5/2-46F7/2", "MEMS_differential_EIT_AT",
    ["MEMS", "miniaturization", "AT_splitting", "differential_detection"], [
    exp("mems", "MEMS 셀 RF 검출", [
        q("silicon_thickness", 6, "mm", "pp. 1-2", "Fabrication thickness; not automatically an effective optical length"),
        q("rf_frequency", 29.75, "GHz", "p. 4, Fig. 3"),
        q("probe_beam_size", 500, "um", "p. 4", "Beam size, convention unspecified"),
        q("coupling_beam_size", 600, "um", "p. 4", "Beam size, convention unspecified"),
    ], [q("minimum_detectable_field", 2.8, "mV/cm", "p. 7, Fig. 7 discussion", "Minimum detected microwave field, not Hz-normalized sensitivity", uncertainty=0.3)]),
], [finding("mems", "고저항 silicon 기반 MEMS 셀과 광학 차분 검출로 소형 RF 검출 장치를 구성.", "pp. 2-6", "miniaturization")], pid="ma_2026_micronano", group="UESTC_MEMS_Ma")

add("offresonant_2026_sensors", "85Rb", "5S-5P-5D-55F7/2", "three_photon_offresonant_Stark_heterodyne",
    ["off_resonant", "three_photon", "sapphire_cell", "optimization", "field_calibration"], [
    exp("ism_40", "40.68 MHz 비공명 carrier", [
        q("rf_frequency", 40.68, "MHz", "p. 9, Table 1"), q("intermediate_frequency", 50, "kHz", "p. 8"),
    ], [q("field_sensitivity", 333, "uV/m/sqrt(Hz)", "p. 9, Table 1", "Table 1 summary estimate", uncertainty=6)],
        notes=["p. 9 Fig. 7의 별도 고정 SIG 평가 352.2±13.7 µV/m/√Hz와 구분."]),
    exp("ism_27", "27.12 MHz 비공명 carrier", [q("rf_frequency", 27.12, "MHz", "p. 9, Table 1")], [
        q("field_sensitivity", 449, "uV/m/sqrt(Hz)", "p. 9, Table 1", "Table 1 summary estimate", uncertainty=8)]),
    exp("ism_13", "13.56 MHz 비공명 carrier", [q("rf_frequency", 13.56, "MHz", "p. 9, Table 1")], [
        q("field_sensitivity", 399, "uV/m/sqrt(Hz)", "p. 9, Table 1", "Table 1 summary estimate", uncertainty=10)]),
    exp("ism_6", "6.78 MHz 비공명 carrier", [q("rf_frequency", 6.78, "MHz", "p. 9, Table 1")], [
        q("field_sensitivity", 125, "uV/m/sqrt(Hz)", "p. 9, Table 1", "Table 1 summary estimate", uncertainty=4)]),
], [finding("optimization", "LO 세기와 coupling detuning을 자동 최적화. Sapphire 셀의 저주파 투과와 외부 RF 장의 보정 조건을 명시.", "pp. 2-8", "optimization"),
    finding("metric", "감도는 photodiode PSD를 전기장 단위로 환산한 값. 최소 검출장·dynamic range와 별도 보고.", "pp. 4, 8-9", "sensitivity_definition")],
    pid="viray_2026_sensors", group="GTRI_Viray_Sawyer")

(ROOT / "reference_policy.json").write_text(json.dumps(POLICY, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print("Curated", len(list((ROOT / "papers").glob("*.json"))), "papers; Gold", POLICY["gold_reference_paper_id"])
