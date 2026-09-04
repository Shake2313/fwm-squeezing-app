"""
SABES front-end — the lab-facing page.

Rendered by `streamlit_app.py` when the URL carries `?app=sabes`, so it deploys
with GABES but presents as its own site: its own name, its own sidebar, its own
vocabulary. Nothing here is reachable from the GABES scheme dropdown; the only
way in is the link at the bottom of the GABES sidebar.

Three compute tiers, which is what makes an optics-heavy UI usable:

  A  optics      the source chain and post-cell geometry. Closed form, sub-
                 millisecond, runs every rerun. Beam sizes, the power budget,
                 spectral purity and the derived-parameter table all live here.
  B  solve       the Bloch equations, seconds, cached on the *derived* GABES
                 parameters. Moving a knob that does not change what reaches the
                 cell is a cache hit, so detection knobs feel instant.
  C  readout     photodiode powers, linearity, clearance. Needs the gains from B
                 but costs nothing, so it re-renders freely.

The tier split is the reason the cache key is `bridge.solve_key(params)` and not
the settings: two different sets of waveplate angles that deliver the same power
to the cell are the same solve.
"""
from dataclasses import fields, replace
from html import escape
from threading import RLock

import numpy as np
import streamlit as st

from gabes import constants
from gabes.core import blas_single_thread
from gabes.schemes import fwm

from sabes import bridge, detection as detection_model, instruments, layout
from sabes.components import canvas as canvas_module
from sabes.beamline import (SetupSettings, solve_seed_polarizer_deg,
                            solve_seed_trim_deg, solve_split_angle_deg)
from sabes.calibration import default_calibration
from sabes.detection import DetectionSettings

SESSION_PREFIX = "_sabes_"
SETTINGS_VERSION = "sabes-eom-noise-v3"

FIDELITIES = (fwm.FIDELITY_FAST, fwm.FIDELITY_ULTRA)

# Matplotlib's font and layout caches are process-global and Streamlit can
# overlap reruns, so figure construction is serialised the same way the
# GABES app does it.
_PLOT_LOCK = RLock()

#: Every knob, keyed by name. One definition serves the optic dialog and the
#: sidebar, so the two cannot drift apart -- which matters now that almost
#: everything lives behind a click on the drawing rather than in a fixed list.
#:
#: (label, unit, min, max, step, help)
CONTROLS = {
    "ecdl_power_mw": ("ECDL output", "mW", 1.0, 200.0, 1.0,
                      "Seed power into the tapered amplifier."),
    "ta_current_a": ("TA current", "A", 0.5, 4.0, 0.05,
                     "Amplifier drive current. Output saturates, so this is "
                     "not a linear knob."),
    "hwp_ta_deg": ("HWP angle", "°", 0.0, 90.0, 0.5,
                   "Rotates the seed onto the amplifier's gain axis."),
    "hwp_split_deg": ("HWP angle", "°", 0.0, 45.0, 0.01,
                      "Sets how the amplifier output divides between pump and "
                      "seed. This is the pump-power knob."),
    "seed_polarizer_deg": ("Polarizer angle", "°", 0.0, 90.0, 0.05,
                           "Holds the fiber EOM below its optical rating. "
                           "Independent of the split, because the split alone "
                           "cannot satisfy both constraints."),
    "eom_offset_mhz": ("Generator offset from ν_HF", "MHz", -200.0, 200.0, 0.1,
                       "Signal-generator frequency relative to the 3.0357324 "
                       "GHz ground hyperfine splitting. This IS the two-photon "
                       "detuning, with the sign flipped and no calibration "
                       "coefficient in between."),
    "eom_rf_dbm": ("RF drive", "dBm", 0.0, 30.0, 0.5,
                   "Sets the modulation index β = π V_peak / V_π and hence how "
                   "much power lands in the wanted −1 sideband."),
    "eom_unwanted_rin_db_per_hz": (
        "Unwanted-mode RIN", "dBc/Hz", -180.0, -60.0, 1.0,
        "One-sided fractional-intensity-noise PSD for each unwanted "
        "post-etalon EOM mode. It is an input, not an RF/EOM calculation. "
        "Measure it at the detector plane before using it as a calibrated value. "
        "The -100 dBc/Hz default is a Sim-scale sensitivity setting, not a "
        "measurement of the installed system."),
    "hwp_eom_deg": ("HWP angle", "°", 0.0, 90.0, 0.5,
                    "Aligns the input to the modulator's polarization axis."),
    "etalon_detune_ghz_1": ("Detuning", "GHz", -2.0, 2.0, 0.01,
                            "Tilt/temperature detuning of this stage from the "
                            "wanted sideband."),
    "etalon_detune_ghz_2": ("Detuning", "GHz", -2.0, 2.0, 0.01,
                            "Tilt/temperature detuning of this stage from the "
                            "wanted sideband."),
    "etalon_detune_ghz_3": ("Detuning", "GHz", -2.0, 2.0, 0.01,
                            "Tilt/temperature detuning of this stage from the "
                            "wanted sideband."),
    "seed_trim_hwp_deg": ("HWP angle", "°", 0.0, 45.0, 0.01,
                          "With the Glan-Taylor, sets the seed power at the "
                          "cell. The seed arm carries far more than the FWM "
                          "needs, so this mostly throws power away."),
    "seed_trim_qwp_deg": ("QWP angle", "°", 0.0, 90.0, 0.5,
                          "Fine polarization control ahead of the "
                          "Glan-Taylor."),
    "seed_gtp_deg": ("Analyser angle", "°", 0.0, 90.0, 0.5,
                     "Glan-Taylor analyser angle for the trim pair."),
    "pump_telescope_f1_mm": ("Focal length", "mm", 25.0, 1000.0, 5.0,
                             "Cell waist = collimator diameter × f2/f1 × the "
                             "shared scale factor."),
    "pump_telescope_f2_mm": ("Focal length", "mm", 25.0, 1000.0, 5.0,
                             "Second lens of the pump telescope."),
    "seed_telescope_f1_mm": ("Focal length", "mm", 25.0, 1000.0, 5.0,
                             "First lens of the seed telescope."),
    "seed_telescope_f2_mm": ("Focal length", "mm", 25.0, 1000.0, 5.0,
                             "Second lens of the seed telescope."),
    "opd_ghz": ("One-photon detuning Δ", "GHz", -3.0, 3.0, 0.05,
                "Where the pump sits. Set by the ECDL lock, not modelled from "
                "its temperature / current / PZT — mode hops make that "
                "unpredictable."),
    "cell_temp_c": ("Temperature", "°C", 60.0, 150.0, 1.0,
                    "Vapour temperature from the heater controller."),
    "dmirror_separation_mm": ("Separation here", "mm", 0.5, 15.0, 0.05,
                              "Alignment is done by overlapping at the cell "
                              "and separating by this much 23 inch "
                              "downstream, so this is the crossing-angle "
                              "knob."),
    "pump_leakage_dbm": ("Pump leakage", "dBm", -120.0, -10.0, 1.0,
                         "OBSERVED pump power reaching the detector — measure "
                         "it, do not derive it. A Gaussian tail claims 1e-9 "
                         "rejection, which is useless and contradicted by the "
                         "paper naming pump scatter as a real limit. It enters "
                         "the spectrum analyser through its calculable shot noise. "
                         "Classical pump RIN remains unapplied until a measured or "
                         "predeclared spectral density is supplied."),
    "probe_lens_focal_mm": ("Focal length", "mm", 25.0, 500.0, 5.0,
                            "Focusing lens onto the probe photodiode."),
    "conjugate_lens_focal_mm": ("Focal length", "mm", 25.0, 500.0, 5.0,
                                "Focusing lens onto the conjugate photodiode."),
    "pd_defocus_mm": ("Defocus", "mm", 0.0, 100.0, 1.0,
                      "Distance of the diode from the focal plane. Spreading "
                      "the spot is the cheapest way out of the linearity "
                      "limit."),
}

#: The etalon stages, in chain order. Each owns its own detuning: the model has
#: always accepted a 3-tuple and only the UI was collapsing it, while in the lab
#: each stage has its own knob and its own temperature.
ETALON_KEYS = ("etalon_detune_ghz_1", "etalon_detune_ghz_2",
               "etalon_detune_ghz_3")

# ----------------------------------------------------------------------
# Tier B — the only cached layer
# ----------------------------------------------------------------------
@st.cache_data(show_spinner=False, max_entries=64)
def _cached_solve(solve_items, cache_version):
    scheme = fwm.FWMScheme()
    params = dict(scheme.defaults())
    params.update(dict(solve_items))
    with blas_single_thread():
        return scheme.compute(params)


@st.cache_data(show_spinner=False, max_entries=64)
def _cached_observables(raw, param_items, cache_version):
    scheme = fwm.FWMScheme()
    return scheme.observables(raw, dict(param_items))


def _key(name):
    return SESSION_PREFIX + name


def _defaults():
    base = SetupSettings()
    detection = DetectionSettings()
    values = {f.name: getattr(base, f.name) for f in fields(base)}
    values.update({f.name: getattr(detection, f.name) for f in fields(detection)})
    # The generator is exposed as an offset from the hyperfine splitting, because
    # "3.0437 GHz" is nine digits of which only the last three ever move.
    values["eom_offset_mhz"] = (base.eom_frequency_hz
                                - constants.NU_GROUND_HF) / 1e6
    for index, name in enumerate(ETALON_KEYS):
        values[name] = float(base.etalon_detune_ghz[index])
    values.pop("etalon_detune_ghz", None)
    values["resolution"] = FIDELITIES[0]
    values["instrument"] = INSTRUMENTS[0]
    values["scope_mode"] = "Scan"
    return values


def _seed_session_state():
    version_key = _key("_version")
    values = _defaults()
    if st.session_state.get(version_key) != SETTINGS_VERSION:
        for name, value in values.items():
            st.session_state[_key(name)] = value
        st.session_state[version_key] = SETTINGS_VERSION
    else:
        for name, value in values.items():
            st.session_state.setdefault(_key(name), value)


def _current():
    """Session state -> (SetupSettings, DetectionSettings)."""
    get = lambda name: st.session_state[_key(name)]
    settings = SetupSettings(
        ecdl_power_mw=get("ecdl_power_mw"),
        ta_current_a=get("ta_current_a"),
        hwp_ta_deg=get("hwp_ta_deg"),
        hwp_split_deg=get("hwp_split_deg"),
        seed_polarizer_deg=get("seed_polarizer_deg"),
        hwp_eom_deg=get("hwp_eom_deg"),
        eom_frequency_hz=(constants.NU_GROUND_HF
                          + get("eom_offset_mhz") * 1e6),
        eom_rf_dbm=get("eom_rf_dbm"),
        etalon_detune_ghz=tuple(float(get(name)) for name in ETALON_KEYS),
        seed_trim_qwp_deg=get("seed_trim_qwp_deg"),
        seed_trim_hwp_deg=get("seed_trim_hwp_deg"),
        seed_gtp_deg=get("seed_gtp_deg"),
        pump_telescope_f1_mm=get("pump_telescope_f1_mm"),
        pump_telescope_f2_mm=get("pump_telescope_f2_mm"),
        seed_telescope_f1_mm=get("seed_telescope_f1_mm"),
        seed_telescope_f2_mm=get("seed_telescope_f2_mm"),
        opd_ghz=get("opd_ghz"),
        cell_temp_c=get("cell_temp_c"),
        dmirror_separation_mm=get("dmirror_separation_mm"),
    )
    detection = DetectionSettings(
        pump_leakage_dbm=get("pump_leakage_dbm"),
        probe_lens_focal_mm=get("probe_lens_focal_mm"),
        conjugate_lens_focal_mm=get("conjugate_lens_focal_mm"),
        pd_defocus_mm=get("pd_defocus_mm"),
        eom_unwanted_rin_db_per_hz=get("eom_unwanted_rin_db_per_hz"),
    )
    return settings, detection


# ----------------------------------------------------------------------
# Sidebar
# ----------------------------------------------------------------------
def _apply_solved_knobs():
    """Re-solve the three power knobs onto the paper operating point."""
    settings, _ = _current()
    calibration = _session_calibration()
    try:
        split = solve_split_angle_deg(0.600, settings, calibration)
        settings = replace(settings, hwp_split_deg=split)
        polarizer = solve_seed_polarizer_deg(7.0e-3, settings, calibration)
        settings = replace(settings, seed_polarizer_deg=polarizer)
        trim = solve_seed_trim_deg(8.0e-6, settings, calibration)
    except ValueError as exc:
        st.session_state[_key("_solver_error")] = str(exc)
        return
    st.session_state.pop(_key("_solver_error"), None)
    st.session_state[_key("hwp_split_deg")] = round(split, 4)
    st.session_state[_key("seed_polarizer_deg")] = round(polarizer, 4)
    st.session_state[_key("seed_trim_hwp_deg")] = round(trim, 4)


def _back_to_gabes():
    st.query_params.clear()


def _reset_defaults():
    for name, value in _defaults().items():
        st.session_state[_key(name)] = value


def _render_sidebar(host):
    st.sidebar.markdown(
        "<div class='sabes-brand'>"
        "<div class='sabes-brand-name'>SABES</div>"
        "<div class='sabes-brand-sub'>Specific Atomic Bloch Equation Solver</div>"
        "<div class='sabes-brand-note'>EOM-based ⁸⁵Rb seeded-FWM "
        "gain diagnostic</div></div>",
        unsafe_allow_html=True,
    )
    st.sidebar.button("← Back to GABES", on_click=_back_to_gabes,
                      use_container_width=True)
    st.sidebar.divider()

    left, right = st.sidebar.columns(2)
    left.button("Solve knobs", on_click=_apply_solved_knobs,
                use_container_width=True,
                help="Set the split, polarizer and trim angles that put 600 mW "
                     "on the pump, 7 mW into the EOM and 8 µW of seed at the "
                     "cell, for the calibration currently in force.")
    right.button("Reset", on_click=_reset_defaults, use_container_width=True)

    error = st.session_state.get(_key("_solver_error"))
    if error:
        st.sidebar.warning(error)

    st.sidebar.selectbox(
        "Model fidelity", FIDELITIES, key=_key("resolution"),
        help="Fast is for exploring. Ultra is the anchored path and the only one "
             "that carries the Gaussian crossing overlap — geometry optimisation "
             "must use it.")

    st.sidebar.divider()
    st.sidebar.caption(
        "Every optic's parameters live on the optic. Click it in the "
        "**Optical table** tab to open its editor.")

    # Safety net, not a feature: anything the drawing does not claim would
    # otherwise be unreachable. A test asserts every knob has an owner, so in a
    # healthy build this expander never appears.
    orphans = sorted(set(CONTROLS) - _owned_parameters())
    if orphans:
        with st.sidebar.expander(f"Unplaced knobs ({len(orphans)})",
                                 expanded=False):
            st.caption("No optic on the drawing owns these yet.")
            for name in orphans:
                _render_control(st, name)


def _session_calibration():
    """The shipped calibration with any per-optic transmission edits applied."""
    calibration = default_calibration()
    updates = {
        name: st.session_state[SESSION_PREFIX + name]
        for name in _editable_transmission_keys()
        if SESSION_PREFIX + name in st.session_state
    }
    return calibration.with_values(updates)


def _editable_transmission_keys():
    return {node.transmission_key for part in layout.LAYOUTS
            for node in part.nodes if node.transmission_key}


def _seed_transmission_state(calibration):
    for name in _editable_transmission_keys():
        st.session_state.setdefault(SESSION_PREFIX + name,
                                    calibration.value(name))


def _owned_parameters():
    """Every parameter some optic on some layout claims."""
    return {name for part in layout.LAYOUTS
            for node in part.nodes for name in node.params}


def _render_control(container, name):
    """One numeric input, from the single shared definition."""
    label, unit, low, high, step, help_text = CONTROLS[name]
    container.number_input(
        f"{label} [{unit}]" if unit else label,
        min_value=float(low), max_value=float(high), step=float(step),
        key=_key(name), help=help_text,
        format="%.4f" if step < 0.05 else "%.2f")


# ----------------------------------------------------------------------
# Main panel
# ----------------------------------------------------------------------
def _twin_states(result):
    """Return post-cell probe and conjugate beams.

    Parametric gain applies only to the selected seed mode. Other EOM modes
    remain in the probe arm and are absent from the conjugate arm.
    """
    if result.raw is None:
        return None
    gain_s, gain_c = result.gains
    default_transmission = result.geometry.optics_transmission
    transmissions = {
        arm.name: arm.optics_transmission
        for arm in getattr(result.geometry, "arms", ())
    }
    probe_transmission = transmissions.get("probe", default_transmission)
    conjugate_transmission = transmissions.get("conjugate",
                                                default_transmission)
    seed = result.chain.seed
    probe_lines = []
    conjugate_lines = []
    for line in seed.lines:
        if abs(line.offset_hz - result.chain.seed_offset_hz) <= 1.0:
            probe_lines.append(line.scaled(gain_s * probe_transmission))
            conjugate_lines.append(replace(
                line,
                offset_hz=-line.offset_hz,
                power_w=line.power_w * gain_c * conjugate_transmission,
                label="conjugate",
            ))
        else:
            probe_lines.append(line.scaled(probe_transmission))
    return {"probe": seed.with_lines(probe_lines),
            "conjugate": seed.with_lines(conjugate_lines)}


def _probe_key(part_key):
    return _key(f"probe_{part_key}")


def _handle_canvas_event(event, part_key):
    """Route a click: an optic selects, empty table drops a probe.

    Returns whether anything changed. The caller reruns on True, because the
    component's value only reaches Python at the *start* of a run -- the panels
    below have already read the old selection by the time this fires, so without
    a rerun the drawing and the readout would trail one click behind.
    """
    if not event:
        return False
    seen = _key(f"seq_{part_key}")
    if st.session_state.get(seen) == event.get("seq"):
        return False
    st.session_state[seen] = event.get("seq")
    if event.get("kind") == canvas_module.KIND_OPTIC:
        st.session_state[_key("selected")] = event.get("id")
        # Clicking an optic opens its editor straight away; the flag survives
        # the rerun the dialog needs, and the dialog clears it on close.
        st.session_state[_key("open_optic")] = event.get("id")
    else:
        st.session_state[_probe_key(part_key)] = [event["x"], event["y"]]
        st.session_state[_key("selected")] = None
    return True


@st.dialog("Optic")
def _optic_dialog(part_key, node_id):
    """Edit one optic's parameters, in a modal on top of the drawing.

    The controls are the same `CONTROLS` entries the sidebar would render, so an
    optic's knob behaves identically wherever it appears. Streamlit writes each
    widget straight into session state, so there is nothing to apply -- closing
    the dialog is enough, and the drawing behind it is already current.
    """
    part = layout.get(part_key)
    try:
        node = part.node(node_id)
    except KeyError:
        st.warning("That optic is not on this part of the table.")
        return

    st.markdown(f"### {node.display_label or node.id}")
    if node.note:
        st.caption(node.note)

    editable = [name for name in node.params if name in CONTROLS]
    for name in editable:
        _render_control(st, name)

    if node.transmission_key:
        st.markdown("**Transmission**")
        st.number_input(
            "Per-pass transmission", min_value=0.0, max_value=1.0, step=0.001,
            key=SESSION_PREFIX + node.transmission_key, format="%.4f",
            help="Fraction of power this optic passes. Editing it marks the "
                 "coefficient as hand-set, which is what it is until somebody "
                 "measures it.")
        if node.lumped:
            st.caption("This optic is **lumped**: a transmission is the whole "
                       "of what the model knows about it.")

    missing = [name for name in node.params if name not in CONTROLS]
    if missing:
        st.caption("Recorded but not modelled here: " + ", ".join(missing))
    if not editable and not missing and not node.transmission_key:
        st.caption("This optic is in the physics but owns no adjustable "
                   "parameter.")

    if st.button("Done", use_container_width=True):
        st.session_state.pop(_key("open_optic"), None)
        st.rerun()


def _render_optic_panel(part, node_id, settings, detection):
    """Summary of the selected optic, with a way into its editor."""
    try:
        node = part.node(node_id)
    except KeyError:
        return
    values = {f.name: getattr(settings, f.name) for f in fields(settings)}
    values.update({f.name: getattr(detection, f.name) for f in fields(detection)})
    values["eom_offset_mhz"] = (settings.eom_frequency_hz
                                - constants.NU_GROUND_HF) / 1e6
    for index, name in enumerate(ETALON_KEYS):
        values[name] = settings.etalon_detune_ghz[index]

    st.markdown(f"**{node.display_label or node.id}**")
    if node.note:
        st.caption(node.note)
    rows = [(CONTROLS[name][0] if name in CONTROLS else name,
             f"{values.get(name, float('nan')):.4g}"
             + (f" {CONTROLS[name][1]}" if name in CONTROLS else ""))
            for name in node.params]
    if rows:
        st.markdown(_markdown_table(("Parameter", "Value"), rows))
    if st.button("Edit…", key=f"edit_{part.key}_{node_id}",
                 use_container_width=True):
        st.session_state[_key("open_optic")] = node_id
        st.rerun()


#: What can be dropped on a point of the table. The last two need a photodiode,
#: exactly as on the bench, so they read through one rather than off the beam.
INSTRUMENTS = ("Power meter", "Beam profiler", "Wavemeter", "Photodiode",
               "Oscilloscope", "Spectrum analyser")


def _detector_from_calibration(calibration):
    c = calibration.value
    return instruments.Photodiode(
        responsivity_a_per_w=c("pd_responsivity_a_per_w"),
        transimpedance_v_per_a=c("bpd_transimpedance_v_per_a"),
        bandwidth_hz=c("bpd_bandwidth_hz"),
        nep_w_per_rthz=c("bpd_nep_w_per_rthz"),
        saturation_w=c("bpd_cw_saturation_w"))


def _read_instrument(choice, reading, result, calibration, detection):
    """Run the chosen instrument against a probed beam."""
    beam = reading.beam
    seed_offset = result.chain.seed_offset_hz
    # Only the seed path has a "wanted" line; asking the pump about one would
    # report zero and mean nothing.
    wanted = {
        "seed": seed_offset,
        "probe": seed_offset,
        "conjugate": -seed_offset,
    }.get(reading.link.beam)

    if choice == "Power meter":
        return instruments.PowerMeter(wanted_offset_hz=wanted).measure(beam)
    if choice == "Beam profiler":
        return instruments.BeamProfiler().measure(beam)
    if choice == "Wavemeter":
        return instruments.Wavemeter(wanted_offset_hz=wanted).measure(beam)

    detector = _detector_from_calibration(calibration)
    if choice == "Photodiode":
        return detector.measure(beam)

    signal = detector.convert(beam)
    if choice == "Oscilloscope":
        return _scope_reading(signal, result, reading)
    if choice == "Spectrum analyser":
        readout = result.readout
        balanced = reading.link.beam in ("probe", "conjugate")
        total = (readout.total_power_w
                 if balanced and readout else beam.total_power_w)
        return instruments.SpectrumAnalyzer().analyze(
            signal,
            result.gain_referred_noise_db if balanced else None,
            total_power_w=total,
            pump_leakage_dbm=(detection.pump_leakage_dbm
                              if balanced else None),
            eom_noise=(readout.eom_noise
                       if balanced and readout else None))
    raise KeyError(choice)


def _scope_reading(signal, result, reading):
    """Return a swept-frequency trace or a synthetic time trace."""
    mode = st.session_state.get(_key("scope_mode"), "Scan")
    if mode == "Scan" and result.raw is not None:
        axis = np.asarray(result.raw["probe_axis_GHz"])
        beam_name = reading.link.beam
        gain_key = {"probe": "G_s", "conjugate": "G_c"}.get(beam_name)
        if gain_key is None:
            power = np.full_like(axis, reading.beam.total_power_w, dtype=float)
        else:
            arm_index = 0 if beam_name == "probe" else 1
            operating_gain = result.gains[arm_index]
            wanted_offset = (
                -result.chain.seed_offset_hz
                if beam_name == "conjugate"
                else result.chain.seed_offset_hz
            )
            wanted = reading.beam.power_at(wanted_offset)
            scale = wanted / operating_gain if operating_gain > 0.0 else 0.0
            background = max(reading.beam.total_power_w - wanted, 0.0)
            power = np.asarray(result.raw[gain_key]) * scale + background
        return instruments.Oscilloscope().scan(
            signal, axis, power, x_label="Probe frequency", x_unit="GHz")
    return instruments.Oscilloscope().timeseries(signal)


def _render_reading(reading):
    """Render an instrument reading and its provenance."""
    if reading.synthesised:
        st.markdown(
            "<div class='sabes-synth'>SYNTHESISED — an illustrative trace built "
            "from the declared model inputs. It is not a measured-shape "
            "reproduction or validation evidence.</div>", unsafe_allow_html=True)
    if reading.note and not reading.synthesised:
        st.caption(reading.note)

    st.markdown(_markdown_table(
        ("Quantity", "Value"),
        [(q.label, q.formatted()) for q in reading.quantities]))
    for warning in reading.warnings:
        st.warning(warning, icon="⚠️")
    if reading.trace is not None:
        _render_trace(reading.trace)


def _render_trace(trace):
    import matplotlib.pyplot as plt

    from gabes.plot_style import apply_gabes_plot_style

    with _PLOT_LOCK:
        figure, axis = plt.subplots(figsize=(7.2, 3.0))
        for name, values in trace.series.items():
            if trace.kind == "stem":
                axis.vlines(trace.x, 0 if min(values) >= 0 else min(values),
                            values, linewidth=2.0)
                axis.plot(trace.x, values, "o", markersize=4, label=name)
            else:
                axis.plot(trace.x, values, linewidth=1.4, label=name)
        axis.set_xlabel(_axis_label(trace.x_label, trace.x_unit))
        axis.set_ylabel(_axis_label(trace.y_label, trace.y_unit))
        if len(trace.series) > 1:
            axis.legend(fontsize=8, loc="best")
        axis.grid(alpha=0.3)
        apply_gabes_plot_style(figure)
        st.pyplot(figure)
        plt.close(figure)


def _axis_label(label, unit):
    """Matplotlib strings stay ASCII -- the mathtext layout lock in the GABES
    app is there because unicode in axis labels has crashed layout before."""
    text = f"{label} [{unit}]" if unit else label
    return text.replace("µ", "u").replace("√", "sqrt").replace("²", "^2")


def _render_probe_panel(reading):
    if reading is None:
        st.caption("No beam there. Click a beam to probe it, or an optic to "
                   "select it.")
        return
    link = reading.link
    st.markdown(f"**{link.beam} beam** · {link.a} → {link.b}")
    rows = [
        ("Power", f"{reading.power_w * 1e3:.4g} mW"),
        ("1/e² radius", f"{reading.radius_m * 1e6:.1f} µm"),
        ("Along this segment", f"{reading.distance_along_m * 1e3:.1f} mm"
                               + ("" if link.path_m is not None
                                  else "  (segment has no modelled length)")),
        ("Spectral lines", str(len(reading.beam.lines))),
    ]
    if link.path_m is not None:
        rows.append(("Segment length",
                     f"{link.path_m * 1e3:.1f} mm  ({link.provenance})"))
    st.markdown(_markdown_table(("Quantity", "Value"), rows))


def _render_table_tab(result, calibration, host, detection):
    settings, detection = _current()
    twin = _twin_states(result)
    selected = st.session_state.get(_key("selected"))

    # A radio rather than nested tabs: Streamlit does not re-lay-out custom
    # components inside an inactive tab, so a canvas that first mounts hidden
    # stays collapsed at zero height. Drawing one part at a time sidesteps that
    # and halves the work per rerun.
    names = [p.title.split("—")[0].strip() for p in layout.LAYOUTS]
    choice = st.radio("Setup", names, key=_key("part"), horizontal=True,
                      label_visibility="collapsed")
    part = layout.LAYOUTS[names.index(choice)]

    st.caption(part.title)
    probe = st.session_state.get(_probe_key(part.key))
    spec = layout.build_spec(part, result.chain, selected=selected,
                             probe=probe, twin_states=twin)
    event = canvas_module.svg_canvas(spec, selected=selected,
                                     key=f"sabes_table_{part.key}")
    if _handle_canvas_event(event, part.key):
        st.rerun()

    pending = st.session_state.get(_key("open_optic"))
    if pending:
        _optic_dialog(part.key, pending)

    swatches = " ".join(
        f"<span class='sabes-key'><i style='background:{colour};"
        f"{'opacity:.75' if dashed else ''}'></i>{escape(name)}</span>"
        for colour, dashed, name in layout.legend_rows())
    st.markdown(
        f"<div class='sabes-legend'>{swatches}"
        "<span class='sabes-key sabes-key--note'>stroke width encodes "
        "power on a log scale — it is not the beam size</span></div>",
        unsafe_allow_html=True)

    left, right = st.columns([1, 1])
    with left:
        st.markdown("###### Selected optic")
        if selected:
            _render_optic_panel(part, selected, settings, detection)
        else:
            st.caption("Click a highlighted optic. Greyed hardware is drawn "
                       "for orientation but owns no parameter.")
    with right:
        st.markdown("###### Instrument")
        choice = st.selectbox(
            "Instrument", INSTRUMENTS, key=_key("instrument"),
            label_visibility="collapsed",
            help="Click anywhere on the table to put its head there. The "
                 "oscilloscope and the spectrum analyser read through a "
                 "photodiode, as they do on the bench.")
        if choice == "Oscilloscope":
            st.radio("Mode", ("Scan", "Time series"), key=_key("scope_mode"),
                     horizontal=True, label_visibility="collapsed",
                     help="Scan is a swept-laser trace and is computed. Time "
                          "series is sampled from the modelled noise PSD.")
        probed = None
        if probe:
            probed = layout.beam_at(part, probe[0], probe[1], result.chain,
                                    twin_states=twin)
        _render_probe_panel(probed)

    if probed is not None:
        st.divider()
        st.markdown(f"###### {choice} on the {probed.link.beam} beam")
        _render_reading(_read_instrument(
            choice, probed, result, calibration, detection))


def _markdown_table(header, rows):
    head = "| " + " | ".join(header) + " |\n"
    rule = "|" + "|".join("---" for _ in header) + "|\n"
    body = "".join("| " + " | ".join(str(c) for c in row) + " |\n" for row in rows)
    return head + rule + body


def _budget_table(chain):
    rows = []
    for path, name, power_w, waist_m in chain.budget_rows():
        power = (f"{power_w * 1e3:.4f} mW" if power_w >= 1e-6
                 else f"{power_w * 1e6:.4f} µW")
        rows.append((path, name, power, f"{waist_m * 1e6:.0f} µm"))
    return _markdown_table(("Path", "Stage", "Power", "Waist w₀"), rows)


def _etalon_transfer_table(chain):
    """Return per-mode transmission through each etalon."""
    transfers = [row for row in chain.etalon_mode_transfers
                 if row.eom_power_fraction >= 1.0e-8]
    headers = ["Mode", "Offset", "At EOM"]
    headers += [f"T{k + 1}" for k in range(len(chain.etalon_filters))]
    headers += ["Net pass", "Filtered", "|Eout/Ein|", "At cell"]

    def power(value_w):
        if value_w >= 1.0e-6:
            return f"{value_w * 1e6:.3g} µW"
        if value_w >= 1.0e-9:
            return f"{value_w * 1e9:.3g} nW"
        return f"{value_w * 1e12:.3g} pW"

    rows = []
    for transfer in transfers:
        role = " wanted" if transfer.order == -1 else ""
        row = [f"n={transfer.order:+d}{role}",
               f"{transfer.offset_hz / 1e9:+.4f} GHz",
               f"{100.0 * transfer.eom_power_fraction:.3g} %"]
        row += [f"{100.0 * value:.3g} %"
                for value in transfer.stage_power_transmissions]
        row += [f"{100.0 * transfer.total_power_transmission:.3g} %",
                f"{100.0 * transfer.filtered_power_fraction:.3g} %",
                f"{100.0 * transfer.field_magnitude_transmission:.3g} %",
                power(transfer.cell_power_w)]
        rows.append(row)
    return _markdown_table(tuple(headers), rows)


def _provenance_panel(calibration):
    unverified = calibration.unverified()
    verified = len(calibration) - len(unverified)
    st.markdown(
        f"**{verified} of {len(calibration)}** coefficients are backed by a "
        f"datasheet, the paper, or a fit. The rest are placeholders — results "
        f"that depend on them are *predicted*, not *reproduced*."
    )
    rows = [(c.name, f"{c.value:g}", c.unit or "—", c.provenance)
            for c in unverified]
    st.markdown(_markdown_table(("Coefficient", "Value", "Unit", "Provenance"),
                                rows))


def _set_browser_title(title):
    """Retitle the browser tab, and keep it that way.

    `st.set_page_config` can only run once per app and GABES already called it,
    so the tab would otherwise still say GABES on a page meant to read as a
    separate site. A zero-height component reaching out to the top document fixes
    that -- but Streamlit re-asserts its configured title on every rerun, so a
    one-shot assignment flips back on the next interaction. An observer re-applies
    it, and disconnects itself the moment the URL stops being SABES so it cannot
    follow the user back to GABES.

    It targets `window.top`, not `window.parent`: Streamlit Community Cloud wraps
    the app in a second same-origin iframe, so reaching only one level up retitles
    a document nobody sees while the real browser tab keeps its old name.
    """
    import streamlit.components.v1 as components
    components.html(
        """
        <script>
        (function () {
          var want = %s;
          // Highest same-origin ancestor: the real tab locally, and the tab
          // rather than the Community Cloud wrapper when deployed.
          var host = window.parent, doc = null;
          try {
            window.top.document.title;
            host = window.top;
          } catch (e) { host = window.parent; }
          try { doc = host.document; } catch (e) { return; }

          if (host.__sabesTitleObserver) { host.__sabesTitleObserver.disconnect(); }
          function inSabes() {
            try { return host.location.search.indexOf('app=sabes') !== -1; }
            catch (e) { return false; }
          }
          function apply() {
            if (!inSabes()) {
              if (host.__sabesTitleObserver) {
                host.__sabesTitleObserver.disconnect();
                host.__sabesTitleObserver = null;
              }
              return;
            }
            if (doc.title !== want) { doc.title = want; }
          }
          var head = doc.querySelector('head');
          var obs = new MutationObserver(apply);
          obs.observe(head, {subtree: true, childList: true, characterData: true});
          host.__sabesTitleObserver = obs;
          apply();
        })();
        </script>
        """ % repr(title),
        height=0,
    )


def render(host=None):
    """Entry point called by `streamlit_app.py`."""
    _set_browser_title("SABES — Specific Atomic Bloch Equation Solver")
    _inject_css()

    _seed_session_state()
    _render_sidebar(host)

    settings, detection = _current()
    calibration = _session_calibration()
    _seed_transmission_state(calibration)
    fidelity = st.session_state[_key("resolution")]

    # ---- Tier A: optics ----
    result = bridge.run(settings, calibration, detection, solve=False,
                        resolution=fidelity)
    params = result.params

    st.markdown(
        "<div class='sabes-header'><h1>SABES</h1>"
        "<p>Lab settings in, mean-field diagnostics out — the EOM-based ⁸⁵Rb "
        "twin-beam setup "
        "of Sim, Kim &amp; Moon, <em>Sci. Rep.</em> <strong>15</strong>, 7727 "
        "(2025), driven by the knobs you actually turn.</p></div>",
        unsafe_allow_html=True,
    )

    # ---- Tier B: the solve ----
    raw = _cached_solve(bridge.solve_key(params), fwm.FWMScheme().cache_version)
    result = replace(result, raw=raw)
    readout = detection_model.readout(
        result.geometry, result.chain, result.gains, settings, detection,
        calibration)
    result = replace(result, readout=readout)

    # ---- Tier C: readout ----
    gate = raw.get("claim_gate", {})
    metrics = [
        {"label": "Validation level",
         "value": gate.get("level", "MEAN_FIELD_DIAGNOSTIC"),
         "help": "Quantitative gain and physical squeezing claims are blocked."},
        {"label": "Gain-referred diagnostic",
         "value": f"{result.gain_referred_noise_db:+.2f} dB",
         "help": "Algebraic diagnostic at the selected two-photon detuning; the "
                 "microscopic atomic noise covariance is unavailable."},
        {"label": "RIN-loaded diagnostic",
         "value": f"{result.eom_rin_loaded_noise_db:+.2f} dB",
         "help": "Gain-only estimate after unpaired-mode shot noise and the "
                 "specified EOM RIN. This is not a squeezing prediction."},
        {"label": "Algebraic detection floor",
         "value": f"{result.gain_referred_detection_floor_db:+.2f} dB",
         "help": "10·log10(1−η) within the ideal-vacuum completion only; not a "
                 "physical squeezing bound for the current atomic model."},
        {"label": "Detector clearance",
         "value": f"{readout.clearance_db:.1f} dB",
         "help": "Shot noise above the amplifier noise floor. This is detector "
                 "headroom only and does not validate the modeled dB diagnostic."},
        {"label": "Carrier : seed",
         "value": (f"{params['eom_residual_carrier_to_wanted_ratio'] * 100:.3g} %"),
         "help": "Residual pump-frequency light in the seed mode after the "
                 "etalon chain, computed from the quantized powers passed to "
                 "GABES. It remains unapplied to the reduced atomic response; "
                 "the balanced-detector layer applies the separately declared "
                 "effective RIN."},
        {"label": "Pump gain G_s", "value": f"{result.gains[0]:.2f}"},
        {"label": "Conjugate gain G_c", "value": f"{result.gains[1]:.2f}"},
    ]
    _render_metrics(metrics)

    for warning in result.warnings:
        st.warning(warning, icon="⚠️")

    tabs = st.tabs(["Optical table", "Spectrum", "Derived parameters",
                    "Power budget", "Detection", "Calibration"])

    with tabs[0]:
        _render_table_tab(result, calibration, host, detection)

    with tabs[1]:
        observables = _cached_observables(
            raw, tuple(sorted(params.items(), key=lambda kv: kv[0])),
            fwm.FWMScheme().cache_version)
        figure = observables.get("figure")
        if figure is not None and host is not None:
            host.render_fig(figure)
        elif figure is not None:
            st.pyplot(figure)

    with tabs[2]:
        st.markdown(
            "Left column is what you set on the table; right column is what "
            "GABES receives. Every SABES number can be checked by driving GABES "
            "by hand with the right column."
        )
        st.markdown(_markdown_table(
            ("Setting", "Value", "GABES quantity", "Derived"),
            bridge.derived_table(result)))

    with tabs[3]:
        st.markdown(_budget_table(result.chain))
        st.markdown("#### EOM mode → etalon transfer")
        st.caption(
            "Each etalon attenuates each EOM line by the listed power "
            "transmission. The table reports √(T) because phase was not "
            "measured. Modes below 10⁻⁸ of EOM power are omitted.")
        st.markdown(_etalon_transfer_table(result.chain))

    with tabs[4]:
        geom = result.geometry
        eom_penalty = result.eom_rin_penalty_db
        eom_loaded_text = f"{result.eom_rin_loaded_noise_db:+.3g} dB"
        if eom_penalty is not None:
            eom_loaded_text += f" ({eom_penalty:+.3g} dB penalty)"
        else:
            eom_shift = (result.eom_rin_loaded_noise_db
                         - result.gain_referred_noise_db)
            eom_loaded_text += (
                f" ({eom_shift:+.3g} dB normalized shift; absolute EOM PSD "
                "is positive)")
        rows = [
            ("Twin radius at D-mirrors",
             f"{geom.twin_radius_at_dmirror_m * 1e6:.0f} µm"),
            ("Pump–twin separation",
             f"{geom.pump_separation_m * 1e3:.2f} mm"),
            ("Separation margin",
             f"{geom.separation_margin:.1f} twin radii"),
            ("Mean post-cell transmission",
             f"{geom.optics_transmission * 100:.2f} %"),
            ("Observed pump leakage",
             f"{detection.pump_leakage_dbm:.0f} dBm"),
            ("Effective unwanted-mode RIN",
             f"{readout.eom_noise.rin_db_per_hz:.1f} dBc/Hz"),
            ("Unwanted-mode fractional RMS",
             f"{100.0 * readout.eom_noise.fractional_intensity_rms:.3g} % "
             f"over {readout.eom_noise.analysis_bandwidth_hz / 1e3:.0f} kHz"),
            ("Unwanted EOM power at detector",
             f"{readout.eom_noise.unwanted_detector_power_w * 1e9:.3g} nW"),
            ("Classical EOM RIN excess",
             f"{readout.eom_noise.classical_rin_excess_sql:.3g} × SQL"),
            ("Shot-noise-only diagnostic",
             f"{result.eom_shot_noise_only_db:+.2f} dB"),
            ("RIN-loaded diagnostic",
             eom_loaded_text),
            ("Total power on the detector",
             f"{readout.total_power_w * 1e6:.1f} µW"),
            ("Margin above electronic noise",
             f"{readout.margin_above_electronic_db(result.gain_referred_noise_db):.1f} dB"),
        ]
        for arm in readout.arms:
            arm_geometry = next(item for item in geom.arms
                                if item.name == arm.name)
            rows += [
                (f"{arm.name}: post-cell transmission",
                 f"{arm_geometry.optics_transmission * 100:.2f} %"),
                (f"{arm.name}: power", f"{arm.power_w * 1e6:.1f} µW"),
                (f"{arm.name}: spot radius", f"{arm.spot_radius_m * 1e6:.0f} µm"),
                (f"{arm.name}: power density",
                 f"{arm.power_density_w_per_cm2:.2f} W/cm²"),
                (f"{arm.name}: residual pump",
                 f"{arm.residual_pump_w * 1e9:.2f} nW"),
            ]
        st.markdown(_markdown_table(("Quantity", "Value"), rows))
        st.caption(
            "The shipped RIN is an assumption. Unwanted modes pass through the "
            "probe arm without atomic gain or a conjugate partner. Treat the "
            "result as calibrated only after measuring detector-plane mode "
            "powers and RIN.")

    with tabs[5]:
        _provenance_panel(calibration)


def _render_metrics(metrics):
    cards = "".join(
        "<div class='sabes-metric' title='{help}'>"
        "<div class='sabes-metric-label'>{label}</div>"
        "<div class='sabes-metric-value'>{value}</div></div>".format(
            help=escape(str(m.get("help", ""))),
            label=escape(str(m["label"])),
            value=escape(str(m["value"])))
        for m in metrics)
    st.markdown(f"<section class='sabes-metrics'>{cards}</section>",
                unsafe_allow_html=True)


def _inject_css():
    st.markdown(
        """
        <style>
        .sabes-brand { padding: 0.35rem 0 0.6rem 0; }
        .sabes-brand-name { font-size: 2rem; font-weight: 800; letter-spacing: .06em;
            line-height: 1; color: #0F766E; }
        .sabes-brand-sub { font-size: .78rem; font-weight: 600; color: #475569;
            margin-top: .2rem; }
        .sabes-brand-note { font-size: .7rem; color: #64748B; margin-top: .25rem; }
        .sabes-header h1 { margin-bottom: .1rem; color: #0F766E; letter-spacing: .04em; }
        .sabes-header p { color: #475569; margin-top: 0; font-size: .92rem; }
        .sabes-metrics { display: grid; gap: .6rem; margin: .8rem 0 1.1rem 0;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); }
        .sabes-metric { border: 1px solid rgba(15,118,110,.22); border-radius: 12px;
            padding: .6rem .75rem; background: rgba(15,118,110,.05); }
        .sabes-metric-label { font-size: .72rem; text-transform: uppercase;
            letter-spacing: .06em; color: #0F766E; font-weight: 700; }
        .sabes-metric-value { font-size: 1.35rem; font-weight: 700; color: #0F172A;
            margin-top: .15rem; }
        .sabes-legend { display: flex; flex-wrap: wrap; gap: .1rem 1.1rem;
            margin: .35rem 0 .9rem 0; font-size: .78rem; color: #475569; }
        .sabes-key { display: inline-flex; align-items: center; gap: .35rem; }
        .sabes-key i { width: 18px; height: 4px; border-radius: 2px;
            display: inline-block; }
        .sabes-key--note { color: #94A3B8; font-style: italic; }
        .sabes-synth { border-left: 3px solid #B45309;
            background: #FEF3C7; color: #78350F; padding: .45rem .7rem;
            border-radius: 4px; font-size: .78rem; margin-bottom: .6rem; }
        </style>
        """,
        unsafe_allow_html=True,
    )
