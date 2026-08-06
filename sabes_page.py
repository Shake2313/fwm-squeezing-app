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

import streamlit as st

from gabes import constants
from gabes.core import blas_single_thread
from gabes.schemes import fwm

from sabes import bridge, detection as detection_model, layout
from sabes.components import canvas as canvas_module
from sabes.beamline import (SetupSettings, build_source_chain,
                            solve_seed_polarizer_deg, solve_seed_trim_deg,
                            solve_split_angle_deg)
from sabes.calibration import default_calibration
from sabes.detection import DetectionSettings

SESSION_PREFIX = "_sabes_"
SETTINGS_VERSION = "sabes-p4-v1"

FIDELITIES = (fwm.FIDELITY_FAST, fwm.FIDELITY_ULTRA)

#: (session key, label, unit, min, max, step, help). Grouped exactly the way the
#: table is laid out, because that is how someone walks the setup.
CONTROL_GROUPS = (
    ("Laser & amplifier", (
        ("ecdl_power_mw", "ECDL output", "mW", 1.0, 200.0, 1.0,
         "Seed power into the tapered amplifier."),
        ("ta_current_a", "TA current", "A", 0.5, 4.0, 0.05,
         "Amplifier drive current. Output saturates, so this is not a linear knob."),
        ("hwp_ta_deg", "HWP before TA", "°", 0.0, 90.0, 0.5,
         "Rotates the seed onto the amplifier's gain axis."),
    )),
    ("Pump / seed split", (
        ("hwp_split_deg", "Split HWP", "°", 0.0, 45.0, 0.01,
         "Sets how the amplifier output divides between pump and seed. This is "
         "the pump-power knob."),
        ("seed_polarizer_deg", "Seed polarizer", "°", 0.0, 90.0, 0.05,
         "Holds the fiber EOM below its optical rating. Independent of the "
         "split, because the split alone cannot satisfy both constraints."),
    )),
    ("Seed modulation", (
        ("eom_offset_mhz", "Generator offset from ν_HF", "MHz", -200.0, 200.0, 0.1,
         "Signal-generator frequency relative to the 3.0357324 GHz ground "
         "hyperfine splitting. This IS the two-photon detuning, with the sign "
         "flipped and no calibration coefficient in between."),
        ("eom_rf_dbm", "EOM RF drive", "dBm", 0.0, 30.0, 0.5,
         "Sets the modulation index β = π V_peak / V_π and hence how much power "
         "lands in the wanted −1 sideband."),
        ("hwp_eom_deg", "HWP before EOM", "°", 0.0, 90.0, 0.5,
         "Aligns the input to the modulator's polarization axis."),
        ("etalon_detune_ghz", "Etalon detuning", "GHz", -2.0, 2.0, 0.01,
         "Common tilt/temperature detuning of all three filter stages from the "
         "wanted sideband."),
    )),
    ("Seed trim", (
        ("seed_trim_hwp_deg", "Trim HWP", "°", 0.0, 45.0, 0.01,
         "With the Glan-Taylor, sets the seed power at the cell. The seed arm "
         "carries far more than the FWM needs, so this mostly throws power away."),
        ("seed_trim_qwp_deg", "Trim QWP", "°", 0.0, 90.0, 0.5,
         "Fine polarization control ahead of the Glan-Taylor."),
        ("seed_gtp_deg", "Glan-Taylor", "°", 0.0, 90.0, 0.5,
         "Analyser angle for the trim pair."),
    )),
    ("Beam shaping", (
        ("pump_telescope_f1_mm", "Pump telescope L1", "mm", 25.0, 1000.0, 5.0,
         "First lens. Cell waist = collimator diameter × f2/f1 × scale factor."),
        ("pump_telescope_f2_mm", "Pump telescope L2", "mm", 25.0, 1000.0, 5.0,
         "Second lens."),
        ("seed_telescope_f1_mm", "Seed telescope L1", "mm", 25.0, 1000.0, 5.0,
         "First lens of the seed telescope."),
        ("seed_telescope_f2_mm", "Seed telescope L2", "mm", 25.0, 1000.0, 5.0,
         "Second lens of the seed telescope."),
    )),
    ("Cell & geometry", (
        ("opd_ghz", "One-photon detuning Δ", "GHz", -3.0, 3.0, 0.05,
         "Where the pump sits. Set by the ECDL lock, not modelled from its "
         "temperature / current / PZT — mode hops make that unpredictable."),
        ("cell_temp_c", "Cell temperature", "°C", 60.0, 150.0, 1.0,
         "Vapour temperature from the heater controller."),
        ("dmirror_separation_mm", "Separation at D-mirrors", "mm", 0.5, 15.0, 0.05,
         "Alignment is done by overlapping at the cell and separating by this "
         "much 23 inch downstream, so this is the crossing-angle knob."),
    )),
    ("Detection", (
        ("iris_radius_mm", "Iris radius", "mm", 0.2, 6.0, 0.05,
         "Spatial filter for pump rejection. Closing it past the twin beam "
         "turns straight into detection loss."),
        ("probe_lens_focal_mm", "Probe lens", "mm", 25.0, 500.0, 5.0,
         "Focusing lens onto the probe photodiode."),
        ("conjugate_lens_focal_mm", "Conjugate lens", "mm", 25.0, 500.0, 5.0,
         "Focusing lens onto the conjugate photodiode."),
        ("pd_defocus_mm", "Photodiode defocus", "mm", 0.0, 100.0, 1.0,
         "Distance of the diode from the focal plane. Spreading the spot is the "
         "cheapest way out of the linearity limit."),
    )),
)

DETECTION_KEYS = {"iris_radius_mm", "probe_lens_focal_mm",
                  "conjugate_lens_focal_mm", "pd_defocus_mm"}


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
    values["etalon_detune_ghz"] = float(base.etalon_detune_ghz[0])
    values["resolution"] = FIDELITIES[0]
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
    detune = float(get("etalon_detune_ghz"))
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
        etalon_detune_ghz=(detune, detune, detune),
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
        iris_radius_mm=get("iris_radius_mm"),
        probe_lens_focal_mm=get("probe_lens_focal_mm"),
        conjugate_lens_focal_mm=get("conjugate_lens_focal_mm"),
        pd_defocus_mm=get("pd_defocus_mm"),
    )
    return settings, detection


# ----------------------------------------------------------------------
# Sidebar
# ----------------------------------------------------------------------
def _apply_solved_knobs():
    """Re-solve the three power knobs onto the paper operating point."""
    settings, _ = _current()
    calibration = default_calibration()
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


def _leave_dev():
    st.query_params.pop("dev", None)


def _reset_defaults():
    for name, value in _defaults().items():
        st.session_state[_key(name)] = value


def _render_sidebar(host):
    st.sidebar.markdown(
        "<div class='sabes-brand'>"
        "<div class='sabes-brand-name'>SABES</div>"
        "<div class='sabes-brand-sub'>Specific Atomic Bloch Equation Solver</div>"
        "<div class='sabes-brand-note'>EOM-based ⁸⁵Rb intensity-difference "
        "squeezing setup</div></div>",
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

    for group, controls in CONTROL_GROUPS:
        with st.sidebar.expander(group, expanded=(group == "Pump / seed split")):
            for name, label, unit, low, high, step, help_text in controls:
                st.number_input(
                    f"{label} [{unit}]" if unit else label,
                    min_value=float(low), max_value=float(high),
                    step=float(step), key=_key(name), help=help_text,
                    format="%.4f" if step < 0.05 else "%.2f")


# ----------------------------------------------------------------------
# Main panel
# ----------------------------------------------------------------------
def _twin_states(result):
    """Probe and conjugate as drawable beams, once the solve knows their gain.

    Scaling the seed state is right for a stroke width and only that: the twins
    do not carry the seed's residual carrier line. Real per-arm numbers come from
    the detection readout, which is what the panels show.
    """
    if result.raw is None:
        return None
    gain_s, gain_c = result.gains
    transmission = result.geometry.optics_transmission
    seed = result.chain.seed
    return {"probe": seed.scaled(gain_s * transmission),
            "conjugate": seed.scaled(gain_c * transmission)}


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
    else:
        st.session_state[_probe_key(part_key)] = [event["x"], event["y"]]
        st.session_state[_key("selected")] = None
    return True


def _render_optic_panel(part, node_id, settings, detection):
    """What the selected optic is and what it owns.

    Read-only for now: turning these into editors is the next step, and showing
    the values first makes it obvious the picture and the model agree.
    """
    try:
        node = part.node(node_id)
    except KeyError:
        return
    values = {f.name: getattr(settings, f.name) for f in fields(settings)}
    values.update({f.name: getattr(detection, f.name) for f in fields(detection)})
    values["eom_offset_mhz"] = (settings.eom_frequency_hz
                                - constants.NU_GROUND_HF) / 1e6
    values["etalon_detune_ghz"] = settings.etalon_detune_ghz[0]

    st.markdown(f"**{node.label or node.id}**")
    if node.note:
        st.caption(node.note)
    rows = [(name, f"{values.get(name, float('nan')):.4g}")
            for name in node.params]
    if rows:
        st.markdown(_markdown_table(("Parameter", "Value"), rows))
        st.caption("Editing these in place is the next stage; for now they live "
                   "in the sidebar.")


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


def _render_table_tab(result, host):
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
        st.markdown("###### Probe")
        reading = None
        if probe:
            reading = layout.beam_at(part, probe[0], probe[1], result.chain,
                                     twin_states=twin)
        _render_probe_panel(reading)


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

    # Development route for the canvas spike (Stage B of the optical-table
    # plan). It rides the deployment so the component can be verified where it
    # will actually run, but stays off the normal page. Remove with the spike.
    if st.query_params.get("dev") == "canvas":
        from sabes.components import spike
        st.sidebar.button("← Back to SABES", on_click=_leave_dev,
                          use_container_width=True)
        spike.render()
        return
    _seed_session_state()
    _render_sidebar(host)

    settings, detection = _current()
    calibration = default_calibration()
    fidelity = st.session_state[_key("resolution")]

    # ---- Tier A: optics ----
    result = bridge.run(settings, calibration, detection, solve=False,
                        resolution=fidelity)
    params = result.params

    st.markdown(
        "<div class='sabes-header'><h1>SABES</h1>"
        "<p>Lab settings in, squeezing out — the EOM-based ⁸⁵Rb twin-beam setup "
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
    metrics = [
        {"label": "Squeezing", "value": f"{result.squeezing_db:+.2f} dB",
         "help": "Intensity-difference squeezing at the two-photon detuning the "
                 "signal generator is set to."},
        {"label": "Loss-limited floor",
         "value": f"{result.squeezing_floor_db:+.2f} dB",
         "help": "10·log10(1−η). No operating point can beat this; detector QE "
                 "is the highest-leverage number in the setup."},
        {"label": "Detector clearance",
         "value": f"{readout.clearance_db:.1f} dB",
         "help": "Shot noise above the amplifier noise floor. A squeezed trace "
                 "sits |S| dB below shot noise, so this must exceed |S| with "
                 "margin or the measurement cannot resolve it."},
        {"label": "Carrier : seed",
         "value": f"{result.chain.carrier_ratio * 100:.4f} %",
         "help": "Residual pump-frequency light in the seed mode after the "
                 "etalon chain. Sim et al. Fig. 5 shows squeezing degrading as "
                 "this rises."},
        {"label": "Pump gain G_s", "value": f"{result.gains[0]:.2f}"},
        {"label": "Conjugate gain G_c", "value": f"{result.gains[1]:.2f}"},
    ]
    _render_metrics(metrics)

    for warning in result.warnings:
        st.warning(warning, icon="⚠️")

    tabs = st.tabs(["Optical table", "Spectrum", "Derived parameters",
                    "Power budget", "Detection", "Calibration"])

    with tabs[0]:
        _render_table_tab(result, host)

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

    with tabs[4]:
        geom = result.geometry
        rows = [
            ("Twin radius at D-mirrors",
             f"{geom.twin_radius_at_iris_m * 1e6:.0f} µm"),
            ("Pump–twin separation",
             f"{geom.pump_separation_m * 1e3:.2f} mm"),
            ("Separation margin",
             f"{geom.separation_margin:.1f} twin radii"),
            ("Iris transmission",
             f"{geom.arms[0].iris_transmission * 100:.3f} %"),
            ("Total power on the detector",
             f"{readout.total_power_w * 1e6:.1f} µW"),
            ("Margin above electronic noise",
             f"{readout.margin_above_electronic_db(result.squeezing_db):.1f} dB"),
        ]
        for arm in readout.arms:
            rows += [
                (f"{arm.name}: power", f"{arm.power_w * 1e6:.1f} µW"),
                (f"{arm.name}: spot radius", f"{arm.spot_radius_m * 1e6:.0f} µm"),
                (f"{arm.name}: power density",
                 f"{arm.power_density_w_per_cm2:.2f} W/cm²"),
                (f"{arm.name}: residual pump",
                 f"{arm.residual_pump_w * 1e9:.2f} nW"),
            ]
        st.markdown(_markdown_table(("Quantity", "Value"), rows))

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
        </style>
        """,
        unsafe_allow_html=True,
    )
