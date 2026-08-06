"""
The Sim et al. table, as two layouts matching `references/IDS setup part *.png`.

Display coordinates are a clean schematic in the drawings' arrangement, not a
pixel trace of them: the topology, the ordering and the rough placement are
faithful, the exact positions are for Stage C to refine against the images. That
is a deliberate limit, not an oversight -- fake pixel precision would suggest the
drawing had been measured when it has not.

Path lengths are the opposite: every one that carries physics is either a design
value (a telescope is afocal, so its lenses sit `f1 + f2` apart by definition) or
a measurement (23 inch from the cell to the D-shaped mirrors), and the ones that
are neither say `nominal` out loud. Segments where the beam does not meaningfully
diverge carry `None` and exist only for the eye.

Corrections applied here against the part 1 drawing, which the earlier linear
model got wrong or omitted:
  - a Glan-Taylor after the three etalons and before the delivery fiber (the
    model had only the one in part 2, so it was incomplete rather than misplaced)
  - a half-wave plate between etalon 1 and etalon 2
  - distinct collimators: F230 on the pump, F220 on the seed
  - the noise eater, the optical isolators, the MOPA's cylindrical lens and the
    Rb-natural locking arm now exist as nodes, drawn but not modelled
"""
from .table import (BEAM_COMMON, BEAM_CONJUGATE, BEAM_PROBE, BEAM_PUMP,
                    BEAM_REFERENCE, BEAM_SEED, Layout, Link, Node,
                    PROVENANCE_DESIGN, PROVENANCE_MEASURED, PROVENANCE_NOMINAL)

# Distances the detection model needs. Kept here so the layout owns the
# geometry and `sabes.detection` stops carrying loose coefficients.
CELL_TO_DMIRROR_M = 0.5842        # 23 inch, measured
DMIRROR_TO_IRIS_M = 0.10          # nominal
IRIS_TO_LENS_M = 0.2158           # nominal; cell->lens then totals 0.90 m


def _optic(node_id, label, x, y, **kw):
    return Node(node_id, label, x, y, **kw)


def _decor(node_id, label, x, y, **kw):
    """Drawn, but the physics does not know about it -- so not clickable."""
    kw.setdefault("modelled", False)
    return Node(node_id, label, x, y, **kw)


def _mirror(node_id, x, y, angle=45.0):
    return _decor(node_id, "", x, y, shape="rect", w=22, h=6, angle=angle,
                  kind="mirror")


# ---------------------------------------------------------------------------
# Part 1 -- source, amplifier, split, seed modulation and filtering
# ---------------------------------------------------------------------------

_PART1_NODES = (
    # --- master laser and the locking arm ---
    _optic("ecdl", "795 nm ECDL", 880, 60, kind="source", w=86, h=52,
           modelled=True, params=("ecdl_power_mw",),
           note="Detuning is set by locking to the reference arm, not by "
                "predicting frequency from temperature / current / PZT."),
    _decor("oi_ecdl", "OI", 880, 140, shape="circle", w=22, h=22,
           note="Optical isolator. Loss only; not modelled."),
    _optic("hwp_lock", "HWP", 880, 186, w=8, h=24, modelled=False),
    _optic("pbs_lock", "PBS", 806, 214, w=26, h=26, modelled=False),
    _decor("lens_lock", "f=300", 806, 168, shape="circle", w=14, h=14),
    _decor("rb_ref", "Rb natural", 560, 56, kind="cell", w=74, h=26,
           note="Saturated-absorption reference. The oscilloscope demo lives "
                "here: this trace is what the detuning is actually set by."),
    _decor("pbs_ref", "PBS", 450, 56, w=22, h=22),
    _decor("pd_ref", "PD", 356, 56, kind="detector", w=24, h=34),
    _mirror("m_lock", 806, 56, angle=-45),

    # --- amplifier chain ---
    _mirror("m_ecdl_1", 806, 260, angle=45),
    _optic("hwp_mopa", "HWP", 742, 260, w=8, h=24, modelled=True,
           params=("hwp_ta_deg",)),
    _mirror("m_ecdl_2", 690, 292, angle=-45),
    _decor("cl", "CL", 500, 330, shape="circle", w=12, h=18,
           note="Cylindrical lens shaping the amplifier input. It is what sets "
                "the amplifier's output beam quality, i.e. ta_m2_out."),
    _optic("mopa", "MOPA", 620, 340, kind="source", w=110, h=54, modelled=True,
           params=("ta_current_a",)),
    _mirror("m_mopa_out", 436, 330, angle=45),
    _decor("oi_mopa", "OI", 350, 480, shape="circle", w=22, h=22),
    _optic("hwp_split", "HWP", 292, 480, w=8, h=24, modelled=True,
           params=("hwp_split_deg",),
           note="The pump-power knob: it sets how the amplifier output divides."),
    _mirror("m_split_1", 236, 480, angle=-45),
    _mirror("m_split_2", 236, 424, angle=45),
    _optic("pbs_split", "PBS", 176, 424, w=26, h=26, modelled=False,
           note="Pump transmits, seed reflects."),

    # --- pump arm out of part 1 ---
    _optic("f230_pump", "F230", 86, 424, kind="fiber", w=30, h=40,
           modelled=False, note="Pump delivery collimator; continues in part 2."),

    # --- seed arm: attenuate, modulate, filter ---
    _optic("pol_seed", "Polarizer", 176, 350, w=26, h=8, modelled=True,
           params=("seed_polarizer_deg",),
           note="Holds the fiber EOM under its ~10 mW rating. The split alone "
                "cannot satisfy both that and the pump power."),
    _decor("noise_eater", "Noise eater", 250, 306, w=62, h=24,
           note="Intensity-noise servo. Not modelled, but it acts on exactly "
                "the noise a squeezing measurement cares about."),
    _mirror("m_seed_1", 176, 258, angle=45),
    _optic("hwp_eom", "HWP", 240, 258, w=8, h=24, modelled=True,
           params=("hwp_eom_deg",)),
    _optic("f220_eom_in", "F220", 306, 258, kind="fiber", w=26, h=34,
           modelled=False),
    _optic("eom", "fiber EOM", 306, 190, w=92, h=22, modelled=True,
           params=("eom_offset_mhz", "eom_rf_dbm"),
           note="EOSPACE LiNbO3 phase modulator. Its drive frequency IS the "
                "two-photon detuning, exactly and with no calibration factor."),
    _optic("f220_eom_out", "F220", 400, 190, kind="fiber", w=26, h=34,
           modelled=False),
    _optic("etalon_1", "Etalon 1", 470, 190, w=48, h=30, modelled=True,
           params=("etalon_detune_ghz",)),
    _mirror("m_et_1", 560, 190, angle=45),
    _optic("hwp_etalon", "HWP", 560, 244, w=8, h=24, modelled=False,
           note="Between etalon 1 and etalon 2 in the drawing."),
    _mirror("m_et_2", 560, 292, angle=-45),
    _optic("etalon_2", "Etalon 2", 484, 292, w=48, h=30, modelled=True,
           params=("etalon_detune_ghz",)),
    _optic("etalon_3", "Etalon 3", 402, 292, w=48, h=30, modelled=True,
           params=("etalon_detune_ghz",)),
    _optic("gt_seed", "GT", 330, 292, w=18, h=28, modelled=True,
           params=("seed_gtp_deg",),
           note="Glan-Taylor after the filter chain. Part 2 has further GTPs on "
                "each arm; the model previously had only one."),
    _mirror("m_seed_out", 262, 292, angle=45),
    _mirror("m_seed_down", 262, 396, angle=-45),
    _optic("f220_seed", "F220", 86, 350, kind="fiber", w=30, h=40,
           modelled=False, note="Seed delivery collimator; continues in part 2."),
)

_PART1_LINKS = (
    Link("ecdl", "oi_ecdl", BEAM_COMMON, carries="ECDL output"),
    Link("oi_ecdl", "hwp_lock", BEAM_COMMON, carries="ECDL output"),
    Link("hwp_lock", "pbs_lock", BEAM_COMMON, carries="ECDL output"),
    # locking arm
    Link("pbs_lock", "lens_lock", BEAM_REFERENCE),
    Link("lens_lock", "m_lock", BEAM_REFERENCE),
    Link("m_lock", "rb_ref", BEAM_REFERENCE),
    Link("rb_ref", "pbs_ref", BEAM_REFERENCE),
    Link("pbs_ref", "pd_ref", BEAM_REFERENCE),
    # to the amplifier
    Link("pbs_lock", "m_ecdl_1", BEAM_COMMON, carries="ECDL output"),
    Link("m_ecdl_1", "hwp_mopa", BEAM_COMMON, carries="to amplifier"),
    Link("hwp_mopa", "m_ecdl_2", BEAM_COMMON, carries="to amplifier"),
    Link("m_ecdl_2", "mopa", BEAM_COMMON, via=((690, 340),),
         carries="to amplifier"),
    Link("cl", "mopa", BEAM_COMMON, carries="to amplifier"),
    # amplifier out, to the split
    Link("mopa", "m_mopa_out", BEAM_COMMON, carries="tapered amplifier"),
    Link("m_mopa_out", "oi_mopa", BEAM_COMMON, via=((436, 480),),
         carries="tapered amplifier"),
    Link("oi_mopa", "hwp_split", BEAM_COMMON, carries="tapered amplifier"),
    Link("hwp_split", "m_split_1", BEAM_COMMON, carries="tapered amplifier"),
    Link("m_split_1", "m_split_2", BEAM_COMMON, carries="tapered amplifier"),
    Link("m_split_2", "pbs_split", BEAM_COMMON, carries="tapered amplifier"),
    # pump arm
    Link("pbs_split", "f230_pump", BEAM_PUMP, carries="PBS transmitted"),
    # seed arm
    Link("pbs_split", "pol_seed", BEAM_SEED, carries="PBS reflected"),
    Link("pol_seed", "m_seed_1", BEAM_SEED, carries="input polarizer"),
    Link("m_seed_1", "hwp_eom", BEAM_SEED, carries="input polarizer"),
    Link("hwp_eom", "f220_eom_in", BEAM_SEED, carries="input polarizer"),
    Link("f220_eom_in", "eom", BEAM_SEED, carries="into fiber EOM"),
    Link("eom", "f220_eom_out", BEAM_SEED, carries="EOM sidebands"),
    Link("f220_eom_out", "etalon_1", BEAM_SEED, carries="EOM sidebands"),
    Link("etalon_1", "m_et_1", BEAM_SEED, carries="etalon 1"),
    Link("m_et_1", "hwp_etalon", BEAM_SEED, carries="etalon 1"),
    Link("hwp_etalon", "m_et_2", BEAM_SEED, carries="etalon 1"),
    Link("m_et_2", "etalon_2", BEAM_SEED, carries="etalon 1"),
    Link("etalon_2", "etalon_3", BEAM_SEED, carries="etalon 2"),
    Link("etalon_3", "gt_seed", BEAM_SEED, carries="etalon 3"),
    Link("gt_seed", "m_seed_out", BEAM_SEED, carries="etalon 3"),
    Link("m_seed_out", "m_seed_down", BEAM_SEED, carries="etalon 3"),
    Link("m_seed_down", "f220_seed", BEAM_SEED, via=((160, 396),),
         carries="etalon 3"),
)

PART1 = Layout(
    key="part1",
    title="Part 1 — source, amplifier, seed modulation and filtering",
    width=960, height=520,
    nodes=_PART1_NODES,
    links=_PART1_LINKS,
)


# ---------------------------------------------------------------------------
# Part 2 -- delivery, the cell, and detection
# ---------------------------------------------------------------------------

_PART2_NODES = (
    # --- seed delivery (red in the drawing) ---
    _optic("col_seed", "5724 · 1.6 mm", 890, 60, kind="fiber", w=34, h=42,
           modelled=False, note="Seed collimator. Its diameter times the "
                                "telescope magnification is the cell waist."),
    _optic("qwp_seed", "QWP", 836, 60, w=8, h=24, modelled=True,
           params=("seed_trim_qwp_deg",)),
    _optic("hwp_seed", "HWP", 812, 60, w=8, h=24, modelled=True,
           params=("seed_trim_hwp_deg",),
           note="With the Glan-Taylor, sets the seed power at the cell."),
    _optic("seed_l1", "f=400", 720, 60, shape="circle", w=16, h=16,
           modelled=True, params=("seed_telescope_f1_mm",)),
    _optic("seed_l2", "f=150", 330, 60, shape="circle", w=16, h=16,
           modelled=True, params=("seed_telescope_f2_mm",)),
    _mirror("m_seed_a", 250, 60, angle=45),
    _optic("gtp_seed", "GTP", 250, 132, w=18, h=28, modelled=False),

    # --- pump delivery (orange in the drawing) ---
    _optic("col_pump", "5725 · 2.0 mm", 890, 152, kind="fiber", w=34, h=42,
           modelled=False),
    _optic("qwp_pump", "QWP", 836, 152, w=8, h=24, modelled=False),
    _optic("hwp_pump", "HWP", 812, 152, w=8, h=24, modelled=False),
    _optic("pump_l1", "f=300", 720, 152, shape="circle", w=16, h=16,
           modelled=True, params=("pump_telescope_f1_mm",)),
    _optic("pump_l2", "f=150", 400, 152, shape="circle", w=16, h=16,
           modelled=True, params=("pump_telescope_f2_mm",)),
    _mirror("m_pump_a", 330, 152, angle=45),
    _optic("gtp_pump", "GTP", 330, 208, w=18, h=28, modelled=False),

    # --- the cell ---
    _optic("pbs_in", "PBS", 250, 256, w=26, h=26, modelled=False,
           note="Combines pump and seed at the crossing angle."),
    _optic("cell", "Rb cell", 372, 256, kind="cell", w=76, h=32, modelled=True,
           params=("cell_temp_c", "opd_ghz"),
           note="Where the four-wave mixing happens. 12.5 mm."),
    _optic("pbs_out", "PBS", 496, 256, w=26, h=26, modelled=True,
           params=(), note="Polarization rejection of the pump, 3000:1."),
    _decor("dump", "Beam block", 496, 322, w=54, h=22),

    # --- twin-beam separation and detection ---
    _optic("dm_probe", "D-mirror", 700, 214, w=24, h=8, angle=45,
           modelled=True, params=("dmirror_separation_mm",),
           note="Alignment is done by overlapping at the cell and separating "
                "by a stated number of mm here, so this is the crossing angle."),
    _optic("dm_conj", "D-mirror", 700, 300, w=24, h=8, angle=-45,
           modelled=True, params=("dmirror_separation_mm",)),
    _optic("iris_probe", "Iris", 790, 214, shape="circle", w=16, h=16,
           modelled=True, params=("iris_radius_mm",)),
    _optic("iris_conj", "Iris", 790, 300, shape="circle", w=16, h=16,
           modelled=True, params=("iris_radius_mm",)),
    _optic("lens_probe", "f=75", 862, 214, shape="circle", w=14, h=14,
           modelled=True, params=("probe_lens_focal_mm",)),
    _optic("lens_conj", "f=100", 862, 300, shape="circle", w=14, h=14,
           modelled=True, params=("conjugate_lens_focal_mm",)),
    _optic("bpd", "Balanced detector", 920, 257, kind="detector", w=54, h=76,
           modelled=True, params=("pd_defocus_mm",),
           note="PDB450A with Hamamatsu S3883 photodiodes fitted."),
)

_PART2_LINKS = (
    # seed delivery
    Link("col_seed", "qwp_seed", BEAM_SEED, carries="delivery fiber + collimator"),
    Link("qwp_seed", "hwp_seed", BEAM_SEED, carries="delivery fiber + collimator"),
    Link("hwp_seed", "seed_l1", BEAM_SEED, path_m=0.400,
         provenance=PROVENANCE_DESIGN, carries="collimator",
         note="Waist sits at L1's front focal plane in the afocal arrangement."),
    Link("seed_l1", "seed_l2", BEAM_SEED, path_m=0.550,
         provenance=PROVENANCE_DESIGN, carries="telescope L1",
         note="f1 + f2 = 400 + 150 mm, by definition of an afocal telescope."),
    Link("seed_l2", "m_seed_a", BEAM_SEED, path_m=0.150,
         provenance=PROVENANCE_DESIGN, carries="trim + telescope"),
    Link("m_seed_a", "gtp_seed", BEAM_SEED, carries="at cell"),
    Link("gtp_seed", "pbs_in", BEAM_SEED, carries="at cell"),

    # pump delivery
    Link("col_pump", "qwp_pump", BEAM_PUMP, carries="delivery fiber"),
    Link("qwp_pump", "hwp_pump", BEAM_PUMP, carries="delivery fiber"),
    Link("hwp_pump", "pump_l1", BEAM_PUMP, path_m=0.300,
         provenance=PROVENANCE_DESIGN, carries="collimator"),
    Link("pump_l1", "pump_l2", BEAM_PUMP, path_m=0.450,
         provenance=PROVENANCE_DESIGN, carries="telescope L1",
         note="f1 + f2 = 300 + 150 mm."),
    Link("pump_l2", "m_pump_a", BEAM_PUMP, path_m=0.150,
         provenance=PROVENANCE_DESIGN, carries="collimator + telescope"),
    Link("m_pump_a", "gtp_pump", BEAM_PUMP, carries="at cell"),
    Link("gtp_pump", "pbs_in", BEAM_PUMP, carries="at cell"),

    # through the cell
    Link("pbs_in", "cell", BEAM_PUMP, path_m=0.05,
         provenance=PROVENANCE_NOMINAL, carries="at cell"),
    Link("cell", "pbs_out", BEAM_PUMP, path_m=0.10,
         provenance=PROVENANCE_NOMINAL, carries="at cell"),
    Link("pbs_out", "dump", BEAM_PUMP, carries="at cell"),

    # twins to the detector
    Link("pbs_out", "dm_probe", BEAM_PROBE, path_m=CELL_TO_DMIRROR_M - 0.10,
         provenance=PROVENANCE_MEASURED, carries="probe",
         note="Cell to the D-shaped mirrors is 23 inch; this is the part after "
              "the rejection PBS."),
    Link("pbs_out", "dm_conj", BEAM_CONJUGATE, path_m=CELL_TO_DMIRROR_M - 0.10,
         provenance=PROVENANCE_MEASURED, carries="conjugate"),
    Link("dm_probe", "iris_probe", BEAM_PROBE, path_m=DMIRROR_TO_IRIS_M,
         provenance=PROVENANCE_NOMINAL, carries="probe"),
    Link("dm_conj", "iris_conj", BEAM_CONJUGATE, path_m=DMIRROR_TO_IRIS_M,
         provenance=PROVENANCE_NOMINAL, carries="conjugate"),
    Link("iris_probe", "lens_probe", BEAM_PROBE, path_m=IRIS_TO_LENS_M,
         provenance=PROVENANCE_NOMINAL, carries="probe"),
    Link("iris_conj", "lens_conj", BEAM_CONJUGATE, path_m=IRIS_TO_LENS_M,
         provenance=PROVENANCE_NOMINAL, carries="conjugate"),
    Link("lens_probe", "bpd", BEAM_PROBE, carries="probe"),
    Link("lens_conj", "bpd", BEAM_CONJUGATE, carries="conjugate"),
)

PART2 = Layout(
    key="part2",
    title="Part 2 — delivery, the vapour cell, and detection",
    width=960, height=400,
    nodes=_PART2_NODES,
    links=_PART2_LINKS,
)


LAYOUTS = (PART1, PART2)


def get(key) -> Layout:
    for layout in LAYOUTS:
        if layout.key == key:
            return layout
    raise KeyError(f"no layout {key!r}; have {[l.key for l in LAYOUTS]}")


#: Node route from the cell to each detection landmark. The detection model asks
#: the layout for these rather than carrying loose distance coefficients, so a
#: table change and a physics change are the same edit.
ROUTE_TO_DMIRROR = ("cell", "pbs_out", "dm_probe")
ROUTE_TO_LENS = ("cell", "pbs_out", "dm_probe", "iris_probe", "lens_probe")


def cell_to_dmirror_m(layout=None):
    """Measured 23 inch, summed over the layout rather than quoted."""
    return (layout or PART2).path_length_m(*ROUTE_TO_DMIRROR)


def cell_to_lens_m(layout=None):
    """Cell to the focusing lens, summed over the layout rather than guessed."""
    return (layout or PART2).path_length_m(*ROUTE_TO_LENS)


def with_detection_distances(cell_to_dmirror=None, dmirror_to_iris=None,
                             iris_to_lens=None, layout=None):
    """PART2 with the post-cell distances replaced.

    Two of these three are nominal today, so this is how a measurement gets
    folded in -- and, in the meantime, how their influence is swept.
    """
    from dataclasses import replace as _replace

    base = layout or PART2
    overrides = {}
    if cell_to_dmirror is not None:
        # The measured 23 inch spans the cell through the rejection PBS, so the
        # override has to land on the segment after it.
        after_pbs = base.link_between("cell", "pbs_out").path_m or 0.0
        overrides[("pbs_out", "dm_probe")] = float(cell_to_dmirror) - after_pbs
        overrides[("pbs_out", "dm_conj")] = float(cell_to_dmirror) - after_pbs
    if dmirror_to_iris is not None:
        overrides[("dm_probe", "iris_probe")] = float(dmirror_to_iris)
        overrides[("dm_conj", "iris_conj")] = float(dmirror_to_iris)
    if iris_to_lens is not None:
        overrides[("iris_probe", "lens_probe")] = float(iris_to_lens)
        overrides[("iris_conj", "lens_conj")] = float(iris_to_lens)

    links = tuple(
        _replace(link, path_m=overrides[(link.a, link.b)],
                 provenance=PROVENANCE_MEASURED)
        if (link.a, link.b) in overrides else link
        for link in base.links
    )
    return _replace(base, links=links)
