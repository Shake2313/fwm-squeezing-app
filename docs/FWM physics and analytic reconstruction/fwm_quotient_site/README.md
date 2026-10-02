# GABES FWM Analytic Atlas

Static deployment of the standalone v3 atlas. No runtime network dependencies,
database, authentication code, or physics computation. Sites controls visitor
access. HTML export preserves the complete scientific text and graph data.

`python build.py` builds `out/index.html` from committed `index.html`.
`python build.py --refresh` updates the source from the parent v3 HTML first.
Before refreshing, run the parent `fwm_quotient_structure_v3_builder.py` after
editing its graph data or `fwm_quotient_structure_v3_view.css` / `.js`.

Local provenance links resolve to the embedded reference catalog in this hosted
edition. The original local HTML retains its optional source-file links.
