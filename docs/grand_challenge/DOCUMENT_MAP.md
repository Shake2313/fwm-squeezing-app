# Grand Challenge document map

[Entry](README.md) · [Next session](NEXT_SESSION.md)

Route by question; read only relevant rows. All 46 preexisting Markdown files
listed once below. Snapshot dates and bound evidence outrank historical “latest”
headers. Current design/decisions belong in the living blueprint; official status
in checklist; NEXT_SESSION is its compact execution handoff. The 46-file inventory
below is the pre-audit inventory, not a count of all current files.
Filenames alone do not identify campaign versions.

## Orientation and policy

- [Blueprint](blueprint.md): living master design, scoped package closures, dependency
  gates, current pilot and mandatory end-turn PI review; original technical design
  retained below.
- [Conventions](conventions.md): frequency, reservoir, channel, normalization contracts.
- [Publication policy](publication_policy.md): scientific ownership, preservation, release rules.

## Atomic, field, and readout foundations

- [Atomic diffusion](derivation.md): explicit reservoirs, ordered spectra, independent QRT.
- [Normalization](normalization_derivation.md): pump/weak-field dipoles, manifold populations.
- [Field coupling](field_derivation.md): photon-flux M/D, segment propagation.
- [Readout](readout_derivation.md): four-sideband covariance, temporal modes, intensity difference.
- [Seed validity](seed_validity_derivation.md): omitted bright-carrier terms, finite-seed diagnostics.
- [Input uncertainty](uncertainty_derivation.md): correlated independent inputs, first-order propagation.

## Periodic, spatial, and kinetic models

- [Periodic noise](periodic_noise_derivation.md): finite-seed microscopic atomic noise.
- [Periodic field](periodic_field_derivation.md): harmonic-to-port coupling, intensity readout.
- [RF bands](spectrum_analysis_derivation.md): squeezing-band extraction and refinement.
- [Kinetics](kinetic_derivation.md): velocity-resolved mean/noise, carrier geometry.
- [Spatial atoms](spatial_derivation.md): noncollinear phases, local microscopic noise.
- [Spatial field](spatial_field_derivation.md): finite-area optical modes, exact phase quotient.
- [Spatial uncertainty](spatial_uncertainty_derivation.md): conditional two-input propagation.

## Transport and ensemble methods

- [Inflow](inflow_derivation.md): boundary flux, chords, common entry phase.
- [Transport](transport_derivation.md): atomic characteristics, boundary/intersegment noise.
- [Moving Rb](rb_transport_derivation.md): reduced Rb characteristics, thermal inflow.
- [Passive bridge](ballistic_linear_channel.md): finite-velocity bosonic channel reference.
- [Segmented QRT](segmented_qrt_reference.md): independent full-density reference.
- [Smooth QRT](smooth_qrt_reference.md): continuous-envelope full-density reference.
- [Smooth transport](smooth_transport_derivation.md): Gaussian characteristics, sparse moments.
- [Adjoint transport](adjoint_transport_derivation.md): independent source and response checks.
- [Exponential transport](exponential_transport_derivation.md): eigenmode integration, CF4.
- [Transport ensemble](transport_ensemble_derivation.md): common lab RF, packet conventions.
- [Ensemble refinement](transport_ensemble_refinement.md): constant-atom integral, stronger checks.
- [Thermal Rb ensemble](rb_thermal_ensemble.md): selected-path reuse and convergence evidence.

## Campaign contracts and dated results

- [Portable campaign](portable_thermal_campaign.md): frozen-source execution contract. Opening “latest p2” summary is historical; use handoff for current state. Initial declaration describes v1.
- [Grid execution](thermal_grid_execution.md): batch, aggregation, nested-reuse contracts and commands.
- [Grid report v1](thermal_campaign_grid_v1.md): v1 campaign p1/seed11 snapshot.
- [Grid report v2](thermal_campaign_grid_v2.md): v1 campaign p2/seed11 snapshot; filename `v2` is report version, **not** campaign v2.
- [Seed211 v1](thermal_campaign_seed211_v1.md): v1 independent-seed and extension snapshot.
- [P2 seed211 v2](thermal_campaign_p2_seed211_v2.md): v2 campaign p2/seed211 snapshot.
- [P2 seed811 v2](thermal_campaign_p2_seed811_v2.md): v2 campaign p2 three-seed snapshot.

## Chronology, checkpoints, and plans

- [October 2 governance/science audit](audit_2026_10_02/README.md): independent
  findings, closed scopes, frozen-source geometry pilot forecast and verification.
- [October 2 code review](code_review_2026_10_02/README.md): storage interruption,
  byte/snapshot binding repairs, clean-candidate regression and commit scope.
- [Research log](research_log.md): append-only chronology; search dates/topics, avoid full startup read.
- [Parallel development](parallel_development.md): dated RF-band/input-uncertainty development record.
- [Gain hotfix handoff](gain_hotfix_handoff.md): calibration caveats, observable/input contracts.
- [September 14 checkpoint](checkpoint_20260914.md): historical resume point.
- [September 18 checkpoint](checkpoint_2026_09_18.md): historical independent review.
- [September 28 progress](progress_2026_09_28/README.md): dated execution, diagnostics, tests, evidence links.
- [Boundary-history plan](progress_2026_09_28/boundary_history_plan.md): proposed pre-entry history experiment; plan, not completed validation.

## Publication preservation records

- [Refresh index](publication_refresh_2026_09_18/README.md): dated publication update and verification.
- [Analytic coverage](publication_refresh_2026_09_18/analytic_coverage.md): preserved theory, restorations, checks.
- [Quotient coverage](publication_refresh_2026_09_18/quotient_coverage.md): logical structure, corrections, evidence.
- [Squeezing coverage](publication_refresh_2026_09_18/squeezing_coverage.md): content inventory and editorial rationale.

## Artifact routing

| Location | Role; read when needed |
|---|---|
| [Current publications](current_publications.json) | Active scientific source/artifact paths; resolve here. |
| Root `*_report*.json`, `*.png` | Recorded scientific snapshots and figures; follow relevant derivation links. |
| [Campaign v1](thermal_campaign_v1/), [campaign v2](thermal_campaign_v2/) | Plans, source manifests, frozen `sources.zip`, archived `controllers/`, path/batch/grid/comparison/ensemble reports. Controller bytes bind historical runs. Read declared plan and specific report first. |
| Campaign `cache/` | Large native packets; targeted keys only. Bulk reads add cost without answering status questions. |
| [Earlier path cache](rb_thermal_path_cache_v1/), [reference jobs](rb_thermal_reference_jobs_v1/) | Earlier selected-path/reference evidence; preserve historical identities. |
| [September 28 artifacts](progress_2026_09_28/) | Sealed diagnostics, independent checks, figures, logs, byte inventories. |
| Progress `*.txt`, root `plot_*.py` | Some `.txt` files are executable Python helpers. Inspect code/arguments before running; use fresh outputs, preserve prior records. |

No bulk cache loading, ZIP extraction, controller execution, or snapshot rewriting
needed for orientation. Internal record seals establish content consistency;
they do not replace native audits or scientific convergence gates.
