# Diagram Brief

## User Goal
- Output: editable draw.io workflow plus SVG and PNG for the GitHub landing page.
- Audience: manuscript readers and analysts reusing BRS-Map.
- Must communicate: 7 datasets/16 contrasts, three independent evidence axes,
  signed UCell, two same-Celltype strict-100 weighted-kNN rounds, direction
  consistency, primary-first selection and final visualization rules.
- Must not do: imply cross-Celltype transfer, per-slice scaling, repeated scaling,
  guessed AP order or silent replacement of failed programs.

## Source Inventory
| id | source | role | priority | notes |
|---|---|---|---|---|
| S1 | `R/*.R` | content and structure | must | executable contracts |
| S2 | `config/*.tsv` | parameters | must | public frozen values |
| S3 | `docs/METHODS.md` | terminology | must | concise public Methods |
| S4 | existing workflow | style continuity | should | restrained pastel stage cards |

## Requirement Traceability
| id | requirement | planned visual encoding |
|---|---|---|
| R1 | three evidence dimensions remain independent | three labelled rows in the evidence card |
| R2 | one anchor and one 41-type transfer per query | identity card |
| R3 | exact-Celltype Round 1 and Round 2, k=20, PCs 1-30, source >=100 | two blue/purple cards |
| R4 | direction reversals removed | explicit red audit card before outputs |
| R5 | final scaling rules | readout card with landscape, Q95 spatial and pooled balance |
| R6 | primary-first sensitivity fallback and fallback marker | audit card with S, *, S* |

## Semantic Model
The diagram is a left-to-right data pipeline. Every connector means "passes an
audited table/object to the next stage". Evidence dimensions annotate program
eligibility but do not merge into one composite quality score.

## Style Contract
- Canvas: 1800 x 560 px, white background.
- Font: Arial; title 26 pt; card heading 15 pt; body 11-12 pt.
- Stroke: 2 px, rounded cards, orthogonal grey arrows.
- Palette: teal input, slate identity, amber scoring, blue Round 1, purple Round
  2, red direction audit, green readouts.
- No decorative icons or raster assets.

## Open Assumptions
None. The public workflow records AP order as `UNRESOLVED_DO_NOT_GUESS`.
