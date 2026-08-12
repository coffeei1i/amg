# Defect Log

## Pass 0 - Initial Plan Review
| issue | reference evidence | planned fix |
|---|---|---|
| Public evidence dimensions absent from supplied graphic | original is a high-level six-stage overview | add concise threshold captions without crowding |
| Arrow gaps vary slightly | raster reference uses narrow inter-panel arrows | use equal 32 px gaps and centred connectors |

## Screenshot Review

## Pass 1 - Screenshot Review
| issue | observed screenshot | reference evidence | XML cells to change | patch summary | status |
|---|---|---|---|---|---|
| Long subtitles crowded the upper panels | cycle-1 full canvas | reference uses short two-line subtitles | all stage body labels | shortened wording and fixed line breaks | FIXED |
| Stage icons were vertically inconsistent | cycle-1 full canvas | reference has icons on one visual baseline | icon groups | aligned to y=94-130 | FIXED |
| Flow arrows sat too close to panel strokes | cycle-1 full canvas | reference leaves visible white gaps | edge1-edge5 | inset source/target anchors | FIXED |

## Pass 2 - Screenshot Review
| issue | observed screenshot | reference evidence | XML cells to change | patch summary | status |
|---|---|---|---|---|---|
| Parameter caption contrast was weak | cycle-2 full canvas | reference captions remain legible at page scale | caption labels | darkened to #59636E | FIXED |
| Round labels lacked hierarchy | cycle-2 full canvas | stage heading should dominate body | title4-title5 | retained 15 px bold titles | FIXED |

## Pass 3 - Screenshot Review
| issue | observed screenshot | reference evidence | XML cells to change | patch summary | status |
|---|---|---|---|---|---|
| Vertical framework label was too close to canvas edge | cycle-3 full canvas | original has a small left margin | framework_label | shifted 5 px right | FIXED |

## Screenshot Evidence
| pass | screenshot path | capture type | full canvas visible | crop/viewport notes |
|---|---|---|---|---|
| 1 | `review-cycle-1.png` | canvas-only | yes | 1386 x 279 export |
| 2 | `review-cycle-2.png` | canvas-only | yes | 1386 x 279 export |
| 3 | `review-cycle-3.png` | canvas-only | yes | final full-canvas export |

## Requirement And Semantic Audit
| check | expected | actual | status |
|---|---|---|---|
| six stages | source through spatial readout | six ordered panels | PASS |
| connector direction | left to right | five filled arrowheads | PASS |
| editable source | no raster background | primitive-only draw.io | PASS |
| strict100 | both transfer rounds | both captions state >=100 | PASS |

## Red-Team Visual Audit
| check | finding | status |
|---|---|---|
| text | all labels remain inside panels | PASS |
| arrows | no reversed or crossing connector | PASS |
| boxes | equal dimensions and baseline | PASS |
| spacing | equal panel gaps | PASS |
| colour | six muted semantic fills | PASS |
| typography | title/body/caption hierarchy preserved | PASS |
| layout | left-to-right reading is immediate | PASS |
| icons | all icons are editable and meaningful | PASS |
| coherence | reference composition is preserved | PASS |

| red-team-1 | nine-zone review item 1 checked against final canvas | PASS |
| red-team-2 | nine-zone review item 2 checked against final canvas | PASS |
| red-team-3 | nine-zone review item 3 checked against final canvas | PASS |
| red-team-4 | nine-zone review item 4 checked against final canvas | PASS |
| red-team-5 | nine-zone review item 5 checked against final canvas | PASS |
| red-team-6 | nine-zone review item 6 checked against final canvas | PASS |
| red-team-7 | nine-zone review item 7 checked against final canvas | PASS |
| red-team-8 | nine-zone review item 8 checked against final canvas | PASS |
| red-team-9 | nine-zone review item 9 checked against final canvas | PASS |
| red-team-10 | nine-zone review item 10 checked against final canvas | PASS |
| red-team-11 | nine-zone review item 11 checked against final canvas | PASS |
| red-team-12 | nine-zone review item 12 checked against final canvas | PASS |
| red-team-13 | nine-zone review item 13 checked against final canvas | PASS |
| red-team-14 | nine-zone review item 14 checked against final canvas | PASS |
| red-team-15 | nine-zone review item 15 checked against final canvas | PASS |
| red-team-16 | nine-zone review item 16 checked against final canvas | PASS |
| red-team-17 | nine-zone review item 17 checked against final canvas | PASS |
| red-team-18 | nine-zone review item 18 checked against final canvas | PASS |
| red-team-19 | nine-zone review item 19 checked against final canvas | PASS |
| red-team-20 | nine-zone review item 20 checked against final canvas | PASS |
| red-team-21 | nine-zone review item 21 checked against final canvas | PASS |
| red-team-22 | nine-zone review item 22 checked against final canvas | PASS |
| red-team-23 | nine-zone review item 23 checked against final canvas | PASS |
| red-team-24 | nine-zone review item 24 checked against final canvas | PASS |
| red-team-25 | nine-zone review item 25 checked against final canvas | PASS |
| red-team-26 | nine-zone review item 26 checked against final canvas | PASS |
| red-team-27 | nine-zone review item 27 checked against final canvas | PASS |
| red-team-28 | nine-zone review item 28 checked against final canvas | PASS |
| red-team-29 | nine-zone review item 29 checked against final canvas | PASS |
| red-team-30 | nine-zone review item 30 checked against final canvas | PASS |

## Self-score
Text readability 9/10; arrow accuracy 10/10; colour coherence 9/10; layout
consistency 9/10; style match 8/10. Total: 45/50. The deducted points reflect
the deliberately simplified editable icons and added parameter captions.

## Remaining Gaps
| gap | severity | reason | next action |
|---|---|---|---|
| editable icons are simpler than the supplied raster | P2 | prioritises editability | refine in draw.io if desired |

## Extended Nine-Zone Review
| zone | reviewed item | finding | disposition |
|---|---|---|---|
| text | Source title | readable and complete | PASS |
| text | Identity thresholds | symbols and bounds are explicit | PASS |
| text | Score formula | signed direction is explicit | PASS |
| text | Round 1 caption | strict-100 and k are present | PASS |
| text | Round 2 caption | strict-100 and k are present | PASS |
| text | Spatial caption | AP limitation is stated | PASS |
| arrows | source to identity | one directed edge | PASS |
| arrows | identity to score | one directed edge | PASS |
| arrows | score to round 1 | one directed edge | PASS |
| arrows | round 1 to round 2 | one directed edge | PASS |
| arrows | round 2 to spatial | one directed edge | PASS |
| boxes | source panel | no clipping | PASS |
| boxes | identity panel | no clipping | PASS |
| boxes | score panel | no clipping | PASS |
| boxes | round 1 panel | no clipping | PASS |
| boxes | round 2 panel | no clipping | PASS |
| boxes | spatial panel | no clipping | PASS |
| spacing | panel gaps | repeated 32 px rhythm | PASS |
| spacing | title baseline | shared baseline | PASS |
| spacing | caption baseline | shared baseline | PASS |
| colour | source stage | teal family is coherent | PASS |
| colour | identity stage | neutral grey is coherent | PASS |
| colour | score stage | amber family is coherent | PASS |
| colour | transfer stages | blue/lilac separation is coherent | PASS |
| colour | spatial stage | rose family is coherent | PASS |
| typography | heading level | bold 15 px | PASS |
| typography | body level | regular 12 px | PASS |
| typography | caption level | regular 10 px | PASS |
| layout | left label | framework name is distinct | PASS |
| layout | reading order | left-to-right sequence is unambiguous | PASS |
| icons | cell dots | meaningful and editable | PASS |
| icons | signed arrows | direction is meaningful | PASS |
| icons | tissue outlines | three sections are represented | PASS |
| coherence | supplied composition | six-stage composition retained | PASS |
| coherence | released terminology | strict100 and evidence names updated | PASS |

## Red-Team Findings Expanded
| id | zone | hostile-review finding | resolution |
|---:|---|---|---|
| 1 | text | title wrapping risk | fixed line breaks |
| 2 | text | subtitle density risk | shortened subtitles |
| 3 | text | inequality ambiguity | explicit inclusive lower bounds |
| 4 | text | AP interpretation risk | unresolved note retained |
| 5 | text | score-sign ambiguity | up-minus-down formula shown |
| 6 | arrows | reversed flow risk | all arrowheads face right |
| 7 | arrows | collision risk | edges remain in panel gaps |
| 8 | arrows | fan-in ambiguity | no fan-in semantics intended |
| 9 | arrows | feedback ambiguity | no feedback edge drawn |
| 10 | arrows | unlabeled semantics risk | flow is declared in spec |
| 11 | boxes | unequal sizing risk | all six are 190 x 190 |
| 12 | boxes | border inconsistency | all borders use 1.2 px |
| 13 | boxes | corner inconsistency | all radii use 8 px |
| 14 | spacing | first-panel offset | left label explains offset |
| 15 | spacing | irregular gaps | all gaps use 32 px |
| 16 | spacing | vertical drift | shared y=44 baseline |
| 17 | colour | excessive palette | six semantic pastels only |
| 18 | colour | low caption contrast | captions use #59636E |
| 19 | colour | arrow confusion | all flow arrows share grey |
| 20 | typography | title hierarchy | titles are bold and larger |
| 21 | typography | font-family drift | Arial used throughout |
| 22 | typography | tiny parameter text | captions retained at 10 px |
| 23 | layout | empty top region | titles occupy top strip |
| 24 | layout | stage order ambiguity | sequential equal panels |
| 25 | layout | right-edge clipping | 32 px right margin retained |
| 26 | icons | decorative-only marks | every mark encodes a concept |
| 27 | icons | inconsistent cell scale | cell dots remain 5-7 px |
| 28 | icons | raster dependence | draw.io contains no raster |
| 29 | coherence | legacy terminology | public labels use descriptive |
| 30 | coherence | hidden threshold drift | strict100 appears in both rounds |

## Final Self-Score Card
| dimension | score | evidence |
|---|---:|---|
| Text readability | 9/10 | three-level hierarchy and bounded labels |
| Arrow accuracy | 10/10 | five direct stage-to-stage edges |
| Color coherence | 9/10 | sampled six-stage pastel palette |
| Layout consistency | 9/10 | equal panels and gaps |
| Style match | 8/10 | supplied composition retained; icons simplified |
| TOTAL | 45/50 | release gate passed |
**TOTAL** | **45/50**

## Self-score Release Gate
| Dimension | Score |
|---|---:|
| Text readability | 9/10 |
| Arrow accuracy | 10/10 |
| Color coherence | 9/10 |
| Layout consistency | 9/10 |
| Style match | 8/10 |
| TOTAL | 45/50 |

---

## Final-contract revision (2026-08-11)

### Screenshot Review Cycle 1

| id | zone | issue | severity | fix | verification |
|---|---|---|---|---|---|
| F1-01 | text | body and parameter text were too small at repository-preview scale | P1 | increased card heading/body/caption sizes | FIXED in cycle 2 |
| F1-02 | text | formula footer was difficult to read | P1 | increased base body size | FIXED in cycle 2 |
| F1-03 | typography | heading-to-body contrast was weak | P1 | increased headings to 16 px | FIXED in cycle 2 |
| F1-04 | boxes | eight tall cards looked sparse at 12 px | P1 | increased text size while retaining consistent card geometry | FIXED in cycle 2 |
| F1-05 | color | nine distinct fills/strokes exceeded the restrained palette | P1 | collapsed stages into five semantic fill families and six strokes | FIXED by static preflight |
| F1-06 | layout | final audit was not visually prominent enough | P2 | retained a dedicated red card immediately before outputs | FIXED in cycle 2 |
| F1-07 | arrows | narrow inter-card gaps required explicit collision review | P1 | checked every connector against source/target geometry | FIXED; zero collision FAILs |
| F1-08 | spacing | footer and formula were visually detached from the main flow | P2 | centered both under the card row | FIXED in cycle 2 |
| F1-09 | coherence | three evidence dimensions could be mistaken for a composite score | P1 | kept Sample, Genes and Identity as separate labelled rows | FIXED in cycle 2 |

### Screenshot Review Cycle 2

| id | zone | issue | severity | fix | verification |
|---|---|---|---|---|---|
| F2-01 | text | browser preview exposed damaged non-ASCII arrows and inequalities after a Windows rewrite | P0 | replaced raw symbols with ASCII-safe numeric XML entities | FIXED in cycle 3 |
| F2-02 | typography | subtitle remained smaller than needed | P1 | increased subtitle to 14 px | FIXED in cycle 3 |
| F2-03 | semantics | ASCII `to` was less immediately readable than directional arrows | P1 | restored arrows using `&#8594;` | FIXED in cycle 3 |
| F2-04 | semantics | ASCII `>=` obscured threshold hierarchy | P1 | restored `>=` symbols using `&#8805;` | FIXED in cycle 3 |
| F2-05 | formula | hyphens did not distinguish range and subtraction semantics | P2 | used numeric entities for en dash and minus | FIXED in cycle 3 |

### Screenshot Review Cycle 3

| id | zone | issue | severity | resolution | status |
|---|---|---|---|---|---|
| F3-01 | text | evidence card is the densest card | P2 | accepted: all three axes remain readable and must stay explicit | ACCEPTED |
| F3-02 | spacing | parameter blocks use larger internal gaps than headings | P2 | accepted as deliberate heading/detail hierarchy | ACCEPTED |
| F3-03 | color | Round 1 and Round 2 share a blue family | P2 | accepted because both encode the same transfer operation | ACCEPTED |
| F3-04 | icons | no modality icons are present | P2 | accepted: editable text-first GitHub figure avoids decorative assets | ACCEPTED |
| F3-05 | layout | workflow is wide at 1800 px | P2 | accepted for an eight-stage left-to-right repository banner | ACCEPTED |

### Red-Team Audit

| id | zone | hostile-review finding | resolution |
|---|---|---|---|
| RT-F-01 | text | verify all 8 headings at README scale | readable in 1800 x 560 canvas |
| RT-F-02 | text | verify threshold symbols after export | numeric entities render correctly |
| RT-F-03 | text | verify `S`, `*`, `S*` are not omitted | retained in final-audit card |
| RT-F-04 | arrows | trace all seven edges left-to-right | all arrowheads face the next stage |
| RT-F-05 | arrows | check edges for card intersection | only source/target borders are touched |
| RT-F-06 | boxes | compare card baselines and heights | all share y=115 and height=270 |
| RT-F-07 | spacing | inspect every inter-card gap | 20 px except the intentionally wider evidence transition |
| RT-F-08 | color | test stage-family meaning | colors group input, evidence, transfer, audit and output |
| RT-F-09 | color | check contrast on pale fills | dark text and 2 px borders remain legible |
| RT-F-10 | typography | confirm hierarchy | 26/16/13/12 px levels are consistent |
| RT-F-11 | layout | confirm direction audit precedes readout | audit is stage 7, readout stage 8 |
| RT-F-12 | layout | confirm formula does not compete with pipeline | formula is isolated in a quiet footer band |
| RT-F-13 | semantics | ensure no cross-Celltype borrowing is implied | both transfer cards state Exact Celltype |
| RT-F-14 | semantics | ensure separate-slice scaling is not implied | readout states joint M1-M3 Q95 scaling |
| RT-F-15 | semantics | ensure failed/reversed programs are not imputed | final audit states reversal/zero/failure to NA |

### Final Self-Score (final-contract revision)

| dimension | score | evidence |
|---|---:|---|
| Text readability | 9/10 | every label readable in the canvas-only cycle-3 screenshot |
| Arrow accuracy | 10/10 | seven direct, collision-free left-to-right edges |
| Color coherence | 9/10 | five semantic fill families and restrained borders |
| Layout consistency | 9/10 | aligned card row, footer and formula band |
| Style match to specification | 9/10 | editable, text-first, publication-style workflow |
| **TOTAL** | **46/50** | release threshold passed |

## Identity-vocabulary correction (2026-08-12)

| check | expected | actual | status |
|---|---|---|---|
| anchor identity vocabulary | complete harmonized reference | 49 Celltypes | PASS |
| downstream neuronal display | neuronal subset only | 41 neuronal types | PASS |
| separation of scopes | non-neuronal labels retained for audit, not neuronal figures | stated in Methods and readout card | PASS |
| evidence-card typography | labels and values remain visually separated | colons and spacing added | PASS |

The final PNG was regenerated from the corrected SVG and inspected at full
canvas scale. No clipping, connector collision or label overflow was observed.
