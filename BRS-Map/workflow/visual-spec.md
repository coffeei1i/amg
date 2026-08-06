# Visual Spec

## Source
- Reference image: `BRS-Map_workflow.png` (1,386 x 279 px)
- Target drawio: `BRS-Map_workflow.drawio`
- Canvas: 1,386 x 279 px, white
- Font policy: Arial/Helvetica, black, publication-safe sans serif

## Global Style
- Background: `#FFFFFF`
- Primary font: Arial
- Stroke style: 1 px, muted stage-specific stroke
- Arrow style: straight dark grey `#59636E`, filled classic arrowhead
- Color palette: teal `#E8F4F0/#238B77`; grey `#EEF1F4/#64727E`;
  amber `#FFF4DA/#C18A24`; blue `#EAF0FA/#3567A8`; lilac
  `#F5EBF8/#8E5AA1`; rose `#FBEDEC/#C67974`.

## Regions
| id | bbox x,y,w,h | role | visual notes |
|---|---|---|---|
| label | 4,30,40,220 | framework name | vertical text |
| source | 54,44,190,190 | source data | pale teal rounded box |
| identity | 276,44,190,190 | identity mapping | pale grey rounded box |
| score | 498,44,190,190 | score construction | pale amber rounded box |
| round1 | 720,44,190,190 | round 1 | pale blue rounded box |
| round2 | 942,44,190,190 | round 2 | pale lilac rounded box |
| spatial | 1164,44,190,190 | spatial readout | pale rose rounded box |

## Text Blocks
Stage titles use 15 px bold; central mechanism labels use 12 px; evidence and
parameter captions use 9-10 px. All text is centre-aligned within its stage.

## Shapes
Six rounded stage boxes are repeated on a shared baseline. Small editable
circles, document/card primitives and directional symbols indicate cells,
scores and spatial sections without embedding raster content.

## Connectors
Five left-to-right connectors join consecutive boxes. Each means that the
output of one stage is the input to the next; there are no feedback edges.

## Semantic Relations And Flow
Source objects receive unified identities, yield cell-type-specific signed
program scores, transfer to snRNA cells, transfer again to Stereo-seq CellBins,
and become spatial response maps.

## Icons And Images
All icons in the draw.io file are editable primitives. The supplied PNG is
preserved as a reference/release image only and is not embedded in the draw.io.
