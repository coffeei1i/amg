# BRS-Map: concise Methods

## Framework and inputs

BRS-Map converts contrast-specific transcriptional responses in public
single-cell RNA-sequencing datasets into signed, cell-type-resolved scores and
transfers these scores through a fixed amygdala snRNA-seq reference to
Stereo-seq CellBins. The released analysis retained seven accessions, 48
objects, 16 prespecified comparisons and 861,176 post-QC external cells.
External objects were required to contain raw RNA counts, condition and
independent `experimental_unit_id` metadata, and labels mapped to a common
41-cell-type vocabulary. Identity mapping itself was performed upstream with
one Seurat anchor set, one transfer of the 41-type label and one reference-UMAP
projection per query. Queries with reference-feature coverage below 0.70 were
rejected.

Three evidence dimensions were kept orthogonal throughout the analysis. The
sample/statistical dimension was labelled `formal`, `low_power` or
`descriptive`; the gene-program dimension was labelled `formal_signature` or
`fallback`; and the cell-identity dimension was labelled `primary` or
`sensitivity`. Primary analyses retained identity assignments with prediction
scores >= 0.50. Sensitivity analyses additionally retained assignments with
scores from 0.30 to < 0.50; scores < 0.30 were treated as unassigned.

## Experimental-unit pseudobulk and sample evidence

Cells were aggregated within each `experimental_unit_id` and exact Celltype,
so cells were not treated as biological replicates. Unit-by-cell-type groups
containing fewer than 20 cells were removed before arm-level replication was
assessed. A comparison was `formal` when both arms contained at least three independent units. Comparisons with at least two independent units in each arm
but fewer than three in either arm were `low_power`. Pooled comparisons without
independent replication were `descriptive` and were not assigned replicate-
based P values or FDR values.

For replicated comparisons, differential expression was fitted separately for
each Celltype using edgeR. Genes were filtered with `filterByExpr`; library
sizes were normalized by trimmed mean of M values; dispersions were estimated
robustly; and robust quasi-likelihood generalized linear models were tested
with `glmQLFTest`. For descriptive comparisons, expression was summarized as
the difference in log2 counts per million using an offset of 0.5.

## Signed response programs and per-cell scores

A `formal_signature` required FDR < 0.05, absolute log2 fold change >= 0.25,
and at least 10 genes in each direction. Significant genes were ordered by FDR
and then signed fold change, with at most 100 genes retained per direction. If
these requirements were not met, a `fallback` program used the 20 genes with
the largest positive finite nonzero fold changes and the 20 genes with the
largest-magnitude negative fold changes. After intersecting each direction
with genes expressed in the target object, at least 10 genes per direction
were still required; otherwise that program-object combination was not scored.

UCell was calculated separately for the upregulated and downregulated gene
sets with `maxRank = 1,500`. The signed response score for each eligible
external cell was defined as `UCell_up - UCell_down`. Program provenance,
sample evidence, identity mode and the exact gene lists were retained with
every score sidecar.

## Two-stage same-type weighted-kNN transfer

Both transfer rounds used reference PCA coordinates, PCs 1-30 and k = 20.
Neighbors were restricted to the exact Celltype; no cross-type borrowing was
allowed. Round 1 transferred external-cell scores to snRNA-seq reference cells after experimental-unit size balancing
and was performed only when at least 100 scored external cells of the matching
type were available. Round 2 transferred the resulting snRNA scores to
Stereo-seq CellBins and required at least 100 scored snRNA reference cells of
the matching type. Stereo-seq identities were read from
`predicted_classes4`.

For a target cell, the Gaussian distance weight of source neighbor *i* was
`g_i = exp[-(d_i/sigma)^2]`, where `sigma` was the distance to the kth neighbor
with a positive-distance fallback for degenerate neighborhoods. To prevent a
large experimental unit from dominating Round 1, `g_i` was divided by the
number of source cells from that experimental unit, and all weights were then
renormalized to sum to one. The transferred score was the weighted mean. The
same weighting rule was used in Round 2, with the snRNA source unit retained
from Round 1.

## Audit and spatial readout

For every target, BRS-Map records the transferred score, weighted standard
deviation, source-cell count, numbers of neighbors and contributing units,
effective neighbor number (`1/sum(w^2)`), maximum unit weight, and nearest and
kth-neighbor distances. Programs failing the identity, gene-coverage or
strict-100 source requirement remain explicit audit records rather than being
silently imputed. Spatial maps display CellBins; they do not represent spot
averages. The anterior-posterior ordering of the three Stereo-seq sections was
not independently established and is therefore stored as
`UNRESOLVED_DO_NOT_GUESS`.

## Scope

BRS-Map begins with harmonized, post-QC, identity-mapped external objects. Raw
download, cell calling, QC, doublet removal and anchor-based cell-identity
mapping are upstream operations and are not recomputed by this module. Formal
inference is restricted to replicated comparisons; low-power and descriptive
results are reported with their evidence labels and are not interpreted as
equivalent to formal tests.
