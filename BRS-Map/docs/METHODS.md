# BRS-Map: concise Methods

## Framework and inputs

BRS-Map converts contrast-specific transcriptional responses in public
single-cell RNA-sequencing datasets into signed, cell-type-resolved scores and
transfers them through a fixed amygdala snRNA-seq reference to Stereo-seq
CellBins. The released analysis retained seven accessions, 48 objects, 16
prespecified comparisons and 861,176 post-QC external cells. External objects
contained raw RNA counts, condition and independent `experimental_unit_id`
metadata, and labels mapped to a common 41-Celltype vocabulary. Identity
mapping was performed upstream with one Seurat anchor set, one transfer of the
41-type label and one reference-UMAP projection per query. Queries with
reference-feature coverage below 0.70 were rejected.

Three evidence dimensions were retained independently. Sample/statistical
evidence was `formal`, `low_power` or `descriptive`; the gene program was
`formal_signature` or `fallback`; and identity mapping was `primary` or
`sensitivity`. Primary assignments required prediction score >= 0.50.
Sensitivity assignments included scores from 0.30 to < 0.50; scores < 0.30
were unassigned. Final presentation followed a primary-first policy: an
eligible sensitivity result was used only when its corresponding primary
result was unavailable and was marked `S`.

## Experimental-unit pseudobulk and sample evidence

Cells were aggregated within each `experimental_unit_id` and exact Celltype;
cells were never treated as biological replicates. Unit-by-Celltype groups
with fewer than 20 cells were removed before arm-level replication was
assessed. A comparison was `formal` when both arms contained at least three independent
units. A comparison was `low_power` when both arms contained at
least two independent units but either arm contained fewer than three.
Comparisons consisting of one pooled unit per arm were `descriptive` and were
not assigned replicate-based P values or FDR values. A Celltype that failed the
minimum unit requirement after the 20-cell filter was unavailable; it was not
silently downgraded to another evidence tier.

For replicated comparisons, cell-type-specific differential expression used
edgeR. Genes were filtered with `filterByExpr`; library sizes were normalized
by trimmed mean of M values; dispersions were estimated robustly; and robust
quasi-likelihood generalized linear models were tested with `glmQLFTest`. For
descriptive comparisons, expression was summarized as the difference in log2
counts per million with an offset of 0.5.

The three GSE256522 comparisons—STFP versus odour-only exposure, STFP versus
home cage, and odour-only exposure versus home cage—were treated as three
independent formal contrasts. No comparison was used as a gene blacklist for
another comparison.

## Signed response programs and per-cell scores

A `formal_signature` required FDR < 0.05, absolute log2 fold change >= 0.25
and at least 10 genes in each direction. Significant genes were ordered by FDR
and signed fold change, with at most 100 genes retained per direction. If the
formal requirements were not met, a ranked-logFC `fallback` retained 20 genes
with the largest positive finite non-zero fold changes and 20 genes with the
largest-magnitude negative fold changes. After intersection with genes expressed in the target
object, at least 10 genes per direction remained mandatory. Fallback programs
were marked `*` in Celltype-level displays.

UCell was calculated separately for the upregulated and downregulated sets
with `maxRank = 1,500`. The signed score of an eligible external cell was
`UCell_up - UCell_down`. Exact gene lists, sample evidence, identity mode and
program provenance accompanied every score sidecar.

## Two-stage same-type weighted-kNN transfer

Both rounds used reference PCA coordinates, PCs 1-30 and k = 20. Neighbors
were restricted to the exact Celltype; cross-type borrowing was prohibited.
Round 1 transferred external-cell scores to snRNA-seq reference cells after
experimental-unit size balancing and required at least 100 scored external
cells of the matching Celltype. Round 2 transferred snRNA scores to Stereo-seq
CellBins and required at least 100 scored snRNA reference cells of the matching
Celltype. Stereo identities were read from `predicted_classes4`.

For target cell j and source neighbor i, the Gaussian distance weight was
`g_ij = exp[-(d_ij / sigma_j)^2]`, where `sigma_j` was the distance to the kth
neighbor, with a positive-distance fallback for degenerate neighborhoods.
During Round 1, `g_ij` was divided by the number of source cells contributed by
that experimental unit, and weights were normalized to sum to one. The
transferred score was the weighted mean. Round 2 used the same rule with the
Round-1 snRNA source unit retained.

## Direction-consistency and identity selection

For each program and exact Celltype, mean valid transferred scores were
calculated independently after Round 1 and Round 2. A program entered final
Celltype-level and spatial summaries only when both means were finite, non-zero
and had the same sign. A direction reversal, a zero mean, a failed transfer or
an unavailable source produced an explicit NA state rather than a plotted
response. Direction filtering preceded identity selection. Primary results
were preferred; when primary was unavailable but sensitivity passed the same
transfer and direction criteria, sensitivity was used and marked `S`.
Sensitivity plus fallback was therefore marked `S*`.

## Final scaling and visualization

For the neuronal response landscape, individual valid Round-1 snRNA scores
were scaled within each behavior-by-Celltype group. Each score was divided by
the maximum absolute score in that group and clipped to [-1,1]; the scaled
cell scores were then averaged. Scaling therefore preceded averaging and was
not applied to already summarized contrast columns.

For each individual Round-2 spatial program, one denominator was calculated
jointly across all valid target neurons in M1, M2 and M3. The denominator was
the 95th percentile of the absolute raw Round-2 score. Raw scores were divided
by this shared denominator and clipped to [-1,1]. Sections were neither scaled
separately nor averaged. Integrated spatial figures reused these frozen
program-level values without additional scaling.

For the subnuclear response landscape, unique valid neurons from M1, M2 and M3
were pooled within each behavior contrast and subnucleus. The pooled neuronal directional balance
was `(positive_n - negative_n) / (positive_n + negative_n)`; zero scores did
not enter the sign denominator. Point size
represented the proportion of neurons with a valid Round-2 score.

## Audit and scope

BRS-Map records each transferred score, weighted standard deviation,
source-cell count, neighbor and contributing-unit counts, effective neighbor
number (`1 / sum(w^2)`), maximum unit weight, nearest- and kth-neighbor
distances, direction status, identity source, program class and scaling
denominator. Programs failing identity, gene coverage or strict-100 criteria
remain explicit audit records. Spatial maps show CellBins, not spot averages.
The anterior-posterior order of the three Stereo-seq sections was not
independently established and remains `UNRESOLVED_DO_NOT_GUESS`.

BRS-Map begins with harmonized, post-QC, identity-mapped external objects. Raw
download, cell calling, QC, doublet removal and anchor-based identity mapping
are upstream operations. Formal inference is restricted to replicated
comparisons; low-power and descriptive results retain their evidence labels
and are not interpreted as equivalent to formal tests.
