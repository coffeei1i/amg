# BRS-Map public release design

## Objective

Add a self-contained, reproducible `BRS-Map/` module to `coffeei1i/amg` that documents and implements behavior-response score construction and two-stage weighted-kNN transfer from harmonized external single-cell RNA-seq data to the study snRNA-seq reference and then to Stereo-seq CellBin data.

The public release must contain a concise English Methods section, all frozen parameters, executable R code, contract tests, environment metadata, and the supplied BRS-Map workflow figure in both publication-ready and editable formats.

## Scope

### Included

1. Input validation for harmonized external single-cell objects with sample metadata and 41-type identity mapping.
2. Three orthogonal evidence dimensions:
   - sample/statistical evidence;
   - gene-program evidence;
   - cell-identity mapping evidence.
3. Experimental-unit pseudobulk aggregation and contrast eligibility.
4. edgeR quasi-likelihood differential expression or descriptive pooled log-fold change.
5. Formal and fallback signed gene-program construction.
6. Per-cell signed UCell scoring.
7. Round 1 same-cell-type weighted-kNN transfer from external cells to the snRNA-seq reference.
8. Round 2 same-cell-type weighted-kNN transfer from snRNA-seq cells to Stereo-seq CellBins.
9. Audit tables, reproducibility manifests, and example plotting entry points.
10. English Methods and workflow assets.

### Excluded

- raw-data download;
- cell calling;
- upstream single-cell QC and doublet removal;
- re-computation of the 41-type Seurat anchor mapping;
- raw or processed biological data;
- large result matrices and project-specific absolute paths.

## Public terminology

Internal legacy values may be accepted as backward-compatible inputs, but all public documentation and output labels use:

- `formal` for replicate-supported inference;
- `low_power` for eligible but under-replicated contrasts;
- `descriptive` for pooled, non-replicated contrasts formerly labelled `pool_only` or `EXPLORATORY_POOL`;
- `formal_signature` and `fallback` for gene-program evidence;
- `primary` and `sensitivity` for cell-identity evidence.

## Frozen analytical contract

### Input cohort

- Seven retained accessions: E-MTAB-12096, GSE103976, GSE231790, GSE246147, GSE254360, GSE256522 and GSE310656.
- Forty-eight harmonized external objects.
- Sixteen pre-specified behavioral comparisons.
- GSE207128 remains excluded.
- The current external-cell total is 861,176 post-QC cells.

### Dimension 1: sample/statistical evidence

- Cells are aggregated by `experimental_unit_id × Celltype`; cells are never treated as biological replicates.
- A pseudobulk unit-cell-type group requires at least 20 cells.
- `formal`: both contrast arms contain at least three independent experimental units.
- `low_power`: both arms contain at least two independent experimental units, but at least one arm contains fewer than three.
- `descriptive`: pooled, non-replicated comparison. Log2 CPM values use an offset of 0.5 and no replicate-based P value or FDR is reported.
- Comparisons below their declared minimum independent-unit requirement are not released for that cell type.
- Replicated contrasts use edgeR: `filterByExpr`, TMM normalization, robust dispersion estimation, robust quasi-likelihood fitting, and a two-group case-versus-control test.

### Dimension 2: gene-program evidence

- Formal genes satisfy `FDR < 0.05` and `abs(logFC) >= 0.25`.
- Each released signed program requires at least 10 upregulated and 10 downregulated genes.
- Formal programs retain at most 100 genes per direction, ranked by FDR and then signed log-fold change.
- If the formal criterion is not met, the fallback ranks all non-zero finite log-fold changes and retains up to the top 20 positive and top 20 negative genes.
- Fallback programs still require at least 10 genes per direction.
- After intersection with a target expression matrix, at least 10 genes per direction must remain.
- UCell is calculated with `maxRank = 1500`.
- The signed response score is `UCell_up - UCell_down`.

### Dimension 3: cell-identity evidence

- All labels belong to the frozen 41-cell-type vocabulary.
- `primary` uses only `mapping_status == accepted`, corresponding to prediction score `>= 0.50`.
- `sensitivity` uses accepted and sensitivity mappings; the sensitivity band is `0.30 <= score < 0.50`.
- Scores below 0.30 are unassigned and excluded.
- Minimum query feature coverage is 0.70.
- Identity evidence is orthogonal to sample/statistical and gene-program evidence; the three dimensions are recorded independently.

### Signed two-stage weighted-kNN

- Both rounds use `k = 20` and frozen reference PCA dimensions 1–30.
- Neighbors must have the exact same unified cell type; no cross-cell-type borrowing is allowed.
- Round 1 requires at least 100 eligible external source cells for the program-cell-type combination.
- Round 2 requires at least 100 scored snRNA-seq reference cells for the program-cell-type combination.
- For neighbor distance `d_i`, the Gaussian component is `exp[-(d_i/sigma)^2]`, where `sigma` is the kth-neighbor distance (with a deterministic positive-distance fallback).
- Each neighbor weight is divided by the size of its source experimental unit and then renormalized to sum to one.
- Outputs include the transferred score, weighted standard deviation, neighbor count, source-unit count, effective sample size, maximum source-unit weight, nearest-neighbor distance and kth-neighbor distance.
- Stereo-seq identity is read from `predicted_classes4`.
- Anterior-posterior slice order remains `UNRESOLVED_DO_NOT_GUESS`.

## Repository layout

```text
BRS-Map/
├── README.md
├── R/
│   ├── contracts.R
│   ├── pseudobulk.R
│   ├── programs.R
│   ├── signed_ucell.R
│   ├── weighted_knn.R
│   └── audit.R
├── scripts/
│   ├── run_brs_map.R
│   └── preflight.R
├── config/
│   ├── parameters.tsv
│   ├── contrasts.tsv
│   └── datasets.tsv
├── tests/
├── docs/
│   └── METHODS.md
├── workflow/
│   ├── BRS-Map_workflow.png
│   ├── BRS-Map_workflow.svg
│   └── BRS-Map_workflow.drawio
└── environment/
    └── session-info.txt
```

## Workflow figure contract

The user-supplied six-stage figure is the content and layout reference:

1. Source Data;
2. Identity Mapping;
3. Score Construction;
4. Round 1 Transfer;
5. Round 2 Transfer;
6. Spatial Readout.

The supplied PNG is retained. An editable draw.io recreation and SVG export are added. Text and parameter labels must match the released code and English Methods.

## Methods writing contract

The English Methods section is concise but complete. It must:

- define BRS-Map before using the acronym;
- state the input boundary;
- describe the three evidence dimensions as orthogonal annotations;
- report every threshold and frozen parameter listed above;
- distinguish inferential (`formal`) from descriptive evidence;
- state that sensitivity identity mappings do not replace primary results;
- state that cell-type matching is exact in both transfer rounds;
- avoid unsupported biological claims and avoid inventing anterior-posterior order.

## Verification and publication

Before push:

1. parse all R files;
2. run contract tests for the three evidence dimensions, signed UCell and both strict-100 transfer gates;
3. validate parameter tables against source constants;
4. validate draw.io XML and visually inspect the rendered figure through at least three repair cycles;
5. confirm the repository excludes data, generated matrices, credentials and absolute local paths;
6. review `git diff`, create one focused feature commit series, and push to the existing `master` branch only after all checks pass.

## Risks and mitigations

- **Terminology drift:** public labels are centralized in one parameter/contract layer; legacy aliases are documented.
- **Threshold drift:** tests bind documentation/configuration to the frozen constants.
- **Data leakage:** `.gitignore` and pre-push checks reject RDS, H5AD, large tables and absolute server paths.
- **Identity overstatement:** primary and sensitivity mappings remain separate, and unassigned cells are excluded.
- **Pseudoreplication:** aggregation and evidence-tier tests enforce independent experimental units.
- **Figure/code mismatch:** workflow labels are checked against the parameter table and Methods text.
