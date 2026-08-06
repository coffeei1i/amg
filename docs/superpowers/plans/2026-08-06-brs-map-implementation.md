# BRS-Map Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Publish a self-contained, tested BRS-Map module in `coffeei1i/amg` with exact evidence thresholds, signed UCell scoring, two-stage strict-100 weighted-kNN transfer, concise English Methods, and editable workflow assets.

**Architecture:** BRS-Map is a configuration-driven R research module that starts from harmonized, identity-mapped external single-cell objects. Focused R files implement contracts, pseudobulk evidence classification, signed program construction, UCell scoring, weighted-kNN transfer, and audit output; one orchestrator sources the modules and validates inputs before running stages. Documentation and workflow assets are generated from the same frozen parameter contract.

**Tech Stack:** R 4.2+, Matrix, data.table, edgeR, UCell, RANN, Seurat/SeuratObject, jsonlite, digest; base-R contract tests; draw.io XML and SVG/PNG workflow assets; Git/GitHub.

---

## File map

- `BRS-Map/R/contracts.R`: frozen parameters, public terminology, schema and input validators.
- `BRS-Map/R/pseudobulk.R`: mapping-mode selection, unit-cell-type aggregation and statistical evidence classification.
- `BRS-Map/R/programs.R`: edgeR/descriptive DE and formal/fallback signed program selection.
- `BRS-Map/R/signed_ucell.R`: per-cell signed UCell scoring and score-sidecar validation.
- `BRS-Map/R/weighted_knn.R`: reusable unit-balanced weighted-kNN and both transfer stages.
- `BRS-Map/R/audit.R`: cross-stage audit and release manifest creation.
- `BRS-Map/scripts/preflight.R`: dependency and configuration validation without producing formal output.
- `BRS-Map/scripts/run_brs_map.R`: explicit stage orchestrator with configurable paths.
- `BRS-Map/config/parameters.tsv`: machine-readable frozen parameter table.
- `BRS-Map/config/datasets.tsv`: seven retained datasets and exclusion record.
- `BRS-Map/config/contrasts.tsv`: 16 pre-specified comparisons and sample-evidence rules.
- `BRS-Map/tests/*.R`: base-R contract tests for every public threshold and formula.
- `BRS-Map/docs/METHODS.md`: concise English manuscript Methods.
- `BRS-Map/workflow/*`: supplied PNG plus editable draw.io and SVG.
- `BRS-Map/environment/session-info.txt`: tested package/runtime declaration.
- `BRS-Map/README.md`: input boundary, quick start, outputs and audit interpretation.
- `.gitignore`: reject data objects, generated outputs and local brainstorm state.
- `README.md`: link the new BRS-Map module from the repository root.

### Task 1: Isolate the implementation and establish repository hygiene

**Files:**
- Modify: `.gitignore`
- Modify: `README.md`
- Create: `BRS-Map/README.md`

- [ ] **Step 1: Create a feature branch and verify a clean tree**

Run:

```bash
git status --short
git switch -c feat/brs-map
```

Expected: no pre-existing user changes; branch `feat/brs-map` is active.

- [ ] **Step 2: Add a failing repository-hygiene test**

Create `BRS-Map/tests/test_repository_hygiene.R` that recursively rejects `.rds`, `.h5ad`, generated `output/` trees, credentials, and `/data*` or `Z:/` absolute paths in public R/config/documentation files.

```r
root <- normalizePath(file.path(getwd(), "BRS-Map"), mustWork = TRUE)
files <- list.files(root, recursive = TRUE, full.names = TRUE,
                    all.files = TRUE, no.. = TRUE)
stopifnot(!any(grepl("\\.(rds|h5ad|h5|loom)$", files, ignore.case = TRUE)))
text_files <- files[grepl("\\.(R|md|tsv|txt|yml|yaml|json)$", files)]
contents <- unlist(lapply(text_files, readLines, warn = FALSE), use.names = FALSE)
stopifnot(!any(grepl("(^|[^A-Za-z])/data[0-9]*/|Z:/data", contents)))
cat("REPOSITORY HYGIENE TEST PASS\n")
```

- [ ] **Step 3: Run the test and verify it fails because the module is absent**

Run: `Rscript BRS-Map/tests/test_repository_hygiene.R`

Expected: FAIL because `BRS-Map` has not yet been created.

- [ ] **Step 4: Add module-level ignore rules and README boundary**

Add patterns for `.Rhistory`, `.RData`, `.Rproj.user/`, `.superpowers/`, `output/`, `output.staging/`, `*.rds`, `*.h5ad`, `*.h5`, `*.loom`, and credential files. Create the BRS-Map README skeleton that defines accepted input and excluded upstream work.

- [ ] **Step 5: Re-run the hygiene test**

Expected: `REPOSITORY HYGIENE TEST PASS`.

- [ ] **Step 6: Commit**

```bash
git add .gitignore README.md BRS-Map/README.md BRS-Map/tests/test_repository_hygiene.R
git commit -m "chore: establish BRS-Map module boundary"
```

### Task 2: Freeze the public parameter and terminology contract

**Files:**
- Create: `BRS-Map/R/contracts.R`
- Create: `BRS-Map/config/parameters.tsv`
- Create: `BRS-Map/config/datasets.tsv`
- Create: `BRS-Map/config/contrasts.tsv`
- Create: `BRS-Map/tests/test_contracts.R`

- [ ] **Step 1: Write failing assertions for all frozen constants**

The test sources `contracts.R` and checks:

```r
p <- brs_parameters()
stopifnot(
  p$minimum_cells_per_unit_celltype == 20L,
  p$formal_minimum_units_per_arm == 3L,
  p$low_power_minimum_units_per_arm == 2L,
  p$formal_fdr < 0.05,
  p$formal_abs_logfc == 0.25,
  p$minimum_genes_per_direction == 10L,
  p$formal_maximum_genes_per_direction == 100L,
  p$fallback_maximum_genes_per_direction == 20L,
  p$ucell_max_rank == 1500L,
  p$identity_primary_score == 0.50,
  p$identity_sensitivity_score == 0.30,
  p$minimum_query_feature_coverage == 0.70,
  p$knn_k == 20L,
  identical(p$pca_dimensions, 1:30),
  p$round1_minimum_cells == 100L,
  p$round2_minimum_cells == 100L
)
stopifnot(identical(brs_public_evidence("EXPLORATORY_POOL"), "descriptive"))
```

- [ ] **Step 2: Run and verify failure**

Run: `Rscript BRS-Map/tests/test_contracts.R`

Expected: FAIL because `brs_parameters()` does not exist.

- [ ] **Step 3: Implement the contract**

`brs_parameters()` returns the constants above. `brs_public_evidence()` maps `FORMAL` to `formal`, `AUXILIARY_LOW_POWER` to `low_power`, and `EXPLORATORY_POOL`/`pool_only` to `descriptive`, rejecting unknown values. Add validators for the 41-type vocabulary, identity scores in `[0,1]`, unique cell IDs, required metadata fields, and contrast-table schemas.

- [ ] **Step 4: Add the three TSV configurations**

`parameters.tsv` contains one row per constant. `datasets.tsv` contains the seven retained datasets plus GSE207128 as `excluded`. `contrasts.tsv` contains the 16 frozen comparisons, case/control labels, role, minimum units, formal minimum units and public evidence label.

- [ ] **Step 5: Run tests and commit**

```bash
Rscript BRS-Map/tests/test_contracts.R
git add BRS-Map/R/contracts.R BRS-Map/config BRS-Map/tests/test_contracts.R
git commit -m "feat: freeze BRS-Map evidence contract"
```

### Task 3: Implement experimental-unit pseudobulk and statistical evidence

**Files:**
- Create: `BRS-Map/R/pseudobulk.R`
- Create: `BRS-Map/tests/test_pseudobulk.R`

- [ ] **Step 1: Write failing tests**

Cover accepted-only primary selection, accepted-plus-sensitivity selection, exclusion of unassigned cells, exact count conservation, the 20-cell unit-cell-type threshold, formal 3-vs-3, low-power 2-vs-2, and descriptive pooled evidence.

```r
formal <- classify_sample_evidence(3L, 4L, minimum = 2L,
                                   formal_minimum = 3L,
                                   pooled = FALSE)
low <- classify_sample_evidence(2L, 3L, minimum = 2L,
                                formal_minimum = 3L,
                                pooled = FALSE)
pool <- classify_sample_evidence(1L, 1L, minimum = 1L,
                                 formal_minimum = 3L,
                                 pooled = TRUE)
stopifnot(formal == "formal", low == "low_power", pool == "descriptive")
```

- [ ] **Step 2: Verify failure**

Run: `Rscript BRS-Map/tests/test_pseudobulk.R`

- [ ] **Step 3: Implement minimal functions**

Implement `select_cells_by_identity_mode()`, `aggregate_pseudobulk_counts()`, `classify_sample_evidence()` and `assess_contrast_eligibility()`. Aggregation key is `experimental_unit_id × Celltype`; groups with fewer than 20 cells are excluded before arm counts are assessed.

- [ ] **Step 4: Verify count conservation and evidence labels**

Expected: `PSEUDOBULK TEST PASS`.

- [ ] **Step 5: Commit**

```bash
git add BRS-Map/R/pseudobulk.R BRS-Map/tests/test_pseudobulk.R
git commit -m "feat: add unit-aware pseudobulk evidence"
```

### Task 4: Implement DE and signed program construction

**Files:**
- Create: `BRS-Map/R/programs.R`
- Create: `BRS-Map/tests/test_programs.R`

- [ ] **Step 1: Write failing tests for formal, fallback and rejection**

```r
formal <- build_signed_program(formal_de, list(program_id = "formal",
                               sample_evidence = "formal"))
stopifnot(formal$signature_method == "formal_signature",
          length(formal$up_genes) == 100L,
          length(formal$down_genes) == 100L)
fallback <- build_signed_program(weak_de, list(program_id = "fallback",
                                 sample_evidence = "low_power"))
stopifnot(fallback$signature_method == "fallback",
          length(fallback$up_genes) == 20L,
          length(fallback$down_genes) == 20L)
```

- [ ] **Step 2: Verify failure**

Run: `Rscript BRS-Map/tests/test_programs.R`

- [ ] **Step 3: Implement DE and ranking**

Implement `run_edger_contrast()` with `filterByExpr`, `calcNormFactors`, robust `estimateDisp`, robust `glmQLFit`, and `glmQLFTest(coef=2)`. Implement descriptive pooled log2 CPM difference with offset 0.5 and `NA` P/FDR. Formal ranking uses `FDR < 0.05`, `abs(logFC) >= 0.25`, FDR then signed-logFC ordering, capped at 100/direction. Fallback uses non-zero finite logFC and top 20/direction. Both require at least 10/direction.

- [ ] **Step 4: Run tests and commit**

```bash
Rscript BRS-Map/tests/test_programs.R
git add BRS-Map/R/programs.R BRS-Map/tests/test_programs.R
git commit -m "feat: construct signed response programs"
```

### Task 5: Implement signed UCell scoring

**Files:**
- Create: `BRS-Map/R/signed_ucell.R`
- Create: `BRS-Map/tests/test_signed_ucell.R`

- [ ] **Step 1: Write failing scorer-injection tests**

Use a deterministic fake scorer to verify identity-mode cell selection, ≥10 mapped genes per direction, `maxRank=1500`, unique `(cell_id, program_id)`, and `signed_response_score == UCell_up - UCell_down` to tolerance `1e-14`.

- [ ] **Step 2: Verify failure**

Run: `Rscript BRS-Map/tests/test_signed_ucell.R`

- [ ] **Step 3: Implement scoring and validation**

Implement `default_ucell_scorer()`, `extract_ucell_vector()`, `score_one_object_programs()` and `validate_score_sidecar()`. Preserve cell order and write no object-level expression data to the repository.

- [ ] **Step 4: Run tests and commit**

```bash
Rscript BRS-Map/tests/test_signed_ucell.R
git add BRS-Map/R/signed_ucell.R BRS-Map/tests/test_signed_ucell.R
git commit -m "feat: add signed UCell scoring"
```

### Task 6: Implement both strict-100 weighted-kNN transfers

**Files:**
- Create: `BRS-Map/R/weighted_knn.R`
- Create: `BRS-Map/tests/test_weighted_knn.R`

- [ ] **Step 1: Write failing mathematical contract tests**

Test 99-source rejection and exact-100 acceptance in both rounds, `k=20`, PCA 1:30, exact cell-type isolation, Gaussian kth-distance scaling, inverse source-unit-size balancing, normalized weights, and finite audit metrics.

```r
fail <- run_round1_one_program(source99, target, "p", "A")
pass <- run_round1_one_program(source100, target, "p", "A")
stopifnot(all(fail$status == "ROUND1_SOURCE_LT_100"),
          all(pass$status == "PASS"),
          all(pass$n_neighbors == 20L))
```

- [ ] **Step 2: Verify failure**

Run: `Rscript BRS-Map/tests/test_weighted_knn.R`

- [ ] **Step 3: Implement the transfer core**

Implement RANN nearest-neighbor search. For each target and neighbor `i`, calculate `g_i = exp(-(d_i/sigma)^2)`, set `sigma` to the kth distance with the documented positive-distance fallback, calculate `w_i = g_i / n(unit_i)`, and normalize `w`. Return transferred mean, weighted SD, number of neighbors/units, `1/sum(w^2)`, maximum summed unit weight, nearest/kth distances and source count.

- [ ] **Step 4: Implement stage wrappers**

Round 1 transfers external `signed_response_score` to same-type snRNA cells and requires ≥100 external cells. Round 2 transfers passed snRNA scores to same-type Stereo-seq CellBins and requires ≥100 scored snRNA cells. Use `predicted_classes4` and retain `UNRESOLVED_DO_NOT_GUESS` for AP order.

- [ ] **Step 5: Run tests and commit**

```bash
Rscript BRS-Map/tests/test_weighted_knn.R
git add BRS-Map/R/weighted_knn.R BRS-Map/tests/test_weighted_knn.R
git commit -m "feat: add strict-100 two-stage weighted kNN"
```

### Task 7: Add audit, preflight and orchestrator

**Files:**
- Create: `BRS-Map/R/audit.R`
- Create: `BRS-Map/scripts/preflight.R`
- Create: `BRS-Map/scripts/run_brs_map.R`
- Create: `BRS-Map/tests/test_orchestrator.R`

- [ ] **Step 1: Write failing orchestration tests**

Verify stage selection (`preflight`, `pseudobulk`, `de`, `programs`, `ucell`, `round1`, `round2`, `audit`, `all`), refusal to overwrite formal output without `--resume`, absence of server-specific defaults, and complete release audit fields.

- [ ] **Step 2: Verify failure**

Run: `Rscript BRS-Map/tests/test_orchestrator.R`

- [ ] **Step 3: Implement audit and CLI configuration**

All input/output paths come from command-line arguments or environment variables. `preflight.R` checks R/package versions, config hashes and schemas without creating formal output. `run_brs_map.R` sources only module-local R files and writes a SHA256-linked manifest.

- [ ] **Step 4: Run all code tests and commit**

```bash
for f in BRS-Map/tests/test_*.R; do Rscript "$f"; done
git add BRS-Map/R/audit.R BRS-Map/scripts BRS-Map/tests/test_orchestrator.R
git commit -m "feat: orchestrate and audit BRS-Map"
```

### Task 8: Write the concise English Methods and usage documentation

**Files:**
- Create: `BRS-Map/docs/METHODS.md`
- Modify: `BRS-Map/README.md`
- Create: `BRS-Map/environment/session-info.txt`
- Create: `BRS-Map/tests/test_documentation_contract.R`

- [ ] **Step 1: Write failing documentation-contract tests**

Require the Methods to contain all three evidence dimensions, 20-cell/2-unit/3-unit thresholds, identity 0.30/0.50 and 0.70 feature coverage, FDR/logFC/10/20/100/maxRank parameters, k/PCA/strict-100 parameters, exact matching, unit balancing, `predicted_classes4`, and unresolved AP order. Reject public occurrences of `pool-only` except in an explicit legacy-alias note.

- [ ] **Step 2: Draft the English Methods**

Use the argument: “BRS-Map converts contrast-specific, cell-type-resolved transcriptional responses into signed per-cell scores and transfers them through a fixed snRNA reference to spatial CellBins while keeping sample, program and identity evidence explicit and orthogonal.”

Write compact subsections for input harmonization, evidence model, pseudobulk/DE, signed programs/UCell, Round 1, Round 2 and audits. Include every frozen parameter without claiming AP order or formal inference for descriptive contrasts.

- [ ] **Step 3: Complete README and environment declaration**

Document RStudio and command-line entry points, required metadata columns, input/output tables, legacy aliases, and tested package versions. Link the workflow figure and Methods.

- [ ] **Step 4: Run tests and commit**

```bash
Rscript BRS-Map/tests/test_documentation_contract.R
git add BRS-Map/docs BRS-Map/README.md BRS-Map/environment BRS-Map/tests/test_documentation_contract.R README.md
git commit -m "docs: document the BRS-Map framework"
```

### Task 9: Add and validate workflow assets

**Files:**
- Create: `BRS-Map/workflow/BRS-Map_workflow.png`
- Create: `BRS-Map/workflow/BRS-Map_workflow.drawio`
- Create: `BRS-Map/workflow/BRS-Map_workflow.svg`
- Create: `BRS-Map/workflow/visual-spec.md`
- Create: `BRS-Map/workflow/layout-grid.md`
- Create: `BRS-Map/workflow/asset-ledger.md`
- Create: `BRS-Map/workflow/defect-log.md`

- [ ] **Step 1: Preserve the supplied PNG and extract its style contract**

Record canvas, six-column geometry, pastel panel colors, sans-serif typography, dark-grey arrows, editable pictograms and exact stage labels.

- [ ] **Step 2: Author editable draw.io XML**

Recreate the six stages using draw.io primitives, not a raster background. Correct the figure text to the released terminology and parameters while preserving the supplied composition.

- [ ] **Step 3: Run static preflight**

```bash
python <drawio-skill>/scripts/validate_visual_quality.py BRS-Map/workflow/BRS-Map_workflow.drawio
```

Expected: zero FAIL findings.

- [ ] **Step 4: Complete three screenshot-review cycles**

For each cycle: render, crop to canvas, inventory all nine visual zones, fix every P0/P1 issue, regenerate, verify, and append the evidence to `defect-log.md`.

- [ ] **Step 5: Export SVG and validate**

Run strict draw.io validation and verify PNG/SVG dimensions, readable labels, connector directions and no embedded external assets.

- [ ] **Step 6: Commit**

```bash
git add BRS-Map/workflow
git commit -m "docs: add editable BRS-Map workflow"
```

### Task 10: Final verification and GitHub publication

**Files:**
- Modify: `BRS-Map/PACKAGE_SHA256SUMS.txt`

- [ ] **Step 1: Parse every R file**

```bash
Rscript -e 'for (f in list.files("BRS-Map", "[.]R$", recursive=TRUE, full.names=TRUE)) parse(f); cat("R PARSE PASS\n")'
```

- [ ] **Step 2: Run the complete contract suite**

Run every `BRS-Map/tests/test_*.R`; expected: all PASS.

- [ ] **Step 3: Generate and verify SHA256 manifest**

Include all public module files except the checksum file itself; verify a clean readback.

- [ ] **Step 4: Review the release diff and repository size**

```bash
git diff master...HEAD --check
git status --short
git ls-files BRS-Map
```

Expected: only intended source, docs, config, tests and small workflow assets.

- [ ] **Step 5: Push the feature branch and merge/push to master**

Use the configured GitHub credentials and temporarily bypass only the stale local proxy for this repository. Confirm remote commit IDs after push.

- [ ] **Step 6: Report publication**

Provide the GitHub URLs for `BRS-Map/README.md`, `docs/METHODS.md`, and all three workflow files, plus the final commit hash and verification summary.
