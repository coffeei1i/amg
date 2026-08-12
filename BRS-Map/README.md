# BRS-Map

BRS-Map is an evidence-aware framework for constructing cell-type-resolved,
signed behavioral response scores from external single-cell RNA-seq datasets
and transferring those scores through a fixed snRNA-seq reference to
Stereo-seq CellBins.

![BRS-Map workflow](workflow/BRS-Map_workflow.png)

## Public input boundary

BRS-Map starts from harmonized external single-cell objects that already
contain:

- raw RNA counts;
- unique cell identifiers;
- accession and object identifiers;
- condition and independent `experimental_unit_id` metadata;
- a unified `Celltype` label drawn from 49 harmonized amygdala cell types;
- `mapping_status` (`accepted`, `sensitivity`, or `unassigned`);
- coordinates in the frozen reference PCA space.

Raw-data download, cell calling, QC, doublet removal, and Seurat anchor
re-computation are outside this module.

Anchor-based identity mapping uses the complete 49-type reference vocabulary.
Downstream neuronal landscapes and spatial summaries are restricted to the 41
neuronal types within that vocabulary; non-neuronal types remain available for
identity and audit outputs but are not included in these neuronal visualizations.

## Framework

1. Aggregate cells by independent experimental unit and cell type.
2. Classify sample/statistical evidence as `formal`, `low_power`, or `descriptive`.
3. Estimate cell-type-specific differential expression.
4. Build formal or ranked-logFC fallback signed gene programs.
5. Calculate `UCell_up - UCell_down` for each eligible external cell.
6. Transfer scores to same-type snRNA-seq cells by unit-balanced weighted kNN.
7. Transfer snRNA scores to same-type Stereo-seq CellBins using the same rule.
8. Remove Celltype programs whose aggregate Round-1 and Round-2 score directions disagree.
9. Prefer primary identity results and use an eligible sensitivity result only when primary is unavailable (`S` marker); mark ranked-logFC fallback programs with `*`.
10. Export score, uncertainty, coverage, distance, direction, scaling and provenance audits.

## Final readout contract

- Neuronal landscape: scale individual valid snRNA scores within each exact
  behavior-by-Celltype group by that group's maximum absolute value, clip to
  `[-1,1]`, and then calculate the mean scaled score. The method never scales
  already-averaged contrast columns.
- Per-program spatial maps: calculate one 95th percentile absolute-value
  denominator jointly across M1, M2 and M3 for the same program, then clip to
  `[-1,1]`. Slices are neither normalized separately nor averaged.
- Integrated spatial maps reuse those frozen program-level values without a
  second scaling operation.
- Subnuclear landscape: pool unique neurons from all three sections and report
  `(positive_n - negative_n) / (positive_n + negative_n)`; point size records
  the proportion of neurons with a valid Round-2 score.

These operations are implemented in [`R/final_readout.R`](R/final_readout.R).

See [the full English Methods](docs/METHODS.md) and the
[workflow diagram](workflow/BRS-Map_workflow.svg). Exact settings are stored in
[`config/parameters.tsv`](config/parameters.tsv).

## Repository policy

The repository contains code, small configuration tables, tests,
documentation, and workflow artwork only. Biological data and generated result
matrices are intentionally excluded.

## Running on a server

Run `Rscript scripts/run_brs_map.R --stage preflight` first. Formal stages require explicit `--input-root`, `--output-root`, and `--adapter` arguments; no private server path is embedded. See [SERVER_QUICKSTART.md](docs/SERVER_QUICKSTART.md) and copy [server_adapter_template.R](examples/server_adapter_template.R) to a private working directory before editing I/O.

## Frozen evidence model

- Sample/statistical evidence: `formal`, `low_power`, or `descriptive`.
- Gene program: `formal_signature` or ranked-logFC `fallback`.
- Cell identity: `primary` or `sensitivity`.

The dimensions are independent and are retained in every downstream audit table.

For the primary presentation, identity selection is primary-first. An eligible
sensitivity result is used only when the matching primary result is unavailable;
it is explicitly marked `S`. A ranked-logFC fallback program is marked `*`, so
`S*` means sensitivity identity plus fallback gene program.
