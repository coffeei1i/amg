source(file.path(dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]))), "helper.R"))
source_module("contracts.R")

p <- brs_parameters()
stopifnot(
  p$minimum_cells_per_unit_celltype == 20L,
  p$formal_minimum_units_per_arm == 3L,
  p$low_power_minimum_units_per_arm == 2L,
  p$formal_fdr == 0.05,
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
stopifnot(
  brs_public_evidence("formal") == "formal",
  brs_public_evidence("low_power") == "low_power",
  brs_public_evidence("descriptive") == "descriptive"
)
for (invalid_label in c("invalid_label", "unknown")) {
  bad <- try(brs_public_evidence(invalid_label), silent = TRUE)
  stopifnot(inherits(bad, "try-error"))
}

cat("CONTRACT TEST PASS\n")
