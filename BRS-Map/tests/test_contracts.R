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
  p$identity_reference_celltypes == 49L,
  p$visualized_neuronal_celltypes == 41L,
  identical(p$identity_anchor_normalization, "SCT"),
  p$identity_anchor_k == 5L,
  p$identity_anchor_max_features == 3000L,
  p$identity_anchor_min_shared_features == 200L,
  p$identity_score_k_cap == 30L,
  p$identity_weight_k_cap == 50L,
  p$identity_mapping_seed == 1729L,
  p$identity_primary_score == 0.50,
  p$identity_sensitivity_score == 0.30,
  p$minimum_query_feature_coverage == 0.70,
  p$knn_k == 20L,
  identical(p$pca_dimensions, 1:30),
  p$round1_minimum_cells == 100L,
  p$round2_minimum_cells == 100L,
  p$spatial_scaling_quantile == 0.95,
  p$landscape_scaling_quantile == 1,
  p$visualization_lower_limit == -1,
  p$visualization_upper_limit == 1
)

contrast_path <- file.path(
  dirname(dirname(normalizePath(sub("^--file=", "", grep(
    "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]])))),
  "config", "contrasts.tsv"
)
contrasts <- read.delim(contrast_path, check.names = FALSE, stringsAsFactors = FALSE)
gse256522 <- contrasts[contrasts$accession == "GSE256522", , drop = FALSE]
stopifnot(
  nrow(gse256522) == 3L,
  setequal(gse256522$contrast_id, c(
    "stfp_vs_odour_only",
    "stfp_vs_home_cage",
    "odour_only_vs_home_cage"
  )),
  all(gse256522$configured_evidence == "formal")
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
