options(stringsAsFactors = FALSE)

brs_parameters <- function() {
  list(
    minimum_cells_per_unit_celltype = 20L,
    formal_minimum_units_per_arm = 3L,
    low_power_minimum_units_per_arm = 2L,
    formal_fdr = 0.05,
    formal_abs_logfc = 0.25,
    minimum_genes_per_direction = 10L,
    formal_maximum_genes_per_direction = 100L,
    fallback_maximum_genes_per_direction = 20L,
    ucell_max_rank = 1500L,
    identity_reference_celltypes = 49L,
    visualized_neuronal_celltypes = 41L,
    identity_anchor_normalization = "SCT",
    identity_anchor_k = 5L,
    identity_anchor_max_features = 3000L,
    identity_anchor_min_shared_features = 200L,
    identity_score_k_cap = 30L,
    identity_weight_k_cap = 50L,
    identity_mapping_seed = 1729L,
    identity_primary_score = 0.50,
    identity_sensitivity_score = 0.30,
    minimum_query_feature_coverage = 0.70,
    knn_k = 20L,
    pca_dimensions = 1:30,
    round1_minimum_cells = 100L,
    round2_minimum_cells = 100L,
    direction_zero_tolerance = 0,
    landscape_scaling_quantile = 1,
    spatial_scaling_quantile = 0.95,
    visualization_lower_limit = -1,
    visualization_upper_limit = 1,
    identity_selection_policy = "primary_then_sensitivity",
    stereo_identity_column = "predicted_classes4",
    ap_order_status = "UNRESOLVED_DO_NOT_GUESS"
  )
}

brs_public_evidence <- function(value) {
  value <- as.character(value)
  map <- c(
    formal = "formal",
    low_power = "low_power",
    descriptive = "descriptive"
  )
  result <- unname(map[value])
  if (anyNA(result)) {
    stop("Unknown sample/statistical evidence label: ",
         paste(unique(value[is.na(result)]), collapse = ", "), call. = FALSE)
  }
  result
}

brs_identity_mode <- function(mapping_status, mode = c("primary", "sensitivity")) {
  mode <- match.arg(mode)
  status <- tolower(as.character(mapping_status))
  if (mode == "primary") status == "accepted" else
    status %in% c("accepted", "sensitivity")
}

validate_identity_metadata <- function(metadata, celltypes = NULL) {
  required <- c(
    "Celltype", "mapping_status", "condition", "experimental_unit_id"
  )
  missing <- setdiff(required, colnames(metadata))
  if (length(missing)) {
    stop("Identity metadata lacks: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (is.null(rownames(metadata)) || any(!nzchar(rownames(metadata))) ||
      anyDuplicated(rownames(metadata))) {
    stop("Cell identifiers must be non-blank and unique", call. = FALSE)
  }
  allowed <- c("accepted", "sensitivity", "unassigned")
  status <- tolower(as.character(metadata$mapping_status))
  if (any(!status %in% allowed)) {
    stop("Unsupported mapping_status", call. = FALSE)
  }
  if (!is.null(celltypes)) {
    value <- as.character(metadata$Celltype)
    assigned <- status != "unassigned"
    if (any(assigned & !value %in% celltypes)) {
      stop("Mapped Celltype is outside the frozen vocabulary", call. = FALSE)
    }
  }
  invisible(metadata)
}

validate_prediction_scores <- function(score, feature_coverage = NULL) {
  score <- as.numeric(score)
  if (any(!is.finite(score)) || any(score < 0 | score > 1)) {
    stop("Identity prediction scores must be finite probabilities", call. = FALSE)
  }
  if (!is.null(feature_coverage)) {
    coverage <- as.numeric(feature_coverage)
    if (any(!is.finite(coverage)) || any(coverage < 0 | coverage > 1)) {
      stop("Query feature coverage must lie in [0,1]", call. = FALSE)
    }
  }
  invisible(TRUE)
}
