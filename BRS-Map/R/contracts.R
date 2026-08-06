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
    identity_primary_score = 0.50,
    identity_sensitivity_score = 0.30,
    minimum_query_feature_coverage = 0.70,
    knn_k = 20L,
    pca_dimensions = 1:30,
    round1_minimum_cells = 100L,
    round2_minimum_cells = 100L,
    stereo_identity_column = "predicted_classes4",
    ap_order_status = "UNRESOLVED_DO_NOT_GUESS"
  )
}

brs_public_evidence <- function(value) {
  value <- as.character(value)
  map <- c(
    FORMAL = "formal",
    formal = "formal",
    AUXILIARY_LOW_POWER = "low_power",
    low_power = "low_power",
    EXPLORATORY_POOL = "descriptive",
    pool_only = "descriptive",
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
