options(stringsAsFactors = FALSE)

select_cells_by_identity_mode <- function(metadata, mode = c("primary", "sensitivity")) {
  mode <- match.arg(mode)
  validate_identity_metadata(metadata)
  known <- function(x) {
    value <- trimws(as.character(x))
    nzchar(value) & !toupper(value) %in% c("NA", "N/A", "UNKNOWN", "UNRESOLVED")
  }
  brs_identity_mode(metadata$mapping_status, mode) &
    known(metadata$Celltype) & known(metadata$condition) &
    known(metadata$experimental_unit_id)
}

aggregate_pseudobulk_counts <- function(counts, metadata, identity_mode,
                                        min_cells = brs_parameters()$minimum_cells_per_unit_celltype) {
  if (!identical(colnames(counts), rownames(metadata))) {
    stop("Counts columns must exactly match metadata row names", call. = FALSE)
  }
  if (length(min_cells) != 1L || !is.finite(min_cells) || min_cells < 1L) {
    stop("min_cells must be a positive integer", call. = FALSE)
  }
  keep_base <- select_cells_by_identity_mode(metadata, identity_mode)
  key <- paste(metadata$experimental_unit_id, metadata$Celltype, sep = "\036")
  group_sizes <- table(key[keep_base])
  eligible <- names(group_sizes[group_sizes >= as.integer(min_cells)])
  keep <- keep_base & key %in% eligible

  audit <- data.frame(
    cell_id = rownames(metadata), keep = keep, identity_mode = identity_mode,
    group_key = key, stringsAsFactors = FALSE
  )
  if (!any(keep)) {
    return(list(counts = counts[, FALSE, drop = FALSE],
                sample_meta = data.frame(), cell_audit = audit))
  }

  levels <- unique(key[keep])
  if (requireNamespace("Matrix", quietly = TRUE)) {
    design <- Matrix::sparse.model.matrix(~ 0 + factor(key[keep], levels = levels))
  } else {
    design <- stats::model.matrix(~ 0 + factor(key[keep], levels = levels))
  }
  aggregated <- counts[, keep, drop = FALSE] %*% design
  colnames(aggregated) <- levels
  parts <- strsplit(levels, "\036", fixed = TRUE)
  sample_meta <- data.frame(
    group_key = levels,
    experimental_unit_id = vapply(parts, `[[`, character(1), 1L),
    Celltype = vapply(parts, `[[`, character(1), 2L),
    stringsAsFactors = FALSE
  )
  unit_meta <- unique(metadata[keep, c("experimental_unit_id", "condition")])
  if (anyDuplicated(unit_meta$experimental_unit_id)) {
    stop("One experimental unit maps to multiple conditions", call. = FALSE)
  }
  sample_meta$condition <- unit_meta$condition[match(
    sample_meta$experimental_unit_id, unit_meta$experimental_unit_id
  )]
  sample_meta$n_cells <- as.integer(group_sizes[sample_meta$group_key])
  sample_meta$identity_mode <- identity_mode
  if (!isTRUE(all.equal(as.numeric(sum(aggregated)),
                        as.numeric(sum(counts[, keep, drop = FALSE])),
                        tolerance = 0))) {
    stop("Pseudobulk count conservation failed", call. = FALSE)
  }
  list(counts = aggregated, sample_meta = sample_meta, cell_audit = audit)
}

classify_sample_evidence <- function(n_case, n_control, minimum = 2L,
                                     formal_minimum = 3L, pooled = FALSE) {
  n_case <- as.integer(n_case)
  n_control <- as.integer(n_control)
  if (isTRUE(pooled)) {
    return(if (n_case >= minimum && n_control >= minimum) "descriptive" else NA_character_)
  }
  if (n_case < minimum || n_control < minimum) return(NA_character_)
  if (n_case >= formal_minimum && n_control >= formal_minimum) {
    "formal"
  } else {
    "low_power"
  }
}

assess_contrast_eligibility <- function(sample_meta, case_label, control_label,
                                        minimum = 2L, formal_minimum = 3L,
                                        pooled = FALSE) {
  split_meta <- split(sample_meta, as.character(sample_meta$Celltype))
  do.call(rbind, lapply(names(split_meta), function(celltype) {
    x <- split_meta[[celltype]]
    n_case <- length(unique(x$experimental_unit_id[x$condition == case_label]))
    n_control <- length(unique(x$experimental_unit_id[x$condition == control_label]))
    evidence <- classify_sample_evidence(
      n_case, n_control, minimum, formal_minimum, pooled
    )
    data.frame(
      Celltype = celltype, n_case_units = n_case, n_control_units = n_control,
      sample_evidence = evidence, eligible = !is.na(evidence),
      formal_inference_allowed = identical(evidence, "formal"),
      stringsAsFactors = FALSE
    )
  }))
}
