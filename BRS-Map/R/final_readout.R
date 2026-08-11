options(stringsAsFactors = FALSE)

require_brs_columns <- function(x, required, label = "table") {
  if (!is.data.frame(x)) stop(label, " must be a data.frame", call. = FALSE)
  missing <- setdiff(required, colnames(x))
  if (length(missing)) stop(label, " lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(x)
}

brs_group_key <- function(x, columns) {
  require_brs_columns(x, columns)
  values <- lapply(x[columns], as.character)
  invalid <- Reduce(`|`, lapply(values, function(v) is.na(v) | !nzchar(v)))
  if (any(invalid)) stop("Grouping columns contain blank values", call. = FALSE)
  do.call(paste, c(values, sep = "\036"))
}

aggregate_program_mean <- function(x, keys, score_col, status_col = "status",
                                   pass_status = "PASS") {
  required <- c(keys, score_col)
  if (!is.null(status_col)) required <- c(required, status_col)
  require_brs_columns(x, required)
  keep <- is.finite(as.numeric(x[[score_col]]))
  if (!is.null(status_col)) keep <- keep & as.character(x[[status_col]]) == pass_status
  x <- x[keep, c(keys, score_col), drop = FALSE]
  if (!nrow(x)) {
    out <- as.data.frame(setNames(replicate(length(keys), character(), simplify = FALSE), keys))
    out$mean_score <- numeric()
    out$n_scored <- integer()
    return(out)
  }
  pieces <- split(seq_len(nrow(x)), brs_group_key(x, keys))
  do.call(rbind, lapply(pieces, function(index) {
    row <- x[index[[1L]], keys, drop = FALSE]
    row$mean_score <- mean(as.numeric(x[[score_col]][index]))
    row$n_scored <- length(index)
    row
  }))
}

direction_status <- function(round1_mean, round2_mean,
                             tolerance = brs_parameters()$direction_zero_tolerance) {
  round1_mean <- as.numeric(round1_mean)
  round2_mean <- as.numeric(round2_mean)
  if (length(round1_mean) != length(round2_mean)) {
    stop("Round-1 and Round-2 vectors must have equal length", call. = FALSE)
  }
  result <- rep("NOT_EVALUATED", length(round1_mean))
  finite <- is.finite(round1_mean) & is.finite(round2_mean)
  zero <- finite & (abs(round1_mean) <= tolerance | abs(round2_mean) <= tolerance)
  result[finite & !zero & sign(round1_mean) == sign(round2_mean)] <- "SAME_DIRECTION"
  result[finite & !zero & sign(round1_mean) != sign(round2_mean)] <- "DIRECTION_REVERSAL"
  result[zero] <- "AMBIGUOUS_ZERO"
  result
}

audit_direction_consistency <- function(
    round1, round2, keys = c("program_id", "Celltype"),
    score_col = "transferred_score", status_col = "status",
    pass_status = "PASS",
    tolerance = brs_parameters()$direction_zero_tolerance) {
  r1 <- aggregate_program_mean(round1, keys, score_col, status_col, pass_status)
  r2 <- aggregate_program_mean(round2, keys, score_col, status_col, pass_status)
  names(r1)[names(r1) == "mean_score"] <- "round1_mean_score"
  names(r1)[names(r1) == "n_scored"] <- "round1_n_scored"
  names(r2)[names(r2) == "mean_score"] <- "round2_mean_score"
  names(r2)[names(r2) == "n_scored"] <- "round2_n_scored"
  out <- merge(r1, r2, by = keys, all = TRUE, sort = FALSE)
  out$direction_status <- direction_status(out$round1_mean_score, out$round2_mean_score, tolerance)
  out$eligible_for_final <- out$direction_status == "SAME_DIRECTION"
  out
}

select_final_identity <- function(
    primary, sensitivity, keys = c("program_id", "Celltype"),
    status_col = "status", pass_status = "PASS",
    direction_col = "direction_status") {
  required <- c(keys, status_col, direction_col, "signature_method", "sample_evidence")
  require_brs_columns(primary, required, "primary")
  require_brs_columns(sensitivity, required, "sensitivity")
  eligible <- function(x) x[
    as.character(x[[status_col]]) == pass_status &
      as.character(x[[direction_col]]) == "SAME_DIRECTION", , drop = FALSE
  ]
  p <- eligible(primary)
  s <- eligible(sensitivity)
  if ((nrow(p) && anyDuplicated(brs_group_key(p, keys))) ||
      (nrow(s) && anyDuplicated(brs_group_key(s, keys)))) {
    stop("Identity candidates must contain one row per key", call. = FALSE)
  }
  p$identity_source <- rep("primary", nrow(p))
  p$identity_marker <- rep("", nrow(p))
  p_key <- if (nrow(p)) brs_group_key(p, keys) else character()
  s_key <- if (nrow(s)) brs_group_key(s, keys) else character()
  s <- s[!s_key %in% p_key, , drop = FALSE]
  s$identity_source <- rep("sensitivity", nrow(s))
  s$identity_marker <- rep("S", nrow(s))
  columns <- union(colnames(p), colnames(s))
  add_missing <- function(x) {
    for (column in setdiff(columns, colnames(x))) x[[column]] <- rep(NA, nrow(x))
    x[columns]
  }
  out <- rbind(add_missing(p), add_missing(s))
  fallback <- tolower(as.character(out$signature_method)) %in%
    c("fallback", "ranked_logfc_fallback", "ranked_logfc_fallback_top20")
  out$program_marker <- ifelse(fallback, "*", "")
  out$point_marker <- paste0(out$identity_marker, out$program_marker)
  rownames(out) <- NULL
  out
}

scale_scores_by_group <- function(
    x, group_cols, score_col, quantile_probability = 1,
    lower = brs_parameters()$visualization_lower_limit,
    upper = brs_parameters()$visualization_upper_limit) {
  require_brs_columns(x, c(group_cols, score_col))
  probability <- as.numeric(quantile_probability)
  if (length(probability) != 1L || !is.finite(probability) || probability <= 0 || probability > 1) {
    stop("quantile_probability must lie in (0,1]", call. = FALSE)
  }
  out <- x
  out$scale_denominator <- NA_real_
  out$scaled_score <- NA_real_
  groups <- split(seq_len(nrow(out)), brs_group_key(out, group_cols))
  for (index in groups) {
    raw <- as.numeric(out[[score_col]][index])
    finite <- is.finite(raw)
    denominator <- if (any(finite)) as.numeric(stats::quantile(
      abs(raw[finite]), probs = probability, names = FALSE, type = 7
    )) else NA_real_
    out$scale_denominator[index] <- denominator
    if (is.finite(denominator) && denominator > 0) {
      out$scaled_score[index[finite]] <- pmax(lower, pmin(upper, raw[finite] / denominator))
    } else if (is.finite(denominator) && denominator == 0) {
      out$scaled_score[index[finite]] <- 0
    }
  }
  out
}

summarise_celltype_landscape <- function(
    scores, group_cols = c("accession", "contrast_id", "Celltype"),
    score_col = "transferred_score") {
  require_brs_columns(scores, c(group_cols, score_col))
  if (!nrow(scores)) {
    out <- scores[FALSE, group_cols, drop = FALSE]
    out$n_cells <- integer()
    out$n_valid_scores <- integer()
    out$valid_score_fraction <- numeric()
    out$scale_denominator <- numeric()
    out$mean_scaled_score <- numeric()
    return(out)
  }
  scaled <- scale_scores_by_group(
    scores, group_cols, score_col,
    quantile_probability = brs_parameters()$landscape_scaling_quantile
  )
  groups <- split(seq_len(nrow(scaled)), brs_group_key(scaled, group_cols))
  do.call(rbind, lapply(groups, function(index) {
    row <- scaled[index[[1L]], group_cols, drop = FALSE]
    valid <- is.finite(scaled$scaled_score[index])
    row$n_cells <- length(index)
    row$n_valid_scores <- sum(valid)
    row$valid_score_fraction <- mean(valid)
    row$scale_denominator <- unique(scaled$scale_denominator[index])[[1L]]
    row$mean_scaled_score <- if (any(valid)) mean(scaled$scaled_score[index][valid]) else NA_real_
    row
  }))
}

scale_spatial_joint_q95 <- function(scores, program_col = "program_id",
                                    score_col = "transferred_score") {
  out <- scale_scores_by_group(
    scores, program_col, score_col,
    quantile_probability = brs_parameters()$spatial_scaling_quantile
  )
  out$scaling_scope <- "PROGRAM_ACROSS_ALL_SLICES"
  out$scaling_method <- "Q95_ABSOLUTE_VALUE_CLIPPED_TO_MINUS1_PLUS1"
  out
}

summarise_subnuclear_directional_balance <- function(
    scores, group_cols = c("accession", "contrast_id", "subnucleus"),
    score_col = "transferred_score", target_col = "target_id") {
  require_brs_columns(scores, c(group_cols, score_col, target_col))
  if (!nrow(scores)) {
    out <- scores[FALSE, group_cols, drop = FALSE]
    out$total_neurons <- integer()
    out$valid_neurons <- integer()
    out$valid_neuron_proportion <- numeric()
    out$positive_n <- integer()
    out$negative_n <- integer()
    out$zero_n <- integer()
    out$directional_balance <- numeric()
    return(out)
  }
  pieces <- split(seq_len(nrow(scores)), brs_group_key(scores, c(group_cols, target_col)))
  unique_target <- do.call(rbind, lapply(pieces, function(index) {
    row <- scores[index[[1L]], c(group_cols, target_col), drop = FALSE]
    value <- as.numeric(scores[[score_col]][index])
    row$pooled_score <- if (any(is.finite(value))) mean(value[is.finite(value)]) else NA_real_
    row
  }))
  groups <- split(seq_len(nrow(unique_target)), brs_group_key(unique_target, group_cols))
  do.call(rbind, lapply(groups, function(index) {
    row <- unique_target[index[[1L]], group_cols, drop = FALSE]
    value <- unique_target$pooled_score[index]
    valid <- is.finite(value)
    positive_n <- sum(value[valid] > 0)
    negative_n <- sum(value[valid] < 0)
    zero_n <- sum(value[valid] == 0)
    denominator <- positive_n + negative_n
    row$total_neurons <- length(index)
    row$valid_neurons <- sum(valid)
    row$valid_neuron_proportion <- mean(valid)
    row$positive_n <- positive_n
    row$negative_n <- negative_n
    row$zero_n <- zero_n
    row$directional_balance <- if (denominator > 0) (positive_n - negative_n) / denominator else NA_real_
    row
  }))
}
