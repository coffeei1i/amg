options(stringsAsFactors = FALSE)

pc_columns <- function(x, dimensions) {
  columns <- paste0("PC_", as.integer(dimensions))
  missing <- setdiff(columns, colnames(x))
  if (length(missing)) stop("Missing PCA columns", call. = FALSE)
  columns
}

empty_knn_result <- function(target_id, status, source_n) {
  data.frame(
    target_id = as.character(target_id), transferred_score = NA_real_,
    weighted_sd = NA_real_, n_neighbors = 0L, n_units = 0L,
    effective_n = NA_real_, max_unit_weight = NA_real_,
    nearest_distance = NA_real_, kth_distance = NA_real_,
    source_n = as.integer(source_n), status = status,
    stringsAsFactors = FALSE
  )
}

nearest_neighbors <- function(source, target, k) {
  if (requireNamespace("RANN", quietly = TRUE)) {
    return(RANN::nn2(source, target, k = k, searchtype = "standard"))
  }
  index <- matrix(NA_integer_, nrow(target), k)
  distance <- matrix(NA_real_, nrow(target), k)
  for (i in seq_len(nrow(target))) {
    d <- sqrt(rowSums((source - matrix(
      target[i, ], nrow(source), ncol(source), byrow = TRUE
    ))^2))
    order_i <- order(d, seq_along(d))[seq_len(k)]
    index[i, ] <- order_i
    distance[i, ] <- d[order_i]
  }
  list(nn.idx = index, nn.dists = distance)
}

weighted_knn_predict <- function(source_matrix, target_matrix, source_score,
                                 source_unit, target_id,
                                 k = brs_parameters()$knn_k) {
  source_matrix <- as.matrix(source_matrix)
  target_matrix <- as.matrix(target_matrix)
  storage.mode(source_matrix) <- "double"
  storage.mode(target_matrix) <- "double"
  source_score <- as.numeric(source_score)
  source_unit <- as.character(source_unit)
  target_id <- as.character(target_id)
  if (ncol(source_matrix) != ncol(target_matrix)) {
    stop("PCA dimension mismatch", call. = FALSE)
  }
  if (nrow(source_matrix) != length(source_score) ||
      nrow(source_matrix) != length(source_unit)) {
    stop("Source input length mismatch", call. = FALSE)
  }
  if (nrow(target_matrix) != length(target_id) || anyDuplicated(target_id)) {
    stop("Target identity contract failed", call. = FALSE)
  }
  if (any(!is.finite(source_matrix)) || any(!is.finite(target_matrix)) ||
      any(!is.finite(source_score)) || any(!nzchar(source_unit))) {
    stop("Non-finite or blank kNN input", call. = FALSE)
  }
  if (!length(target_id)) {
    return(empty_knn_result(character(), "NO_MATCHING_TARGETS", nrow(source_matrix)))
  }
  k_use <- min(as.integer(k), nrow(source_matrix))
  nearest <- nearest_neighbors(source_matrix, target_matrix, k_use)
  unit_size <- table(source_unit)
  result <- vector("list", nrow(target_matrix))
  for (i in seq_len(nrow(target_matrix))) {
    index <- nearest$nn.idx[i, ]
    distance <- nearest$nn.dists[i, ]
    unit <- source_unit[index]
    sigma <- distance[[length(distance)]]
    if (!is.finite(sigma) || sigma <= 0) {
      positive <- distance[distance > 0]
      sigma <- if (length(positive)) max(positive) else 1
    }
    weight <- exp(-(distance / sigma)^2) / as.numeric(unit_size[unit])
    if (!is.finite(sum(weight)) || sum(weight) <= 0) weight[] <- 1
    weight <- weight / sum(weight)
    value <- source_score[index]
    estimate <- sum(weight * value)
    unit_weight <- rowsum(weight, group = unit, reorder = FALSE)
    result[[i]] <- data.frame(
      target_id = target_id[[i]], transferred_score = estimate,
      weighted_sd = sqrt(max(0, sum(weight * (value - estimate)^2))),
      n_neighbors = length(index), n_units = nrow(unit_weight),
      effective_n = 1 / sum(weight^2), max_unit_weight = max(unit_weight),
      nearest_distance = distance[[1L]],
      kth_distance = distance[[length(distance)]],
      source_n = nrow(source_matrix), status = "PASS",
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, result)
}

run_round1_one_program <- function(source, target, program_id, celltype,
                                   k = brs_parameters()$knn_k,
                                   dimensions = brs_parameters()$pca_dimensions,
                                   minimum_source_cells = brs_parameters()$round1_minimum_cells) {
  source <- source[source$Celltype == celltype, , drop = FALSE]
  target <- target[target$Celltype == celltype, , drop = FALSE]
  source_n <- nrow(source)
  if (!nrow(target)) {
    out <- empty_knn_result(character(), "NO_MATCHING_TARGETS", source_n)
  } else if (source_n < minimum_source_cells) {
    out <- empty_knn_result(target$cell_id, "ROUND1_SOURCE_LT_100", source_n)
  } else {
    out <- weighted_knn_predict(
      source[, pc_columns(source, dimensions), drop = FALSE],
      target[, pc_columns(target, dimensions), drop = FALSE],
      source$signed_response_score, source$experimental_unit_id,
      target$cell_id, k
    )
  }
  out$program_id <- program_id
  out$Celltype <- celltype
  out$round1_minimum_source_cells <- minimum_source_cells
  out
}

run_round2_one_program <- function(source, stereo, program_id, celltype,
                                   k = brs_parameters()$knn_k,
                                   dimensions = brs_parameters()$pca_dimensions,
                                   minimum_reference_cells = brs_parameters()$round2_minimum_cells,
                                   identity_column = brs_parameters()$stereo_identity_column) {
  source <- source[
    source$Celltype == celltype & source$status == "PASS" &
      is.finite(source$transferred_score), , drop = FALSE
  ]
  target <- stereo[stereo$Celltype == celltype, , drop = FALSE]
  source_n <- nrow(source)
  if (!nrow(target)) {
    out <- empty_knn_result(character(), "NO_MATCHING_TARGETS", source_n)
  } else if (source_n < minimum_reference_cells) {
    out <- empty_knn_result(target$target_id, "ROUND2_REFERENCE_LT_100", source_n)
  } else {
    out <- weighted_knn_predict(
      source[, pc_columns(source, dimensions), drop = FALSE],
      target[, pc_columns(target, dimensions), drop = FALSE],
      source$transferred_score, source$reference_unit,
      target$target_id, k
    )
  }
  out$program_id <- program_id
  out$Celltype <- celltype
  out$round2_minimum_reference_cells <- minimum_reference_cells
  out$spatial_identity_column <- identity_column
  out$ap_order_status <- brs_parameters()$ap_order_status
  out
}
