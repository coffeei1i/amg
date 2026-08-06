options(stringsAsFactors = FALSE)

default_ucell_scorer <- function(expression, features, maxRank) {
  if (!requireNamespace("UCell", quietly = TRUE)) {
    stop("Package UCell is required", call. = FALSE)
  }
  UCell::ScoreSignatures_UCell(
    expression, features = features, maxRank = as.integer(maxRank)
  )
}

extract_ucell_vector <- function(scores, prefix, cells) {
  if (is.null(dim(scores)) || length(dim(scores)) != 2L) {
    stop("UCell returned a non-tabular result", call. = FALSE)
  }
  column <- grep(paste0("^", prefix, "(_UCell)?$"),
                 colnames(scores), value = TRUE)
  if (length(column) != 1L) stop("Unexpected UCell score columns", call. = FALSE)
  if (!is.null(rownames(scores))) {
    index <- match(cells, rownames(scores))
    if (anyNA(index)) stop("UCell cell IDs do not match", call. = FALSE)
    scores <- scores[index, , drop = FALSE]
  }
  value <- as.numeric(scores[, column])
  if (length(value) != length(cells) || any(!is.finite(value))) {
    stop("Invalid UCell score vector", call. = FALSE)
  }
  value
}

validate_score_sidecar <- function(scores) {
  required <- c(
    "cell_id", "program_id", "Celltype", "identity_mode",
    "sample_evidence", "signature_method", "n_up_genes_mapped",
    "n_down_genes_mapped", "UCell_up", "UCell_down",
    "signed_response_score"
  )
  missing <- setdiff(required, colnames(scores))
  if (length(missing)) stop("Score sidecar lacks required columns", call. = FALSE)
  if (nrow(scores)) {
    if (anyDuplicated(scores[c("cell_id", "program_id")])) {
      stop("Duplicate cell-program scores", call. = FALSE)
    }
    if (any(scores$n_up_genes_mapped < 10L) ||
        any(scores$n_down_genes_mapped < 10L)) {
      stop("Mapped signed program has fewer than 10 genes per direction",
           call. = FALSE)
    }
    expected <- scores$UCell_up - scores$UCell_down
    if (!isTRUE(all.equal(scores$signed_response_score, expected,
                          tolerance = 1e-14))) {
      stop("Signed score must equal UCell_up minus UCell_down", call. = FALSE)
    }
  }
  invisible(scores)
}

score_one_object_programs <- function(counts, metadata, programs, gene_sets,
                                      scorer = default_ucell_scorer,
                                      accession, object_id,
                                      max_rank = brs_parameters()$ucell_max_rank,
                                      minimum_mapped = brs_parameters()$minimum_genes_per_direction) {
  if (!identical(colnames(counts), rownames(metadata))) {
    stop("Counts and metadata order mismatch", call. = FALSE)
  }
  selected <- programs[
    programs$accession == accession & programs$status == "released", , drop = FALSE
  ]
  output <- list()
  for (i in seq_len(nrow(selected))) {
    program <- selected[i, , drop = FALSE]
    cells <- which(
      select_cells_by_identity_mode(metadata, program$identity_mode) &
        metadata$Celltype == program$Celltype
    )
    if (!length(cells)) next
    genes <- gene_sets[[program$program_id]]
    up <- intersect(genes$up, rownames(counts))
    down <- intersect(genes$down, rownames(counts))
    if (length(up) < minimum_mapped || length(down) < minimum_mapped) {
      stop("Fewer than 10 mapped genes per direction", call. = FALSE)
    }
    cell_id <- colnames(counts)[cells]
    result <- scorer(
      counts[, cells, drop = FALSE], features = list(up = up, down = down),
      maxRank = max_rank
    )
    up_score <- extract_ucell_vector(result, "up", cell_id)
    down_score <- extract_ucell_vector(result, "down", cell_id)
    output[[length(output) + 1L]] <- data.frame(
      cell_id = cell_id, accession = accession, object_id = object_id,
      Celltype = metadata$Celltype[cells], condition = metadata$condition[cells],
      experimental_unit_id = metadata$experimental_unit_id[cells],
      mapping_status = metadata$mapping_status[cells],
      program_id = program$program_id, contrast_id = program$contrast_id,
      identity_mode = program$identity_mode,
      sample_evidence = program$sample_evidence,
      signature_method = program$signature_method,
      n_up_genes_mapped = length(up), n_down_genes_mapped = length(down),
      UCell_up = up_score, UCell_down = down_score,
      signed_response_score = up_score - down_score,
      stringsAsFactors = FALSE
    )
  }
  scores <- if (length(output)) do.call(rbind, output) else data.frame(
    cell_id = character(), program_id = character(), Celltype = character(),
    identity_mode = character(), sample_evidence = character(),
    signature_method = character(), n_up_genes_mapped = integer(),
    n_down_genes_mapped = integer(), UCell_up = numeric(), UCell_down = numeric(),
    signed_response_score = numeric(), stringsAsFactors = FALSE
  )
  validate_score_sidecar(scores)
  scores
}
