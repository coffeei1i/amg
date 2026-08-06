options(stringsAsFactors = FALSE)

validate_de_table <- function(de) {
  required <- c("gene", "logFC", "FDR")
  if (!all(required %in% colnames(de))) {
    stop("DE table must contain gene, logFC, and FDR", call. = FALSE)
  }
  de$gene <- as.character(de$gene)
  de$logFC <- suppressWarnings(as.numeric(de$logFC))
  de$FDR <- suppressWarnings(as.numeric(de$FDR))
  if (any(!nzchar(de$gene)) || anyDuplicated(de$gene)) {
    stop("DE gene names must be non-blank and unique", call. = FALSE)
  }
  de
}

remove_program_blacklist <- function(de, blacklist = character()) {
  blacklist <- unique(as.character(blacklist))
  de[!de$gene %in% blacklist[nzchar(blacklist)], , drop = FALSE]
}

rank_formal_genes <- function(de, blacklist = character(),
                              fdr_cutoff = brs_parameters()$formal_fdr,
                              abs_logfc_min = brs_parameters()$formal_abs_logfc,
                              maximum = brs_parameters()$formal_maximum_genes_per_direction) {
  de <- remove_program_blacklist(validate_de_table(de), blacklist)
  keep <- is.finite(de$FDR) & is.finite(de$logFC) &
    de$FDR < fdr_cutoff & abs(de$logFC) >= abs_logfc_min
  up <- de[keep & de$logFC > 0, , drop = FALSE]
  down <- de[keep & de$logFC < 0, , drop = FALSE]
  up <- up[order(up$FDR, -up$logFC, up$gene), , drop = FALSE]
  down <- down[order(down$FDR, down$logFC, down$gene), , drop = FALSE]
  list(up = head(up, maximum), down = head(down, maximum))
}

rank_fallback_genes <- function(de, blacklist = character(),
                                maximum = brs_parameters()$fallback_maximum_genes_per_direction) {
  de <- remove_program_blacklist(validate_de_table(de), blacklist)
  keep <- is.finite(de$logFC) & de$logFC != 0
  up <- de[keep & de$logFC > 0, , drop = FALSE]
  down <- de[keep & de$logFC < 0, , drop = FALSE]
  up <- up[order(-up$logFC, up$gene), , drop = FALSE]
  down <- down[order(down$logFC, down$gene), , drop = FALSE]
  list(up = head(up, maximum), down = head(down, maximum))
}

build_signed_program <- function(de, metadata, blacklist = character(),
                                 minimum = brs_parameters()$minimum_genes_per_direction) {
  sample_evidence <- brs_public_evidence(metadata$sample_evidence)
  formal <- rank_formal_genes(de, blacklist)
  formal_ok <- sample_evidence != "descriptive" &&
    nrow(formal$up) >= minimum && nrow(formal$down) >= minimum
  if (formal_ok) {
    selected <- formal
    method <- "formal_signature"
  } else {
    selected <- rank_fallback_genes(de, blacklist)
    method <- "fallback"
  }
  released <- nrow(selected$up) >= minimum && nrow(selected$down) >= minimum
  c(metadata, list(
    sample_evidence = sample_evidence,
    status = if (released) "released" else "rejected_insufficient_genes",
    signature_method = if (released) method else "none",
    up_genes = as.character(selected$up$gene),
    down_genes = as.character(selected$down$gene),
    n_up = nrow(selected$up), n_down = nrow(selected$down),
    blacklist_size = length(unique(blacklist))
  ))
}

run_descriptive_contrast <- function(counts, sample_meta, case_label,
                                     control_label, offset = 0.5) {
  library_size <- colSums(counts)
  if (any(library_size <= 0)) stop("Zero library size", call. = FALSE)
  log_cpm <- log2(t(t(as.matrix(counts)) / library_size) * 1e6 + offset)
  case <- sample_meta$condition == case_label
  control <- sample_meta$condition == control_label
  if (!any(case) || !any(control)) stop("Contrast arm is empty", call. = FALSE)
  data.frame(
    gene = rownames(counts),
    logFC = rowMeans(log_cpm[, case, drop = FALSE]) -
      rowMeans(log_cpm[, control, drop = FALSE]),
    logCPM = rowMeans(log_cpm), PValue = NA_real_, FDR = NA_real_,
    inference_status = "DESCRIPTIVE_NO_REPLICATE_FDR",
    stringsAsFactors = FALSE
  )
}

run_edger_contrast <- function(counts, sample_meta, case_label, control_label,
                               sample_evidence) {
  evidence <- brs_public_evidence(sample_evidence)
  keep <- sample_meta$condition %in% c(case_label, control_label)
  x <- counts[, keep, drop = FALSE]
  m <- sample_meta[keep, , drop = FALSE]
  if (evidence == "descriptive") {
    return(run_descriptive_contrast(x, m, case_label, control_label))
  }
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package edgeR is required", call. = FALSE)
  }
  group <- factor(m$condition, levels = c(control_label, case_label))
  design <- stats::model.matrix(~ group)
  y <- edgeR::DGEList(counts = x, group = group)
  expressed <- edgeR::filterByExpr(y, design = design)
  if (sum(expressed) < 10L) stop("Too few expressed genes", call. = FALSE)
  y <- y[expressed, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  test <- edgeR::glmQLFTest(fit, coef = 2L)
  tab <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
  data.frame(
    gene = rownames(tab), logFC = tab$logFC, logCPM = tab$logCPM,
    PValue = tab$PValue, FDR = tab$FDR, inference_status = "EDGER_QLF",
    stringsAsFactors = FALSE
  )
}
