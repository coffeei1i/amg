source(file.path(dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]))), "helper.R"))
source_module("contracts.R")
source_module("pseudobulk.R")

stopifnot(
  classify_sample_evidence(3L, 4L, 2L, 3L, FALSE) == "formal",
  classify_sample_evidence(2L, 3L, 2L, 3L, FALSE) == "low_power",
  classify_sample_evidence(1L, 1L, 1L, 3L, TRUE) == "descriptive",
  is.na(classify_sample_evidence(1L, 2L, 2L, 3L, FALSE))
)

meta <- data.frame(
  Celltype = "A", condition = rep(c("case", "control"), each = 40),
  experimental_unit_id = rep(c("u1", "u2", "u3", "u4"), each = 20),
  mapping_status = rep(c("accepted", "sensitivity"), 40),
  row.names = paste0("c", 1:80), stringsAsFactors = FALSE
)
counts <- matrix(1:240, nrow = 3, dimnames = list(
  paste0("g", 1:3), rownames(meta)
))
primary <- aggregate_pseudobulk_counts(counts, meta, "primary", min_cells = 10L)
sensitivity <- aggregate_pseudobulk_counts(
  counts, meta, "sensitivity", min_cells = 20L
)
stopifnot(sum(primary$counts) == sum(counts[, seq(1, 80, by = 2)]))
stopifnot(sum(sensitivity$counts) == sum(counts))
stopifnot(all(sensitivity$sample_meta$n_cells == 20L))

cat("PSEUDOBULK TEST PASS\n")
