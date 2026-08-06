source(file.path(dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]))), "helper.R"))
source_module("contracts.R")
source_module("pseudobulk.R")
source_module("signed_ucell.R")

counts <- matrix(1, nrow = 24, ncol = 4, dimnames = list(
  paste0("g", 1:24), paste0("c", 1:4)
))
meta <- data.frame(
  accession = "GSE1", object_id = "o1", Celltype = "A",
  condition = "case", experimental_unit_id = "u1",
  mapping_status = c("accepted", "sensitivity", "unassigned", "accepted"),
  row.names = colnames(counts), stringsAsFactors = FALSE
)
programs <- data.frame(
  program_id = "p1", accession = "GSE1", contrast_id = "x",
  Celltype = "A", identity_mode = "primary", sample_evidence = "formal",
  signature_method = "formal_signature", status = "released",
  stringsAsFactors = FALSE
)
sets <- list(p1 = list(up = paste0("g", 1:12), down = paste0("g", 13:24)))
fake <- function(expression, features, maxRank) {
  stopifnot(maxRank == 1500L)
  data.frame(up_UCell = seq_len(ncol(expression)) / 10,
             down_UCell = seq_len(ncol(expression)) / 20,
             row.names = colnames(expression))
}
scores <- score_one_object_programs(
  counts, meta, programs, sets, scorer = fake,
  accession = "GSE1", object_id = "o1"
)
stopifnot(nrow(scores) == 2L)
stopifnot(all(scores$mapping_status == "accepted"))
stopifnot(isTRUE(all.equal(
  scores$signed_response_score, scores$UCell_up - scores$UCell_down,
  tolerance = 1e-14
)))

cat("SIGNED UCELL TEST PASS\n")
