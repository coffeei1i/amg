script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
source(file.path(dirname(script_path), "helper.R"))
source_module("contracts.R")
source_module("final_readout.R")

round1 <- data.frame(
  program_id = c("same", "flip", "zero"),
  Celltype = "A",
  status = "PASS",
  transferred_score = c(0.4, 0.3, 0),
  stringsAsFactors = FALSE
)
round2 <- data.frame(
  program_id = c("same", "flip", "zero"),
  Celltype = "A",
  status = "PASS",
  transferred_score = c(0.2, -0.1, 0.5),
  stringsAsFactors = FALSE
)
direction <- audit_direction_consistency(round1, round2)
status <- setNames(direction$direction_status, direction$program_id)
stopifnot(
  status[["same"]] == "SAME_DIRECTION",
  status[["flip"]] == "DIRECTION_REVERSAL",
  status[["zero"]] == "AMBIGUOUS_ZERO",
  direction$eligible_for_final[direction$program_id == "same"],
  !direction$eligible_for_final[direction$program_id == "flip"]
)

primary <- data.frame(
  program_id = c("p1", "p2"), Celltype = "A",
  status = c("PASS", "ROUND1_SOURCE_LT_100"),
  direction_status = c("SAME_DIRECTION", "NOT_EVALUATED"),
  signature_method = c("formal_signature", "fallback"),
  sample_evidence = c("formal", "formal"), stringsAsFactors = FALSE
)
sensitivity <- data.frame(
  program_id = c("p1", "p2"), Celltype = "A", status = "PASS",
  direction_status = "SAME_DIRECTION",
  signature_method = c("formal_signature", "fallback"),
  sample_evidence = "formal", stringsAsFactors = FALSE
)
selected <- select_final_identity(primary, sensitivity)
stopifnot(
  selected$identity_source[selected$program_id == "p1"] == "primary",
  selected$identity_source[selected$program_id == "p2"] == "sensitivity",
  selected$point_marker[selected$program_id == "p1"] == "",
  selected$point_marker[selected$program_id == "p2"] == "S*"
)

cell_scores <- data.frame(
  accession = "X", contrast_id = "c", Celltype = "A",
  cell_id = paste0("c", 1:4), transferred_score = c(-2, 1, 4, 0),
  stringsAsFactors = FALSE
)
landscape <- summarise_celltype_landscape(cell_scores)
stopifnot(
  landscape$scale_denominator == 4,
  isTRUE(all.equal(landscape$mean_scaled_score, 0.1875)),
  landscape$valid_score_fraction == 1
)

spatial <- data.frame(
  program_id = "p", slice_id = c("M1", "M1", "M2", "M3"),
  target_id = paste0("t", 1:4), transferred_score = c(-2, 1, 4, 100),
  stringsAsFactors = FALSE
)
joint <- scale_spatial_joint_q95(spatial)
stopifnot(
  length(unique(joint$scale_denominator)) == 1L,
  unique(joint$scaling_scope) == "PROGRAM_ACROSS_ALL_SLICES",
  max(joint$scaled_score, na.rm = TRUE) == 1,
  min(joint$scaled_score, na.rm = TRUE) >= -1
)

subnuclear <- data.frame(
  accession = "X", contrast_id = "c", subnucleus = "BA",
  target_id = paste0("n", 1:5), transferred_score = c(1, 2, -1, -2, 0),
  stringsAsFactors = FALSE
)
balance <- summarise_subnuclear_directional_balance(subnuclear)
stopifnot(
  balance$positive_n == 2L, balance$negative_n == 2L,
  balance$zero_n == 1L, balance$directional_balance == 0
)

cat("FINAL READOUT TEST PASS\n")
