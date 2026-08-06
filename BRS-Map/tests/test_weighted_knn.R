source(file.path(dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]))), "helper.R"))
source_module("contracts.R")
source_module("weighted_knn.R")

make_source <- function(n, celltype = "A") {
  x <- data.frame(
    cell_id = paste0("s", seq_len(n)), Celltype = celltype,
    signed_response_score = seq(-1, 1, length.out = n),
    experimental_unit_id = rep(paste0("u", 1:5), length.out = n),
    stringsAsFactors = FALSE
  )
  for (i in 1:30) x[[paste0("PC_", i)]] <- seq_len(n) / n + i / 100
  x
}
target <- data.frame(cell_id = c("t1", "t2"), Celltype = "A")
for (i in 1:30) target[[paste0("PC_", i)]] <- c(0.2, 0.8) + i / 100

fail <- run_round1_one_program(make_source(99), target, "p", "A")
pass <- run_round1_one_program(make_source(100), target, "p", "A")
stopifnot(all(fail$status == "ROUND1_SOURCE_LT_100"))
stopifnot(all(pass$status == "PASS"), all(pass$n_neighbors == 20L))
stopifnot(all(is.finite(pass$effective_n)), all(pass$max_unit_weight <= 1))

round2_source <- transform(
  make_source(100), target_id = paste0("r", 1:100), status = "PASS",
  transferred_score = signed_response_score, reference_unit = experimental_unit_id
)
stereo <- data.frame(target_id = c("b1", "b2"), Celltype = "A")
for (i in 1:30) stereo[[paste0("PC_", i)]] <- c(0.3, 0.7) + i / 100
r2 <- run_round2_one_program(round2_source, stereo, "p", "A")
stopifnot(all(r2$status == "PASS"), all(r2$n_neighbors == 20L))
stopifnot(all(r2$spatial_identity_column == "predicted_classes4"))

cat("WEIGHTED KNN TEST PASS\n")
