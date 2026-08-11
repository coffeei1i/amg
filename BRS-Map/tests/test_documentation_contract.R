options(stringsAsFactors = FALSE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
path <- file.path(root, "docs", "METHODS.md")
stopifnot(file.exists(path))
text <- paste(readLines(path, warn = FALSE), collapse = "\n")
required <- c(
  "formal", "low_power", "descriptive", "formal_signature", "fallback",
  "primary", "sensitivity", "20 cells", "three independent", "two independent",
  "0.30", "0.50", "0.70", "FDR < 0.05", "0.25", "10 genes",
  "20 genes", "100 genes", "maxRank = 1,500", "k = 20", "PCs 1-30",
  "at least 100", "exact Celltype", "experimental-unit size",
  "direction reversal", "95th percentile", "primary-first",
  "pooled neuronal directional balance",
  "predicted_classes4", "UNRESOLVED_DO_NOT_GUESS"
)
stopifnot(all(vapply(required, grepl, logical(1L), x = text, fixed = TRUE)))
cat("DOCUMENTATION CONTRACT TEST PASS\n")
