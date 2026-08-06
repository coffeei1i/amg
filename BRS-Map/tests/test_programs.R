source(file.path(dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]))), "helper.R"))
source_module("contracts.R")
source_module("programs.R")

formal_de <- data.frame(
  gene = paste0("g", 1:220),
  logFC = c(seq(3, 0.26, length.out = 110),
            seq(-3, -0.26, length.out = 110)),
  FDR = 0.001, stringsAsFactors = FALSE
)
formal <- build_signed_program(
  formal_de, list(program_id = "formal", sample_evidence = "formal")
)
stopifnot(
  formal$status == "released", formal$signature_method == "formal_signature",
  length(formal$up_genes) == 100L, length(formal$down_genes) == 100L
)

weak_de <- transform(formal_de, FDR = 0.2)
fallback <- build_signed_program(
  weak_de, list(program_id = "fallback", sample_evidence = "low_power")
)
stopifnot(
  fallback$status == "released", fallback$signature_method == "fallback",
  length(fallback$up_genes) == 20L, length(fallback$down_genes) == 20L
)

small <- weak_de[c(1:9, 111:119), ]
rejected <- build_signed_program(
  small, list(program_id = "small", sample_evidence = "low_power")
)
stopifnot(rejected$status == "rejected_insufficient_genes")

cat("PROGRAM TEST PASS\n")
