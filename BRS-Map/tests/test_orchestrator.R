options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

source(file.path(root, "R", "audit.R"))
source(file.path(root, "scripts", "preflight.R"))

stages <- brs_valid_stages()
stopifnot(identical(
  stages,
  c("preflight", "pseudobulk", "de", "programs", "ucell",
    "round1", "round2", "audit", "all")
))

temporary <- tempfile("brs-map-audit-")
dir.create(temporary)
writeLines("x", file.path(temporary, "input.txt"))
manifest <- build_release_manifest(
  stage = "round2",
  input_files = file.path(temporary, "input.txt"),
  output_files = character(),
  parameters = list(k = 20L, minimum = 100L)
)
stopifnot(
  manifest$schema_version == "1.0",
  manifest$status == "COMPLETE",
  manifest$stage == "round2",
  nzchar(manifest$inputs[[1L]]$sha256),
  manifest$parameters$k == 20L
)

occupied <- file.path(temporary, "occupied")
dir.create(occupied)
writeLines("keep", file.path(occupied, "existing.txt"))
stopifnot(inherits(try(assert_output_destination(occupied, resume = FALSE),
                       silent = TRUE), "try-error"))
stopifnot(isTRUE(assert_output_destination(occupied, resume = TRUE)))

scripts <- c(
  file.path(root, "scripts", "preflight.R"),
  file.path(root, "scripts", "run_brs_map.R")
)
text <- unlist(lapply(scripts, readLines, warn = FALSE), use.names = FALSE)
stopifnot(!any(grepl(paste0("/", "data[0-9]*/|Z:", "/data"), text)))

cat("ORCHESTRATOR TEST PASS\n")
