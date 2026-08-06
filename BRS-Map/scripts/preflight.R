options(stringsAsFactors = FALSE)

brs_valid_stages <- function() {
  c("preflight", "pseudobulk", "de", "programs", "ucell",
    "round1", "round2", "audit", "all")
}

read_brs_config <- function(module_root) {
  config_root <- file.path(module_root, "config")
  paths <- file.path(config_root, c("parameters.tsv", "datasets.tsv", "contrasts.tsv"))
  if (any(!file.exists(paths))) stop("BRS-Map configuration is incomplete", call. = FALSE)
  list(
    parameters = read.delim(paths[[1L]], check.names = FALSE),
    datasets = read.delim(paths[[2L]], check.names = FALSE),
    contrasts = read.delim(paths[[3L]], check.names = FALSE),
    paths = paths
  )
}

run_brs_preflight <- function(module_root, require_runtime_packages = TRUE) {
  config <- read_brs_config(module_root)
  stopifnot(
    sum(config$datasets$status == "retained") == 7L,
    nrow(config$contrasts) == 16L,
    all(config$contrasts$accession %in%
          config$datasets$accession[config$datasets$status == "retained"])
  )
  required <- c("Matrix", "data.table", "edgeR", "UCell", "RANN",
                "jsonlite", "digest")
  if (isTRUE(require_runtime_packages)) {
    missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
    if (length(missing)) {
      stop("Missing R package(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
  }
  list(status = "PASS", n_datasets = 7L, n_contrasts = 16L,
       packages = required, config = config)
}
