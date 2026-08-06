#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

parse_cli <- function(args) {
  result <- list(stage = "preflight", module_root = NULL, input_root = NULL,
                 output_root = NULL, adapter = NULL, resume = FALSE)
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--resume") {
      result$resume <- TRUE
      i <- i + 1L
    } else {
      if (i == length(args)) stop("Missing value for ", key, call. = FALSE)
      value <- args[[i + 1L]]
      name <- switch(key, "--stage" = "stage", "--module-root" = "module_root",
                     "--input-root" = "input_root", "--output-root" = "output_root",
                     "--adapter" = "adapter", stop("Unknown argument: ", key, call. = FALSE))
      result[[name]] <- value
      i <- i + 2L
    }
  }
  result
}

find_module_root <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."),
                mustWork = TRUE)
}

run_brs_map <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  cli <- parse_cli(arguments)
  module_root <- if (is.null(cli$module_root)) find_module_root() else
    normalizePath(cli$module_root, mustWork = TRUE)
  for (file in list.files(file.path(module_root, "R"), "[.]R$", full.names = TRUE)) source(file)
  source(file.path(module_root, "scripts", "preflight.R"))
  if (!cli$stage %in% brs_valid_stages()) stop("Unsupported stage", call. = FALSE)
  preflight <- run_brs_preflight(module_root, require_runtime_packages = cli$stage != "preflight")
  if (cli$stage == "preflight") {
    message("BRS-MAP PREFLIGHT PASS: datasets=", preflight$n_datasets,
            " contrasts=", preflight$n_contrasts)
    return(invisible(preflight))
  }
  if (is.null(cli$input_root) || is.null(cli$output_root) || is.null(cli$adapter)) {
    stop("Formal stages require --input-root, --output-root, and --adapter", call. = FALSE)
  }
  input_root <- normalizePath(cli$input_root, mustWork = TRUE)
  output_root <- normalizePath(cli$output_root, mustWork = FALSE)
  adapter <- normalizePath(cli$adapter, mustWork = TRUE)
  assert_output_destination(output_root, cli$resume)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  adapter_env <- new.env(parent = globalenv())
  sys.source(adapter, envir = adapter_env)
  requested <- if (cli$stage == "all") brs_valid_stages()[2:8] else cli$stage
  produced <- character()
  for (stage in requested) {
    function_name <- paste0("brs_stage_", stage)
    if (!exists(function_name, envir = adapter_env, mode = "function", inherits = FALSE)) {
      stop("Adapter lacks function ", function_name, call. = FALSE)
    }
    message("[BRS-Map] ", stage)
    value <- get(function_name, envir = adapter_env)(
      input_root = input_root, output_root = output_root,
      module_root = module_root, resume = cli$resume
    )
    produced <- c(produced, as.character(value))
  }
  manifest <- build_release_manifest(
    stage = cli$stage,
    input_files = c(preflight$config$paths, adapter),
    output_files = produced[file.exists(produced)],
    parameters = as.list(setNames(preflight$config$parameters$value,
                                  preflight$config$parameters$parameter))
  )
  write_release_manifest(manifest, file.path(output_root, "BRS_MAP_COMPLETE.json"))
  invisible(manifest)
}

if (sys.nframe() == 0L) run_brs_map()
