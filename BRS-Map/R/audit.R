options(stringsAsFactors = FALSE)

sha256_file <- function(path) {
  if (!file.exists(path)) stop("Missing audit input: ", path, call. = FALSE)
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(file = path, algo = "sha256"))
  }
  value <- tools::md5sum(path)
  unname(value[[1L]])
}

file_record <- function(path) {
  info <- file.info(path)
  list(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    bytes = unname(info$size),
    sha256 = sha256_file(path)
  )
}

build_release_manifest <- function(stage, input_files, output_files,
                                   parameters, status = "COMPLETE") {
  list(
    schema_version = "1.0",
    status = status,
    stage = as.character(stage),
    generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    parameters = parameters,
    inputs = lapply(as.character(input_files), file_record),
    outputs = lapply(as.character(output_files), file_record),
    r_version = R.version.string
  )
}

write_release_manifest <- function(manifest, path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package jsonlite is required to write the audit manifest", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(manifest, path, auto_unbox = TRUE, pretty = TRUE,
                       null = "null", na = "null")
  invisible(path)
}

assert_output_destination <- function(path, resume = FALSE) {
  if (dir.exists(path)) {
    contents <- list.files(path, all.files = TRUE, no.. = TRUE)
    if (length(contents) && !isTRUE(resume)) {
      stop("Output destination is not empty; use --resume: ", path,
           call. = FALSE)
    }
  }
  TRUE
}
