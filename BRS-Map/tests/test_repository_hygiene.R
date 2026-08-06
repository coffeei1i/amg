options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
module_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

files <- list.files(
  module_root, recursive = TRUE, full.names = TRUE,
  all.files = TRUE, no.. = TRUE
)
stopifnot(!any(grepl("\\.(rds|h5ad|h5|loom)$", files, ignore.case = TRUE)))
stopifnot(!any(grepl("(^|[/\\\\])output([/\\\\]|$)", files, ignore.case = TRUE)))

text_files <- files[grepl(
  "\\.(R|md|tsv|txt|yml|yaml|json)$", files, ignore.case = TRUE
)]
contents <- unlist(lapply(text_files, readLines, warn = FALSE), use.names = FALSE)
unix_data_pattern <- paste0("(^|[^A-Za-z])", "/", "data[0-9]*/")
windows_data_pattern <- paste0("Z:", "/", "data")
stopifnot(!any(grepl(
  paste(unix_data_pattern, windows_data_pattern, sep = "|"),
  contents
)))
stopifnot(!any(grepl("(token|password|secret)[[:space:]]*=", contents,
                     ignore.case = TRUE)))

cat("REPOSITORY HYGIENE TEST PASS\n")
