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
absolute_path_pattern <- paste0(
  "(^|[[:space:]='\"`(])(",
  "[A-Za-z]:[/\\\\]|",
  "/(data[0-9]*|home|Users|mnt|your)/",
  ")"
)
stopifnot(!any(grepl(absolute_path_pattern, contents, perl = TRUE)))
markdown_files <- text_files[grepl("[.]md$", text_files, ignore.case = TRUE)]
for (markdown_path in markdown_files) {
  markdown <- paste(readLines(markdown_path, warn = FALSE), collapse = "\n")
  matches <- regmatches(markdown, gregexpr("!?\\[[^]]*\\]\\(([^)]+)\\)", markdown, perl = TRUE))[[1L]]
  if (!length(matches) || identical(matches, character(0))) next
  targets <- sub("^.*\\(([^)]+)\\)$", "\\1", matches)
  local_targets <- targets[!grepl("^(https?://|#|mailto:)", targets)]
  stopifnot(!any(grepl("^(/|[A-Za-z]:[/\\\\]|file:)", local_targets, perl = TRUE)))
  for (target in local_targets) {
    target_path <- file.path(dirname(markdown_path), URLdecode(target))
    stopifnot(file.exists(target_path))
  }
}
stopifnot(!any(grepl("(token|password|secret)[[:space:]]*=", contents,
                     ignore.case = TRUE)))

cat("REPOSITORY HYGIENE TEST PASS\n")
