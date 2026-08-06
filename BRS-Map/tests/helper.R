test_module_root <- function() {
  frame_files <- vapply(sys.frames(), function(x) {
    if (!is.null(x$ofile)) as.character(x$ofile) else ""
  }, character(1))
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files)) {
    return(normalizePath(file.path(dirname(tail(frame_files, 1L)), ".."),
                         mustWork = TRUE))
  }
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  normalizePath(file.path(dirname(sub("^--file=", "", script_arg[[1L]])), ".."),
                mustWork = TRUE)
}

source_module <- function(file) {
  source(file.path(test_module_root(), "R", file), local = FALSE)
}
