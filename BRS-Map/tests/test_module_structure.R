options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
module_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
repository_root <- dirname(module_root)

required <- c(
  file.path(module_root, "README.md"),
  file.path(module_root, "R", "final_readout.R"),
  file.path(module_root, "workflow", "BRS-Map_workflow.svg"),
  file.path(module_root, "workflow", "BRS-Map_workflow.drawio"),
  file.path(repository_root, ".gitignore")
)
stopifnot(all(file.exists(required)))

cat("MODULE STRUCTURE TEST PASS\n")
