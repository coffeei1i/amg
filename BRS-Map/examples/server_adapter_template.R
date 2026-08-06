# Copy this file outside the Git repository and edit only the I/O adapters.
# The BRS-Map mathematical functions are already sourced by run_brs_map.R.

brs_stage_pseudobulk <- function(input_root, output_root, module_root, resume) {
  stop("Implement project-specific object loading, then call aggregate_pseudobulk_counts()")
}

brs_stage_de <- function(input_root, output_root, module_root, resume) {
  stop("Load pseudobulk checkpoints, then call run_edger_contrast() or run_descriptive_contrast()")
}

brs_stage_programs <- function(input_root, output_root, module_root, resume) {
  stop("Load DE tables, then call build_signed_program()")
}

brs_stage_ucell <- function(input_root, output_root, module_root, resume) {
  stop("Load an expression object and programs, then call score_one_object_programs()")
}

brs_stage_round1 <- function(input_root, output_root, module_root, resume) {
  stop("Load external and snRNA PCA tables, then call run_round1_one_program()")
}

brs_stage_round2 <- function(input_root, output_root, module_root, resume) {
  stop("Load Round-1 and Stereo PCA tables, then call run_round2_one_program()")
}

brs_stage_audit <- function(input_root, output_root, module_root, resume) {
  character()
}
