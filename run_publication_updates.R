if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]])) {
 x<-suppressWarnings(Sys.setlocale("LC_CTYPE","English_United States.utf8"))
 if(!nzchar(x))stop("Use a UTF-8 Windows R environment")
}
# Rscript --vanilla run_publication_updates.R RUN_ID [MODULE]
local({
  entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (!grepl("^([A-Za-z]:|/)", entry)) entry <- file.path(getwd(), entry)
  root <- dirname(entry)
  source(file.path(root, "R/functions/project_paths.R"), local = environment())
  activate_existing_project_library(root)
  for (name in c("assertions", "model_helpers", "clinical_correlations", "pooled_pathway_maaslin", "group_display", "source_resources",
                 "data_import", "diversity_models", "abundant_taxa_models", "publication_display_data", "medication_audit",
                 "metadata_audit", "checkpoints", "figure_exports", "figure_panels",
                 "pipeline_runner")) {
    source(file.path(root, "R/functions", paste0(name, ".R")), local = environment())
  }
  args <- commandArgs(TRUE)
  stopifnot(length(args) %in% 1:2, grepl("^Run_[0-9]{8}_[0-9]{6}$", args[1]))
  modules <- c("clinical_correlations", "pooled_pathways", "sample_flow", "alpha_details", "raw_beta", "medication_sensitivity", "abundant_taxa", "publication_display")
  if (length(args) == 2L) {
    stopifnot(args[2] %in% modules)
    modules <- args[2]
  }
  source(file.path(root,"tests/release_access_guard.R"),local=environment())
  rc3_install_guard(root,args[1])
  result <- run_project_pipeline(root, selected = modules, resume = args[1])
  quit(status = result$status)
})
