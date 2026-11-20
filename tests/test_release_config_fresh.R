#!/usr/bin/env Rscript
local({
 if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
 entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));root<-normalizePath(dirname(dirname(entry)),winslash="/",mustWork=TRUE)
 source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
 source(file.path(root,"R/functions/publication_bootstrap.R"),local=environment());load_publication_support(root,environment())
 source(file.path(root,"tests/release_validation_config.R"),local=environment())
 # Isolate the sensitivity validator in a fresh R process and inspect its config
 # at the exact publication_data boundary; do not fit or draw anything.
 publication_data<-function(run,config){stopifnot(config$publication$pathway_candidate_limit$genus==20,config$publication$pathway_candidate_limit$species==30);stop("CONFIG_COMPLETE_AT_PUBLICATION_DATA",call.=FALSE)}
 source(file.path(root,"tests/validate_sensitivity_panel_saved.R"),local=environment())
 err<-tryCatch({validate_sensitivity_panel_saved(root,list(),tempdir());NULL},error=function(e)e)
 stopifnot(inherits(err,"error"),identical(conditionMessage(err),"CONFIG_COMPLETE_AT_PUBLICATION_DATA"))
 cat("PASS fresh-process complete publication config before publication_data\n")
})
