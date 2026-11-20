#!/usr/bin/env Rscript
# Environment regression only. Does not call diagnose_reml_warnings or fit a model.
local({
 if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
 entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));root<-normalizePath(dirname(dirname(entry)),winslash="/",mustWork=TRUE)
 source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
 source(file.path(root,"tests/reml_diagnostic_environment.R"),local=environment())
 diagnostic_assert_raw_version("4.8-0","4.8-0","metafor")
 for(actual in c("4.8-1","4.7-0","4.8.0")){
   e<-tryCatch({diagnostic_assert_raw_version(actual,"4.8-0","metafor");NULL},error=function(e)e)
   stopifnot(inherits(e,"error"),grepl("expected DESCRIPTION Version=4.8-0",conditionMessage(e),fixed=TRUE))
 }
 for(actual in list(NULL,NA_character_,"",c("4.8-0","4.8-1")))stopifnot(inherits(tryCatch({diagnostic_assert_raw_version(actual,"4.8-0","metafor");NULL},error=function(e)e),"error"))
 v<-diagnostic_metafor_version_check(root);stopifnot(isTRUE(v$pass))
 parse(file.path(root,"diagnose_reml_warnings.R"),encoding="UTF-8")
 cat("PASS same raw version accepted; normalized display cannot cause false rejection\nPASS actual version mismatch and invalid descriptors rejected\nPASS loaded metafor namespace DESCRIPTION matches frozen manifest\nPASS diagnostic script syntax\n")
 print(v[,setdiff(names(v),"loaded_package_path"),drop=FALSE])
})
