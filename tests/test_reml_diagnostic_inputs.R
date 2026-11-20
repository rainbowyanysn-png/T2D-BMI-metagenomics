#!/usr/bin/env Rscript
local({
 if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
 entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));root<-normalizePath(dirname(dirname(entry)),winslash="/",mustWork=TRUE)
 source(file.path(root,"tests/reml_diagnostic_inputs.R"),local=environment())
 tmp<-tempfile("diagnostic_dirs_");dir.create(tmp);on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
 stopifnot(identical(diagnostic_result_dir(c(inputs=tmp),"inputs"),tmp),identical(diagnostic_result_dir(list(inputs=tmp),"inputs"),tmp))
 for(x in list(c(other=tmp),c(inputs=tmp,inputs=tmp),c(inputs="missing_directory"),unname(tmp))){err<-tryCatch({diagnostic_result_dir(x,"inputs");NULL},error=function(e)e);stopifnot(inherits(err,"error"))}
 source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
 target<-diagnostic_fit_file(tmp,"pathways","BMI25_Total_T2D_vs_HC","a long pathway / feature name")
 stopifnot(nchar(substring(target,nchar(tmp)+2L))<50L,!identical(target,diagnostic_fit_file(tmp,"pathways","BMI25_OW_Ctrl_vs_NonOW_Ctrl","a long pathway / feature name")))
 dir.create(dirname(target),recursive=TRUE);fixture<-list(test_fixture=TRUE);saveRDS(fixture,target);stopifnot(identical(readRDS(target),fixture))
 parse(file.path(root,"diagnose_reml_warnings.R"),encoding="UTF-8")
 cat("PASS named-character-vector and named-list module paths\nPASS missing/duplicate/invalid module paths rejected\nPASS short artifact path and temporary fixture serialization\nPASS diagnostic syntax\n")
})
