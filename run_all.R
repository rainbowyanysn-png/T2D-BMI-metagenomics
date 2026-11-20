#!/usr/bin/env Rscript
if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]])) {
 x<-suppressWarnings(Sys.setlocale("LC_CTYPE","English_United States.utf8"))
 if(!nzchar(x))stop("Use a UTF-8 Windows R environment")
}
local({
 entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));stopifnot(length(entry)==1L)
 root<-normalizePath(dirname(entry),winslash="/",mustWork=TRUE)
 source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
 args<-commandArgs(TRUE);admission<-yaml::read_yaml(file.path(root,"config/module_admission.yml"))
 required<-c("inputs",utils::read.csv(file.path(root,"config/pipeline_modules.csv"))$module)
 ready<-isTRUE(admission$production_ready) && identical(admission$clean_rc4_run_pending,FALSE) && all(required%in%unlist(admission$admitted_modules)) && all(vapply(c("clean_run_pass","independent_reference_comparison_pass","manual_review_complete","citation_and_license_complete"),function(k)isTRUE(admission[[k]]),logical(1)))
 candidate<-identical(args,"--candidate")
 if(!candidate && !(length(args)==0L && ready))
  stop("Release admission is closed. Use --candidate for author-run clean acceptance; see RELEASE_CHECKLIST.md",call.=FALSE)
 exe<-file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe"else"Rscript")
 code<-system2(exe,c("--vanilla",shQuote(file.path(root,"run_publication_candidate.R"))))
 quit(status=code)
})
