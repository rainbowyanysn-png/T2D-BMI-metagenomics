#!/usr/bin/env Rscript
if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]])) {
 x<-suppressWarnings(Sys.setlocale("LC_CTYPE","English_United States.utf8"))
 if(!nzchar(x))stop("Use a UTF-8 Windows R environment")
}
local({
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
  if(length(entry)!=1L)stop("Use Rscript with an explicit entry script")
  if(!grepl("^([A-Za-z]:|/)",entry))entry<-file.path(getwd(),entry)
  root<-normalizePath(dirname(entry),winslash="/",mustWork=TRUE)
  source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
  source(file.path(root,"R/functions/analysis_bootstrap.R"),local=environment());load_analysis_support(root,environment())
  args<-parse_pipeline_arguments(commandArgs(TRUE))
  if(!is.null(args$resume)){
   state<-latest_pipeline_state(file.path(root,"output",args$resume))
   assert_true(!is.null(state),"Only a Run created inside this RC4 checkout can be resumed")
  }
  source(file.path(root,"tests/release_access_guard.R"),local=environment());rc3_install_guard(root,args$resume)
  Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
  result<-do.call(run_project_pipeline,c(list(root=root),args));quit(status=result$status)
})
