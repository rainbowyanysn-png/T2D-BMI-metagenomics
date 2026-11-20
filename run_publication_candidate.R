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
  contract<-utils::read.csv(file.path(root,"config/statistical_implementation_sha256.csv"))
  assert_true(all(vapply(file.path(root,contract$file),sha256_file,character(1))==contract$sha256),"Statistical implementation hash differs")
  assert_true(!length(commandArgs(TRUE)),"Clean candidate always starts a NEW Run")
  source(file.path(root,"tests/release_access_guard.R"),local=environment());rc3_install_guard(root)
  Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
  result<-run_project_pipeline(root,"publication_candidate")
  if(result$status!=0L)quit(status=result$status)
  # run_project_pipeline persists the access audit in its invocation directory.
  safe_write_csv(get("v3_file_access_report",.GlobalEnv)(),file.path(result$run$dir,"file_access_guard.csv"))
  exe<-file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe"else"Rscript")
  status<-system2(exe,c("--vanilla",shQuote(file.path(root,"export_publication.R")),result$run$id,"all"))
  quit(status=status)
})
