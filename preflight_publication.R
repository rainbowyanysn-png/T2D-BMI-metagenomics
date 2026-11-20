#!/usr/bin/env Rscript
local({
  if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
  root<-normalizePath(dirname(entry),winslash="/",mustWork=TRUE)
  source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
  source(file.path(root,"R/functions/publication_bootstrap.R"),local=environment());load_publication_support(root,environment())
  args<-commandArgs(TRUE);assert_true(length(args)%in%1:2,"Use preflight_publication.R RUN_ID [all|tables|figures]")
  source(file.path(root,"tests/release_access_guard.R"),local=environment());rc3_install_guard(root,args[1])
  x<-publication_preflight(root,args[1],if(length(args)==2L)args[2] else "all")
  quit(status=if(x$pass)0L else 1L)
})
