if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]])) {
 x<-suppressWarnings(Sys.setlocale("LC_CTYPE","English_United States.utf8"))
 if(!nzchar(x))stop("Use a UTF-8 Windows R environment")
}
#!/usr/bin/env Rscript
local({
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
  stopifnot(length(entry)==1L)
  if(!grepl("^([A-Za-z]:|/)",entry))entry<-file.path(getwd(),entry)
  root<-dirname(entry)
  source(file.path(root,"R/functions/project_paths.R"),local=environment())
  root<-project_root(entry);activate_existing_project_library(root)
  source(file.path(root,"R/functions/publication_bootstrap.R"),local=environment())
  load_publication_support(root, environment())
  args<-commandArgs(TRUE)
  assert_true(length(args)%in%c(1L,2L),"Usage: Rscript --vanilla export_publication.R RUN_ID [all|figures|tables]")
  source(file.path(root,"tests/release_access_guard.R"),local=environment());rc3_install_guard(root,args[1])
  result<-publication_export(root,args[1],if(length(args)==2L)args[2] else "all")
  quit(status=result$status)
})
