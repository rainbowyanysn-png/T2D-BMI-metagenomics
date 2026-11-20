#!/usr/bin/env Rscript
local({
  if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));stopifnot(length(entry)==1L)
  root<-normalizePath(dirname(entry),winslash="/",mustWork=TRUE)
  source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
  stopifnot(requireNamespace("yaml",quietly=TRUE),requireNamespace("digest",quietly=TRUE))
  source(file.path(root,"R/functions/public_release_files.R"),local=environment())
  args<-commandArgs(TRUE)
  if(identical(args,"--manifest")) {
    x<-public_release_manifests(root);cat("Public payload files: ",nrow(x),"\nManifest: ",file.path(root,"public_file_manifest.csv"),"\n",sep="")
  } else if(length(args)==2L&&args[1]=="--copy") {
    public_release_copy(root,args[2]);cat("Copied final public package: ",args[2],"\n",sep="")
  } else stop("Use --manifest, or --copy NEW_DESTINATION after release admission. Neither option runs statistics.",call.=FALSE)
})
