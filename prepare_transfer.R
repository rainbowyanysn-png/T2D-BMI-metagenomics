#!/usr/bin/env Rscript
# Author-invoked candidate handoff only, not production publication approval.
# Copies the public payload with empty output/data_derived; never runs analysis/export.
local({
  if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));stopifnot(length(entry)==1L)
  root<-normalizePath(dirname(entry),winslash="/",mustWork=TRUE)
  source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
  stopifnot(requireNamespace("yaml",quietly=TRUE),requireNamespace("digest",quietly=TRUE))
  source(file.path(root,"R/functions/public_release_files.R"),local=environment())
  args<-commandArgs(TRUE);zip_requested<-FALSE;destination<-NULL;i<-1L
  while(i<=length(args)) {
    if(args[i]=="--zip"){stopifnot(!zip_requested);zip_requested<-TRUE;i<-i+1L}
    else if(args[i]=="--dest"&&i<length(args)){stopifnot(is.null(destination));destination<-args[i+1L];i<-i+2L}
    else stop("Use prepare_transfer.R [--zip] [--dest NEW_DIRECTORY]",call.=FALSE)
  }
  if(zip_requested&&!requireNamespace("zip",quietly=TRUE))stop("Restore renv.lock first; zip is not installed automatically.",call.=FALSE)
  stamp<-paste0(format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid())
  if(is.null(destination))destination<-file.path(dirname(root),paste0("T2D_BMI_RC4_CodeCandidate_",stamp))
  destination<-normalizePath(destination,winslash="/",mustWork=FALSE)
  if(file.exists(destination)||dir.exists(destination))stop("Transfer destination already exists; nothing overwritten.",call.=FALSE)
  canonical_root<-tolower(normalizePath(root,winslash="/",mustWork=TRUE))
  if(startsWith(tolower(destination),paste0(canonical_root,"/")))stop("Transfer destination must be outside this working project.",call.=FALSE)
  zipfile<-paste0(destination,".zip")
  if(zip_requested&&file.exists(zipfile))stop("ZIP already exists; nothing overwritten.",call.=FALSE)
  x<-public_release_selection(root)
  expected<-utils::read.csv(file.path(root,"public_file_manifest.csv"),stringsAsFactors=FALSE)
  if(!identical(x$file,expected$file)||!identical(x$sha256,expected$sha256))stop("Public manifest is stale. Run build_public_package.R --manifest first.",call.=FALSE)
  stopifnot(dir.create(destination,recursive=TRUE))
  copied<-logical(nrow(x))
  for(j in seq_len(nrow(x))) {
    src<-file.path(root,x$file[j]);dst<-file.path(destination,x$file[j]);dir.create(dirname(dst),recursive=TRUE,showWarnings=FALSE)
    stopifnot(file.copy(src,dst,overwrite=FALSE));copied[j]<-identical(digest::digest(file=dst,algo="sha256"),x$sha256[j]);stopifnot(copied[j])
    stopifnot(identical(digest::digest(file=src,algo="sha256"),x$sha256[j]))
  }
  # Metadata copied alongside the payload does not recursively hash itself.
  utils::write.csv(x,file.path(destination,"public_file_manifest.csv"),row.names=FALSE,fileEncoding="UTF-8")
  utils::write.csv(x[c("file","sha256")],file.path(destination,"public_file_sha256.csv"),row.names=FALSE,fileEncoding="UTF-8")
  for(n in c("output","data_derived")){dir.create(file.path(destination,n));file.create(file.path(destination,n,".gitkeep"))}
  allowed<-c(x$file,"public_file_manifest.csv","public_file_sha256.csv","output/.gitkeep","data_derived/.gitkeep")
  actual<-list.files(destination,recursive=TRUE,all.files=TRUE,no..=TRUE)
  stopifnot(setequal(actual,allowed))
  if(zip_requested) {
    zip::zipr(zipfile=zipfile,files=basename(destination),root=dirname(destination),recurse=TRUE,include_directories=TRUE,mode="mirror")
    expected_members<-paste(basename(destination),allowed,sep="/")
    contents<-zip::zip_list(zipfile)$filename
    stopifnot(setequal(contents[!endsWith(contents,"/")],expected_members))
  }
  audit_root<-Sys.getenv("T2DBMI_AUDIT_DIR",unset=file.path(dirname(root),paste0(basename(root),"_Audit")))
  dir.create(audit_root,recursive=TRUE,showWarnings=FALSE)
  receipt<-data.frame(file=x$file,expected_sha256=x$sha256,copy_verified=copied)
  utils::write.csv(receipt,file.path(audit_root,paste0("Transfer_",stamp,".csv")),row.names=FALSE)
  cat("CANDIDATE HANDOFF ONLY; production admission unchanged.\nFolder: ",destination,"\nPayload files: ",nrow(x),"\n",sep="")
  if(zip_requested)cat("ZIP: ",zipfile,"\nZIP SHA-256: ",digest::digest(file=zipfile,algo="sha256"),"\n",sep="")
  cat("No existing outputs, caches or audit folders were copied. No statistics or figures were run.\n")
})
