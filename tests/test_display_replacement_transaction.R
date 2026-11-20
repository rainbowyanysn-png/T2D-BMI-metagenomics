#!/usr/bin/env Rscript
# Real file I/O tests in temporary fixtures only. Never renders or touches project outputs.
local({
 if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
 entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));root<-normalizePath(dirname(dirname(entry)),winslash="/",mustWork=TRUE)
 source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
 stopifnot(requireNamespace("digest",quietly=TRUE));sha256_file<-function(f)digest::digest(file=f,algo="sha256")
 source(file.path(root,"tests/display_replacement_transaction.R"),local=environment())
 sandbox<-tempfile("replacement_tests_");dir.create(sandbox);on.exit(unlink(sandbox,recursive=TRUE),add=TRUE)
 fixture<-function(id){
  d<-file.path(sandbox,id);for(n in c("old","stage","backup"))dir.create(file.path(d,n),recursive=TRUE)
  names<-c("panel.pdf","panel.tiff","panel.png","manifest.csv");dst<-file.path(d,"old",names);src<-file.path(d,"stage",names)
  for(i in seq_along(names)){writeBin(as.raw(c(0,1,i,255)),dst[i]);writeBin(as.raw(c(10,20,i,0,200)),src[i])}
  protected<-file.path(d,"unrelated_statistic.csv");writeLines("unchanged",protected)
  list(dst=dst,src=src,targets=setNames(src,dst),backup=file.path(d,"backup"),old_hashes=unname(vapply(dst,sha256_file,character(1))),protected=setNames(sha256_file(protected),protected))
 }
 original<-function(x)identical(unname(vapply(x$dst,sha256_file,character(1))),x$old_hashes)
 failure<-function(expr){e<-tryCatch({force(expr);NULL},error=function(e)e);stopifnot(inherits(e,"error"));e}
 x<-fixture("success");receipt<-replacement_commit(x$targets,x$backup,x$protected)
 stopifnot(identical(unname(vapply(x$dst,sha256_file,character(1))),unname(vapply(x$src,sha256_file,character(1)))),identical(receipt$before_sha256,x$old_hashes),all(file.exists(receipt$backup)))
 cat("PASS named-hash backup verification and successful image/metadata fixture commit\n")
 x<-fixture("backup_failure")
 copy<-function(from,to,overwrite=FALSE){if(startsWith(basename(to),"02_"))return(FALSE);file.copy(from,to,overwrite=overwrite)}
 e<-failure(replacement_commit(x$targets,x$backup,x$protected,copy_file=copy));stopifnot(original(x),grepl("Backup copy failed",conditionMessage(e),fixed=TRUE))
 cat("PASS backup failure leaves every target unchanged\n")
 x<-fixture("commit_failure");failed<-FALSE
 copy<-function(from,to,overwrite=FALSE){if(overwrite&&identical(to,x$dst[[2]])&&!failed){failed<<-TRUE;return(FALSE)};file.copy(from,to,overwrite=overwrite)}
 e<-failure(replacement_commit(x$targets,x$backup,x$protected,copy_file=copy));stopifnot(original(x),grepl("All touched targets restored",conditionMessage(e),fixed=TRUE))
 cat("PASS mid-commit failure rolls back all touched targets and verifies original hashes\n")
 x<-fixture("locked_unchanged_target")
 copy<-function(from,to,overwrite=FALSE){if(overwrite&&identical(to,x$dst[[2]]))return(FALSE);file.copy(from,to,overwrite=overwrite)}
 e<-failure(replacement_commit(x$targets,x$backup,x$protected,copy_file=copy));stopifnot(original(x),grepl("All touched targets restored",conditionMessage(e),fixed=TRUE))
 cat("PASS unchanged locked destination recognized by original hash during rollback\n")
 x<-fixture("hash_mismatch");damaged<-FALSE
 copy<-function(from,to,overwrite=FALSE){ok<-file.copy(from,to,overwrite=overwrite);if(overwrite&&identical(to,x$dst[[3]])&&!damaged){damaged<<-TRUE;writeBin(as.raw(99),to)};ok}
 e<-failure(replacement_commit(x$targets,x$backup,x$protected,copy_file=copy));stopifnot(original(x),grepl("byte-hash mismatch",conditionMessage(e),fixed=TRUE))
 cat("PASS post-copy hash mismatch is detected and rolled back\n")
 x<-fixture("protected_change");changed<-FALSE
 copy<-function(from,to,overwrite=FALSE){ok<-file.copy(from,to,overwrite=overwrite);if(overwrite&&!changed){changed<<-TRUE;writeLines("external fixture edit",names(x$protected))};ok}
 e<-failure(replacement_commit(x$targets,x$backup,x$protected,copy_file=copy));stopifnot(original(x),grepl("protected file hashes changed",conditionMessage(e),fixed=TRUE))
 cat("PASS unrelated-file mutation detected; owned replacement targets rolled back\n")
 parse(file.path(root,"replace_sfig4c_annotation.R"),encoding="UTF-8")
 cat("PASS replacement entry syntax. No project image/table/statistic changed.\n")
})
