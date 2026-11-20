#!/usr/bin/env Rscript
# Author-invoked archive/review workflow only; not part of analysis/export.
local({
 if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
 entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));root<-normalizePath(dirname(entry),winslash="/",mustWork=TRUE)
 source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
 source(file.path(root,"R/functions/publication_bootstrap.R"),local=environment());load_publication_support(root,environment())
 source(file.path(root,"tests/release_validation_config.R"),local=environment());source(file.path(root,"tests/release_validation_methods.R"),local=environment());source(file.path(root,"tests/rc4_release_review_methods.R"),local=environment())
 args<-commandArgs(TRUE);assert_true(length(args)>=1L,"Use prepare, freeze-clr or archive-reml plus explicit paths")
 mode<-args[1];opts<-list();i<-2L
 while(i<=length(args)){assert_true(i<length(args)&&startsWith(args[i],"--"),"Options require name/value pairs");name<-sub("^--","",args[i]);assert_true(is.null(opts[[name]]),paste("Duplicate option:",name));opts[[name]]<-args[i+1L];i<-i+2L}
 need<-function(names){assert_true(setequal(names(opts),names),paste("Required exact options:",paste(paste0("--",names),collapse=" ")))}
 if(mode=="prepare"){need(c("run","export","diagnostics","baseline-attempt","out"));review_prepare(root,opts$run,opts$export,opts$diagnostics,opts[["baseline-attempt"]],opts$out)}
 else if(mode=="freeze-clr"){need(c("review","out"));review_freeze_clr(root,opts$review,opts$out)}
 else if(mode=="archive-reml"){need(c("review","out"));review_archive_reml(root,opts$review,opts$out)}
 else stop("Unknown review operation",call.=FALSE)
})
