#!/usr/bin/env Rscript
# Pure display-string test only; no rendering or statistical calculation.
local({
 if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
 entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));root<-normalizePath(dirname(dirname(entry)),winslash="/",mustWork=TRUE)
 source(file.path(root,"R/functions/assertions.R"),local=environment());source(file.path(root,"R/functions/publication_refinement_adapter.R"),local=environment())
 x<-data.frame(comparison="Test comparison",label="n = 10\nConcordant = 80.0%\nr = 0.56\nShared FDR n = 2");original<-x
 z<-final_sensitivity_strip_labels(x);stopifnot(identical(x,original),length(strsplit(z,"\n")[[1]])==3L)
 for(line in strsplit(x$label,"\n")[[1]])stopifnot(grepl(line,z,fixed=TRUE))
 x$label<-"incomplete";err<-tryCatch({final_sensitivity_strip_labels(x);NULL},error=function(e)e);stopifnot(inherits(err,"error"))
 parse(file.path(root,"replace_sfig4c_annotation.R"),encoding="UTF-8")
 cat("PASS saved annotation strings preserved in outside-panel strip labels; invalid format rejected; replacement syntax valid.\n")
})
