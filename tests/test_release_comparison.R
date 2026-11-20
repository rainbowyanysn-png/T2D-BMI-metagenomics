#!/usr/bin/env Rscript
local({
 entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));root<-dirname(dirname(entry))
 source(file.path(root,"R/functions/assertions.R"),local=environment())
 source(file.path(root,"tests/release_validation_methods.R"),local=environment())
 a<-data.frame(feature=c("a","b"),coef=c(.2,-.4),pval=c(.001,.3),formal_q=c(.01,.5))
 folder<-tempfile("release_fixture_");dir.create(folder)
 release_compare_frames(a,a,"equal",folder,1e-8)
 fails<-function(b,label){dir<-file.path(folder,label);dir.create(dir);inherits(tryCatch({release_compare_frames(a,b,label,dir,1e-8);NULL},error=function(e)e),"error")}
 b<-a;b$coef[1]<--b$coef[1];stopifnot(fails(b,"direction"))
 b<-a[2:1,];stopifnot(fails(b,"ID_order"))
 b<-a;b$pval[1]<-.1;stopifnot(fails(b,"P_difference"))
 x<-data.frame(feature="a",coef=1,pval=.01,formal_q=.05-1e-10);y<-x;y$formal_q<-.05+1e-10
 dir<-file.path(folder,"boundary");dir.create(dir)
 stopifnot(inherits(tryCatch({release_compare_frames(x,y,"threshold",dir,1e-8);NULL},error=function(e)e),"error"))
 cat("PASS comparison fixtures: numeric, ID order, direction and exact significance boundary; no models fitted\n")
})
