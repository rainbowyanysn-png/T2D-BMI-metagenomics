#!/usr/bin/env Rscript
# Provenance regression only. No prepare/archive/freeze/model/export operation.
local({
 if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
 entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));project<-normalizePath(dirname(dirname(entry)),winslash="/",mustWork=TRUE)
 source(file.path(project,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(project)
 source(file.path(project,"R/functions/publication_bootstrap.R"),local=environment());load_publication_support(project,environment())
 source(file.path(project,"tests/rc4_release_review_methods.R"),local=environment())
 tmp<-tempfile("export_binding_");dir.create(tmp);on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
 run<-list(root=tmp,id="Run_20000101_000000",dir=file.path(tmp,"output/Run_20000101_000000"));export<-file.path(tmp,"export");dir.create(file.path(export,"manifests"),recursive=TRUE)
 make<-function(id){base<-file.path(run$dir,"modules/clr",id);folder<-file.path(base,"tables/08_CLR_sensitivity");dir.create(folder,recursive=TRUE,showWarnings=FALSE);for(n in c("result.csv","metadata.csv"))utils::write.csv(data.frame(feature="fixture_species",coef=1,pval=.01,qval=.04),file.path(folder,n),row.names=FALSE);base}
 old<-make("a002");current<-make("a003")
 checkpoint<-function(base,signature="same"){
  files<-list.files(file.path(base,"tables"),recursive=TRUE,full.names=TRUE);index<-file.path(base,"checkpoint_artifacts.csv")
  utils::write.csv(data.frame(project_file=substring(files,nchar(tmp)+2L),sha256=vapply(files,sha256_file,character(1))),index,row.names=FALSE)
  yaml::write_yaml(list(status="COMPUTATION_PASS",module="clr",signature=signature,artifact_index=substring(index,nchar(run$dir)+2L),artifact_index_sha256=sha256_file(index)),file.path(base,"COMPLETED.yml"))
 }
 checkpoint(old);checkpoint(current)
 files<-list.files(file.path(old,"tables/08_CLR_sensitivity"),full.names=TRUE)
 manifest<-data.frame(input_path=substring(files,nchar(tmp)+2L),sha256=unname(vapply(files,sha256_file,character(1))))
 write_manifest<-function(x)utils::write.csv(x,file.path(export,"manifests/input_sha256.csv"),row.names=FALSE)
 write_manifest(manifest)
 z<-review_export_clr_binding(tmp,run,export,current);stopifnot(z$relation=="BYTE_IDENTICAL_COMPLETE_ATTEMPTS_WITH_IDENTICAL_SIGNATURE",nrow(z$evidence)==2L)
 stopifnot(review_export_clr_binding(tmp,run,export,old)$relation=="SAME_ATTEMPT")
 reject<-function(){e<-tryCatch({review_export_clr_binding(tmp,run,export,current);NULL},error=function(e)e);stopifnot(inherits(e,"error"))}
 checkpoint(current,"different");reject();checkpoint(current)
 utils::write.csv(data.frame(feature="fixture_species",coef=2,pval=.01,qval=.04),file.path(current,"tables/08_CLR_sensitivity/result.csv"),row.names=FALSE);checkpoint(current);reject();utils::write.csv(data.frame(feature="fixture_species",coef=1,pval=.01,qval=.04),file.path(current,"tables/08_CLR_sensitivity/result.csv"),row.names=FALSE);checkpoint(current)
 write_manifest(manifest[-1,,drop=FALSE]);reject();write_manifest(rbind(manifest,manifest[1,]));reject();bad<-manifest;bad$sha256[1]<-"wrong";write_manifest(bad);reject();write_manifest(manifest)
 cat("PASS identical attempt and fully identical rerun accepted\nPASS signature/content (same significant count, changed coefficient)/missing/duplicate/manifest-hash defects rejected\n")
 args<-commandArgs(TRUE)
 if(length(args)){
  stopifnot(length(args)==2L);actual<-review_run(project,args[1]);binding<-review_export_clr_binding(project,actual,file.path(project,"output",args[2]),actual$result_dirs[["clr"]])
  cat("PASS READ-ONLY ACTUAL BINDING: ",binding$export_attempt," -> ",basename(actual$result_dirs[["clr"]]),"; files=",nrow(binding$evidence),"; ",binding$relation,"\n",sep="")
 }
 cat("No review folder, archive, reference freeze or model generated.\n")
})
