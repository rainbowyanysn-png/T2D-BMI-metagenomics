#!/usr/bin/env Rscript
# Saved text, source definitions and public inventory checks only; no model/export/device.
local({
  if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));stopifnot(length(entry)==1L)
  root<-normalizePath(dirname(dirname(entry)),winslash="/",mustWork=TRUE)
  source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
  source(file.path(root,"R/functions/publication_bootstrap.R"),local=environment());load_publication_support(root,environment())
  source(file.path(root,"R/functions/public_release_files.R"),local=environment())
  config<-read_config(root);config$publication<-yaml::read_yaml(file.path(root,"config/publication_panel_config.yml"))
  z<-pub_source_review(root,config,"confounder");clean<-pub_reader_table(z)
  drug<-z$variable%in%c("medication","metformin")
  stopifnot(nrow(z)==75L,sum(drug)==10L,identical(z$reason_not_adjusted[drug],clean$reason_not_adjusted[drug]),
    identical(z$use_in_current_model,clean$use_in_current_model),all(z$use_in_current_model%in%c("yes","no")))
  recorded<-z$cohort%in%c("KarlssonFH_2013","MetaCardis_2020_a","QinJ_2012")
  stopifnot(all(z$use_in_current_model[z$variable=="medication"&recorded]=="yes"),
    all(z$use_in_current_model[z$variable=="medication"&!recorded]=="no"),all(z$use_in_current_model[z$variable=="metformin"]=="no"))
  # True source limitation is not removed by text cleaning.
  i<-z$cohort=="QinJ_2012"&z$variable=="metformin"
  stopifnot(grepl("cannot be assumed metformin-free",z$source_verification_note[i],fixed=TRUE),
    grepl("cannot be assumed metformin-free",clean$source_verification_note[i],fixed=TRUE))
  x<-public_release_selection(root)
  stopifnot(!"preview_publication.R"%in%x$file,!file.exists(file.path(root,"preview_publication.R")),
    !any(grepl("(^|/)(audit|output|data_derived)/|(^|/)(preview|resize_|replace_)",x$file)),
    all(c("run_all.R","run_publication_candidate.R","run_all_development.R","export_publication.R","validate_release.R")%in%x$file))
  # No selected source dependency tree has been dropped.
  for(tree in c("R","config","tests","data_raw")) {
    all<-file.path(tree,list.files(file.path(root,tree),recursive=TRUE,full.names=FALSE,all.files=TRUE,no..=TRUE))
    stopifnot(all(all%in%x$file))
  }
  saved<-utils::read.csv(file.path(root,"public_file_manifest.csv"),stringsAsFactors=FALSE)
  stopifnot(identical(x$file,saved$file),identical(x$sha256,saved$sha256),identical(x$bytes,as.numeric(saved$bytes)))
  contract<-utils::read.csv(file.path(root,"config/statistical_implementation_sha256.csv"))
  stopifnot(all(vapply(file.path(root,contract$file),function(f)digest::digest(file=f,algo="sha256"),character(1))==contract$sha256))
  allr<-unique(c(list.files(file.path(root,"R"),pattern="[.]R$",recursive=TRUE,full.names=TRUE),list.files(file.path(root,"tests"),pattern="[.]R$",recursive=TRUE,full.names=TRUE),list.files(root,pattern="[.]R$",full.names=TRUE)))
  for(f in allr)parse(f,encoding="UTF-8")
  cat("PASS S15 schema, scoped model-use flags, CSV -> source reader -> text cleaning\nPASS True metformin-free source limitation retained\nPASS Development preview absent from public candidate\nPASS Five formal entries and complete R/config/tests/data-source trees\nPASS Public inventory hashes\nPASS Statistical implementation hashes\nPASS All R syntax\nPublic payload files: ",nrow(x),"\n",sep="")
})
