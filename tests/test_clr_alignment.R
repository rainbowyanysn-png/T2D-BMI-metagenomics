#!/usr/bin/env Rscript
# Model specifications/designs and export metadata only. No statistical fitting.
local({
  if(.Platform$OS.type=="windows")Sys.setlocale("LC_CTYPE","English_United States.utf8")
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));stopifnot(length(entry)==1L)
  root<-normalizePath(dirname(dirname(entry)),winslash="/",mustWork=TRUE)
  support<-Sys.getenv("T2DBMI_TEST_SUPPORT_ROOT",unset=root)
  source(file.path(support,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(support)
  source(file.path(support,"R/functions/assertions.R"),local=environment())
  source(file.path(support,"R/functions/medication_audit.R"),local=environment())
  source(file.path(root,"R/pipeline/08_CLR_sensitivity.R"),local=environment())
  check<-function(n,f){f();message("PASS ",n)}
  failure<-function(f)inherits(tryCatch({f();NULL},error=function(e)e),"error")
  studies<-c("KarlssonFH_2013","MetaCardis_2020_a","QinJ_2012","SankaranarayananK_2015","YuJ_2015")
  bmi<-c(rep(22,261),rep(24,46),rep(27,119),rep(22,105),rep(24,33),rep(27,622))
  n<-length(bmi);md<-data.frame(sample_id=paste0("fixture_",seq_len(n)),study_name=rep(studies,length.out=n),disease=c(rep("control",426),rep("T2D",760)),BMI=bmi,age=40+seq_len(n)%%43,sex=ifelse(seq_len(n)%%3==0,"male","female"),treatment=rep(c("metformin;insulin","insulin","dppiv",NA),length.out=n))
  md$study_name[md$BMI==24]<-"QinJ_2012";md$BMI23_threshold<-ifelse(md$study_name%in%c("QinJ_2012","YuJ_2015"),23,25)
  group<-function(cut)ifelse(md$disease=="T2D",ifelse(md$BMI<cut,"Non_Overweight_T2D","Overweight_T2D"),ifelse(md$BMI<cut,"Non_Overweight_Control","Overweight_Control"))
  md$analysis_group<-group(25);md$BMI23_group<-group(md$BMI23_threshold)
  check("cohort-adapted thresholds and primary group identities",function(){z<-clr_prepared_metadata(md);stopifnot(identical(z$BMI23_binary,as.integer(md$BMI>=md$BMI23_threshold)),sum(z$T2D_binary&!z$BMI23_binary)==105L)})
  check("uniform-23 or inconsistent source labels fail",function(){bad<-md;bad$BMI23_threshold<-23;stopifnot(failure(function()clr_prepared_metadata(bad)));bad<-md;bad$BMI23_group[1]<-"Overweight_Control";stopifnot(failure(function()clr_prepared_metadata(bad)))})
  check("T2D species covariates/reference match primary drug mapping",function(){z<-clr_prepared_metadata(md);s<-clr_pair_model_spec("Overweight_vs_Lean_T2D","Taxonomy",z,c("Non_Overweight_T2D","Overweight_T2D"));stopifnot(s$medication_adjusted,nrow(s$metadata)==760L,"antidiabetes_group"%in%all.vars(as.formula(s$formula)),s$medication_reference==medication_model_reference(medication_model_covariate(md[md$disease=="T2D",])));d<-model.matrix(as.formula(s$formula),s$metadata);stopifnot(nrow(d)==760L,qr(d)$rank==ncol(d))})
  check("pathways and non-internal species do not gain drug adjustment",function(){z<-clr_prepared_metadata(md);for(s in list(clr_pair_model_spec("Overweight_vs_Lean_T2D","Pathway",z,c("Non_Overweight_T2D","Overweight_T2D")),clr_pair_model_spec("Lean_T2D_vs_Lean_Control","Taxonomy",z,c("Non_Overweight_Control","Non_Overweight_T2D"))))stopifnot(!s$medication_adjusted,!"antidiabetes_group"%in%all.vars(as.formula(s$formula)))})
  check("combination therapy and unknown records preserve original classification",function(){z<-clr_prepared_metadata(md);stopifnot(all(as.character(z$antidiabetes_group[md$treatment%in%"metformin;insulin"])=="Metformin"),all(as.character(z$antidiabetes_group[is.na(md$treatment)&md$disease=="T2D"])=="Other_Anti-diabetic_Drugs"))})
  args<-commandArgs(TRUE)
  if(length(args)){
    stopifnot(length(args)==2L,args[1]=="--metadata-run",grepl("^Run_[0-9]{8}_[0-9]{6}$",args[2]))
    files<-list.files(file.path(support,"output",args[2],"modules","inputs"),pattern="^01_sample_manifest[.]csv$",recursive=TRUE,full.names=TRUE);stopifnot(length(files)==1L)
    actual<-read.csv(files,check.names=FALSE);z<-clr_prepared_metadata(actual)
    check("real saved sample metadata and full-rank adjusted CLR design",function(){s<-clr_pair_model_spec("Overweight_vs_Lean_T2D","Taxonomy",z,c("Non_Overweight_T2D","Overweight_T2D"));d<-model.matrix(as.formula(s$formula),s$metadata);stopifnot(nrow(z)==1186L,nrow(d)==760L,qr(d)$rank==ncol(d),identical(as.integer(table(s$metadata$Group_binary)),c(138L,622L)))})
  }
  message("PASS CLR alignment design checks. No model fitted; no saved result replaced.")
})
