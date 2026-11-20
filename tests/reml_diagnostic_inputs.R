# Shared real-input preparation, used by preflight and diagnostic refits alike.
# Defines functions only. No fitting, correction, plotting or source-time I/O.
diagnostic_result_dir <- function(result_dirs,module) {
  if(is.null(names(result_dirs))||sum(names(result_dirs)==module)!=1L)stop("Diagnostic result module must exist exactly once: ",module,call.=FALSE)
  path<-result_dirs[[module]]
  if(!is.character(path)||length(path)!=1L||is.na(path)||!dir.exists(path))stop("Invalid diagnostic result directory for module: ",module,call.=FALSE)
  path
}
diagnostic_prepare_inputs <- function(root,run,export,read) {
  inputdir<-diagnostic_result_dir(run$result_dirs,"inputs")
  md<-read(file.path(inputdir,"tables/01_sample_manifest.csv"))
  assert_columns(md,c("sample_id","study_name","disease","BMI25_group","BMI23_group"),"Diagnostic sample manifest")
  assert_true(nrow(md)==1186L&&!anyDuplicated(md$sample_id)&&sum(md$disease=="T2D")==760L&&sum(md$disease=="control")==426L,"Diagnostic sample composition differs")
  studies<-unique(as.character(md$study_name));assert_true(length(studies)==5L,"Diagnostic cohort count differs")
  context<-read(file.path(export,"manifests/execution_context.csv"))
  assert_true(nrow(context)==1L&&identical(context$run_id,run$id),"Diagnostic export belongs to a different Run")
  forest_features<-list()
  for(panel in c("C1","C2")){
    f<-read(file.path(export,"plot_data",paste0("Fig6",panel,"_REML_representative_forest.csv")))
    assert_columns(f,c("pathway","study","pooled"),paste("Diagnostic forest",panel))
    assert_true(nrow(f)==6L&&sum(f$pooled)==1L&&length(unique(f$pathway))==1L,"Diagnostic forest must contain one feature and five cohorts plus pooled")
    forest_features[[panel]]<-unique(f$pathway)
  }
  specs<-list(
    pathways=list(BMI25_Total_T2D_vs_HC=c("T2D","control"),BMI25_NonOW_T2D_vs_NonOW_Ctrl=c("Non_Overweight_T2D","Non_Overweight_Control"),BMI25_OW_T2D_vs_OW_Ctrl=c("Overweight_T2D","Overweight_Control"),BMI25_OW_Ctrl_vs_NonOW_Ctrl=c("Overweight_Control","Non_Overweight_Control"),BMI25_OW_T2D_vs_NonOW_T2D=c("Overweight_T2D","Non_Overweight_T2D")),
    bmi23_pathways=list(BMI23_Lean_T2D_vs_Lean_Control=c("Non_Overweight_T2D","Non_Overweight_Control"),BMI23_Overweight_T2D_vs_Overweight_Control=c("Overweight_T2D","Overweight_Control"),BMI23_Overweight_vs_Lean_Controls=c("Overweight_Control","Non_Overweight_Control"),BMI23_Overweight_vs_Lean_T2D=c("Overweight_T2D","Non_Overweight_T2D")))
  jobs<-audit<-list()
  for(m in names(specs)){
    base<-diagnostic_result_dir(run$result_dirs,m)
    pooledfiles<-list.files(file.path(base,"tables"),pattern="_all[.]csv$",recursive=TRUE,full.names=TRUE)
    actual<-sub("_all[.]csv$","",basename(pooledfiles))
    assert_true(!anyDuplicated(actual)&&setequal(actual,names(specs[[m]])),paste("Diagnostic complete comparison file set differs:",m))
    model_folder<-file.path(base,"models",if(m=="pathways")"03_primary_BMI25_pathways"else"07_BMI23_pathways")
    for(file in pooledfiles){
      formal<-read(file);comparison<-sub("_all[.]csv$","",basename(file));cohorts<-list()
      assert_columns(formal,c("pathway","combined_log2FC","p_value","padj",if(m=="pathways")"stderr"else"se"),paste("Saved formal",comparison))
      assert_true(!anyNA(formal$pathway)&&!anyDuplicated(formal$pathway),paste("Duplicate/missing formal feature:",comparison))
      for(study in studies){
        f<-file.path(model_folder,paste0(comparison,"_",study,"_DESeq2.csv"))
        group_field<-if(comparison=="BMI25_Total_T2D_vs_HC")"disease"else if(m=="pathways")"BMI25_group"else"BMI23_group"
        expected<-all(specs[[m]][[comparison]]%in%as.character(md[[group_field]][md$study_name==study]))
        if(!file.exists(f)){assert_true(!expected,paste("Missing expected saved cohort input:",comparison,study));next}
        d<-read(f);study_column<-if(m=="pathways")"study"else"study_name"
        assert_columns(d,c("pathway",study_column,"log2FoldChange","lfcSE"),paste("Saved cohort",comparison,study))
        assert_true(!anyNA(d$pathway)&&!anyDuplicated(d$pathway)&&all(as.character(d[[study_column]])==study),paste("Cohort input IDs/study mismatch:",comparison,study))
        assert_true(is.numeric(d$log2FoldChange)&&is.numeric(d$lfcSE),paste("Cohort effect/SE must be numeric:",comparison,study))
        cohorts[[study]]<-d
      }
      assert_true(length(cohorts)>=2L,paste("Insufficient diagnostic cohort inputs:",comparison))
      combined<-dplyr::bind_rows(cohorts)
      paths<-if(m=="pathways")combined |> dplyr::count(pathway) |> dplyr::filter(n>1) |> dplyr::pull(pathway)else names(split(combined,combined$pathway))
      # Exactly the original module-specific validity rule; record, do not alter it.
      valid<-if(m=="pathways")!is.na(combined$lfcSE)&combined$lfcSE>0&is.finite(combined$log2FoldChange)else is.finite(combined$lfcSE)&combined$lfcSE>0&is.finite(combined$log2FoldChange)
      counts<-table(combined$pathway[valid]);ready<-names(counts)[counts>=2L]
      assert_true(all(formal$pathway%in%ready),paste("Formal saved features cannot be reconstructed from cohort inputs:",comparison))
      key<-paste(m,comparison);jobs[[key]]<-list(module=m,formal=formal,comparison=comparison,combined=combined,paths=paths)
      audit[[key]]<-data.frame(module=m,comparison=comparison,cohorts_read=length(cohorts),saved_formal_features=nrow(formal),candidate_paths=length(paths),fit_ready_features=sum(paths%in%ready),status="INPUT_PREFLIGHT_PASS_NO_REFITS")
    }
  }
  assert_true(length(jobs)==9L,"Diagnostic requires five BMI25 plus four BMI23 comparisons")
  list(jobs=jobs,forest_features=forest_features,audit=do.call(rbind,audit))
}

# Short diagnostic artifact paths avoid Windows path-length failures; identities stay in CSV.
diagnostic_fit_file <- function(out,module,comparison,feature) {
  family<-substr(digest::digest(paste(module,comparison,sep="|"),algo="sha256",serialize=FALSE),1,8)
  id<-substr(digest::digest(feature,algo="sha256",serialize=FALSE),1,20)
  file.path(out,"fits",family,paste0(id,".rds"))
}
