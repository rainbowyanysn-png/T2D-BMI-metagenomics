# Read-only saved-evidence review. No fitted model, optimizer or statistic is changed.
release_reml_saved_review <- function(root,run,out,warnings,export=NULL) {
  dir<-file.path(out,"reml_warning_review");dir.create(dir,recursive=TRUE,showWarnings=FALSE)
  x<-warnings[grepl("Fisher scoring|local maximum|Setting tau",warnings$warning,ignore.case=TRUE),,drop=FALSE]
  x$warning_id<-seq_len(nrow(x));x$feature<-NA_character_;x$comparison<-NA_character_
  x$k<-x$coef<-x$stderr<-x$pval<-x$padj<-x$tau2<-x$I2<-NA_real_
  x$formal_significant<-x$entered_representative_forest<-NA
  x$discussion_status<-"NOT_ASSESSED_MANUSCRIPT_NOT_PROVIDED"
  x$localization_status<-"UNRESOLVED: aggregate warning log has no feature/comparison; no REML fitted object saved"
  safe_write_csv(x,file.path(dir,"unlocalized_warning_inventory.csv"))
  universe<-list();inventory<-list()
  for(m in intersect(c("pathways","bmi23_pathways"),names(run$result_dirs))){
    folder<-run$result_dirs[[m]];files<-list.files(file.path(folder,"tables"),pattern="_all[.]csv$",recursive=TRUE,full.names=TRUE)
    inventory[[m]]<-data.frame(module=m,saved_fitted_REML_objects=length(list.files(folder,pattern="[.]rds$",recursive=TRUE)),tau2_saved=m=="bmi23_pathways",warning_context="module and message only; no pathway/contrast")
    for(f in files){
      d<-utils::read.csv(f,check.names=FALSE,stringsAsFactors=FALSE);comparison<-sub("_all[.]csv$","",basename(f))
      z<-data.frame(module=m,comparison=comparison,feature=d$pathway,k=if("n_studies"%in%names(d))d$n_studies else d$k,
        coef=d$combined_log2FC,stderr=if("stderr"%in%names(d))d$stderr else d$se,pval=d$p_value,padj=d$padj,
        tau2=if("tau2"%in%names(d))d$tau2 else NA_real_,I2=if("I2"%in%names(d))d$I2 else d$Heterogeneity_I2,
        formal_significant=d$padj<.05,warning_link="NOT_IDENTIFIABLE_FROM_SAVED_LOGS",entered_representative_forest=NA,
        discussion_status="NOT_ASSESSED_MANUSCRIPT_NOT_PROVIDED")
      if(!is.null(export)){
        z$entered_representative_forest<-FALSE
        for(panel in c("C1","C2")){
          relevant<-m=="pathways"&&comparison==if(panel=="C1")"BMI25_Total_T2D_vs_HC"else"BMI25_NonOW_T2D_vs_NonOW_Ctrl"
          if(relevant){g<-file.path(export,"plot_data",paste0("Fig6",panel,"_REML_representative_forest.csv"));assert_file(g);forest<-utils::read.csv(g,check.names=FALSE);z$entered_representative_forest<-z$feature%in%forest$pathway}
        }
      }
      universe[[length(universe)+1L]]<-z
    }
  }
  safe_write_csv(do.call(rbind,universe),file.path(dir,"saved_pathway_result_universe_NOT_warning_attribution.csv"))
  safe_write_csv(do.call(rbind,inventory),file.path(dir,"saved_evidence_inventory.csv"))
  safe_write_lines(c("UNRESOLVED REML NUMERICAL WARNINGS: DO NOT TREAT AS HARMLESS.",
    "The module-level warning messages do not identify contrast/pathway. BMI25 saved pooled tables omit tau2; neither module saved REML fitted objects.",
    "Saved result universe contains actual k/effect/SE/P/padj/I2 and tau2 only where saved. These rows are NOT claimed to be the warning features. Zero I2/tau2 alone cannot attribute a warning.",
    "Forest membership is readable, but whether one of the 31 warnings affected those features remains unknown. Manuscript discussion was not supplied.",
    "Author may run diagnose_reml_warnings.R --run RUN_ID --export EXPORT_ID. It refits REML only in a new diagnostic directory, never substitutes formal estimates.",
    "Serialized cohort CSV effect/SE precision can differ from the original in-memory objects. Diagnostic warnings/counts cannot retrospectively prove exact identity with the historical warnings."),file.path(dir,"README_UNRESOLVED.txt"))
  invisible(nrow(x))
}
