#!/usr/bin/env Rscript
local({
 if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
 entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));root<-normalizePath(dirname(dirname(entry)),winslash="/",mustWork=TRUE)
 source(file.path(root,"R/functions/project_paths.R"));activate_existing_project_library(root)
 source(file.path(root,"R/functions/publication_bootstrap.R"));load_publication_support(root,environment())
 reject<-function(f){e<-tryCatch({f();NULL},error=function(e)e);stopifnot(inherits(e,"error"))}
 ids<-submission_table_manifest()$table_id;index<-data.frame(table_id=ids)
 submission_validate_table_identity(index)
 reject(function()submission_validate_table_identity(rbind(index,data.frame(table_id="Supplementary_Table_S15"))))
 bad_identity<-data.frame(table_id=ids,content_id=submission_table_manifest()$content_id);bad_identity$content_id[16]<-"medication_adjustment_sensitivity";reject(function()submission_validate_table_identity(bad_identity))
 wrong<-list(Models=data.frame(x=1),Paired_coefficients=data.frame(x=1),Notes=data.frame(x=1))
 reject(function()submission_validate_table_identity(index,wrong))
 d<-data.frame(feature=paste0("feature",1:5945),coef=2,qval=.001,visualization_candidate=TRUE,
  taxonomic_level="genus",genus="Excluded",taxonomic_contributor="Excluded",base_pathway_ID="PWY",pathway_description="Description",
  pval=.00001,stderr=.1,selected_for_genus_panel=FALSE,selected_for_species_panel=FALSE)
 species<-rep(c("Bacteroides fragilis","Escherichia coli","Eubacterium eligens","Eubacterium eligens CAG 72","Faecalibacterium prausnitzii","Bacteroides fragilis CAG 47"),c(19,12,12,8,7,6))
 d$taxonomic_level[1:64]<-"species";d$taxonomic_contributor[1:64]<-species
 d$genus[1:64]<-sub(" .*","",species);d$base_pathway_ID[1:64]<-paste0("PWY",rep(1:42,length.out=64));d$coef[38:64]<- -2
 view<-submission_pathway_view(root,d);submission_validate_s7(root,d,view)
 stopifnot(sum(view$displayed_in_Fig5B)==64L,identical(view$coef,d$coef),identical(view$pval,d$pval),identical(view$qval,d$qval))
 bad<-view;bad$selected_for_genus_panel<-seq_len(nrow(bad))<=20;reject(function()submission_validate_s7(root,d,bad))
 bad<-view;bad$selected_for_species_panel<-seq_len(nrow(bad))<=30;reject(function()submission_validate_s7(root,d,bad))
 bad<-view[-1,];reject(function()submission_validate_s7(root,d,bad))
 bad<-rbind(view,view[1,]);reject(function()submission_validate_s7(root,d,bad))
 bad<-view;bad$feature[1]<-"UNKNOWN";reject(function()submission_validate_s7(root,d,bad))
 bad<-view;bad$coef[1]<-3;reject(function()submission_validate_s7(root,d,bad))
 z<-data.frame(variable=c("medication","metformin","diet"),reason_not_adjusted=c(rep("Categorical antidiabetic treatment was included in the T2D-internal species models, with omission of this covariate assessed in sensitivity analyses. A reliably defined metformin-free analysis was not available.",2),"Unchanged evidence"),evidence_grade="A")
 v<-submission_source_view(z);stopifnot(identical(v$evidence_grade,z$evidence_grade),identical(v$reason_not_adjusted[3],z$reason_not_adjusted[3]))
 cat("PASS16 table identities; retired/misnamed medication workbooks rejected\nPASS64 full-ID display membership; stale top30/missing/duplicate/mismatched features or changed numbers rejected\nPASS scoped source text mapping preserves evidence\n")
})
