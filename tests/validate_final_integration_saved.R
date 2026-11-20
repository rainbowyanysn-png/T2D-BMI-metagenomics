# Saved current-run export acceptance only; never used by a model or feature filter.
validate_final_integration_saved <- function(root,run,export) {
  pub_clr_alignment_contract(run)
  source(file.path(root,"tests/release_validation_config.R"),local=environment())
  read<-function(relative)utils::read.csv(file.path(export,relative),check.names=FALSE,stringsAsFactors=FALSE)
  sheets<-function(id)unique(read(paste0("tables/csv/",id,"__data_dictionary.csv"))$sheet)
  tables<-read("table_index.csv");submission_validate_table_identity(tables)
  submission_validate_saved(root,run,export)
  assert_true(identical(sheets("Supplementary_Table_S3"),c("Cohort_PERMANOVA","Marginal_PERMANOVA","PERMDISP","Interactions","BMI23","Descriptive_ComBat","Raw_Aitchison","Notes")),"S3 eight-sheet contract differs")
  s4<-sheets("Supplementary_Table_S4")
  required<-c("LOCO_primary_significant","LOCO_feature_summary","LOCO_all_results","LOCO_model_summary","LOCO_exclusions","LOCO_data_dictionary")
  assert_true(all(required%in%s4)&&"All_five_comparisons"%in%s4,"S4 original/LOCO sheets missing")
  counts<-c(LOCO_primary_significant=120L,LOCO_feature_summary=79L,LOCO_all_results=5000L,LOCO_model_summary=25L)
  for(n in names(counts))assert_true(nrow(read(paste0("tables/csv/Supplementary_Table_S4__",n,".csv")))==counts[[n]],paste("S4 row count differs:",n))
  pairwise<-read("tables/csv/Supplementary_Table_S2__Pairwise.csv")
  assert_columns(pairwise,"metric","S2 pairwise metric identity")
  assert_true(nrow(pairwise)==12L&&setequal(pairwise$metric,c("Shannon","Richness"))&&all(table(pairwise$metric)==6L)&&!anyDuplicated(pairwise[c("metric","group1","group2")]),"S2 pairwise metrics missing or ambiguous")
  exported_dispersion<-read("tables/csv/Supplementary_Table_S3__PERMDISP.csv")
  actual<-final_beta_sources(run)$dispersion
  assert_true(identical(exported_dispersion$comparison,actual$comparison)&&isTRUE(all.equal(exported_dispersion$F,actual$F,tolerance=1e-12))&&all(exported_dispersion[["Pr(>F)"]]==.001),"S3 dispersion differs from the current Run")
  index<-read("figure_index.csv");config<-release_publication_config(root)
  pub_validate_panel_counts(index,config)
  source(file.path(root,"tests/validate_sensitivity_panel_saved.R"),local=environment())
  validate_sensitivity_panel_saved(root,run,export)
  b<-publication_data(run,config)
  for(panel in c("C1","C2")) {
    id<-paste0("Fig6",panel,"_REML_representative_forest");z<-final_pathway_forest(b,panel)
    values<-read(paste0("plot_data/",id,".csv"));ranking<-read(paste0("plot_data/",id,"_candidate_ranking.csv"))
    assert_true(nrow(values)==6L&&sum(values$pooled)==1L&&all(values$pathway==z$values$pathway[1]),"Forest displayed feature/cohort count differs")
    assert_true(identical(ranking$pathway[ranking$selected],z$ranking$pathway[z$ranking$selected]),"Forest ranking selection changed")
  }
  heat<-index[index$figure_id=="Fig6B_REML_pathway_effect_heatmap",]
  assert_true(nrow(heat)==1L&&heat$width_mm==140&&heat$height_mm==180,"Latest Fig6B proportion missing")
  invisible(TRUE)
}
