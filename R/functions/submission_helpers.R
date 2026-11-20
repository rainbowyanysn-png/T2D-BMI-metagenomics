# Submission views only. Frozen statistics, sources and plot implementations are immutable.
submission_table_manifest <- function() data.frame(
  content_id=c("participant_characteristics","analysis_summary",paste0("supplement_",1:13),"confounder_availability"),
  table_id=c("Main_Table_1","Analysis_Summary_Table",paste0("Supplementary_Table_S",1:14)),
  previous_table_id=c("Main_Table_1","Analysis_Summary_Table",paste0("Supplementary_Table_S",1:13),"Supplementary_Table_S15"),
  stringsAsFactors=FALSE)

submission_source_view <- function(x) {
  before<-x
  rows<-x$variable%in%c("medication","metformin")
  text<-", with omission of this covariate assessed in sensitivity analyses"
  x$reason_not_adjusted[rows]<-gsub(text,"",x$reason_not_adjusted[rows],fixed=TRUE)
  evidence<-setdiff(names(x),"reason_not_adjusted")
  assert_true(identical(x[evidence],before[evidence]),"Source evidence changed in submission view")
  assert_true(!any(grepl("omission of this covariate assessed",x$reason_not_adjusted[rows],fixed=TRUE)),"Obsolete medication scope remains")
  x
}

submission_pathway_scope <- function(root,d) {
  env<-new.env(parent=baseenv());env$r3_need<-assert_true
  source(file.path(root,"R/functions/refinement/pathway_pure_helpers.R"),local=env)
  list(audit=env$r3_candidate_filter_audit(d),selected=env$r3_expected_candidate_data(d))
}

submission_pathway_view <- function(root,d) {
  original<-d;scope<-submission_pathway_scope(root,d);selected<-scope$selected$feature
  assert_true(!anyDuplicated(d$feature),"Duplicate display feature IDs")
  flag<-d$feature%in%selected
  for(n in c("selected_for_genus_panel","selected_for_species_panel","selected_for_panel"))if(n%in%names(d))d[[n]]<-flag
  d$displayed_in_Fig5B<-d$displayed_in_Fig5C<-d$displayed_in_SFig3<-flag
  d$display_selection_rank<-match(d$feature,selected)
  d$display_selection_rule<-"All eligible species-stratified contributions: saved visualization_candidate; full-model qval <0.01; abs(coef)>1; named genus/species/pathway; exclude UNINTEGRATED/UNMAPPED; no top-N truncation"
  d$display_exclusion_reason<-vapply(seq_len(nrow(d)),function(i){
    if(flag[i])return("")
    predicates<-setdiff(names(scope$audit),c("feature","coef","qval","included"))
    paste(predicates[!as.logical(scope$audit[i,predicates])],collapse="; ")
  },character(1))
  science<-setdiff(names(original),c("selected_for_genus_panel","selected_for_species_panel","selected_for_panel"))
  assert_true(identical(d[science],original[science]),"Submission display view changed source scientific columns")
  d
}

submission_validate_s7 <- function(root,full,s7) {
  expected<-submission_pathway_scope(root,full)$selected$feature
  assert_true(!anyDuplicated(s7$feature),"S7 duplicate feature")
  assert_true(all(s7$feature%in%full$feature),"S7 unknown feature")
  for(n in c("displayed_in_Fig5B","displayed_in_Fig5C","displayed_in_SFig3","selected_for_genus_panel","selected_for_species_panel")){
    assert_columns(s7,n,"S7 display identity")
    assert_true(!anyNA(s7[[n]])&&setequal(s7$feature[as.logical(s7[[n]])],expected),paste("S7 display set differs:",n))
  }
  original<-full[match(s7$feature,full$feature),,drop=FALSE]
  for(n in intersect(c("coef","stderr","pval","qval","N","N.not.zero"),names(s7)))
    assert_true(isTRUE(all.equal(unname(s7[[n]]),unname(original[[n]]),tolerance=1e-12)),paste("S7 scientific column changed:",n))
  assert_true(identical(s7$qval<.05,original$qval<.05)&&identical(s7$qval<.01,original$qval<.01)&&identical(sign(s7$coef),sign(original$coef)),"S7 significance state or coefficient direction changed")
  invisible(TRUE)
}

submission_validate_table_identity <- function(index,sheets=NULL) {
  expected<-submission_table_manifest()$table_id
  assert_true(nrow(index)==length(expected)&&!anyDuplicated(index$table_id)&&setequal(index$table_id,expected),"Submission table identity set must be two main tables plus S1-S14")
  if("content_id"%in%names(index))assert_true(identical(as.character(index$content_id),submission_table_manifest()$content_id[match(index$table_id,expected)]),"Submission content identity does not match table number")
  if(!is.null(sheets)) {
    assert_true(identical(names(sheets),c("Source_and_database","Database_completeness","Notes")),"S14 must contain clinical/confounder availability sheets, not medication sensitivity")
    assert_columns(sheets$Source_and_database,c("cohort","variable","source_study_collected","evidence_grade"),"S14 source evidence")
    assert_columns(sheets$Database_completeness,c("variable_id"),"S14 database completeness")
    text<-paste(unlist(sheets),collapse=" ")
    assert_true(grepl("Categorical antidiabetic treatment was included",text,fixed=TRUE)&&
      grepl("reliably defined metformin-free analysis was not available",text,fixed=TRUE)&&
      !grepl("with omission of this covariate assessed in sensitivity analyses",text,fixed=TRUE)&&
      grepl("Metformin was not a separate binary covariate",text,fixed=TRUE)&&grepl("not verified nonuse",text,fixed=TRUE)&&
      grepl("Pathway meta-analysis and formal interaction models were not medication-adjusted",text,fixed=TRUE),"S14 medication scope text differs")
  }
  invisible(TRUE)
}

submission_validate_saved <- function(root,run,export,figure_export=export) {
  read<-function(p)utils::read.csv(file.path(export,p),check.names=FALSE,stringsAsFactors=FALSE)
  index<-read("table_index.csv");assert_columns(index,"content_id","Submission identity index");submission_validate_table_identity(index)
  files<-list.files(file.path(export,"tables/supplementary"),pattern="[.]xlsx$")
  assert_true(setequal(files,paste0("Supplementary_Table_S",1:14,".xlsx")),"Unexpected or retired workbook in formal directory")
  formal_csv<-list.files(file.path(export,"tables/csv"))
  assert_true(!any(grepl("^Supplementary_Table_S15|^Supplementary_Table_S14__(Models|Paired_coefficients|Concordance)[.]csv$",formal_csv)),"Retired supplement CSV remains in formal submission directory")
  dict<-read("tables/csv/Supplementary_Table_S14__data_dictionary.csv")
  sheets<-setNames(lapply(unique(dict$sheet),function(s)read(paste0("tables/csv/Supplementary_Table_S14__",s,".csv"))),unique(dict$sheet))
  submission_validate_table_identity(index,sheets)
  config<-read_config(root);config$project_root<-root;config$publication<-yaml::read_yaml(file.path(root,"config/publication_panel_config.yml"))
  b<-publication_data(run,config);s7<-read("tables/csv/Supplementary_Table_S7__Candidates_and_selection.csv")
  submission_validate_s7(root,b$candidates,s7)
  expected<-submission_pathway_scope(root,b$candidates)$selected$feature
  reports<-list()
  for(id in c("Fig5B_MaAsLin2_genus_pathway_features","Fig5C_MaAsLin2_species_pathway_features","SFig3_taxon_pathway_assignment_Sankey")) {
    f<-file.path(figure_export,"plot_data",paste0(id,"_display_candidate_selection_audit.csv"));assert_file(f,"Panel source-selection audit")
    x<-utils::read.csv(f,check.names=FALSE)
    assert_true(!anyDuplicated(x$feature)&&setequal(x$feature[x$included],expected),paste("S7/panel source identity mismatch:",id))
    reports[[id]]<-data.frame(panel=id,n_source_records=length(expected),same_feature_set=TRUE)
  }
  invisible(pub_bind(reports))
}
