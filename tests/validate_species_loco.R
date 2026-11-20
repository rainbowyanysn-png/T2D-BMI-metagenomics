validate_species_loco <- function(root,run) {
  state<-latest_pipeline_state(run$dir)
  assert_true("species_loco"%in%names(state$completed),"LOCO checkpoint missing")
  verify_module_checkpoint(run,checkpoint_record(run,state$completed$species_loco))
  directory<-file.path(run$result_dirs[["species_loco"]],"tables/02a_species_LOCO")
  read<-function(name)utils::read.csv(file.path(directory,name),check.names=FALSE,stringsAsFactors=FALSE)
  x<-read("LOCO_all_results.csv");m<-read("LOCO_model_summary.csv");summary<-read("LOCO_primary_significant.csv");features<-read("LOCO_feature_summary.csv")
  assert_true(nrow(x)==5000L&&!anyDuplicated(x[c("model_id","feature")])&&all(table(x$model_id)==200L),"Complete 200x25 result grid required")
  assert_true(nrow(m)==25L&&all(m$status=="PASS")&&all(table(m$comparison_id)==5L)&&all(table(m$omitted_study)==5L),"25 LOCO fits required")
  assert_true(all(m$full_rank)&&all(m$design_rank==m$design_columns)&&all(m$omitted_present==0L)&&all(m$n_cohorts==4L),"LOCO cohort/design mismatch")
  assert_true(all(m$medication_covariate==(m$comparison_id=="OT2D_vs_NT2D_BMI25"))&&all(m$medication_levels[m$medication_covariate]>1L),"Medication-adjusted LOCO identity changed")
  assert_true(all(m$standardize)&all(m$normalization=="TSS")&all(m$transform=="LOG")&all(m$min_prevalence==.10),"LOCO parameters changed")
  tested<-x$tested_after_prevalence
  assert_true(all(is.na(as.matrix(x[!tested,c("coef","stderr","pval","qval","qval_target_bh")]))),"Excluded estimates must remain missing")
  assert_true(all(is.finite(x$coef[tested]))&&all(x$qval[tested]>=0 & x$qval[tested]<=1),"Invalid fitted estimates/q")
  assert_true(isTRUE(all.equal(summary,species_loco_summarise_all(read("LOCO_primary_significant_input.csv"),x),tolerance=1e-12,check.attributes=FALSE)),"LOCO summary cannot be reconstructed")
  species_loco_acceptance(summary,features)
  for(id in m$model_id) {
    file<-file.path(run$result_dirs[["species_loco"]],"models",id,"all_results.tsv")
    fit<-utils::read.delim(file,check.names=FALSE)
    assert_true(isTRUE(all.equal(fit$qval,p.adjust(fit$pval,"BH"),tolerance=1e-12,check.attributes=FALSE)),"MaAsLin2 full-model BH mismatch")
    z<-x[x$model_id==id & tested,,drop=FALSE]
    assert_true(isTRUE(all.equal(z$qval_target_bh,p.adjust(z$pval,"BH"),tolerance=1e-12,check.attributes=FALSE)),"Audit target BH mismatch")
  }
  invisible(TRUE)
}
