# BMI25 species LOCO: all features are fitted before the prespecified report is extracted.
species_loco_specs <- function() list(
  Total_T2D_vs_Control = list(reference="control", comparison="T2D", variable="disease", column="disease", fixed=c("disease","BMI","age","sex","study_name")),
  NT2D_vs_NC_BMI25 = list(reference="Non_Overweight_Control", comparison="Non_Overweight_T2D", variable="Group", column="BMI25_group", fixed=c("Group","age","sex","study_name")),
  OT2D_vs_OC_BMI25 = list(reference="Overweight_Control", comparison="Overweight_T2D", variable="Group", column="BMI25_group", fixed=c("Group","age","sex","study_name")),
  OC_vs_NC_BMI25 = list(reference="Non_Overweight_Control", comparison="Overweight_Control", variable="Group", column="BMI25_group", fixed=c("Group","age","sex","study_name")),
  OT2D_vs_NT2D_BMI25 = list(reference="Non_Overweight_T2D", comparison="Overweight_T2D", variable="Group", column="BMI25_group", fixed=c("Group","age","sex","study_name","antidiabetes_group")))

species_loco_metadata <- function(data, spec, omitted) {
  m <- data$metadata
  m <- m[as.character(m$study_name)!=omitted & as.character(m[[spec$column]]) %in% c(spec$reference,spec$comparison),,drop=FALSE]
  m[[spec$variable]] <- factor(as.character(m[[spec$column]]),levels=c(spec$reference,spec$comparison))
  drug <- "antidiabetes_group" %in% spec$fixed
  if(drug) m <- medication_model_covariate(m)
  m <- droplevels(as.data.frame(m))
  assert_true(!anyNA(m[,c("sample_id",spec$fixed),drop=FALSE]),"Missing LOCO sample/covariate value")
  assert_true(!anyDuplicated(m$sample_id) && all(table(m[[spec$variable]])>0),"Missing comparison group or duplicate sample")
  assert_true(!omitted %in% as.character(m$study_name),"Omitted cohort is still present")
  refs <- c(paste(spec$variable,spec$reference,sep=","),categorical_references(m,spec$fixed))
  if(drug) {
    assert_true(nlevels(m$antidiabetes_group)>1L,"Medication is not estimable; covariate must not be dropped")
    refs <- c(refs,paste("antidiabetes_group",medication_model_reference(m),sep=","))
  }
  design_data <- m
  for(ref in refs) {
    pair <- strsplit(ref,",",fixed=TRUE)[[1]]
    design_data[[pair[1]]] <- stats::relevel(factor(design_data[[pair[1]]]),ref=pair[2])
  }
  design <- stats::model.matrix(stats::reformulate(spec$fixed),design_data)
  assert_true(nrow(design)==nrow(m) && qr(design)$rank==ncol(design),"LOCO design is rank deficient; no automatic covariate removal")
  list(metadata=m,references=refs,design=design,medication=drug,
       medication_reference=if(drug) medication_model_reference(m) else NA_character_)
}

species_loco_primary_gate <- function(run, expected) {
  assert_true("primary" %in% names(run$result_dirs),"Same-run full primary results are required before LOCO")
  p <- run$result_dirs[["primary"]]
  assert_true(startsWith(p,paste0(run$base_dir,"/modules/primary/")),"Primary results must come from this new run")
  summary <- utils::read.csv(file.path(p,"tables/02_primary_BMI25_species_summary.csv"))
  ids <- names(species_loco_specs())
  assert_true(identical(summary$comparison_id,ids),"Primary comparison IDs/order changed")
  assert_true(identical(as.integer(summary$n_tested),c(200L,192L,187L,193L,190L)),"Full-data primary tested counts differ")
  assert_true(identical(as.integer(summary$n_formal_q_lt_005),c(69L,32L,10L,4L,5L)),"Full-data primary significant counts differ; LOCO is prohibited")
  fits <- setNames(lapply(ids,function(id) {
    x <- utils::read.csv(file.path(p,"tables/02_primary_BMI25_species",paste0(id,"_all.csv")),check.names=FALSE)
    s <- species_loco_specs()[[id]]
    assert_true(!anyDuplicated(x$feature) && all(x$metadata==s$variable) && all(x$value==s$comparison),"Primary target/reference identity changed")
    assert_true(all(x$reference_group==s$reference) && all(x$comparison_group==s$comparison),"Primary effect direction changed")
    assert_true(all(x$formal_q==x$qval) && all(x$q_family=="MaAsLin2 full model qval"),"Primary formal q definition changed")
    assert_true(all(x$medication_covariate==("antidiabetes_group" %in% s$fixed)),"Primary medication identity changed")
    x
  }),ids)
  fits
}

species_loco_primary_significant <- function(primary) {
  do.call(rbind,lapply(names(primary),function(id) {
    x<-primary[[id]];x<-x[x$qval<.05,,drop=FALSE]
    data.frame(comparison_id=id,feature=x$feature,
      display_name=gsub("[.]"," ",sub("^species[.]","",x$feature)),full_coef=x$coef,full_qval=x$qval,stringsAsFactors=FALSE)
  }))
}

species_loco_summarise_all <- function(significant,results) {
  out<-lapply(seq_len(nrow(significant)),function(i) {
    f<-significant[i,,drop=FALSE]
    x<-results[results$comparison_id==f$comparison_id & results$feature==f$feature,,drop=FALSE]
    assert_true(nrow(x)==5L&&!anyDuplicated(x$omitted_study),"Five unique LOCO omissions required")
    tested<-x$tested_after_prevalence;v<-x$coef[tested];n<-sum(tested)
    same<-sum(sign(v)==sign(f$full_coef));delta<-abs(x$coef-f$full_coef)
    cbind(f,data.frame(n_LOCO_expected=5L,n_LOCO_testable=n,n_direction_concordant=same,
      direction_concordance=if(n)same/n else NA_real_,n_LOCO_q_lt_0.05=sum(x$qval[tested]<.05),
      maximum_absolute_coefficient_change=if(n)max(delta[tested])else NA_real_,
      largest_change_omitted_study=if(n)x$omitted_study[which.max(replace(delta,!tested,-Inf))]else NA_character_,
      any_direction_reversal=if(n)any(v*f$full_coef<0)else NA,
      all_testable_directions_concordant=n>0L&&same==n,
      prevalence_excluded_count=sum(!tested & x$exclusion_reason=="prevalence",na.rm=TRUE),
      variance_excluded_count=sum(!tested & x$exclusion_reason=="zero_variance",na.rm=TRUE),stringsAsFactors=FALSE))
  })
  z<-do.call(rbind,out);rownames(z)<-NULL;z
}

species_loco_feature_summary <- function(summary) {
  z<-lapply(unique(summary$feature),function(f) {
    x<-summary[summary$feature==f,,drop=FALSE]
    data.frame(feature=f,display_name=x$display_name[1],n_primary_significant_contrasts=nrow(x),
      n_contrasts_without_reversal=sum(!x$any_direction_reversal,na.rm=TRUE),
      all_significant_contrasts_direction_stable=all(x$all_testable_directions_concordant),
      minimum_testable_LOCO_count=min(x$n_LOCO_testable),
      maximum_absolute_coefficient_change=if(any(is.finite(x$maximum_absolute_coefficient_change)))max(x$maximum_absolute_coefficient_change,na.rm=TRUE)else NA_real_,stringsAsFactors=FALSE)
  })
  x<-do.call(rbind,z);rownames(x)<-NULL;x
}

species_loco_acceptance <- function(summary,features) {
  observed<-c(records=nrow(summary),species=nrow(features),stable=sum(summary$all_testable_directions_concordant),
    no_reversal_species=sum(features$n_contrasts_without_reversal==features$n_primary_significant_contrasts),
    q_four=sum(summary$n_LOCO_q_lt_0.05>=4L),q_five=sum(summary$n_LOCO_q_lt_0.05==5L))
  assert_true(identical(as.integer(observed),c(120L,79L,115L,74L,53L,8L)),
    paste("LOCO frozen acceptance mismatch; do not change q/thresholds:",paste(names(observed),observed,collapse="; ")))
  invisible(observed)
}

