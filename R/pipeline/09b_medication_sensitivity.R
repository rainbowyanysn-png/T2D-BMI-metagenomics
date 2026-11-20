# Compares saved current-run primary and sensitivity species estimates; no fitting.
run_medication_sensitivity <- function(data, run, config) {
  read <- function(x) read_current_result(run, file.path("tables",x))
  assert_true(identical(config$models$primary_maaslin_standardize,TRUE) &&
    identical(config$models$bmi23_maaslin_standardize,TRUE) &&
    identical(config$models$no_medication_maaslin_standardize,TRUE), "Medication species models require standardize TRUE")
  results <- list()
  for (threshold in c("BMI25","BMI23")) {
    a <- read(if(threshold=="BMI25") "02_primary_BMI25_species/OT2D_vs_NT2D_BMI25_all.csv" else "07_BMI23_sensitivity/OT2D_vs_NT2D_BMI23_all.csv")
    u <- read(paste0("09_no_medication_covariate/",threshold,"_all.csv"))
    assert_true(nrow(a)==190L && nrow(u)==190L && !anyDuplicated(a$feature) && !anyDuplicated(u$feature) && setequal(a$feature,u$feature), "Medication comparison feature panel differs")
    u <- u[match(a$feature,u$feature),,drop=FALSE]
    assert_true(all(a$N==760L) && all(u$N==760L), "Medication model N differs")
    paired <- data.frame(feature=a$feature,BMI_definition=threshold,
      coef_adjusted=a$coef,coef_unadjusted=u$coef,stderr_adjusted=a$stderr,stderr_unadjusted=u$stderr,
      pval_adjusted=a$pval,pval_unadjusted=u$pval,qval_adjusted=a$qval,qval_unadjusted=u$qval,
      adjusted_significant=a$qval<.05,unadjusted_significant=u$qval<.05,
      direction_agrees=sign(a$coef)==sign(u$coef))
    paired$significance_set <- ifelse(paired$adjusted_significant & paired$unadjusted_significant,"Both",
      ifelse(paired$adjusted_significant,"Adjusted-only",ifelse(paired$unadjusted_significant,"Unadjusted-only","Neither")))
    summary <- data.frame(BMI_definition=threshold,matched_features=nrow(paired),
      adjusted_significant=sum(paired$adjusted_significant),unadjusted_significant=sum(paired$unadjusted_significant),
      shared_significant=sum(paired$significance_set=="Both"),adjusted_only=sum(paired$significance_set=="Adjusted-only"),
      unadjusted_only=sum(paired$significance_set=="Unadjusted-only"),Pearson_r=cor(a$coef,u$coef),
      direction_agreement=mean(paired$direction_agrees),direction_agreement_percent=100*mean(paired$direction_agrees),
      adjusted_standardize=TRUE,unadjusted_standardize=TRUE,
      interpretation="Medication-adjusted primary model versus model without the medication covariate. Not a metformin-naive analysis.")
    expected <- if(threshold=="BMI25") c(5L,4L,4L,1L) else c(3L,3L,3L,0L)
    assert_true(identical(as.integer(unlist(summary[1,c("adjusted_significant","unadjusted_significant","shared_significant","adjusted_only")],use.names=FALSE)),expected), "Medication significant-set counts differ")
    assert_true(abs(summary$Pearson_r-if(threshold=="BMI25") .9680455 else .9768742)<1e-6, "Medication coefficient correlation differs")
    assert_true(abs(summary$direction_agreement_percent-if(threshold=="BMI25")92.1053 else 94.2105)<1e-4,"Medication direction agreement differs")
    if(threshold=="BMI25") {
      bif <- paired[sub("^species[.]","",paired$feature)=="Bifidobacterium.adolescentis",,drop=FALSE]
      assert_true(nrow(bif)==1L && abs(bif$qval_adjusted-.0459046256)<1e-8 && abs(bif$qval_unadjusted-.0838221802)<1e-8,"B. adolescentis medication sensitivity qval differs")
    }
    tables <- list(adjusted_all=a,unadjusted_all=u,coefficients_paired=paired,concordance=summary,
      adjusted_significant=a[a$qval<.05,,drop=FALSE],unadjusted_significant=u[u$qval<.05,,drop=FALSE],
      shared_significant=paired[paired$significance_set=="Both",,drop=FALSE],
      adjusted_only=paired[paired$significance_set=="Adjusted-only",,drop=FALSE])
    for(name in names(tables)) safe_write_csv(tables[[name]],path_in_run(run,"tables/09b_medication_sensitivity",paste0(threshold,"_",name,".csv")))
    results[[threshold]] <- tables
  }
  invisible(results)
}
