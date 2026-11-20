clr_prepared_metadata <- function(metadata) {
  metadata <- as.data.frame(medication_model_covariate(metadata))
  assert_columns(metadata,c("sample_id","BMI23_group","BMI23_threshold"),"CLR input metadata")
  metadata$BMI_numeric <- as.numeric(as.character(metadata$BMI)); metadata$age_numeric <- as.numeric(as.character(metadata$age))
  metadata$sex <- factor(trimws(as.character(metadata$sex))); metadata$study_name <- factor(trimws(as.character(metadata$study_name)))
  metadata$T2D_binary <- as.integer(metadata$disease == "T2D")
  metadata$BMI25_binary <- as.integer(metadata$BMI_numeric >= 25)
  metadata$BMI23_binary <- as.integer(as.character(metadata$BMI23_group) %in% c("Overweight_Control","Overweight_T2D"))
  assert_true(all(metadata$BMI23_binary==as.integer(metadata$BMI_numeric>=metadata$BMI23_threshold)),"CLR BMI23 must reuse harmonized cohort-adapted grouping")
  assert_true(all(metadata$BMI23_threshold==ifelse(as.character(metadata$study_name)%in%c("QinJ_2012","YuJ_2015"),23,25)),"CLR BMI23 cohort thresholds differ")
  assert_true(identical(as.integer(table(factor(as.character(metadata$BMI23_group),levels=c("Non_Overweight_Control","Overweight_Control","Non_Overweight_T2D","Overweight_T2D")))),c(261L,165L,105L,655L)),"CLR adapted BMI23 group counts differ")
  metadata$T2D_x_BMI25 <- metadata$T2D_binary * metadata$BMI25_binary
  metadata$T2D_x_BMI23 <- metadata$T2D_binary * metadata$BMI23_binary
  metadata$BMI_centered <- metadata$BMI_numeric - mean(metadata$BMI_numeric, na.rm = TRUE)
  metadata$T2D_x_BMI_continuous <- metadata$T2D_binary * metadata$BMI_centered
  metadata
}

clr_pair_model_spec <- function(id,domain,metadata,pair) {
  md <- metadata[metadata$analysis_group %in% pair,,drop=FALSE]
  md$Group_binary <- as.integer(as.character(md$analysis_group)==pair[2])
  adjusted <- domain=="Taxonomy" && id=="Overweight_vs_Lean_T2D"
  reference <- ""; formula <- "~ Group_binary + age_numeric + sex + study_name"
  if(adjusted) {
    md$antidiabetes_group <- droplevels(factor(as.character(md$antidiabetes_group)))
    reference <- medication_model_reference(md)
    md$antidiabetes_group <- stats::relevel(md$antidiabetes_group,ref=reference)
    assert_true(!anyNA(md$antidiabetes_group),"Missing CLR medication category")
    formula <- paste(formula,"+ antidiabetes_group")
  }
  list(metadata=md,formula=formula,medication_adjusted=adjusted,medication_reference=reference)
}

run_clr_sensitivity <- function(data, run, config) {
  require_project_packages(c("zCompositions", "compositions", "limma", "dplyr"))
  metadata <- clr_prepared_metadata(data$metadata)
  version <- "RC4_CLR_ALIGNED_20261119"
  manifest_columns <- c("sample_id","study_name","disease","BMI_numeric","age_numeric","sex","analysis_group","BMI23_group","BMI23_threshold","antidiabetes_group","T2D_binary","BMI25_binary","BMI23_binary","BMI_centered")
  safe_write_csv(metadata[,manifest_columns],path_in_run(run,"tables/08_CLR_sensitivity/CLR_sample_metadata.csv"))
  matrices <- list(Taxonomy = data$species, Pathway = data$pathways_raw[rowMeans(data$pathways_raw > 0) >= .10, , drop = FALSE])
  output <- contracts <- list()
  for (domain in names(matrices)) {
    raw <- matrices[[domain]]
    assert_equal(nrow(raw), if (domain == "Taxonomy") 200L else 5945L, paste(domain, "CLR global feature panel"))
    replaced <- zCompositions::cmultRepl(t(raw), label = 0, method = "CZM", output = "p-counts", z.warning = .80, z.delete = FALSE, suppress.print = TRUE)
    assert_true(identical(dim(replaced), dim(t(raw))), "CZM changed the formal CLR sample/feature panel")
    clr <- t(as.matrix(compositions::clr(compositions::acomp(as.matrix(replaced)))))
    dimnames(clr) <- dimnames(raw)
    assert_true(all(is.finite(clr)), "Non-finite formal CLR values")
    fit_one <- function(id, md, formula, term, subset_filter = FALSE, medication_adjusted=FALSE, medication_reference="") {
      md$sex <- droplevels(md$sex); md$study_name <- droplevels(md$study_name)
      matrix <- clr[, md$sample_id, drop = FALSE]
      if (subset_filter) matrix <- matrix[rowMeans(raw[, md$sample_id, drop = FALSE] > 0) >= .10, , drop = FALSE]
      design <- stats::model.matrix(stats::as.formula(formula), data = md)
      assert_true(nrow(design)==nrow(md),paste("CLR silently omitted samples:",id))
      assert_true(qr(design)$rank == ncol(design) && term %in% colnames(design), paste("CLR design failed:", id))
      fit <- limma::eBayes(limma::lmFit(matrix, design), trend = FALSE, robust = FALSE)
      index <- match(term, colnames(design))
      result <- data.frame(feature = rownames(matrix), analysis_domain = domain, analysis_id = id, model_formula = formula, target_term = term, coef_CLR = fit$coefficients[, index], t_statistic = fit$t[, index], pval_CLR = fit$p.value[, index], stringsAsFactors = FALSE)
      result$qval_CLR_BH <- stats::p.adjust(result$pval_CLR, method = "BH")
      result$direction_CLR <- ifelse(result$coef_CLR > 0, "Positive", ifelse(result$coef_CLR < 0, "Negative", "Zero"))
      result$n_samples <- nrow(md); result$n_features_tested <- nrow(matrix)
      result <- result[order(result$qval_CLR_BH, result$pval_CLR), , drop = FALSE]
      result$stderr_moderated <- fit$stdev.unscaled[match(result$feature, rownames(matrix)), index] * sqrt(fit$s2.post[match(result$feature, rownames(matrix))])
      result$q_family <- paste(domain, id, "target term BH")
      result$q_family_size <- nrow(matrix)
      result$BMI_definition <- if (startsWith(id, "BMI23")) "Cohort-adapted BMI23: QinJ_2012/YuJ_2015 >=23; other cohorts >=25" else if (startsWith(id, "Continuous")) "BMI centred across all 1186 samples" else "BMI25: BMI <25 versus BMI >=25"
      result$clr_contract_version <- version
      result$medication_adjusted <- medication_adjusted
      result$medication_reference <- medication_reference
      threshold <- if (startsWith(id, "BMI23")) "BMI23_adapted" else if (startsWith(id, "Continuous")) "continuous_BMI" else "BMI25"
      file <- paste0(threshold,"_",domain,"_",id,"_all.csv")
      safe_write_csv(result, path_in_run(run, "tables/08_CLR_sensitivity", file))
      contracts[[paste(domain,id)]] <<- data.frame(result_file=file,analysis_domain=domain,analysis_id=id,model_formula=formula,target_term=term,n_samples=nrow(md),n_features_tested=nrow(matrix),n_significant=sum(result$qval_CLR_BH<.05,na.rm=TRUE),BMI_definition=result$BMI_definition[1],medication_adjusted=medication_adjusted,medication_reference=medication_reference,design_rank=qr(design)$rank,design_columns=ncol(design),clr_contract_version=version)
      result
    }
    results <- list()
    results$Total_T2D_vs_control <- fit_one("Total_T2D_vs_control", metadata, "~ T2D_binary + BMI_numeric + age_numeric + sex + study_name", "T2D_binary", TRUE)
    pairs <- list(Lean_T2D_vs_Lean_Control = c("Non_Overweight_Control", "Non_Overweight_T2D"), Overweight_T2D_vs_Overweight_Control = c("Overweight_Control", "Overweight_T2D"), Overweight_vs_Lean_Controls = c("Non_Overweight_Control", "Overweight_Control"), Overweight_vs_Lean_T2D = c("Non_Overweight_T2D", "Overweight_T2D"))
    for (id in names(pairs)) {
      spec <- clr_pair_model_spec(id,domain,metadata,pairs[[id]])
      results[[id]] <- fit_one(id,spec$metadata,spec$formula,"Group_binary",TRUE,spec$medication_adjusted,spec$medication_reference)
    }
    results$BMI25_Formal_Interaction <- fit_one("BMI25_Formal_Interaction", metadata, "~ T2D_binary + BMI25_binary + T2D_x_BMI25 + age_numeric + sex + study_name", "T2D_x_BMI25")
    results$ContinuousBMI_Formal_Interaction <- fit_one("ContinuousBMI_Formal_Interaction", metadata, "~ T2D_binary + BMI_centered + T2D_x_BMI_continuous + age_numeric + sex + study_name", "T2D_x_BMI_continuous")
    results$BMI23_Formal_Interaction <- fit_one("BMI23_Formal_Interaction", metadata, "~ T2D_binary + BMI23_binary + T2D_x_BMI23 + age_numeric + sex + study_name", "T2D_x_BMI23")
    output[[domain]] <- results
  }
  safe_write_csv(dplyr::bind_rows(contracts),path_in_run(run,"tables/08_CLR_sensitivity/CLR_model_contracts.csv"))
  output
}
