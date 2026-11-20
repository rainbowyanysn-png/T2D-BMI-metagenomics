run_pooled_pathway_maaslin <- function(data, run, config) {
  require_project_packages(c("Maaslin2", "tibble"))
  check_frozen_versions(config)
  cfg <- config$pooled_pathway_maaslin
  md <- as.data.frame(data$metadata)
  stopifnot(nrow(md) == 1186L, sum(md$disease == "T2D") == 760L,
            sum(md$disease == "control") == 426L, nrow(data$pathways) == 5945L,
            identical(colnames(data$pathways), md$sample_id), !anyDuplicated(md$sample_id),
            isTRUE(cfg$standardize), cfg$normalization == "TSS", cfg$transform == "LOG",
            identical(unlist(cfg$fixed_effects), c("disease", "BMI", "age", "sex", "study_name")),
            length(cfg$random_effects) == 0L, cfg$min_prevalence == .10)
  md$disease <- factor(as.character(md$disease), levels = c("control", "T2D"))
  md$sex <- factor(as.character(md$sex))
  md$study_name <- factor(as.character(md$study_name))
  stopifnot(!anyNA(md[, c("disease", "BMI", "age", "sex", "study_name")]))
  references <- c("disease,control", categorical_references(md, unlist(cfg$fixed_effects)))
  set.seed(cfg$seed)
  fit <- audited_maaslin2(
    input_data = data$pathways, input_metadata = tibble::column_to_rownames(md, "sample_id"),
    output = dirname(path_in_run(run, "models/03a_pooled_MaAsLin2_pathways/MaAsLin2/.keep")),
    fixed_effects = c("disease", "BMI", "age", "sex", "study_name"), random_effects = NULL,
    reference = references, normalization = "TSS", transform = "LOG",
    min_prevalence = .10, standardize = TRUE, analysis_method = "LM",
    plot_heatmap = FALSE, plot_scatter = FALSE)
  x <- fit$results[fit$results$metadata == "disease" & fit$results$value == "T2D", , drop = FALSE]
  x$formal_q <- x$qval
  x$qval_disease_term_BH <- stats::p.adjust(x$pval, "BH")
  x$audit_target_family_size <- nrow(x)
  x$q_family <- "MaAsLin2 full-model feature-by-covariate coefficients"
  x$q_family_size <- sum(!is.na(fit$results$pval))
  x$formal_significant <- x$qval < .05
  x$model_id <- "Pooled_MaAsLin2_pathways_T2D_BMI_age_sex_fixed_study"
  x$model_formula <- "pathway ~ disease + BMI + age + sex + study_name"
  x$study_effect <- "Fixed"
  x$comparison_direction <- "T2D versus control"
  stopifnot(nrow(x) == 5945L, !anyDuplicated(x$feature), all(x$N == 1186L),
            all(is.finite(x$coef)), all(is.finite(x$pval)), all(is.finite(x$qval)))
  x <- x[order(x$qval, x$pval, x$feature), , drop = FALSE]
  keep <- rowMeans(data$pathways_raw > 0) >= config$filters$global_prevalence
  ids <- rownames(data$pathways_raw)[keep]
  stopifnot(identical(make.names(ids, unique = TRUE), rownames(data$pathways)))
  x$pathway <- ids[match(x$feature, rownames(data$pathways))]
  stopifnot(!anyNA(x$pathway))
  checks <- data.frame(check = "tested_5945", observed = nrow(x), expected = cfg$expected_n_tested)
  if(!is.null(cfg$expected_full_model_significant))checks<-rbind(checks,data.frame(
    check="new_fixed_study_full_model_q",observed=sum(x$qval<.05),expected=cfg$expected_full_model_significant))
  if(!is.null(cfg$expected_target_BH_audit_significant))checks<-rbind(checks,data.frame(
    check="new_fixed_study_target_BH_audit",observed=sum(x$qval_disease_term_BH<.05),expected=cfg$expected_target_BH_audit_significant))
  checks$pass <- checks$observed == checks$expected
  safe_write_csv(checks, path_in_run(run, "tables/03a_pooled_MaAsLin2_pathways/numeric_checks.csv"))
  safe_write_csv(x, path_in_run(run, "tables/03a_pooled_MaAsLin2_pathways/Pooled_pathways_all.csv"))
  safe_write_csv(data.frame(model_id = unique(x$model_id), n_tested = nrow(x),
    formal_full_model_q_significant = sum(x$qval < .05),
    target_BH_audit_significant = sum(x$qval_disease_term_BH < .05),
    reference_definition = paste(references, collapse = ";"),
    acceptance_status = "USER_AUTHORIZED_FIXED_STUDY_MODEL_REQUIRES_NEW_RESULT_REVIEW"),
    path_in_run(run, "tables/03a_pooled_MaAsLin2_pathways/Model_result_summary.csv"))
  assert_true(all(checks$pass), "Pooled pathway numeric assertions failed; publication export is blocked.")
}
