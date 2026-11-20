run_additive_t2d_bmi <- function(data, run, config, expected) {
  require_project_packages(c("Maaslin2", "tibble"))
  check_frozen_versions(config)
  set.seed(config$models$additive_seed)
  metadata <- data$metadata
  metadata$BMI_per5 <- (metadata$BMI - mean(metadata$BMI)) / config$models$additive_bmi_scale
  assert_true(abs(mean(metadata$BMI) - config$models$additive_bmi_center) <= config$project$numeric_tolerance, "Frozen additive BMI centring mean was not reproduced.")
  assert_true(identical(colnames(data$species), metadata$sample_id), "Additive species/metadata alignment failed.")
  model_dir <- path_in_run(run, "models", "06_additive_T2D_BMI", "MaAsLin2", ".keep")
  model_dir <- dirname(model_dir)
  fixed <- c("disease", "BMI_per5", "age", "sex", "study_name")
  fit <- audited_maaslin2(input_data = data$species, input_metadata = tibble::column_to_rownames(as.data.frame(metadata), "sample_id"), output = model_dir, fixed_effects = fixed, reference = c("disease,control", categorical_references(metadata, fixed)), normalization = config$models$maaslin_normalization, transform = config$models$maaslin_transform, standardize = config$models$additive_maaslin_standardize, min_prevalence = 0.10, plot_heatmap = FALSE, plot_scatter = FALSE)
  all <- fit$results
  if ("N.not.zero" %in% names(all)) names(all)[names(all) == "N.not.zero"] <- "N.not.0"
  disease <- all[all$metadata == "disease", , drop = FALSE]
  bmi <- all[all$metadata == "BMI_per5", , drop = FALSE]
  assert_equal(nrow(disease), 200L, "additive disease coefficient count")
  assert_equal(nrow(bmi), 200L, "additive BMI coefficient count")
  decorate <- function(x, exposure) {
    x$exposure <- exposure
    x$ci_lower <- x$coef - 1.96 * x$stderr
    x$ci_upper <- x$coef + 1.96 * x$stderr
    x$q_maaslin_full_model <- x$qval
    x$q_within_exposure_BH <- stats::p.adjust(x$pval, method = "BH")
    x$formal_q <- x$q_within_exposure_BH
    x$q_family <- paste("All 200 additive", exposure, "coefficients")
    x$q_family_size <- nrow(x)
    x$model_id <- "Additive_disease_BMI_per5_age_sex_study_standardize_FALSE"
    x$comparison_direction <- if (exposure == "T2D_vs_Control") "T2D minus control" else "BMI increase of 5 kg/m2"
    x$direction <- ifelse(x$coef > 0, "Positive", ifelse(x$coef < 0, "Negative", "Zero"))
    x
  }
  disease <- decorate(disease, "T2D_vs_Control")
  bmi <- decorate(bmi, "BMI_per_5kg_m2")
  joint <- stats::p.adjust(c(disease$pval, bmi$pval), method = "BH")
  disease$q_joint_disease_BMI_BH <- joint[seq_len(nrow(disease))]
  bmi$q_joint_disease_BMI_BH <- joint[nrow(disease) + seq_len(nrow(bmi))]
  disease <- disease[order(disease$q_within_exposure_BH, disease$pval), , drop = FALSE]
  bmi <- bmi[order(bmi$q_within_exposure_BH, bmi$pval), , drop = FALSE]
  long <- rbind(disease, bmi)
  columns <- c("feature", "coef", "stderr", "ci_lower", "ci_upper", "pval", "qval", "q_within_exposure_BH", "q_joint_disease_BMI_BH", "direction")
  wide <- merge(disease[, columns], bmi[, columns], by = "feature", suffixes = c("_T2D", "_BMI_per5"))
  wide$association_class <- ifelse(wide$q_within_exposure_BH_T2D < config$project$alpha & wide$q_within_exposure_BH_BMI_per5 < config$project$alpha, "Both", ifelse(wide$q_within_exposure_BH_T2D < config$project$alpha, "T2D-only", ifelse(wide$q_within_exposure_BH_BMI_per5 < config$project$alpha, "BMI-only", "Neither")))
  summary <- data.frame(exposure = c("T2D_vs_Control", "BMI_per_5kg_m2"), n_tested = c(nrow(disease), nrow(bmi)), n_maaslin_q_lt_005 = c(sum(disease$qval < .05), sum(bmi$qval < .05)), n_term_BH_q_lt_005 = c(sum(disease$q_within_exposure_BH < .05), sum(bmi$q_within_exposure_BH < .05)), n_positive_term_BH = c(sum(disease$q_within_exposure_BH < .05 & disease$coef > 0), sum(bmi$q_within_exposure_BH < .05 & bmi$coef > 0)), n_negative_term_BH = c(sum(disease$q_within_exposure_BH < .05 & disease$coef < 0), sum(bmi$q_within_exposure_BH < .05 & bmi$coef < 0)), q_family = c("200 disease coefficients", "200 BMI coefficients"), stringsAsFactors = FALSE)
  assert_true(identical(summary$n_term_BH_q_lt_005, c(expected$additive$t2d_term_bh, expected$additive$bmi_term_bh)), "Additive term-specific BH count does not match frozen gate.")
  assert_true(summary$n_positive_term_BH[2] == expected$additive$bmi_positive && summary$n_negative_term_BH[2] == expected$additive$bmi_negative, "Additive BMI direction count does not match frozen gate.")
  assert_true(sum(wide$association_class == "Both") == expected$additive$both && sum(wide$association_class == "BMI-only") == expected$additive$bmi_only, "Additive overlap counts do not match frozen gate.")
  safe_write_csv(long, path_in_run(run, "tables", "06_additive_all_target_coefficients.csv"))
  safe_write_csv(wide, path_in_run(run, "tables", "06_additive_coefficients_wide.csv"))
  safe_write_csv(summary, path_in_run(run, "tables", "06_additive_summary.csv"))
  list(long = long, wide = wide, summary = summary)
}
