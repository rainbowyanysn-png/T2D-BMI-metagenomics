run_no_medication_covariate <- function(data, run, config) {
  check_frozen_versions(config)
  set.seed(config$models$primary_seed)
  results <- list()
  for (threshold in c("BMI25", "BMI23")) {
    metadata <- as.data.frame(data$metadata[data$metadata$disease == "T2D", , drop = FALSE])
    group <- metadata[[paste0(threshold, "_group")]]
    metadata$Group <- factor(as.character(group), levels = c("Non_Overweight_T2D", "Overweight_T2D"))
    assert_equal(nrow(metadata), 760L, paste(threshold, "T2D samples"))
    expected <- if (threshold == "BMI25") c(138L, 622L) else c(105L, 655L)
    assert_true(identical(as.integer(table(metadata$Group)), expected), paste(threshold, "T2D four-group classification changed"))
    fixed <- c("Group", "age", "sex", "study_name")
    references <- c("Group,Non_Overweight_T2D", categorical_references(metadata, fixed))
    folder <- dirname(path_in_run(run, "models/09_no_medication_covariate", threshold, "MaAsLin2", ".keep"))
    fit <- audited_maaslin2(input_data = data$species[, metadata$sample_id, drop = FALSE], input_metadata = tibble::column_to_rownames(metadata, "sample_id"), output = folder, fixed_effects = fixed, random_effects = NULL, reference = references, normalization = "TSS", transform = "LOG", min_prevalence = .10, standardize = config$models$no_medication_maaslin_standardize, plot_heatmap = FALSE, plot_scatter = FALSE)
    target <- fit$results[fit$results$metadata == "Group" & fit$results$value == "Overweight_T2D", , drop = FALSE]
    target$qval_target_bh <- stats::p.adjust(target$pval, method = "BH")
    target <- target[order(target$qval_target_bh, target$pval), , drop = FALSE]
    target$formal_q <- target$qval
    target$q_family <- "MaAsLin2 original full-model qval"
    target$q_family_size <- sum(!is.na(fit$results$pval))
    target$ci_lower <- target$coef - 1.96 * target$stderr; target$ci_upper <- target$coef + 1.96 * target$stderr
    target$model_id <- paste0(threshold, "_T2D_internal_no_medication_standardize_TRUE")
    target$comparison_direction <- "Overweight_T2D vs Non_Overweight_T2D"
    safe_write_csv(target, path_in_run(run, "tables/09_no_medication_covariate", paste0(threshold, "_all.csv")))
    safe_write_csv(target[target$qval < .05, , drop = FALSE], path_in_run(run, "tables/09_no_medication_covariate", paste0(threshold, "_formal_q_lt_005.csv")))
    assert_equal(nrow(target), 190L, paste(threshold, "tested taxa"))
    assert_equal(sum(target$qval < .05), if (threshold == "BMI25") 4L else 3L, paste(threshold, "formal discoveries"))
    results[[threshold]] <- target
  }
  results
}
