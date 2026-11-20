interaction_metadata <- function(metadata, type = c("BMI25", "continuous")) {
  type <- match.arg(type)
  x <- metadata
  assert_equal(nrow(x), 1186L, "Interaction population before BMI centering")
  x$BMI_numeric <- as.numeric(as.character(x$BMI))
  x$age_numeric <- as.numeric(as.character(x$age))
  x$T2D_binary <- as.integer(as.character(x$disease) == "T2D")
  if (type == "BMI25") {
    x$BMI25_binary <- as.integer(x$BMI_numeric >= 25)
    x$T2D_x_BMI25 <- x$T2D_binary * x$BMI25_binary
    list(metadata = x, bmi_term = "BMI25_binary", target = "T2D_x_BMI25", label = "BMI25")
  } else {
    x$BMI_centered <- x$BMI_numeric - mean(x$BMI_numeric, na.rm = TRUE)
    x$T2D_x_BMI_continuous <- x$T2D_binary * x$BMI_centered
    list(metadata = x, bmi_term = "BMI_centered", target = "T2D_x_BMI_continuous", label = "continuous_BMI")
  }
}

run_one_formal_interaction <- function(matrix, metadata, type, domain, run, config) {
  spec <- interaction_metadata(metadata, type)
  x <- spec$metadata
  assert_true(identical(colnames(matrix), x$sample_id), paste(domain, type, "sample order mismatch"))
  dir <- path_in_run(run, "models", "05_formal_interactions", paste(spec$label, domain, sep = "_"), "MaAsLin2", ".keep")
  dir <- dirname(dir)
  fixed <- c("T2D_binary", spec$bmi_term, spec$target, "age_numeric", "sex", "study_name")
  set.seed(if (type == "BMI25") config$models$bmi25_interaction_seed else config$models$continuous_interaction_seed)
  fit <- audited_maaslin2(input_data = matrix, input_metadata = tibble::column_to_rownames(as.data.frame(x), "sample_id"), output = dir, fixed_effects = fixed, reference = categorical_references(x, fixed, largest_only = TRUE), normalization = config$models$maaslin_normalization, transform = config$models$maaslin_transform, standardize = config$models$interaction_maaslin_standardize, min_prevalence = 0, plot_heatmap = FALSE, plot_scatter = FALSE)
  target <- fit$results[fit$results$metadata == spec$target, , drop = FALSE]
  assert_equal(nrow(target), nrow(matrix), paste(domain, type, "target coefficient count"))
  target$ci_lower <- target$coef - 1.96 * target$stderr
  target$ci_upper <- target$coef + 1.96 * target$stderr
  if (type == "continuous") {
    target$coef_per_5_BMI <- target$coef * 5
    target$stderr_per_5_BMI <- target$stderr * 5
    target$ci_lower_per_5_BMI <- target$ci_lower * 5
    target$ci_upper_per_5_BMI <- target$ci_upper * 5
  }
  target$q_maaslin_full_model <- target$qval
  target$formal_q <- stats::p.adjust(target$pval, method = "BH")
  target$q_family <- paste0("Target ", spec$target, " P values within ", domain, " / ", spec$label)
  target$q_family_size <- nrow(target)
  target$comparison_id <- paste(spec$label, domain, sep = "_")
  target$direction <- ifelse(target$coef > 0, "Positive", ifelse(target$coef < 0, "Negative", "Zero"))
  target$formal_significant <- target$formal_q < config$project$alpha
  target <- target[order(target$formal_q, target$pval), , drop = FALSE]
  safe_write_csv(target, path_in_run(run, "tables", "05_formal_interactions", paste0(spec$label, "_", domain, "_all.csv")))
  target
}

run_formal_interactions <- function(data, run, config, expected) {
  require_project_packages(c("Maaslin2", "tibble"))
  jobs <- list(BMI25_species = list(matrix = data$species, type = "BMI25", domain = "species", expected = expected$interactions$bmi25_species), BMI25_pathways = list(matrix = data$pathways, type = "BMI25", domain = "pathways", expected = expected$interactions$bmi25_pathways), continuous_species = list(matrix = data$species, type = "continuous", domain = "species", expected = expected$interactions$continuous_species), continuous_pathways = list(matrix = data$pathways, type = "continuous", domain = "pathways", expected = expected$interactions$continuous_pathways))
  results <- lapply(names(jobs), function(id) { job <- jobs[[id]]; ans <- run_one_formal_interaction(job$matrix, data$metadata, job$type, job$domain, run, config); assert_equal(sum(ans$formal_significant), job$expected, paste(id, "formal interaction discoveries")); ans })
  names(results) <- names(jobs)
  safe_write_csv(data.frame(analysis = names(results), n_tested = vapply(results, nrow, integer(1)), n_formal_q_lt_005 = vapply(results, function(x) sum(x$formal_significant), integer(1))), path_in_run(run, "tables", "05_formal_interactions_summary.csv"))
  results
}
