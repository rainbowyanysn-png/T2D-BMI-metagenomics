# Replay only the statistical RNG effects of the five historical jitter plots.
# ggplot2 PositionJitter samples one seed and restores the outer RNG state
# around jitter generation. This is verified against actual builds below.
replay_primary_jitter_rng <- function(seed) {
  require_project_packages("ggplot2")
  set.seed(seed)
  example <- data.frame(x = factor(c("a", "a", "b", "b")), y = c(1, 2, 3, 4))
  for (i in 1:5) invisible(ggplot2::ggplot_build(ggplot2::ggplot(example, ggplot2::aes(x, y)) + ggplot2::geom_jitter(width = c(.15, .15, .1, .1, .2)[i])))
  built <- .Random.seed
  set.seed(seed)
  draws <- vapply(1:5, function(i) sample.int(.Machine$integer.max, 1L), integer(1))
  assert_equal(.Random.seed, built, "Five historical jitter RNG effects")
  draws
}

diversity_data <- function(data) {
  x <- as.data.frame(data$metadata)
  x$shannon <- vegan::diversity(t(data$species), index = "shannon")
  x$richness <- colSums(data$species > 0)
  x
}

historical_beta_clr <- function(data) {
  replaced <- zCompositions::cmultRepl(t(data$species), label = 0, method = "CZM", output = "p-counts", z.warning = .8, z.delete = TRUE)
  clr <- as.matrix(compositions::clr(compositions::acomp(replaced)))
  metadata <- as.data.frame(data$metadata[match(rownames(clr), data$metadata$sample_id), , drop = FALSE])
  assert_equal(nrow(clr), 1182L, "Historical beta retained sample count")
  assert_equal(ncol(clr), 138L, "Historical beta retained feature count")
  assert_true(identical(rownames(clr), metadata$sample_id), "Beta CLR sample alignment")
  list(clr = clr, metadata = metadata)
}

run_diversity_interactions <- function(data, run, config) {
  require_project_packages(c("lmerTest", "vegan", "permute", "zCompositions", "compositions"))
  set.seed(config$models$continuous_interaction_seed)
  metadata <- diversity_data(data)
  metadata$Shannon <- metadata$shannon; metadata$Richness <- metadata$richness
  metadata$age_numeric <- suppressWarnings(as.numeric(as.character(metadata$age)))
  metadata$BMI_numeric <- suppressWarnings(as.numeric(as.character(metadata$BMI)))
  metadata$disease <- factor(as.character(metadata$disease), levels = c("control", "T2D"))
  metadata$sex <- factor(trimws(as.character(metadata$sex)))
  metadata$study_name <- factor(trimws(as.character(metadata$study_name)))
  metadata$BMI25 <- factor(ifelse(metadata$BMI_numeric >= 25, "Overweight", "Non_Overweight"), levels = c("Non_Overweight", "Overweight"))
  metadata$BMI_centered <- metadata$BMI_numeric - mean(metadata$BMI_numeric, na.rm = TRUE)
  alpha <- list()
  for (bmi in c("BMI25", "BMI_centered")) {
    rows <- lapply(c("Shannon", "Richness"), function(metric) {
      formula <- stats::as.formula(paste(metric, "~ age_numeric + sex + disease *", bmi, "+ (1 | study_name)"))
      fit <- lmerTest::lmer(formula, data = metadata, REML = FALSE)
      coefficients <- coef(summary(fit))
      term <- if (bmi == "BMI25") "diseaseT2D:BMI25Overweight" else "diseaseT2D:BMI_centered"
      assert_true(term %in% rownames(coefficients), "Missing alpha interaction term")
      data.frame(metric = metric, interaction_term = term, coefficient = coefficients[term, "Estimate"], stderr = coefficients[term, "Std. Error"], df = coefficients[term, "df"], p_interaction = coefficients[term, "Pr(>|t|)"], n_samples = nrow(metadata))
    })
    result <- do.call(rbind, rows)
    result$ci_lower <- result$coefficient - 1.96 * result$stderr
    result$ci_upper <- result$coefficient + 1.96 * result$stderr
    result$q_interaction_BH <- stats::p.adjust(result$p_interaction, method = "BH")
    if (bmi != "BMI25") {
      result$coefficient_per_5_BMI <- result$coefficient * 5
      result$ci_lower_per_5_BMI <- result$ci_lower * 5
      result$ci_upper_per_5_BMI <- result$ci_upper * 5
    }
    result <- result[order(result$q_interaction_BH, result$p_interaction), , drop = FALSE]
    safe_write_csv(result, path_in_run(run, "tables/05_diversity_interactions", paste0("alpha_", bmi, "_all.csv")))
    alpha[[bmi]] <- result
  }
  beta <- historical_beta_clr(data)
  beta_metadata <- metadata[match(beta$metadata$sample_id, metadata$sample_id), , drop = FALSE]
  distance <- stats::dist(beta$clr)
  control <- permute::how(blocks = beta_metadata$study_name, nperm = config$models$diversity_permutations)
  category <- vegan::adonis2(distance ~ age_numeric + sex + study_name + disease * BMI25, data = beta_metadata, permutations = control, by = "terms")
  continuous <- vegan::adonis2(distance ~ age_numeric + sex + study_name + disease * BMI_centered, data = beta_metadata, permutations = control, by = "terms")
  tables <- list(BMI25 = category, BMI_centered = continuous)
  rows <- lapply(names(tables), function(bmi) {
    x <- as.data.frame(tables[[bmi]])
    term <- paste0("disease:", bmi)
    safe_write_csv(data.frame(term = rownames(x), x, check.names = FALSE), path_in_run(run, "tables/05_diversity_interactions", paste0("beta_", bmi, "_all_terms.csv")))
    data.frame(interaction_term = term, df = x[term, "Df"], sum_of_squares = x[term, "SumOfSqs"], R2 = x[term, "R2"], F = x[term, "F"], p_interaction = x[term, "Pr(>F)"], n_samples = nrow(beta_metadata))
  })
  summary <- do.call(rbind, rows)
  summary$q_interaction_BH <- stats::p.adjust(summary$p_interaction, method = "BH")
  safe_write_csv(summary, path_in_run(run, "tables/05_diversity_interactions", "beta_interaction_summary.csv"))
  safe_write_csv(data.frame(sample_id = rownames(beta$clr)), path_in_run(run, "tables/05_diversity_interactions", "beta_retained_sample_ids.csv"))
  safe_write_csv(data.frame(feature = colnames(beta$clr)), path_in_run(run, "tables/05_diversity_interactions", "beta_retained_feature_ids.csv"))
  list(alpha = alpha, beta = summary)
}
