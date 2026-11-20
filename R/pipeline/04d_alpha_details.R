# Recomputes supplementary evidence from the passed current-run original inputs.
run_alpha_details <- function(data, run, config) {
  root <- run$root
  read <- function(x) read_current_result(run, file.path("tables", x))
  saved <- read("01_sample_manifest.csv")
  expected <- read_expected(root)$inputs
  tolerance <- config$project$numeric_tolerance
  compare <- function(x, y, label) {
    assert_true(isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE)), label)
  }
  calculate <- function() {
    require_project_packages(c("vegan", "lmerTest", "lme4", "emmeans", "performance"))
    directory <- "04_primary_diversity"
    md <- diversity_data(data)
    saved_alpha <- read(file.path(directory, "BMI25_alpha_per_sample.csv"))
    saved_summary <- read(file.path(directory, "BMI25_alpha_LMM_summary.csv"))
    stopifnot(identical(md$sample_id, as.character(saved$sample_id)),
              identical(md$sample_id, as.character(saved_alpha$sample_id)),
              identical(as.character(md$analysis_group), saved_alpha$analysis_group))
    compare(md$shannon, saved_alpha$shannon, "Saved Shannon values differ.")
    stopifnot(all(md$richness == saved_alpha$richness))
    set.seed(config$models$primary_seed)
    coefficients <- list()
    variance <- list()
    for (metric in c("shannon", "richness")) {
      total <- lmerTest::lmer(stats::as.formula(paste(metric, "~ disease + BMI + age + sex + (1 | study_name)")), data = md)
      four <- lmerTest::lmer(stats::as.formula(paste(metric, "~ analysis_group + age + sex + (1 | study_name)")), data = md)
      display <- if (metric == "shannon") "Shannon" else "Richness"
      old <- saved_summary[saved_summary$Metric == display, , drop = FALSE]
      compare(
        c(round(coef(summary(total))["diseaseT2D", "Pr(>|t|)"], 4),
          round(anova(four)$`Pr(>F)`[1], 4), round(performance::icc(four)$ICC_adjusted, 3)),
        c(old$Overall_T2D_vs_Control_p, old$Four_Group_p, old$ICC), paste(metric, "LMM summary differs.")
      )
      pair <- as.data.frame(emmeans::emmeans(four, pairwise ~ analysis_group)$contrasts)
      old_pair <- read(file.path(directory, paste0("BMI25_alpha_", metric, "_pairwise.csv")))
      stopifnot(identical(as.character(pair$contrast), paste(old_pair$group1, "-", old_pair$group2)))
      for (column in c("estimate", "SE", "df", "t.ratio", "p.value")) {
        compare(pair[[column]], old_pair[[column]], paste(metric, "pairwise mismatch:", column))
      }
      compare(stats::p.adjust(pair$p.value, method = "BH"), old_pair$p.adj, paste(metric, "pairwise BH mismatch."))
      for (model in c("total", "four")) {
        fit <- if (model == "total") total else four
        key <- paste(metric, model, sep = "_")
        a <- as.data.frame(coef(summary(fit)), check.names = FALSE)
        coefficients[[key]] <- data.frame(metric = display, model = model,
          formula = paste(deparse(stats::formula(fit)), collapse = " "), term = rownames(a),
          n_samples = stats::nobs(fit), a, row.names = NULL, check.names = FALSE)
        variance[[key]] <- data.frame(metric = display, model = model,
          as.data.frame(lme4::VarCorr(fit)), check.names = FALSE)
      }
    }
    list(BMI25_alpha_LMM_fixed_coefficients = do.call(rbind, coefficients),
         BMI25_alpha_LMM_variance_components = do.call(rbind, variance))
  }

  result <- calculate()
  for (file in names(result)) safe_write_csv(result[[file]],
    path_in_run(run, "tables/04d_alpha_details", paste0(file, ".csv")))
  invisible(result)
}
