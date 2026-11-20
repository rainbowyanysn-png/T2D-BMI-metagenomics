run_primary_diversity <- function(data, run, config) {
  require_project_packages(c("lmerTest", "emmeans", "performance", "vegan", "sva"))
  set.seed(config$models$primary_seed)
  metadata <- diversity_data(data)
  alpha <- lapply(c("shannon", "richness"), function(metric) {
    total <- lmerTest::lmer(stats::as.formula(paste(metric, "~ disease + BMI + age + sex + (1 | study_name)")), data = metadata)
    four <- lmerTest::lmer(stats::as.formula(paste(metric, "~ analysis_group + age + sex + (1 | study_name)")), data = metadata)
    pair <- as.data.frame(emmeans::emmeans(four, pairwise ~ analysis_group)$contrasts)
    pair$p.adj <- stats::p.adjust(pair$p.value, method = "BH")
    groups <- strsplit(as.character(pair$contrast), " - ", fixed = TRUE)
    pair$group1 <- vapply(groups, `[`, character(1), 1L); pair$group2 <- vapply(groups, `[`, character(1), 2L)
    pair <- pair[, c("group1", "group2", "estimate", "SE", "df", "t.ratio", "p.value", "p.adj")]
    safe_write_csv(pair, path_in_run(run, "tables/04_primary_diversity", paste0("BMI25_alpha_", metric, "_pairwise.csv")))
    data.frame(Metric = if (metric == "shannon") "Shannon" else "Richness", Overall_T2D_vs_Control_p = round(coef(summary(total))["diseaseT2D", "Pr(>|t|)"], 4), Four_Group_p = round(anova(four)$`Pr(>F)`[1], 4), ICC = round(performance::icc(four)$ICC_adjusted, 3))
  })
  safe_write_csv(do.call(rbind, alpha), path_in_run(run, "tables/04_primary_diversity", "BMI25_alpha_LMM_summary.csv"))
  safe_write_csv(metadata[, c("sample_id", "study_name", "disease", "BMI", "age", "sex", "analysis_group", "shannon", "richness")], path_in_run(run, "tables/04_primary_diversity", "BMI25_alpha_per_sample.csv"))
  draws <- replay_primary_jitter_rng(config$models$primary_seed)
  safe_write_csv(data.frame(historical_jitter_call = 1:5, seed = draws), path_in_run(run, "logs/primary_beta_jitter_rng_replay.csv"))
  beta <- historical_beta_clr(data)
  corrected <- sva::ComBat(dat = t(beta$clr), batch = beta$metadata$study_name, mod = stats::model.matrix(~ analysis_group + age + sex, data = beta$metadata), par.prior = TRUE, prior.plots = FALSE)
  distance <- stats::dist(t(corrected))
  md <- beta$metadata
  raw_distance <- stats::dist(beta$clr)
  # Source lines 416,417,584,605,606,630 and both cohort loops precede
  # the final exported tests. Preserve the complete original RNG sequence.
  before <- vegan::adonis2(raw_distance ~ study_name, data = md, permutations = config$models$diversity_permutations)
  after <- vegan::adonis2(distance ~ study_name, data = md, permutations = config$models$diversity_permutations)
  invisible(vegan::adonis2(raw_distance ~ study_name + analysis_group, data = md, permutations = config$models$diversity_permutations))
  invisible(vegan::adonis2(distance ~ analysis_group + age + sex, data = md, permutations = config$models$diversity_permutations))
  invisible(vegan::permutest(vegan::betadisper(distance, md$analysis_group), permutations = config$models$diversity_permutations))
  invisible(vegan::adonis2(distance ~ disease + age + sex, data = md, permutations = config$models$diversity_permutations))
  cohort_results <- list()
  for (repeat_loop in 1:2) {
    rows <- lapply(sort(unique(as.character(md$study_name))), function(study) {
      one <- md[md$study_name == study, , drop = FALSE]
      assert_true(nrow(one) >= 10L && length(unique(one$analysis_group)) >= 2L, paste("Historical beta cohort QC", study))
      sub_distance <- stats::as.dist(as.matrix(distance)[one$sample_id, one$sample_id])
      fit <- vegan::adonis2(sub_distance ~ analysis_group + age + sex, data = one, permutations = config$models$diversity_permutations)
      data.frame(study_name = study, R2 = fit$R2[1], F_statistic = fit$F[1], p_value = fit$`Pr(>F)`[1], n_samples = nrow(one))
    })
    cohort_results[[repeat_loop]] <- do.call(rbind, rows)
  }
  safe_write_csv(cohort_results[[1]][, c("study_name", "R2", "p_value", "n_samples")], path_in_run(run, "tables/04_primary_diversity/BMI25_PERMANOVA_stratified_by_study.csv"))
  table <- cohort_results[[2]]
  table$Comparison <- "Four-Group Phenotype Stratification"
  table$R2 <- round(table$R2, 4); table$F_statistic <- round(table$F_statistic, 3)
  safe_write_csv(table, path_in_run(run, "tables/04_primary_diversity/BMI25_Table_PERMANOVA_Stratified_Robustness.csv"))
  safe_write_csv(data.frame(Stage = c("Before_Correction", "After_Correction"), Tested_Variable = "study_name", R2 = c(before$R2[1], after$R2[1]), F_value = c(before$F[1], after$F[1]), P_value = c(before$`Pr(>F)`[1], after$`Pr(>F)`[1])), path_in_run(run, "tables/04_primary_diversity/BMI25_batch_effect_permanova.csv"))
  four <- vegan::adonis2(distance ~ analysis_group + age + sex, data = md, permutations = config$models$diversity_permutations)
  total <- vegan::adonis2(distance ~ disease + age + sex, data = md, permutations = config$models$diversity_permutations)
  dispersion <- vegan::permutest(vegan::betadisper(distance, md$analysis_group), permutations = config$models$diversity_permutations)
  result <- data.frame(Analysis = c("4_Group_Comparison", "T2D_vs_Control", "Dispersion_Test"), Variable = c("analysis_group", "disease", "group_dispersion"), R2 = c(four$R2[1], total$R2[1], NA), F_value = c(four$F[1], total$F[1], dispersion$tab$F[1]), P_value = c(four$`Pr(>F)`[1], total$`Pr(>F)`[1], dispersion$tab$`Pr(>F)`[1]))
  safe_write_csv(result, path_in_run(run, "tables/04_primary_diversity", "BMI25_beta_permanova_global.csv"))
  coordinates <- stats::cmdscale(distance, k = 2, eig = TRUE)
  safe_write_csv(data.frame(sample_id = rownames(coordinates$points), PC1 = coordinates$points[,1], PC2 = coordinates$points[,2]), path_in_run(run, "tables/04_primary_diversity", "BMI25_beta_PCoA_coordinates.csv"))
  safe_write_csv(md, path_in_run(run, "tables/04_primary_diversity/BMI25_beta_sample_manifest.csv"))
  safe_write_csv(data.frame(feature = colnames(beta$clr)), path_in_run(run, "tables/04_primary_diversity/BMI25_beta_feature_manifest.csv"))
  safe_write_lines("Primary beta uses historical ComBat CLR and original jitter/permutation call sequence; formal interaction beta uses raw CLR. Only independent frozen-result validation can certify reproduction.", path_in_run(run, "logs", "primary_beta_provenance.txt"))
  invisible(result)
}
