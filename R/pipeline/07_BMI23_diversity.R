# Source: original BMI23 run_beta plus BMI23.R alpha correction, not new methods.
run_bmi23_diversity <- function(data, run, config) {
  require_project_packages(c("lme4", "emmeans", "vegan", "permute", "sva"))
  set.seed(config$models$bmi23_diversity_seed)
  md <- diversity_data(data)
  md$Shannon <- md$shannon; md$Richness <- md$richness
  md$analysis_group_BMI23_adapted <- factor(as.character(md$BMI23_group), levels = unlist(config$models$bmi23_pathway_levels))
  md$age_numeric <- as.numeric(as.character(md$age))
  rownames(md) <- md$sample_id
  write <- function(x, filename) safe_write_csv(x, path_in_run(run, "tables/07_BMI23_diversity", filename))
  alpha <- lapply(c("Shannon", "Richness"), function(metric) {
    full <- lme4::lmer(stats::as.formula(paste(metric, "~ analysis_group_BMI23_adapted + age_numeric + sex + (1 | study_name)")), data = md)
    null <- lme4::lmer(stats::as.formula(paste(metric, "~ age_numeric + sex + (1 | study_name)")), data = md)
    lrt <- as.data.frame(stats::anova(null, full)); lrt$metric <- metric
    pair <- as.data.frame(emmeans::contrast(emmeans::emmeans(full, ~ analysis_group_BMI23_adapted), method = "pairwise", adjust = "BH")); pair$metric <- metric
    coef <- as.data.frame(summary(full)$coefficients); coef$term <- rownames(coef); coef$metric <- metric; rownames(coef) <- NULL
    list(lrt = lrt, pair = pair, coefficients = coef)
  })
  write(md[, c("sample_id", "study_name", "disease", "BMI", "analysis_group_BMI23_adapted", "Shannon", "Richness")], "BMI23_Alpha_SampleLevel_Data.csv")
  for (key in c("lrt", "pair", "coefficients")) {
    filename <- c(lrt = "BMI23_Alpha_Group_Overall_LRT.csv", pair = "BMI23_Alpha_Pairwise_BH.csv", coefficients = "BMI23_Alpha_LMM_Coefficients.csv")[[key]]
    write(do.call(rbind, lapply(alpha, `[[`, key)), filename)
  }
  description <- do.call(rbind, lapply(c("Shannon", "Richness"), function(metric) stats::aggregate(md[[metric]], by = list(group = md$analysis_group_BMI23_adapted), FUN = function(x) c(n = length(x), mean = mean(x), sd = stats::sd(x), median = stats::median(x), IQR = stats::IQR(x)))))
  write(description, "BMI23_Alpha_Descriptive_Statistics.csv")
  beta <- historical_beta_clr(data)
  bmd <- md[match(beta$metadata$sample_id, md$sample_id), , drop = FALSE]
  assert_equal(rownames(beta$clr), rownames(bmd), "BMI23 beta sample order")
  d <- vegan::vegdist(beta$clr, method = "euclidean")
  control <- permute::how(nperm = config$models$diversity_permutations, blocks = factor(bmd$study_name))
  fit <- vegan::adonis2(d ~ analysis_group_BMI23_adapted + age_numeric + sex, data = bmd, permutations = control)
  dispersion <- vegan::betadisper(d, bmd$analysis_group_BMI23_adapted)
  disp_anova <- as.data.frame(anova(dispersion))
  disp_perm <- as.data.frame(vegan::permutest(dispersion, permutations = config$models$diversity_permutations)$tab)
  # Preserve historical row labels in an ordinary explicit column.
  add_terms <- function(x) { x <- as.data.frame(x); data.frame(term = rownames(x), x, check.names = FALSE, row.names = NULL) }
  write(add_terms(fit), "BMI23_Beta_Aitchison_PERMANOVA_Blocked.csv")
  write(add_terms(disp_anova), "BMI23_Beta_Aitchison_PERMDISP.csv")
  write(add_terms(disp_perm), "BMI23_Beta_Aitchison_PERMDISP_Permutation.csv")
  corrected <- sva::ComBat(dat = t(beta$clr), batch = bmd$study_name, mod = stats::model.matrix(~ analysis_group_BMI23_adapted + age_numeric + sex, data = bmd), par.prior = TRUE, prior.plots = FALSE)
  corrected_distance <- vegan::vegdist(t(corrected), method = "euclidean")
  exploratory <- as.data.frame(vegan::adonis2(corrected_distance ~ analysis_group_BMI23_adapted + age_numeric + sex, data = bmd, permutations = control))
  write(exploratory, "BMI23_Beta_ComBat_Exploratory_PERMANOVA.csv")
  excluded <- md[!md$sample_id %in% bmd$sample_id, c("sample_id", "study_name", "disease", "BMI", "analysis_group_BMI23_adapted"), drop = FALSE]
  excluded$exclusion_reason <- "Not retained by original cmultRepl/CLR compositional-data workflow"
  write(excluded, "BMI23_Beta_CLR_Excluded_Samples_Audit.csv")
  write(data.frame(input_samples = nrow(md), retained_CLR_samples = nrow(bmd), input_features = nrow(data$species), retained_CLR_features = ncol(beta$clr)), "BMI23_Beta_CLR_Retention_Audit.csv")
  invisible(NULL)
}
