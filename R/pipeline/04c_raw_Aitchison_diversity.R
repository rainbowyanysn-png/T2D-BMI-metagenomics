run_raw_aitchison_diversity <- function(data, run, config) {
  require_project_packages(c("vegan", "permute", "zCompositions", "compositions"))
  specification <- config$raw_aitchison
  stopifnot(isTRUE(specification$approved), specification$by == "margin",
            specification$permutations == 999L)
  beta <- historical_beta_clr(data)
  md <- beta$metadata
  md$analysis_group <- factor(as.character(md$BMI25_group), levels = analysis_group_levels)
  md$disease <- factor(as.character(md$disease), levels = c("control", "T2D"))
  md$study_name <- factor(md$study_name)
  md$sex <- factor(md$sex)
  stopifnot(nrow(md) == 1182L, ncol(beta$clr) == 138L,
            identical(rownames(beta$clr), md$sample_id),
            !anyNA(md[, c("age", "sex", "BMI", "analysis_group", "disease", "study_name")]))
  distance <- stats::dist(beta$clr)
  write <- function(x, name) safe_write_csv(x,
    path_in_run(run, "tables/04c_raw_Aitchison", paste0("BMI25_raw_Aitchison_", name, ".csv")))
  provenance <- function(x, formula, restriction, seed, target) {
    data.frame(term = rownames(x), as.data.frame(x), model_formula = formula, target_term = target,
      distance = "Uncorrected CZM/CLR Aitchison", permutation_restriction = restriction,
      permutation_seed = seed, n_permutations = specification$permutations,
      n_samples = nrow(md), test_method = "Marginal PERMANOVA", check.names = FALSE, row.names = NULL)
  }
  write(md, "sample_manifest")
  write(data.frame(feature = colnames(beta$clr)), "feature_manifest")
  write(data$metadata[!data$metadata$sample_id %in% md$sample_id, , drop = FALSE], "excluded_samples")
  ordination <- stats::cmdscale(distance, k = 2, eig = TRUE)
  fraction <- ordination$eig[1:2] / sum(ordination$eig[ordination$eig > 0])
  coordinates <- data.frame(sample_id = rownames(ordination$points),
    PC1 = ordination$points[, 1], PC2 = ordination$points[, 2],
    PC1_fraction = fraction[1], PC2_fraction = fraction[2])
  write(coordinates, "PCoA_coordinates")
  set.seed(specification$seed)
  control <- permute::how(blocks = md$study_name, nperm = specification$permutations)
  total <- vegan::adonis2(stats::as.formula(specification$total_formula), data = md,
    permutations = control, by = "margin")
  write(provenance(total, specification$total_formula, "Within study_name", specification$seed, "disease"),
        "PERMANOVA_total_blocked")
  set.seed(specification$seed + 1L)
  four <- vegan::adonis2(stats::as.formula(specification$four_formula), data = md,
    permutations = control, by = "margin")
  write(provenance(four, specification$four_formula, "Within study_name", specification$seed + 1L, "analysis_group"),
        "PERMANOVA_four_groups_blocked")
  # Cohort-specific tests have no study term or study blocking within one cohort.
  cohort <- lapply(seq_along(levels(md$study_name)), function(i) {
    study <- levels(md$study_name)[i]
    one <- md[md$study_name == study, , drop = FALSE]
    subset_distance <- stats::as.dist(as.matrix(distance)[one$sample_id, one$sample_id])
    nuisance <- c("age", "sex")
    estimable <- nuisance[vapply(one[nuisance], function(x) length(unique(x)) > 1L, logical(1))]
    formula <- stats::as.formula(paste("subset_distance ~", paste(c(estimable, "analysis_group"), collapse = " + ")))
    set.seed(specification$seed + 10L + i)
    fit <- vegan::adonis2(formula, data = one, permutations = specification$permutations, by = "margin")
    stopifnot("analysis_group" %in% rownames(fit))
    data.frame(study_name = study, R2 = fit["analysis_group", "R2"], F = fit["analysis_group", "F"],
      p_value = fit["analysis_group", "Pr(>F)"], n_samples = nrow(one),
      model_formula = paste(deparse(formula), collapse = " "),
      nonvarying_covariates = paste(setdiff(nuisance, estimable), collapse = ";"),
      permutation_seed = specification$seed + 10L + i, n_permutations = specification$permutations)
  })
  write(do.call(rbind, cohort), "PERMANOVA_by_cohort")
  for (group in c("disease", "analysis_group")) {
    set.seed(specification$seed + if (group == "disease") 100L else 101L)
    dispersion <- vegan::betadisper(distance, md[[group]], type = "median")
    perm <- as.data.frame(vegan::permutest(dispersion, permutations = specification$permutations)$tab)
    output <- data.frame(term = rownames(perm), perm, group_variable = group,
      n_samples = nrow(md), type = "Spatial median", permutation_restriction = "Unrestricted",
      n_permutations = specification$permutations, check.names = FALSE, row.names = NULL)
    write(output, paste0("PERMDISP_", if (group == "disease") "disease" else "four_groups"))
  }
}
