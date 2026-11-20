run_primary_bmi25_species <- function(data, run, config, expected) {
  require_project_packages(c("Maaslin2", "tibble"))
  check_frozen_versions(config)
  set.seed(config$models$primary_seed)
  specifications <- list(
    Total_T2D_vs_Control = list(reference = "control", comparison = "T2D", variable = "disease", group_column = "disease", fixed = c("disease", "BMI", "age", "sex", "study_name"), expected = expected$primary_bmi25_species$total_t2d_vs_control),
    NT2D_vs_NC_BMI25 = list(reference = "Non_Overweight_Control", comparison = "Non_Overweight_T2D", variable = "Group", group_column = "BMI25_group", fixed = c("Group", "age", "sex", "study_name"), expected = expected$primary_bmi25_species$nt2d_vs_nc),
    OT2D_vs_OC_BMI25 = list(reference = "Overweight_Control", comparison = "Overweight_T2D", variable = "Group", group_column = "BMI25_group", fixed = c("Group", "age", "sex", "study_name"), expected = expected$primary_bmi25_species$ot2d_vs_oc),
    OC_vs_NC_BMI25 = list(reference = "Non_Overweight_Control", comparison = "Overweight_Control", variable = "Group", group_column = "BMI25_group", fixed = c("Group", "age", "sex", "study_name"), expected = expected$primary_bmi25_species$oc_vs_nc),
    OT2D_vs_NT2D_BMI25 = list(reference = "Non_Overweight_T2D", comparison = "Overweight_T2D", variable = "Group", group_column = "BMI25_group", fixed = c("Group", "age", "sex", "study_name", "antidiabetes_group"), expected = expected$primary_bmi25_species$ot2d_vs_nt2d)
  )
  results <- list()
  for (id in names(specifications)) {
    spec <- specifications[[id]]
    metadata <- data$metadata
    if (spec$variable == "disease") {
      metadata <- metadata[as.character(metadata$disease) %in% c(spec$reference, spec$comparison), , drop = FALSE]
      metadata$disease <- factor(as.character(metadata$disease), levels = c(spec$reference, spec$comparison))
    } else {
      metadata <- metadata[as.character(metadata[[spec$group_column]]) %in% c(spec$reference, spec$comparison), , drop = FALSE]
      metadata$Group <- factor(as.character(metadata[[spec$group_column]]), levels = c(spec$reference, spec$comparison))
    }
    if (id == "OT2D_vs_NT2D_BMI25") {
      metadata <- medication_model_covariate(metadata)
      metadata$antidiabetes_group <- droplevels(metadata$antidiabetes_group)
      assert_true(nrow(metadata) == 760L && identical(as.integer(table(metadata$Group)), c(138L,622L)), "BMI25 T2D sample groups changed")
    }
    references <- c(paste(spec$variable, spec$reference, sep=","), categorical_references(metadata, spec$fixed))
    if (id == "OT2D_vs_NT2D_BMI25") references <- c(references,
      paste("antidiabetes_group", medication_model_reference(metadata), sep=","))
    assert_true(identical(colnames(data$species[, metadata$sample_id, drop = FALSE]), metadata$sample_id), paste(id, "sample alignment failed."))
    model_dir <- path_in_run(run, "models", "02_primary_BMI25_species", id, "MaAsLin2", ".keep")
    model_dir <- dirname(model_dir)
    fit <- audited_maaslin2(input_data = data$species[, metadata$sample_id, drop = FALSE], input_metadata = tibble::column_to_rownames(as.data.frame(metadata), "sample_id"), output = model_dir, fixed_effects = spec$fixed, reference = references, normalization = config$models$maaslin_normalization, transform = config$models$maaslin_transform, standardize = config$models$primary_maaslin_standardize, min_prevalence = config$filters$global_prevalence, plot_heatmap = FALSE, plot_scatter = FALSE)
    target <- fit$results[fit$results$metadata == spec$variable, , drop = FALSE]
    assert_equal(nrow(target), expected$primary_bmi25_tested[[id]], paste(id, "tested species"))
    target$qval_target_bh <- stats::p.adjust(target$pval, method="BH")
    target$formal_q <- target$qval
    target$q_family <- "MaAsLin2 full model qval"
    target$q_family_size <- sum(!is.na(fit$results$pval))
    target$n_target_tested <- nrow(target)
    target$model_id <- paste0("BMI25_", id, "_standardize_TRUE")
    target$reference_group <- spec$reference
    target$comparison_group <- spec$comparison
    target$medication_covariate <- id == "OT2D_vs_NT2D_BMI25"
    target$medication_reference <- if(id == "OT2D_vs_NT2D_BMI25") medication_model_reference(metadata) else NA_character_
    target$model_formula <- paste("abundance ~", paste(spec$fixed, collapse = " + "))
    target$comparison_id <- id
    target$comparison_direction <- paste(spec$comparison, "vs", spec$reference)
    target$direction <- ifelse(target$coef > 0, "Positive", ifelse(target$coef < 0, "Negative", "Zero"))
    target$formal_significant <- target$formal_q < config$project$alpha
    target <- target[if (id == "OT2D_vs_NT2D_BMI25") order(target$formal_q, target$pval) else order(target$formal_q), , drop = FALSE]
    assert_equal(sum(target$formal_significant), spec$expected, paste(id, "formal discoveries"))
    if (id == "OT2D_vs_NT2D_BMI25") {
      assert_true(all(target$N == 760L), "BMI25 model sample count")
      expected_taxa <- c("Roseburia.intestinalis", "Sellimonas.intestinalis", "Eisenbergiella.tayi", "Allisonella.histaminiformans", "Bifidobacterium.adolescentis")
      selected <- sub("^species[.]", "", target$feature[target$qval < .05])
      assert_true(setequal(selected, expected_taxa), "BMI25 adjusted significant-species identities differ")
      bif <- target[sub("^species[.]", "", target$feature) == "Bifidobacterium.adolescentis",,drop=FALSE]
      assert_true(nrow(bif)==1L && abs(bif$coef-1.8639516478)<1e-8 && abs(bif$qval-.0459046256)<1e-8, "B. adolescentis adjusted coefficient/qval differ")
      safe_write_csv(as.data.frame(table(metadata$Group, metadata$antidiabetes_group)),
        path_in_run(run,"tables/02_primary_BMI25_species","BMI25_T2D_medication_group_counts.csv"))
    }
    safe_write_csv(target, path_in_run(run, "tables", "02_primary_BMI25_species", paste0(id, "_all.csv")))
    results[[id]] <- target
  }
  summary <- data.frame(comparison_id = names(results), n_tested = vapply(results, nrow, integer(1)), n_formal_q_lt_005 = vapply(results, function(x) sum(x$formal_significant), integer(1)), stringsAsFactors = FALSE)
  safe_write_csv(summary, path_in_run(run, "tables", "02_primary_BMI25_species_summary.csv"))
  results
}
