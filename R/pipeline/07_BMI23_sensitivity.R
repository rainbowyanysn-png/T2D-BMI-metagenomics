# Cohort-adapted BMI23 species primary models.
run_bmi23_species <- function(data, run, config) {
  require_project_packages(c("Maaslin2", "tibble", "dplyr")); check_frozen_versions(config)
  specifications <- list(
    NT2D_vs_NC_BMI23 = list(pair = c("Non_Overweight_Control", "Non_Overweight_T2D"), tested = 192L, significant = 28L, samples = 366L),
    OT2D_vs_OC_BMI23 = list(pair = c("Overweight_Control", "Overweight_T2D"), tested = 190L, significant = 14L, samples = 820L),
    OC_vs_NC_BMI23 = list(pair = c("Non_Overweight_Control", "Overweight_Control"), tested = 193L, significant = 7L, samples = 426L),
    OT2D_vs_NT2D_BMI23 = list(pair = c("Non_Overweight_T2D", "Overweight_T2D"), tested = 190L, significant = 3L, samples = 760L))
  original <- medication_model_covariate(data$metadata)
  results <- list()
  for (id in names(specifications)) {
    specification <- specifications[[id]]; pair <- specification$pair
    metadata <- as.data.frame(original[original$BMI23_group %in% pair, , drop = FALSE])
    metadata$analysis_group_BMI23_adapted <- factor(as.character(metadata$BMI23_group), levels = pair)
    metadata$age_numeric <- as.numeric(as.character(metadata$age))
    metadata$sex <- factor(as.character(metadata$sex)); metadata$study_name <- factor(as.character(metadata$study_name))
    fixed <- c("analysis_group_BMI23_adapted", "age_numeric", "sex", "study_name")
    references <- c(paste("analysis_group_BMI23_adapted", pair[1], sep = ","), categorical_references(metadata, fixed))
    if (id == "OT2D_vs_NT2D_BMI23") {
      metadata$antidiabetes_group <- as.character(metadata$antidiabetes_group)
      assert_true(length(unique(na.omit(metadata$antidiabetes_group))) > 1L, "BMI23 medication categories missing")
      fixed <- c(fixed, "antidiabetes_group")
      reference <- medication_model_reference(metadata)
      references <- c(references, paste("antidiabetes_group", reference, sep = ","))
    }
    assert_equal(nrow(metadata), specification$samples, paste(id, "samples"))
    if (id == "OT2D_vs_NT2D_BMI23") assert_true(identical(as.integer(table(metadata$analysis_group_BMI23_adapted)), c(105L,655L)), "BMI23 T2D sample groups changed")
    folder <- dirname(path_in_run(run, "models/07_BMI23_sensitivity", id, "MaAsLin2", ".keep"))
    fit <- audited_maaslin2(input_data = t(data$species[, metadata$sample_id, drop = FALSE]), input_metadata = tibble::column_to_rownames(metadata, "sample_id"), output = folder, fixed_effects = fixed, random_effects = NULL, reference = references, normalization = "TSS", transform = "LOG", min_prevalence = .10, standardize = config$models$bmi23_maaslin_standardize, plot_heatmap = FALSE, plot_scatter = FALSE)
    target <- fit$results[fit$results$metadata == "analysis_group_BMI23_adapted" & fit$results$value == pair[2], , drop = FALSE]
    target$qval_BH_target <- stats::p.adjust(target$pval, method = "BH")
    target$formal_q <- target$qval; target$q_family <- "MaAsLin2 full model qval"; target$q_family_size <- sum(!is.na(fit$results$pval))
    target$comparison_direction <- paste(pair[2], "vs", pair[1])
    target$model_id <- paste(id, if (id == "OT2D_vs_NT2D_BMI23") "medication_adjusted" else "primary", sep = "_")
    safe_write_csv(target, path_in_run(run, "tables/07_BMI23_sensitivity", paste0(id, "_all.csv")))
    assert_equal(nrow(target), specification$tested, paste(id, "tested taxa"))
    assert_equal(sum(target$qval < .05), specification$significant, paste(id, "formal discoveries"))
    results[[id]] <- target
  }
  results
}

run_bmi23_pathways <- function(data, run, config) {
  require_project_packages(c("DESeq2", "metafor", "dplyr"))
  original_levels <- c("Non_Overweight_Control", "Non_Overweight_T2D", "Overweight_Control", "Overweight_T2D")
  assert_equal(unlist(config$models$bmi23_pathway_levels), original_levels, "Approved original BMI23 DESeq2 levels")
  specifications <- list(
    Lean_T2D_vs_Lean_Control = c("Non_Overweight_T2D", "Non_Overweight_Control"),
    Overweight_T2D_vs_Overweight_Control = c("Overweight_T2D", "Overweight_Control"),
    Overweight_vs_Lean_Controls = c("Overweight_Control", "Non_Overweight_Control"),
    Overweight_vs_Lean_T2D = c("Overweight_T2D", "Non_Overweight_T2D"))
  n_tested <- c(5571L, 5573L, 5573L, 5571L); discoveries <- c(6L, 33L, 0L, 30L)
  cohort_results <- setNames(vector("list", length(specifications)), names(specifications))
  for (study in unique(as.character(data$metadata$study_name))) {
    input <- data$pathway_cohorts[[study]]
    assert_true(!is.null(input), paste("Missing original cohort pathway input", study))
    ids <- rownames(input$metadata)
    md <- as.data.frame(data$metadata[match(ids, data$metadata$sample_id), , drop = FALSE])
    assert_equal(md$sample_id, colnames(input$matrix), paste(study, "BMI23 pathway sample order"))
    md$analysis_group_BMI23_adapted <- droplevels(factor(as.character(md$BMI23_group), levels = original_levels))
    rownames(md) <- md$sample_id
    yy <- input$matrix; yy[is.na(yy)] <- 0
    yy <- if (max(yy, na.rm = TRUE) <= 1) round(yy * 1e6) else round(yy)
    keep <- rowSums(yy) > 0 & rowSums(yy > 0) >= ceiling(ncol(yy) * config$models$pathway_cohort_prevalence) & apply(yy, 1, stats::var) > 0
    yy <- yy[keep, , drop = FALSE]
    assert_true(nrow(yy) >= config$models$pathway_min_features, paste("BMI23 pathway QC", study))
    safe_write_csv(data.frame(feature = rownames(yy), study_name = study), path_in_run(run, "models/07_BMI23_pathways", paste0(study, "_retained_features.csv")))
    safe_write_csv(data.frame(sample_id = md$sample_id, group = md$analysis_group_BMI23_adapted), path_in_run(run, "models/07_BMI23_pathways", paste0(study, "_sample_groups.csv")))
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = yy, colData = md, design = ~ analysis_group_BMI23_adapted)
    dds <- DESeq2::estimateSizeFactors(dds, type = config$models$pathway_size_factors)
    dds <- DESeq2::estimateDispersions(dds, fitType = config$models$pathway_dispersion_fit, quiet = TRUE)
    dds <- DESeq2::nbinomWaldTest(dds, quiet = TRUE)
    for (id in names(specifications)) {
      pair <- specifications[[id]]
      if (!all(pair %in% levels(md$analysis_group_BMI23_adapted))) next
      rr <- as.data.frame(DESeq2::results(dds, contrast = c("analysis_group_BMI23_adapted", pair), cooksCutoff = FALSE))
      rr$pathway <- rownames(rr); rr$study_name <- study
      rr <- rr[, c("pathway", "study_name", "log2FoldChange", "lfcSE", "pvalue", "padj")]
      cohort_results[[id]][[study]] <- rr
      safe_write_csv(rr, path_in_run(run, "models/07_BMI23_pathways", paste0("BMI23_", id, "_", study, "_DESeq2.csv")))
    }
  }
  output <- list()
  for (i in seq_along(specifications)) {
    id <- names(specifications)[i]; all <- dplyr::bind_rows(cohort_results[[id]])
    failures <- list()
    pooled <- lapply(split(all, all$pathway), function(d) {
      d <- d[is.finite(d$log2FoldChange) & is.finite(d$lfcSE) & d$lfcSE > 0, , drop = FALSE]
      if (nrow(d) < 2L) return(NULL)
      fit <- tryCatch(metafor::rma(yi = log2FoldChange, sei = lfcSE, method = config$models$pathway_meta_method, data = d), error = function(e) {failures[[d$pathway[1]]] <<- data.frame(pathway = d$pathway[1], error = conditionMessage(e)); NULL})
      if (is.null(fit)) return(NULL)
      data.frame(pathway = d$pathway[1], k = fit$k, combined_log2FC = unname(fit$b), se = fit$se, z = fit$zval, p_value = fit$pval, tau2 = fit$tau2, I2 = fit$I2)
    })
    result <- dplyr::bind_rows(pooled); result$padj <- stats::p.adjust(result$p_value, method = "BH")
    result <- result[stringi::stri_order(result$pathway, locale = config$output$historical_collation), , drop = FALSE]
    result$ci_lower <- result$combined_log2FC - 1.96 * result$se; result$ci_upper <- result$combined_log2FC + 1.96 * result$se
    result$q_family <- paste("All successful BMI23 REML pathways", id); result$q_family_size <- nrow(result)
    result$comparison_direction <- paste(specifications[[id]][1], "vs", specifications[[id]][2])
    result$model_id <- paste0("BMI23_adapted_DESeq2_REML_", id)
    safe_write_csv(result, path_in_run(run, "tables/07_BMI23_pathways", paste0("BMI23_", id, "_all.csv")))
    if (length(failures)) safe_write_csv(dplyr::bind_rows(failures), path_in_run(run, "logs", paste0("BMI23_", id, "_REML_failures.csv")))
    assert_equal(nrow(result), n_tested[i], paste(id, "BMI23 tested pathways"))
    assert_equal(sum(result$padj < config$project$alpha), discoveries[i], paste(id, "BMI23 discoveries"))
    output[[id]] <- result
  }
  output
}

run_bmi23_interactions <- function(data, run, config) {
  check_frozen_versions(config)
  metadata <- as.data.frame(data$metadata)
  metadata$sex <- factor(as.character(metadata$sex)); metadata$study_name <- factor(as.character(metadata$study_name))
  metadata$age_numeric <- as.numeric(as.character(metadata$age))
  metadata$T2D_binary <- as.integer(metadata$disease == "T2D")
  metadata$BMI_adapted_binary <- as.integer(metadata$BMI >= metadata$BMI23_threshold)
  metadata$T2D_x_BMI_adapted <- metadata$T2D_binary * metadata$BMI_adapted_binary
  fixed <- c("T2D_binary", "BMI_adapted_binary", "T2D_x_BMI_adapted", "age_numeric", "sex", "study_name")
  results <- list()
  for (domain in c("species", "pathways")) {
    folder <- dirname(path_in_run(run, "models/07_BMI23_interactions", domain, "MaAsLin2", ".keep"))
    fit <- audited_maaslin2(input_data = t(data[[domain]]), input_metadata = tibble::column_to_rownames(metadata, "sample_id"), output = folder, fixed_effects = fixed, random_effects = NULL, reference = categorical_references(metadata, fixed), normalization = "TSS", transform = "LOG", min_prevalence = 0, standardize = config$models$interaction_maaslin_standardize, plot_heatmap = FALSE, plot_scatter = FALSE)
    # The original source's value == "1" filter was corrected in the frozen extraction.
    target <- fit$results[fit$results$metadata == "T2D_x_BMI_adapted", , drop = FALSE]
    target$qval_BH_target <- stats::p.adjust(target$pval, method = "BH")
    target$formal_q <- target$qval_BH_target
    target$q_family <- paste("All BMI23 target interaction terms within", domain)
    target$q_family_size <- nrow(target)
    safe_write_csv(target, path_in_run(run, "tables/07_BMI23_interactions", paste0("BMI23_", domain, "_all.csv")))
    assert_equal(nrow(target), if (domain == "species") 200L else 5945L, paste("BMI23 interaction", domain, "tested features"))
    assert_equal(sum(target$formal_q < .05), 0L, paste("BMI23 interaction", domain, "discoveries"))
    results[[domain]] <- target
  }
  results
}
