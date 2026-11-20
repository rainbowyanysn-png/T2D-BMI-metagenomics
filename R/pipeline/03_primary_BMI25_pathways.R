# Reconstruct source Parts 6.2 and 7 from unfiltered original pathway abundances.
run_primary_bmi25_pathways <- function(data, run, config, expected) {
  require_project_packages(c("DESeq2", "metafor", "dplyr"))
  assert_true(!is.null(data$pathways_raw), "Unfiltered original pathway matrix is required.")
  specifications <- list(
    Total_T2D_vs_HC = c("disease", "T2D", "control"),
    NonOW_T2D_vs_NonOW_Ctrl = c("analysis_group", "Non_Overweight_T2D", "Non_Overweight_Control"),
    OW_T2D_vs_OW_Ctrl = c("analysis_group", "Overweight_T2D", "Overweight_Control"),
    OW_Ctrl_vs_NonOW_Ctrl = c("analysis_group", "Overweight_Control", "Non_Overweight_Control"),
    OW_T2D_vs_NonOW_T2D = c("analysis_group", "Overweight_T2D", "Non_Overweight_T2D"))
  studies <- unique(as.character(data$metadata$study_name))
  cohort_results <- list()
  for (study in studies) {
    assert_true(!is.null(data$pathway_cohorts[[study]]), paste("Original per-cohort pathway input missing:", study))
    metadata <- data$pathway_cohorts[[study]]$metadata
    metadata$disease <- factor(metadata$disease, levels = c("control", "T2D"))
    assert_true(nrow(metadata) >= 10 && length(unique(metadata$analysis_group)) >= 2, paste("Historical cohort QC changed:", study))
    counts <- data$pathway_cohorts[[study]]$matrix
    counts[is.na(counts)] <- 0
    counts <- if (max(counts, na.rm = TRUE) <= 1) round(counts * 1e6) else round(counts)
    counts <- counts[rowSums(counts) > 0, , drop = FALSE]
    counts <- counts[rowSums(counts > 0) >= ceiling(ncol(counts) * config$models$pathway_cohort_prevalence), , drop = FALSE]
    counts <- counts[apply(counts, 1, stats::var) > 0, , drop = FALSE]
    assert_true(nrow(counts) >= config$models$pathway_min_features, paste("Insufficient pathways:", study))
    safe_write_csv(data.frame(pathway = rownames(counts), study_name = study), path_in_run(run, "models", "03_primary_BMI25_pathways", paste0(study, "_retained_features.csv")))
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ analysis_group)
    dds <- DESeq2::estimateSizeFactors(dds, type = config$models$pathway_size_factors)
    dds <- DESeq2::estimateDispersions(dds, fitType = config$models$pathway_dispersion_fit)
    dds <- DESeq2::nbinomWaldTest(dds)
    for (id in names(specifications)) {
      spec <- specifications[[id]]
      if (!all(spec[2:3] %in% as.character(metadata[[spec[1]]]))) next
      model <- dds
      if (spec[1] == "disease") {
        DESeq2::design(model) <- ~ disease
        model <- DESeq2::nbinomWaldTest(model)
      }
      result <- as.data.frame(DESeq2::results(model, contrast = spec, cooksCutoff = FALSE))
      result$pathway <- rownames(result); result$study <- study
      cohort_results[[id]][[study]] <- result
      safe_write_csv(result, path_in_run(run, "models", "03_primary_BMI25_pathways", paste0("BMI25_", id, "_", study, "_DESeq2.csv")))
    }
  }
  assert_equal(length(studies), 5L, "Pathway input cohort count")
  output <- list()
  for (id in names(specifications)) {
    combined <- dplyr::bind_rows(cohort_results[[id]])
    paths <- combined |> dplyr::count(pathway) |> dplyr::filter(n > 1) |> dplyr::pull(pathway)
    pooled <- list(); failures <- list()
    for (pathway in paths) {
      one <- combined[combined$pathway == pathway & !is.na(combined$lfcSE) & combined$lfcSE > 0 & is.finite(combined$log2FoldChange), , drop = FALSE]
      if (nrow(one) < 2) next
      fit <- tryCatch(metafor::rma(yi = log2FoldChange, sei = lfcSE, data = one, method = config$models$pathway_meta_method), error = function(e) {failures[[pathway]] <<- data.frame(pathway = pathway, error = conditionMessage(e)); NULL})
      if (is.null(fit)) next
      pooled[[pathway]] <- data.frame(pathway = pathway, combined_log2FC = as.numeric(fit$beta), p_value = as.numeric(fit$pval), n_studies = as.numeric(fit$k), Heterogeneity_I2 = as.numeric(fit$I2), stderr = as.numeric(fit$se), ci_lower = as.numeric(fit$ci.lb), ci_upper = as.numeric(fit$ci.ub))
    }
    result <- dplyr::bind_rows(pooled)
    result$padj <- stats::p.adjust(result$p_value, method = "BH")
    result <- result[order(result$padj), , drop = FALSE]
    result$q_family <- paste("All successful REML pathways within BMI25", id)
    result$q_family_size <- nrow(result)
    result$comparison_direction <- paste(specifications[[id]][2], "vs", specifications[[id]][3])
    result$model_id <- paste0("BMI25_DESeq2_REML_", id)
    safe_write_csv(result, path_in_run(run, "tables", "03_primary_BMI25_pathways", paste0("BMI25_", id, "_all.csv")))
    if (length(failures)) safe_write_csv(dplyr::bind_rows(failures), path_in_run(run, "logs", paste0("BMI25_", id, "_REML_failures.csv")))
    assert_equal(nrow(result), expected$primary_bmi25_pathways[[id]]$n_tested, paste(id, "pooled pathway family size"))
    assert_equal(sum(result$padj < config$project$alpha), expected$primary_bmi25_pathways[[id]]$n_significant, paste(id, "pathway discoveries"))
    output[[id]] <- result
  }
  output
}
