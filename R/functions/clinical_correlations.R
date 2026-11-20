clinical_feature_display <- function(feature) {
  x <- sub("^(species|genus)[.]", "", feature)
  x <- trimws(gsub("[._]+", " ", x))
  gsub(" sp ", " sp. ", x, fixed = TRUE)
}

clinical_correlation_data <- function(data, main, config) {
  cfg <- config$clinical_correlations
  stopifnot(cfg$n_features == 20L, cfg$family_size == 60L,
            cfg$method == "spearman", identical(cfg$exact, FALSE))
  md <- as.data.frame(data$metadata)
  stopifnot(nrow(md) == 1186L, nrow(data$species) == 200L,
            !anyDuplicated(md$sample_id), identical(colnames(data$species), md$sample_id),
            sum(md$disease == "T2D") == 760L, sum(md$disease == "control") == 426L,
            setequal(as.character(md$study_name), unlist(config$cohorts$names)))
  raw_sex <- as.character(md$sex)
  sex <- tolower(trimws(raw_sex))
  assert_true(!anyNA(sex) && setequal(unique(sex), c("female", "male")),
              "Sex must contain exactly recognized female and male levels with no missing values; see sex coding audit.")
  clinical <- data.frame(BMI = md$BMI, age = md$age, sex = as.integer(sex == "male"))
  stopifnot(is.numeric(clinical$BMI), is.numeric(clinical$age),
            all(is.finite(clinical$BMI)), all(is.finite(clinical$age)))
  assert_columns(main, c("feature", "coef", "pval", "qval"), "Current BMI25 total MaAsLin2 result")
  stopifnot(nrow(main) == 200L, !anyDuplicated(main$feature),
            sum(main$qval < config$project$alpha) == 69L,
            all(is.finite(main$qval)), all(is.finite(main$pval)))
  if ("formal_q" %in% names(main)) stopifnot(identical(main$formal_q, main$qval))
  main <- main[order(main$qval, main$pval, main$feature, method = "radix"), , drop = FALSE]
  selected <- main[seq_len(cfg$n_features), c("feature", "coef", "pval", "qval")]
  names(selected) <- c("feature", "main_model_coef", "main_model_pval", "main_model_qval")
  selected$feature_rank <- seq_len(nrow(selected))
  selected$display_name <- clinical_feature_display(selected$feature)
  selected <- selected[, c("feature_rank", "feature", "display_name",
                           "main_model_coef", "main_model_pval", "main_model_qval")]
  stopifnot("genus.Streptococcus" %in% selected$feature,
            all(selected$feature %in% rownames(data$species)))
  tests <- lapply(seq_len(nrow(selected)), function(i) {
    x <- as.numeric(data$species[selected$feature[i], ])
    stopifnot(all(is.finite(x)))
    do.call(rbind, lapply(c("BMI", "age", "sex"), function(trait) {
      y <- clinical[[trait]]
      complete <- stats::complete.cases(x, y)
      stopifnot(sum(complete) >= 3L, length(unique(x[complete])) > 1L,
                length(unique(y[complete])) > 1L)
      test <- stats::cor.test(x[complete], y[complete], method = "spearman", exact = FALSE)
      data.frame(selected[i, ], clinical_trait = trait, n_complete = sum(complete),
        spearman_rho = unname(test$estimate), p_value = test$p.value, row.names = NULL)
    }))
  })
  result <- do.call(rbind, tests)
  stopifnot(nrow(result) == 60L, !anyDuplicated(result[, c("feature", "clinical_trait")]),
            all(result$n_complete[result$clinical_trait %in% c("BMI", "age")] == 1186L),
            all(is.finite(result$p_value)))
  result$q_bh_60 <- stats::p.adjust(result$p_value, method = "BH")
  result$fdr_significant <- result$q_bh_60 < config$project$alpha
  result$direction <- ifelse(result$spearman_rho > 0, "Positive",
                             ifelse(result$spearman_rho < 0, "Negative", "Zero"))
  result$analysis_scope <- "Exploratory correlations between selected T2D-associated microbial features and clinical variables."
  result$correlation_method <- "Spearman, exact FALSE; pairwise complete observations; no covariate adjustment"
  result$q_family <- "All 20 features x BMI, age and sex"
  result$q_family_size <- 60L
  stopifnot(identical(result$q_bh_60, stats::p.adjust(result$p_value, "BH")))
  # YAML stores each benchmark as a record; construct a conventional table.
  benchmarks <- do.call(rbind, lapply(cfg$benchmarks, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
  rows <- match(paste(benchmarks$feature, benchmarks$trait),
                paste(result$feature, result$clinical_trait))
  stopifnot(!anyNA(rows))
  observed <- result[rows, ]
  checks <- list()
  for (field in c("rho", "p", "q")) {
    actual <- observed[[switch(field, rho = "spearman_rho", p = "p_value", q = "q_bh_60")]]
    expected <- benchmarks[[field]]
    bound <- if (field == "rho") rep(cfg$rho_absolute_tolerance, length(actual)) else
      abs(expected) * cfg$probability_relative_tolerance
    checks[[field]] <- data.frame(feature = benchmarks$feature, clinical_trait = benchmarks$trait,
      field = field, expected = expected, observed = actual, absolute_difference = abs(actual - expected),
      allowed_difference = bound, pass = abs(actual - expected) <= bound)
  }
  list(selected = selected, all = result, benchmark_checks = do.call(rbind, checks))
}

clinical_export_workbook <- function(run, config, directory = "tables/02b_exploratory_clinical_correlations") {
  check_export_package_versions(run$root, "openxlsx")
  x <- read_current_result(run, "tables/02b_exploratory_clinical_correlations/Clinical_correlations_all_60_tests.csv")
  stopifnot(nrow(x) == 60L)
  display <- x[, c("feature_rank", "display_name", "clinical_trait", "n_complete",
                   "spearman_rho", "p_value", "q_bh_60", "fdr_significant")]
  names(display) <- c("Rank", "Species", "Clinical trait", "N", "Spearman rho",
                      "P value", "BH-adjusted q value", "FDR significant")
  display[["Clinical trait"]] <- c(BMI = "BMI", age = "Age", sex = "Sex")[display[["Clinical trait"]]]
  style <- read_presentation_config(run$root)
  wb <- openxlsx::createWorkbook(creator = "T2D BMI exploratory correlation pipeline")
  openxlsx::addWorksheet(wb, "Clinical correlations", gridLines = FALSE)
  openxlsx::addWorksheet(wb, "Notes", gridLines = FALSE)
  head_style <- openxlsx::createStyle(fontName = style$font_family, fontSize = 11,
    textDecoration = "bold", fgFill = "#F2F2F2", wrapText = TRUE)
  openxlsx::writeData(wb, "Clinical correlations", display, withFilter = TRUE, headerStyle = head_style)
  openxlsx::addStyle(wb, "Clinical correlations",
    openxlsx::createStyle(fontName = style$font_family, fontSize = 11),
    rows = 2:61, cols = 1:8, gridExpand = TRUE)
  openxlsx::addStyle(wb, "Clinical correlations", openxlsx::createStyle(textDecoration = "italic"),
                     rows = 2:61, cols = 2, gridExpand = TRUE, stack = TRUE)
  openxlsx::addStyle(wb, "Clinical correlations", openxlsx::createStyle(numFmt = "0.000"),
                     rows = 2:61, cols = 5, gridExpand = TRUE, stack = TRUE)
  for (column in c(6L, 7L)) for (i in seq_len(nrow(display))) {
    value <- display[[column]][i]
    format <- if (value < .001 && value > 0) "0.00E+00" else {
      decimals <- if (value == 0) 0L else max(0L, 2L - floor(log10(abs(value))))
      if (decimals == 0L) "0" else paste0("0.", paste(rep("0", decimals), collapse = ""))
    }
    openxlsx::addStyle(wb, "Clinical correlations", openxlsx::createStyle(numFmt = format),
                       rows = i + 1L, cols = column, stack = TRUE)
  }
  openxlsx::freezePane(wb, "Clinical correlations", firstRow = TRUE)
  openxlsx::setColWidths(wb, "Clinical correlations", cols = 1:8,
                         widths = c(8, 38, 16, 10, 18, 18, 23, 18))
  openxlsx::setRowHeights(wb, "Clinical correlations", rows = 1, heights = 32)
  notes <- c(
    "Exploratory Spearman correlations; no covariate adjustment.",
    "Features are the first 20 current BMI25 total T2D MaAsLin2 results sorted by original full-model qval, pval and feature ID. Genus features are retained.",
    "Original saved relative abundance is used without further normalization or transformation.",
    "Sex is Female=0 and Male=1. Positive rho indicates higher ranked abundance in males; negative rho indicates higher ranked abundance in females.",
    "All 60 tests are corrected together using BH. Formal significance uses q_bh_60 <0.05.",
    "N is the number of pairwise complete observations. No significant correlation or absence of correlation establishes independence from host factors or causality.",
    paste("Saved statistical input SHA-256:", sha256_file(current_result_path(run,
      "tables/02b_exploratory_clinical_correlations/Clinical_correlations_all_60_tests.csv")))
  )
  openxlsx::writeData(wb, "Notes", data.frame(Note = notes), headerStyle = head_style)
  openxlsx::setColWidths(wb, "Notes", cols = 1, widths = 110)
  openxlsx::addStyle(wb, "Notes", openxlsx::createStyle(fontName = style$font_family, wrapText = TRUE),
                     rows = 2:8, cols = 1, gridExpand = TRUE)
  openxlsx::setRowHeights(wb, "Notes", rows = 2:8, heights = 45)
  target <- path_in_run(run, directory, "Supplementary_Table_S5_Clinical_correlations.xlsx")
  assert_true(!file.exists(target), "Refusing to overwrite clinical correlation workbook")
  openxlsx::saveWorkbook(wb, target, overwrite = FALSE)
  back <- openxlsx::read.xlsx(target, sheet = "Clinical correlations",
    check.names = FALSE, sep.names = " ", skipEmptyRows = FALSE)
  stopifnot(nrow(back) == 60L, identical(names(back), names(display)),
            identical(as.character(back$Species), as.character(display$Species)),
            all(back$N == display$N))
  for (field in c("Spearman rho", "P value", "BH-adjusted q value")) {
    stopifnot(all(abs(back[[field]] - display[[field]]) <= 1e-12))
  }
  stopifnot(identical(back[["BH-adjusted q value"]] < .05, x$q_bh_60 < .05))
  invisible(target)
}

clinical_export_figure <- function(run, config) {
  require_project_packages("ggplot2")
  style <- read_presentation_config(run$root)
  config$output <- utils::modifyList(config$output, style$figure)
  x <- read_current_result(run, "tables/02b_exploratory_clinical_correlations/Clinical_correlations_all_60_tests.csv")
  stopifnot(nrow(x) == 60L, length(unique(x$feature)) == 20L,
            identical(x$q_bh_60 < .05, x$fdr_significant))
  selected <- x[!duplicated(x$feature), c("feature", "feature_rank", "display_name")]
  selected <- selected[order(selected$feature_rank), ]
  x$plot_label <- factor(x$display_name, levels = rev(selected$display_name))
  x$trait_display <- factor(c(BMI = "BMI", age = "Age", sex = "Sex")[x$clinical_trait],
                            levels = c("BMI", "Age", "Sex"))
  x$stars <- ifelse(x$q_bh_60 < .001, "***", ifelse(x$q_bh_60 < .01, "**",
                          ifelse(x$q_bh_60 < .05, "*", "")))
  limit <- max(.25, max(abs(x$spearman_rho)))
  plot <- ggplot2::ggplot(x, ggplot2::aes(trait_display, plot_label, fill = spearman_rho)) +
    ggplot2::geom_tile(colour = "white", linewidth = .35) +
    ggplot2::geom_text(ggplot2::aes(label = stars), family = style$font_family, size = 4) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#FFFFFF", high = "#B2182B",
      midpoint = 0, limits = c(-limit, limit), name = expression(Spearman~rho)) +
    ggplot2::labs(x = NULL, y = NULL,
      caption = "Significance based on BH-adjusted q values across 60 tests") +
    publication_theme(config) + ggplot2::theme(
      axis.text.y = ggplot2::element_text(face = "italic", size = 10),
      axis.text.x = ggplot2::element_text(face = "bold"),
      axis.line = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 11),
      plot.caption = ggplot2::element_text(size = 9, hjust = 0))
  id <- "02b_exploratory_clinical_correlations/Figure4F_species_clinical_correlation_heatmap"
  safe_write_csv(x, path_in_run(run, "plot_data", paste0(id, ".csv")))
  input <- current_result_path(run,
    "tables/02b_exploratory_clinical_correlations/Clinical_correlations_all_60_tests.csv")
  safe_write_csv(data.frame(input_file = input, sha256 = sha256_file(input),
    colour_column = "spearman_rho", significance_column = "q_bh_60",
    q_family_size = 60L, symmetric_colour_limit = limit),
    path_in_run(run, "plot_data", paste0(id, "_input_manifest.csv")))
  caption <- paste(
    "Exploratory correlations between selected T2D-associated microbial features and clinical variables.",
    "Spearman correlations use the saved relative abundance and pairwise complete observations, with no covariate adjustment.",
    "Rows follow current total-T2D MaAsLin2 full-model qval, pval and feature order.",
    "Stars use BH across all 60 tests. Sex: female=0 and male=1.",
    "This exploratory analysis cannot establish independence from host factors or causal associations.")
  save_publication_figure(plot, id, run, config, 7.5, 8.2, caption)
}
