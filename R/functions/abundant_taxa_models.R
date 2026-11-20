# Secondary/exploratory abundant-taxa analysis; independent of formal MaAsLin2 inference.
abundant_taxa_prepare <- function(data, run, config) {
  require_project_packages(c("curatedMetagenomicData", "SummarizedExperiment"))
  md <- as.data.frame(data$metadata)
  stopifnot(nrow(md) == 1186L, sum(md$disease == "T2D") == 760L,
            sum(md$disease == "control") == 426L, !anyDuplicated(md$sample_id),
            identical(colnames(data$species), md$sample_id), nrow(data$species) == 200L,
            identical(levels(md$BMI25_group), analysis_group_levels),
            setequal(as.character(md$study_name), unlist(config$cohorts$names)))
  retrieve <- prepare_frozen_resources(run, config)
  long <- retrieve(md, dataType = "relative_abundance", rownames = "long")
  short <- retrieve(md, dataType = "relative_abundance", rownames = "short")
  a <- SummarizedExperiment::assay(long)
  b <- SummarizedExperiment::assay(short)
  stopifnot(setequal(colnames(a), md$sample_id), setequal(colnames(b), md$sample_id))
  # mergeData builds sample/cohort-specific rowData and may include rows absent
  # from the package-wide lookup. Use the official annotations returned with
  # this exact assay rather than rematching against a separate global object.
  official <- as.data.frame(SummarizedExperiment::rowData(long))
  assert_columns(official, c("species", "genus", "family", "phylum"),
                 "Provider rowData for the returned long assay")
  j <- match(rownames(a), rownames(official))
  assert_true(!anyNA(j) && !anyDuplicated(rownames(a)) &&
                !anyDuplicated(rownames(official)),
              "Returned long assay and its official rowData do not align uniquely")
  provider <- asNamespace("curatedMetagenomicData")
  empty_fun <- get("taxonomyRankEmpty", provider, inherits = TRUE)
  label_fun <- get("getTaxonomyLabels", provider, inherits = TRUE)
  a <- a[, md$sample_id, drop = FALSE]
  b <- b[, md$sample_id, drop = FALSE]
  tolerance <- config$project$numeric_tolerance
  comparisons <- list()
  maps <- list()
  compare_values <- function(left, right, label) {
    stopifnot(identical(dim(left), dim(right)), identical(rownames(left), rownames(right)),
              identical(colnames(left), colnames(right)),
              all(is.finite(left)), all(is.finite(right)))
    delta <- abs(left - right)
    maximum <- max(delta)
    comparisons[[label]] <<- data.frame(comparison = label, n_features = nrow(left),
      n_samples = ncol(left), left_storage = typeof(left), right_storage = typeof(right),
      max_absolute_difference = maximum, tolerance = tolerance, pass = maximum <= tolerance)
    if (maximum > tolerance) {
      bad <- which(delta > tolerance, arr.ind = TRUE)
      bad <- head(bad, 100L)
      report <- data.frame(feature = rownames(left)[bad[, 1]], sample_id = colnames(left)[bad[, 2]],
        long_value = left[bad], short_value = right[bad], absolute_difference = delta[bad])
      safe_write_csv(report, path_in_run(run, "tables/04b_abundant_taxa",
        paste0("alignment_failure_", label, ".csv")))
      safe_write_csv(do.call(rbind, comparisons),
        path_in_run(run, "tables/04b_abundant_taxa", "abundance_alignment_checks.csv"))
      stop(paste("Abundance mismatch after cohort-specific official mapping:", label,
                 "maximum absolute difference", format(maximum, digits = 16)), call. = FALSE)
    }
  }
  # The provider assigns unique labels before merging cohorts. Reproduce that
  # exact scope: uniqueness suffixes generated on a merged object can differ.
  for (study in unique(as.character(md$study_name))) {
    ids <- md$sample_id[as.character(md$study_name) == study]
    one <- retrieve(md[as.character(md$study_name) == study, , drop = FALSE],
                    dataType = "relative_abundance", rownames = "long")
    raw <- SummarizedExperiment::assay(one)
    rowdata <- as.data.frame(SummarizedExperiment::rowData(one))
    stopifnot(identical(rownames(raw), rownames(rowdata)), setequal(colnames(raw), ids))
    empty <- empty_fun(one)
    labels <- as.character(label_fun(one[!empty, ]))
    match_short <- match(labels, rownames(b))
    stopifnot(length(labels) == sum(!empty), !anyDuplicated(labels), !anyNA(match_short))
    left <- raw[!empty, ids, drop = FALSE]
    rownames(left) <- labels
    compare_values(left, b[match_short, ids, drop = FALSE], study)
    map <- data.frame(study_name = study, taxonomic_string = rownames(raw),
      rowdata[, c("species", "genus", "family", "phylum")],
      provider_short_name = NA_character_, retained_in_short = !empty, row.names = NULL)
    map$provider_short_name[!empty] <- labels
    map$V3_short_name <- ifelse(empty, NA_character_, make.names(map$provider_short_name))
    stopifnot(all(grepl("(^|[|])p__", map$taxonomic_string)))
    map$phylum_string <- sub("[|;].*$", "", sub("^.*p__", "", map$taxonomic_string))
    stopifnot(!anyNA(map$phylum_string), all(nzchar(map$phylum_string)),
              !any(grepl("[[:cntrl:]]", map$phylum_string)))
    # The authoritative historical aggregation used the explicit p__ field.
    # Preserve rowData annotations separately, including disagreement/missingness.
    map$phylum_rowData <- map$phylum
    map$phylum_annotation_status <- ifelse(
      is.na(map$phylum_rowData) | !nzchar(trimws(map$phylum_rowData)),
      "rowData_missing",
      ifelse(map$phylum_rowData == map$phylum_string, "identical", "different"))
    map$phylum_source <- "Explicit p__ field of provider long taxonomy; original analysis rule"
    map$phylum <- map$phylum_string
    safe_write_csv(map[map$phylum_annotation_status != "identical", , drop = FALSE],
      path_in_run(run, "tables/04b_abundant_taxa",
        paste0("phylum_annotation_audit_", study, ".csv")))
    maps[[study]] <- map
  }
  taxonomy <- do.call(rbind, maps)
  # Phyla use every long-form feature, as in the original plotting analysis.
  stopifnot(all(grepl("(^|[|])p__", rownames(a))))
  phylum <- sub("[|;].*$", "", sub("^.*p__", "", rownames(a)))
  stopifnot(!anyNA(phylum), all(nzchar(phylum)),
            !any(grepl("[[:cntrl:]]", phylum)))
  short_names <- make.names(rownames(b))
  stopifnot(!anyDuplicated(short_names))
  ix <- match(rownames(data$species), short_names)
  stopifnot(!anyNA(ix))
  rownames(b) <- short_names
  compare_values(b[ix, , drop = FALSE], data$species, "global_200_panel")
  safe_write_csv(do.call(rbind, comparisons),
    path_in_run(run, "tables/04b_abundant_taxa", "abundance_alignment_checks.csv"))
  resources <- read.csv(file.path(run$root, config$source_data$resource_manifest))
  resources <- resources[resources$data_type == "relative_abundance" &
                           resources$study_name %in% unlist(config$cohorts$names), ]
  stopifnot(nrow(resources) == 5L)
  stopifnot(!anyDuplicated(resources$study_name))
  resource_row <- match(taxonomy$study_name, resources$study_name)
  stopifnot(!anyNA(resource_row))
  taxonomy$resource_id <- resources$resource_id[resource_row]
  taxonomy$data_version <- resources$data_version[resource_row]
  safe_write_csv(taxonomy, path_in_run(run, "tables/04b_abundant_taxa", "taxonomy_mapping.csv"))
  safe_write_csv(resources, path_in_run(run, "tables/04b_abundant_taxa", "taxonomy_resource_manifest.csv"))
  list(metadata = md, species = data$species,
       phylum = rowsum(a, phylum, reorder = TRUE))
}

abundant_taxa_fit <- function(input, config) {
  require_project_packages(c("lmerTest", "emmeans", "lme4"))
  cfg <- config$abundant_taxa
  stopifnot(cfg$pairwise_adjust == "none", cfg$formal_adjust == "BH",
            cfg$family_scope == "taxonomic_level", isTRUE(cfg$REML))
  formula <- stats::as.formula(cfg$formula)
  md <- input$metadata
  md$analysis_group <- factor(as.character(md$BMI25_group), levels = analysis_group_levels)
  md$sex <- factor(md$sex)
  md$study_name <- factor(md$study_name)
  selected <- list()
  samples <- list()
  coefficients <- list()
  comparisons <- list()
  diagnostics <- list()
  variance <- list()
  residuals <- list()
  for (rank in c("phylum", "species")) {
    matrix <- input[[rank]]
    means <- rowMeans(matrix)
    # Original top_n selection includes all ties; never silently break a cutoff tie.
    number <- if (rank == "phylum") 4L else 8L
    threshold <- sort(means, decreasing = TRUE)[number]
    ids <- sort(names(means)[means >= threshold], method = "radix")
    stopifnot(length(ids) == number, setequal(ids, unlist(cfg[[paste0("expected_", rank)]])))
    selected[[rank]] <- data.frame(feature = ids, taxonomic_level = rank,
      ranking_metric = "Overall arithmetic mean relative abundance", ranking_value = unname(means[ids]),
      ranking_scope = if (rank == "species") "200-feature global 10% prevalence panel" else "All original taxonomic features aggregated by phylum",
      in_global_200_panel = if (rank == "species") TRUE else NA,
      selection_rule = paste0("top_n(", number, ", mean), retaining ties; exact feature-count assertion"),
      display_order = seq_along(ids))
    for (id in ids) {
      x <- md
      x$Abundance <- as.numeric(matrix[id, ])
      stopifnot(all(is.finite(x$Abundance)), !anyNA(x[, c("analysis_group", "age", "sex", "study_name")]))
      warnings <- character()
      fit <- withCallingHandlers(lmerTest::lmer(formula, data = x, REML = TRUE),
        warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
      stopifnot(stats::nobs(fit) == 1186L)
      key <- paste(rank, id)
      context <- data.frame(feature = id, taxonomic_level = rank,
        analysis_scope = "Secondary/exploratory abundant-taxa analysis",
        model_formula = cfg$formula, transformation = "None; original percentage relative abundance",
        n_samples = 1186L, n_features = number, reference_group = analysis_group_levels[1],
        p_adjust_method = "BH on unadjusted pairwise P values",
        test_family = paste0(rank, ": all selected features x six pairwise comparisons"),
        family_size = number * 6L, REML = TRUE)
      repeat_context <- function(n) context[rep(1L, n), , drop = FALSE]
      samples[[key]] <- cbind(x[, c("sample_id", "study_name", "analysis_group", "age", "sex")],
        repeat_context(nrow(x)), raw_abundance = x$Abundance, model_abundance = x$Abundance)
      cf <- as.data.frame(coef(summary(fit)), check.names = FALSE)
      ci <- stats::qt(.975, cf$df) * cf[["Std. Error"]]
      coefficients[[key]] <- cbind(repeat_context(nrow(cf)), term = rownames(cf), cf,
        ci_lower = cf$Estimate - ci, ci_upper = cf$Estimate + ci,
        comparison_direction = "Fixed-effect coefficient versus stated reference",
        formal_significant = NA)
      marginal <- emmeans::emmeans(fit, ~ analysis_group)
      pair <- emmeans::contrast(marginal, method = "pairwise", adjust = "none")
      raw <- as.data.frame(summary(pair, adjust = "none"))
      tukey <- as.data.frame(summary(pair, adjust = "tukey"))
      stopifnot(nrow(raw) == 6L, identical(raw$contrast, tukey$contrast))
      parts <- strsplit(as.character(raw$contrast), " - ", fixed = TRUE)
      comparisons[[key]] <- cbind(repeat_context(6L), raw,
        group1 = vapply(parts, "[", character(1), 1L),
        group2 = vapply(parts, "[", character(1), 2L),
        p_raw = raw$p.value, p_Tukey_historical = tukey$p.value,
        p_Tukey_then_BH_historical = stats::p.adjust(tukey$p.value, method = "BH"),
        comparison_direction = as.character(raw$contrast))
      res <- stats::residuals(fit)
      convergence <- unlist(fit@optinfo$conv$lme4$messages)
      diagnostics[[key]] <- cbind(context, singular = lme4::isSingular(fit),
        optimizer_code = paste(fit@optinfo$conv$opt, collapse = ";"),
        convergence_messages = paste(convergence, collapse = ";"),
        warnings = paste(unique(warnings), collapse = ";"),
        residual_mean = mean(res), residual_sd = stats::sd(res),
        residual_min = min(res), residual_max = max(res),
        overall_group_P = anova(fit)[["Pr(>F)"]][1])
      residuals[[key]] <- cbind(repeat_context(length(res)), sample_id = x$sample_id,
        fitted = stats::fitted(fit), residual = res)
      vc <- as.data.frame(lme4::VarCorr(fit))
      variance[[key]] <- cbind(repeat_context(nrow(vc)), vc)
    }
  }
  pairs <- do.call(rbind, comparisons)
  pairs$p_adj_BH <- NA_real_
  for (rank in c("phylum", "species")) {
    ix <- which(pairs$taxonomic_level == rank)
    stopifnot(length(ix) == if (rank == "phylum") 24L else 48L,
              all(is.finite(pairs$p_raw[ix])))
    pairs$p_adj_BH[ix] <- stats::p.adjust(pairs$p_raw[ix], method = "BH")
  }
  pairs$formal_significant <- pairs$p_adj_BH < config$project$alpha
  pairs$historical_significant <- pairs$p_Tukey_then_BH_historical < config$project$alpha
  pairs$bracket_change <- ifelse(pairs$formal_significant,
    ifelse(pairs$historical_significant, "retained", "added"),
    ifelse(pairs$historical_significant, "removed", "absent_in_both"))
  results <- list(selected_top4_phyla = selected$phylum, selected_top8_species = selected$species,
    abundant_taxa_sample_level_data = do.call(rbind, samples),
    abundant_taxa_LMM_coefficients = do.call(rbind, coefficients),
    abundant_taxa_emmeans_pairwise_all = pairs,
    abundant_taxa_emmeans_significant = pairs[pairs$formal_significant, , drop = FALSE],
    abundant_taxa_bracket_changes = pairs,
    abundant_taxa_model_diagnostics = do.call(rbind, diagnostics),
    abundant_taxa_variance_components = do.call(rbind, variance),
    abundant_taxa_model_residuals = do.call(rbind, residuals),
    abundant_taxa_analysis_summary = data.frame(taxonomic_level = c("phylum", "species"),
      n_features = c(4L, 8L), n_comparisons = c(24L, 48L),
      significant = c(sum(pairs$formal_significant[pairs$taxonomic_level == "phylum"]),
                      sum(pairs$formal_significant[pairs$taxonomic_level == "species"])),
      formal_q = "p_adj_BH; unadjusted P values pooled within taxonomic level",
      scope = "Secondary/exploratory abundant-taxa analysis"))
  for (name in names(results)) {
    x <- results[[name]]
    n <- nrow(x)
    fields <- list(model_formula = cfg$formula, transformation = "None; original percentage relative abundance",
      n_samples = 1186L, reference_group = analysis_group_levels[1],
      comparison_direction = "group1 minus group2; fixed coefficients use the stated reference",
      p_adjust_method = "BH on unadjusted emmeans P values within taxonomic level",
      formal_significance_column = "p_adj_BH", analysis_scope = "Secondary/exploratory abundant-taxa analysis")
    for (field in setdiff(names(fields), names(x))) x[[field]] <- rep(fields[[field]], n)
    if (!"family_size" %in% names(x)) x$family_size <- ifelse(x$taxonomic_level == "phylum", 24L, 48L)
    if (!"n_features" %in% names(x)) x$n_features <- ifelse(x$taxonomic_level == "phylum", 4L, 8L)
    if (!"formal_significant" %in% names(x)) x$formal_significant <- rep(NA, n)
    if (!"test_family" %in% names(x)) x$test_family <- paste(x$taxonomic_level, "all selected features x six comparisons")
    results[[name]] <- x
  }
  results
}

abundant_taxa_figures <- function(run, config) {
  require_project_packages(c("ggplot2", "patchwork"))
  style <- read_presentation_config(run$root)
  config$output <- utils::modifyList(config$output, style$figure)
  d <- read_current_result(run, "tables/04b_abundant_taxa/abundant_taxa_sample_level_data.csv")
  all_pairs <- read_current_result(run, "tables/04b_abundant_taxa/abundant_taxa_emmeans_pairwise_all.csv")
  diagnostics <- read_current_result(run, "tables/04b_abundant_taxa/abundant_taxa_model_diagnostics.csv")
  catalogue <- list()
  for (rank in c("phylum", "species")) {
    selected <- read_current_result(run, paste0("tables/04b_abundant_taxa/",
      if (rank == "phylum") "selected_top4_phyla.csv" else "selected_top8_species.csv"))
    panels <- lapply(selected$feature, function(id) {
      x <- d[d$feature == id & d$taxonomic_level == rank, , drop = FALSE]
      x$display_group <- factor(group_display[x$analysis_group], levels = unlist(style$groups$order))
      hit <- all_pairs[all_pairs$feature == id & all_pairs$taxonomic_level == rank &
                         all_pairs$formal_significant, , drop = FALSE]
      maximum <- max(x$raw_abundance)
      hit$y <- if (nrow(hit)) seq(maximum * 1.1, maximum * 1.4, length.out = nrow(hit)) else numeric()
      hit$x1 <- match(group_display[hit$group1], levels(x$display_group))
      hit$x2 <- match(group_display[hit$group2], levels(x$display_group))
      p <- ggplot2::ggplot(x, ggplot2::aes(display_group, raw_abundance)) +
        ggplot2::geom_boxplot(ggplot2::aes(fill = display_group, colour = display_group),
          outlier.shape = NA, width = .4, alpha = .6, linewidth = .8) +
        ggplot2::geom_point(ggplot2::aes(colour = display_group), alpha = .7, size = .5,
          position = ggplot2::position_jitter(width = .1, height = 0, seed = config$output$label_seed)) +
        ggplot2::scale_fill_manual(values = unlist(style$groups$fill)) +
        ggplot2::scale_colour_manual(values = unlist(style$groups$dark)) +
        ggplot2::labs(title = if (rank == "species") parse(text = species_plot_label(id))[[1]] else id,
          subtitle = paste0("LMM group P = ", formatC(diagnostics$overall_group_P[
            diagnostics$feature == id & diagnostics$taxonomic_level == rank], digits = 3, format = "g")),
          x = NULL, y = "Relative abundance (%)") +
        publication_theme(config) +
        ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
        ggplot2::coord_cartesian(ylim = c(0, maximum * 1.5))
      if (nrow(hit)) p <- p +
        ggplot2::geom_segment(data = hit, ggplot2::aes(x = x1, xend = x1, y = y, yend = y - maximum * .02),
          inherit.aes = FALSE, linewidth = .4) +
        ggplot2::geom_segment(data = hit, ggplot2::aes(x = x2, xend = x2, y = y, yend = y - maximum * .02),
          inherit.aes = FALSE, linewidth = .4) +
        ggplot2::geom_segment(data = hit, ggplot2::aes(x = x1, xend = x2, y = y, yend = y),
          inherit.aes = FALSE, linewidth = .4) +
        ggplot2::geom_text(data = hit, ggplot2::aes(x = (x1 + x2) / 2, y = y,
          label = paste0("BH P = ", formatC(p_adj_BH, digits = 2, format = "g"))),
          inherit.aes = FALSE, vjust = -.4, family = style$font_family, size = 3.2)
      p
    })
    id <- paste0("SuppFig2_BMI25_", rank, "_exploratory_LMM")
    safe_write_csv(d[d$taxonomic_level == rank, ], path_in_run(run, "plot_data", paste0(id, ".csv")))
    safe_write_csv(all_pairs[all_pairs$taxonomic_level == rank, ],
                   path_in_run(run, "plot_data", paste0(id, "_comparisons.csv")))
    caption <- paste(
      "Secondary/exploratory abundant-taxa analysis. BMI25 groups use <25 and >=25 kg/m2.",
      "LMM: abundance ~ analysis_group + age + sex + (1 | study_name), REML.",
      "Unadjusted emmeans P values are BH-adjusted across all 24 phylum or 48 species contrasts.",
      "Only p_adj_BH <0.05 is annotated. Historical Tukey-then-BH is an audit column.",
      "Formal species inference is based on the prevalence-filtered MaAsLin2 models.")
    catalogue[[rank]] <- save_publication_figure(patchwork::wrap_plots(panels, ncol = 4),
      id, run, config, 18, if (rank == "phylum") 5 else 10, caption)
  }
  do.call(rbind, catalogue)
}
