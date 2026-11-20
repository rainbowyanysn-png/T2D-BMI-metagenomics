saved_p_label <- function(p) {
  ifelse(is.na(p), "= NA", ifelse(p < .001, "< 0.001", paste0("= ", formatC(p, format = "f", digits = 3))))
}

alpha_reference_panel <- function(alpha, summary, metric, grouping, style, config, tag = NULL) {
  data <- alpha
  data$value <- data[[if (metric == "Shannon") "shannon" else "richness"]]
  total <- grouping == "total"
  data$display_group <- if (total) factor(ifelse(data$disease == "control", "Ctrl", "T2D"), levels = c("Ctrl", "T2D")) else factor(group_display[as.character(data$analysis_group)], levels = unlist(style$groups$order))
  colours <- unlist(if (total) style$groups$total else style$groups$fill)
  dark <- unlist(if (total) style$groups$total_dark else style$groups$dark)
  row <- summary[summary$Metric == metric, , drop = FALSE]
  assert_equal(nrow(row), 1L, "Saved alpha summary metric")
  subtitle <- if (total) paste0("Global LMM P ", saved_p_label(row$Overall_T2D_vs_Control_p)) else paste0("Four-group LMM P ", saved_p_label(row$Four_Group_p), " | ICC = ", row$ICC)
  ggplot2::ggplot(data, ggplot2::aes(display_group, value)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = display_group, colour = display_group), outlier.shape = NA, width = .4, alpha = .6, linewidth = .65) +
    ggplot2::geom_point(ggplot2::aes(colour = display_group), position = ggplot2::position_jitter(width = .1, height = 0, seed = config$output$label_seed), alpha = .55, size = .85) +
    ggplot2::scale_fill_manual(values = colours) + ggplot2::scale_colour_manual(values = dark) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.04, .12))) +
    ggplot2::labs(title = if (total) paste("Global", metric, "diversity") else paste(metric, "across BMI25 groups"), subtitle = subtitle, x = NULL, y = if (metric == "Shannon") "Shannon index" else "Observed species richness", tag = tag) +
    publication_theme(config) + ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(face = "bold"))
}

pcoa_reference_panel <- function(coordinates, style, config, cohorts, total = FALSE, tag = NULL) {
  data <- coordinates
  data$phenotype <- if (total) factor(data$disease_display, levels = c("Ctrl", "T2D")) else data$display_group
  data$study_name <- factor(data$study_name, levels = cohorts)
  palette <- unlist(if (total) style$groups$total else style$groups$dark)
  shapes <- setNames(as.integer(unlist(style$beta$cohort_shapes)), cohorts)
  ggplot2::ggplot(data, ggplot2::aes(PC1, PC2)) +
    ggplot2::geom_point(ggplot2::aes(colour = phenotype, shape = study_name), alpha = style$beta$point_alpha, size = style$beta$point_size) +
    ggplot2::stat_ellipse(ggplot2::aes(group = phenotype, colour = phenotype), type = style$beta$ellipse_type, level = style$beta$ellipse_level, linetype = 2, linewidth = style$beta$ellipse_linewidth, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = palette, drop = FALSE) + ggplot2::scale_shape_manual(values = shapes, drop = FALSE) +
    ggplot2::labs(x = "PCoA 1", y = "PCoA 2", tag = tag) + publication_theme(config) +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 16), axis.text = ggplot2::element_text(size = 13), legend.text = ggplot2::element_text(size = 11))
}

permanova_reference_panel <- function(data, config, tag = NULL) {
  assert_columns(data, c("study_name", "R2", "p_value", "n_samples"), "Saved cohort PERMANOVA")
  data$cohort <- factor(data$study_name, levels = data$study_name[order(data$R2)])
  data$label <- paste0("P ", saved_p_label(data$p_value), " | n = ", data$n_samples)
  ggplot2::ggplot(data, ggplot2::aes(R2, cohort)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = R2, yend = cohort), colour = "grey80", linewidth = .55) +
    ggplot2::geom_point(colour = "#2166AC", size = 2.5) +
    ggplot2::geom_text(ggplot2::aes(label = label), hjust = -.12, size = 3, family = config$output$font_family) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(.02, .6))) +
    ggplot2::labs(x = "PERMANOVA model R-squared", y = NULL, tag = tag) + publication_theme(config)
}

species_heatmap_data <- function(species, comparison_order) {
  hits <- unique(species$feature[species$formal_q < .05])
  grid <- expand.grid(feature = sort(hits), contrast_display = comparison_order, stringsAsFactors = FALSE)
  key <- paste(species$feature, species$contrast_display, sep = "\r")
  i <- match(paste(grid$feature, grid$contrast_display, sep = "\r"), key)
  grid$coef <- species$coef[i]; grid$formal_q <- species$formal_q[i]
  grid$tested <- !is.na(i)
  grid$label <- ifelse(is.na(grid$formal_q), "", ifelse(grid$formal_q < .001, "***", ifelse(grid$formal_q < .01, "**", ifelse(grid$formal_q < .05, "*", ""))))
  grid$contrast_display <- factor(grid$contrast_display, levels = comparison_order)
  grid$species_label <- factor(species_plot_label(grid$feature), levels = rev(species_plot_label(sort(hits))))
  grid
}
raw_aitchison_figures <- function(run, config) {
  style <- read_presentation_config(run$root)
  config$output <- utils::modifyList(config$output, style$figure)
  read <- function(name) read_current_result(run,
    paste0("tables/04c_raw_Aitchison/BMI25_raw_Aitchison_", name, ".csv"))
  coordinates <- read("PCoA_coordinates")
  md <- read("sample_manifest")
  j <- match(coordinates$sample_id, md$sample_id)
  stopifnot(!anyNA(j), !anyDuplicated(coordinates$sample_id), nrow(coordinates) == 1182L)
  coordinates$display_group <- factor(group_display[md$analysis_group[j]], levels = unlist(style$groups$order))
  coordinates$disease_display <- ifelse(md$disease[j] == "control", "Ctrl", "T2D")
  coordinates$study_name <- md$study_name[j]
  catalogue <- list()
  save <- function(plot, id, data, width, height, caption) {
    safe_write_csv(data, path_in_run(run, "plot_data", paste0(id, ".csv")))
    catalogue[[id]] <<- save_publication_figure(plot, id, run, config, width, height, caption)
  }
  for (total in c(TRUE, FALSE)) {
    id <- paste0("Fig2_BMI25_raw_Aitchison_PCoA_", if (total) "total" else "four_groups")
    plot <- pcoa_reference_panel(coordinates, style, config, unlist(config$cohorts$names), total = total) +
      ggplot2::labs(x = sprintf("PCoA 1 (%.1f%%)", 100 * coordinates$PC1_fraction[1]),
                    y = sprintf("PCoA 2 (%.1f%%)", 100 * coordinates$PC2_fraction[1]))
    save(plot, id, coordinates, 9, 6.5,
      "BMI25, uncorrected CZM/CLR Aitchison distances, 1182 samples and 138 features. Ordination and dashed 95% t data ellipses are descriptive. Point shapes identify cohorts. Dispersion results accompany PERMANOVA.")
    geometry <- ggplot2::ggplot_build(plot)$data[[2]]
    safe_write_csv(geometry, path_in_run(run, "plot_data", paste0(id, "_ellipse_geometry.csv")))
  }
  cohort <- read("PERMANOVA_by_cohort")
  save(permanova_reference_panel(cohort, config), "Fig2_BMI25_raw_Aitchison_cohort_PERMANOVA",
    cohort, 9, 4.8,
    "Raw Aitchison four-group marginal PERMANOVA within each cohort, 999 permutations. Nonvarying nuisance covariates are recorded in plot-data. Group dispersion may contribute to observed differences.")
  for (scope in c("disease", "four_groups")) {
    d <- read(paste0("PERMDISP_", scope))
    plotted <- d[d$term == "Groups", , drop = FALSE]
    plotted$label <- sprintf("Permutation P = %.3f", plotted[["Pr(>F)"]])
    plot <- ggplot2::ggplot(plotted, ggplot2::aes(F, group_variable)) +
      ggplot2::geom_point(colour = "#2166AC", size = 3) +
      ggplot2::geom_text(ggplot2::aes(label = label), hjust = -.1, family = style$font_family, size = 3.6) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(.1, .8))) +
      ggplot2::labs(x = "PERMDISP F statistic", y = NULL) + publication_theme(config)
    save(plot, if (scope == "four_groups") "SuppFig_BMI25_raw_Aitchison_PERMDISP" else
      "SuppFig_BMI25_raw_Aitchison_PERMDISP_disease", d, 8, 3.1,
      "Raw Aitchison distances, spatial-median dispersion, 999 unrestricted permutations. Significant dispersion indicates that PERMANOVA differences may involve both group location and within-group spread.")
  }
  # Resolve the source directory without interpreting filenames as directories.
  folder <- dirname(current_result_path(run, "tables/04c_raw_Aitchison/BMI25_raw_Aitchison_PCoA_coordinates.csv"))
  files <- list.files(folder, pattern = "[.]csv$", full.names = TRUE)
  safe_write_csv(data.frame(input_file = basename(files), sha256 = vapply(files, sha256_file, character(1))),
                 path_in_run(run, "plot_data/raw_Aitchison_input_sha256.csv"))
  do.call(rbind, catalogue)
}

