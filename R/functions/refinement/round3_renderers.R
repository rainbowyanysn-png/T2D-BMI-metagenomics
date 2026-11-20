r3_columns <- function(d, columns, id) {
  missing <- setdiff(columns, names(d))
  r3_need(!length(missing), paste(id, "is missing saved columns:", paste(missing, collapse = ", ")))
}

r3_clean_pathway_taxon <- function(x) {
  x <- r3_html_symbols(as.character(x))
  base <- sub("\\|.*$", "", x)
  taxon <- ifelse(grepl("s__", x), sub(".*s__", "", x), ifelse(grepl("g__", x), sub(".*g__", "", x), ""))
  taxon <- gsub("_", " ", taxon, fixed = TRUE)
  ifelse(nzchar(taxon), paste0(base, " (", taxon, ")"), base)
}

r3_no_legend_frame <- function(p) p + ggplot2::theme(
  legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
  legend.box.background = ggplot2::element_rect(fill = "transparent", colour = NA),
  legend.key = ggplot2::element_rect(fill = "transparent", colour = NA))

r3_layout_policy <- function(p, settings = list()) {
  elements <- c("legend.background", "legend.box.background", "legend.key")
  removed <- all(vapply(elements, function(n) {
    e <- p$theme[[n]]
    !is.null(e) && identical(e$fill, "transparent") && length(e$colour) == 1L && is.na(e$colour)
  }, logical(1)))
  values <- c(list(legend_outer_frames_removed = removed), settings)
  data.frame(setting = names(values), value = vapply(values, function(x) paste(as.character(x), collapse = "|"), character(1)), stringsAsFactors = FALSE)
}

r3_taxon_pathway_candidates <- function(d) {
  audit <- r3_candidate_filter_audit(d)
  d$coef <- as.numeric(d$coef); d$qval <- as.numeric(d$qval)
  # The eligibility audit contains every predicate; legacy panel flags/ranks
  # in the source table are retained for audit but never drive selection.
  selected <- dplyr::filter(d, audit$included) |>
    dplyr::arrange(qval, dplyr::desc(abs(coef)), feature)
  r3_assert_candidate_counts(selected)
  selected$display_selection_rank <- seq_len(nrow(selected))
  selected
}

r3_plot_pathway_volcano <- function(d, panel, style) {
  id <- panel$figure_id
  r3_columns(d, c("feature", "coef", "qval", "label_selected", "label_rank", "label_text"), id)
  d$coef <- as.numeric(d$coef); d$qval <- as.numeric(d$qval)
  d$label_selected <- r3_as_logical(d$label_selected)
  r3_need(nrow(d) == 5945L && sum(d$qval < .05) == 1641L && sum(d$label_selected) == 12L,
          "Fig5A frozen 5945/1641/12 counts differ")
  r3_need(!anyDuplicated(d$feature) && all(is.finite(d$coef)) && all(is.finite(d$qval)) && all(d$qval > 0 & d$qval <= 1),
          "Fig5A feature/coefficient/q-value validity failed")
  d$display_class <- ifelse(d$qval < .05 & d$coef > 0, "Positive",
                            ifelse(d$qval < .05 & d$coef < 0, "Negative", "NS"))
  d$minus_log10_q <- -log10(d$qval)
  d$label_display <- ifelse(d$label_selected, r3_html_symbols(d$label_text), NA_character_)
  selected <- d[d$label_selected, , drop = FALSE]
  selected <- selected[order(selected$label_rank), , drop = FALSE]
  neg <- selected[selected$coef < 0, , drop = FALSE]
  pos <- selected[selected$coef > 0, , drop = FALSE]
  nudge <- max(abs(d$coef), na.rm = TRUE) * 0.10
  cols <- c(Positive = style$direction$positive, Negative = style$direction$negative, NS = style$direction$nonsignificant)
  p <- ggplot2::ggplot(d, ggplot2::aes(coef, minus_log10_q)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.35) +
    ggplot2::geom_hline(yintercept = -log10(.05), colour = "grey60", linewidth = 0.35, linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(colour = display_class), size = 0.85, alpha = 0.42) +
    ggplot2::scale_colour_manual(values = cols, breaks = c("Positive", "Negative", "NS"),
                                 labels = c("Positive", "Negative", "Not significant"), name = "Association") +
    ggrepel::geom_text_repel(
      data = neg, ggplot2::aes(label = label_display), nudge_x = -nudge,
      family = style$font$family, size = 2.8, fontface = "plain", lineheight = 0.88,
      direction = "both", force = 2.5, box.padding = 0.45, point.padding = 0.20,
      min.segment.length = 0, segment.color = "grey45", segment.size = 0.25,
      max.overlaps = Inf, seed = 20261116, na.rm = TRUE
    ) +
    ggrepel::geom_text_repel(
      data = pos, ggplot2::aes(label = label_display), nudge_x = nudge,
      family = style$font$family, size = 2.8, fontface = "plain", lineheight = 0.88,
      direction = "both", force = 2.5, box.padding = 0.45, point.padding = 0.20,
      min.segment.length = 0, segment.color = "grey45", segment.size = 0.25,
      max.overlaps = Inf, seed = 20261116, na.rm = TRUE
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(.28, .28))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.02, .18))) +
    ggplot2::labs(title = panel$title, x = "MaAsLin2 coefficient", y = expression(-log[10]("MaAsLin2 q value"))) +
    r3_theme(style, 9) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.35),
      legend.position = "bottom", plot.margin = ggplot2::margin(7, 10, 5, 10)
    )
  p <- r3_no_legend_frame(p)
  audit <- selected[, c("feature", "coef", "qval", "label_rank", "label_text", "label_display")]
  list(plot = p, display_data = d, audits = list(label_audit = audit, layout_policy = r3_layout_policy(p)),
       caption = "All 5,945 frozen disease-term coefficients are shown. Formal significance uses the saved MaAsLin2 full-model qval; 1,641 features have qval<0.05. The 12 saved labels are unchanged.")
}

r3_plot_pathway_taxonomy_bar <- function(d, panel, style) {
  full <- d; d <- r3_taxon_pathway_candidates(full)
  z <- dplyr::count(d, genus, taxonomic_contributor, name = "n")
  r3_need(nrow(z) == 6L && sum(z$n) == 64L, "Fig5B must show all64 assignments across six species")
  genus <- sort(unique(z$genus), method = "radix")
  fill <- stats::setNames(unlist(style$taxonomy_palette), genus)
  border <- colorspace::darken(fill, amount = 0.30)
  order <- z$taxonomic_contributor[order(z$n, z$taxonomic_contributor)]
  z$taxonomic_contributor <- factor(z$taxonomic_contributor, levels = order)
  p <- ggplot2::ggplot(z, ggplot2::aes(taxonomic_contributor, n, fill = genus, colour = genus)) +
    ggplot2::geom_col(width = 0.40, linewidth = 0.50) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = fill, name = "Genus") +
    ggplot2::scale_colour_manual(values = border, name = "Genus") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .08)), breaks = scales::breaks_width(2)) +
    ggplot2::labs(title = panel$title, x = "Species", y = "Number of selected pathway features") +
    r3_theme(style, 9) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "italic", size = 8.5),
      legend.text = ggplot2::element_text(face = "italic"), legend.position = "right"
    )
  p <- r3_no_legend_frame(p)
  list(plot = p, display_data = d, source_features = d$feature,
       audits = list(counts = z, candidate_selection_audit = r3_candidate_filter_audit(full), layout_policy = r3_layout_policy(p)),
       caption = "All 64 eligible taxon-stratified pathway features were included. Bar heights count feature assignments for six species in four genera; selection uses saved visualization_candidate, species-level annotations, full-model qval <0.01 and absolute coefficient >1, with named contributors and UNINTEGRATED/UNMAPPED excluded. No top-N truncation or statistical refit was applied.")
}

r3_plot_pathway_bubble <- function(d, panel, style) {
  full <- d; d <- r3_taxon_pathway_candidates(full)
  r3_need(!anyDuplicated(d[c("base_pathway_ID", "taxonomic_contributor")]), "Multiple feature assignments overlap at one bubble location; do not merge silently")
  path_order <- unique(d$base_pathway_ID)
  d$pathway_display <- d$base_pathway_ID
  species_order <- d |>
    dplyr::count(taxonomic_contributor, sort = TRUE) |>
    dplyr::arrange(n, taxonomic_contributor) |>
    dplyr::pull(taxonomic_contributor)
  d$pathway_display <- factor(d$pathway_display, levels = unique(path_order))
  d$taxonomic_contributor <- factor(d$taxonomic_contributor, levels = unique(species_order))
  d$direction <- ifelse(d$coef > 0, "Positive", "Negative")
  cols <- c(Positive = style$direction$positive, Negative = style$direction$negative)
  p <- ggplot2::ggplot(d, ggplot2::aes(pathway_display, taxonomic_contributor)) +
    ggplot2::geom_point(ggplot2::aes(size = abs(coef), colour = direction), alpha = .82) +
    ggplot2::scale_colour_manual(values = cols, name = "Direction") +
    ggplot2::scale_size_continuous(range = c(2.5, 7), name = "Absolute MaAsLin2 coefficient") +
    ggplot2::labs(title = panel$title, x = "MetaCyc pathway feature", y = "Species") +
    r3_theme(style, 8.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 6.6),
      axis.text.y = ggplot2::element_text(face = "italic", size = 8.5),
      panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = .25),
      legend.position = "right", plot.margin = ggplot2::margin(6, 5, 12, 6)
    )
  p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1, size = 5.8),
    axis.text.y = ggplot2::element_text(face = "italic", size = 7.0, colour = "black"))
  p <- r3_no_legend_frame(p)
  mapping <- d[, c("feature", "base_pathway_ID", "pathway_description", "taxonomic_contributor", "genus", "coef", "qval")]
  list(plot = p, display_data = d, source_features = d$feature,
       audits = list(axis_order = data.frame(pathway_display = path_order, order = seq_along(path_order)),
         pathway_mapping = mapping, candidate_selection_audit = r3_candidate_filter_audit(full),
         layout_policy = r3_layout_policy(p, list(species_label_size = p$theme$axis.text.y$size,
           pathway_label_size = p$theme$axis.text.x$size, pathway_label_angle = p$theme$axis.text.x$angle,
           pathway_labels = "base_pathway_ID", n_bubbles = nrow(d),
           size_legend_title = p$scales$get_scales("size")$name))),
       caption = "All 64 eligible taxon-stratified pathway features were included. The selection contains six species,42 MetaCyc pathway IDs,37 positive and27 negative coefficients. Axes display pathway IDs; full names remain in the CSV and mapping audit. Point size is the absolute MaAsLin2 coefficient and colour is its direction. No top-N truncation, statistical refit or P/q recalculation was applied.")
}

r3_circle <- function(cx, cy, radius = 1, n = 300, set) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + radius * cos(theta), y = cy + radius * sin(theta), set = set)
}

r3_plot_pathway_overlap <- function(d, panel, style) {
  r3_columns(d, c("set", "n"), panel$figure_id)
  d$n <- as.integer(d$n)
  expected <- c("Overall only" = 9L, "Shared" = 4L, "NT2D vs NC only" = 129L)
  got <- stats::setNames(d$n, d$set)
  r3_need(identical(unname(got[names(expected)]), unname(expected)), "Fig6A frozen 9/4/129 overlap differs")
  circles <- rbind(r3_circle(-.62, 0, set = "T2D vs Control"), r3_circle(.62, 0, set = "NT2D vs NC"))
  fills <- c("T2D vs Control" = style$overlap$overall, "NT2D vs NC" = style$overlap$nt2d)
  labels <- data.frame(x = c(-.88, 0, .88), y = 0, n = c(9, 4, 129))
  set_labels <- data.frame(x = c(-.90, .90), y = 1.18, label = c("T2D vs Control", "NT2D vs NC"))
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = circles, ggplot2::aes(x, y, group = set, fill = set), alpha = .60, colour = NA) +
    ggplot2::geom_text(data = labels, ggplot2::aes(x, y, label = n), family = style$font$family, fontface = "bold", size = 5) +
    ggplot2::geom_text(data = set_labels, ggplot2::aes(x, y, label = label), family = style$font$family, fontface = "bold", size = 3.7) +
    ggplot2::scale_fill_manual(values = fills, guide = "none") +
    ggplot2::coord_equal(xlim = c(-1.75, 1.75), ylim = c(-1.18, 1.45), clip = "off") +
    ggplot2::labs(title = panel$title) +
    ggplot2::theme_void(base_family = style$font$family, base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = .5, size = 11))
  list(plot = p, display_data = d, audits = list(overlap_counts = d,
       layout_policy = r3_layout_policy(p, list(subtitle_absent = is.null(p$labels$subtitle)))),
       caption = "Complete frozen significant sets are summarized: 9 overall-only, 4 shared, and 129 NT2D-vs-NC-only pathway features. Circle areas are not proportional to counts.")
}

r3_q_stars <- function(q) ifelse(is.na(q), "", ifelse(q < .001, "***", ifelse(q < .01, "**", ifelse(q < .05, "*", ""))))

# Rectangular left-side tree coordinates from the existing hclust object.
# No new clustering, labels, feature selection or external tree package.
r3_hclust_segments <- function(hc) {
  r3_need(inherits(hc, "hclust") && is.matrix(hc$merge), "Invalid hclust object for Round3 dendrogram")
  n <- length(hc$order)
  r3_need(n >= 2L && nrow(hc$merge) == n - 1L && ncol(hc$merge) == 2L &&
    length(hc$height) == n - 1L && identical(sort(hc$order), seq_len(n)) &&
    all(is.finite(hc$height) & hc$height >= 0), "Invalid hclust object for Round3 dendrogram")
  leaf_position <- integer(n); leaf_position[hc$order] <- seq_len(n)
  cluster_position <- numeric(n - 1L)
  segments <- vector("list", 3L * (n - 1L)); segment_index <- 0L
  child_position <- function(child) if (child < 0L) leaf_position[-child] else cluster_position[child]
  child_height <- function(child) if (child < 0L) 0 else hc$height[child]
  for (i in seq_len(n - 1L)) {
    children <- hc$merge[i, ]
    r3_need(all(children != 0L & ifelse(children < 0L, -children <= n, children < i)), "Invalid dendrogram child reference")
    left <- children[1]; right <- children[2]
    lp <- child_position(left); rp <- child_position(right)
    lh <- child_height(left); rh <- child_height(right); ph <- hc$height[i]
    cluster_position[i] <- (lp + rp) / 2
    segment_index <- segment_index + 1L
    segments[[segment_index]] <- data.frame(x = -lh, y = lp, xend = -ph, yend = lp)
    segment_index <- segment_index + 1L
    segments[[segment_index]] <- data.frame(x = -rh, y = rp, xend = -ph, yend = rp)
    segment_index <- segment_index + 1L
    segments[[segment_index]] <- data.frame(x = -ph, y = lp, xend = -ph, yend = rp)
  }
  do.call(rbind, segments)
}

r3_plot_pathway_heatmap <- function(d, panel, style) {
  id <- panel$figure_id
  r3_columns(d, c("pathway", "comparison", "coef", "formal_q", "tested", "label"), id)
  d$coef <- as.numeric(d$coef); d$formal_q <- as.numeric(d$formal_q); d$tested <- r3_as_logical(d$tested)
  comparisons <- c("T2D vs Control", "NT2D vs NC", "OT2D vs OC", "OC vs NC", "OT2D vs NT2D")
  r3_need(nrow(d) == 145L && length(unique(d$pathway)) == 29L && setequal(unique(d$comparison), comparisons), "Fig6B frozen 29x5 matrix differs")
  r3_need(!anyDuplicated(d[, c("pathway", "comparison")]), "Fig6B has duplicate pathway-comparison cells")
  wide <- tidyr::pivot_wider(d[, c("pathway", "comparison", "coef")], names_from = comparison, values_from = coef)
  wide <- wide[, c("pathway", comparisons)]
  mat <- as.matrix(wide[, comparisons]); rownames(mat) <- wide$pathway
  r3_need(all(is.na(mat) | is.finite(mat)), "Invalid infinite heatmap coefficient; do not replace it with zero")
  cluster_mat <- mat; cluster_mat[is.na(cluster_mat)] <- 0
  hc <- stats::hclust(stats::dist(cluster_mat), method = "complete")
  row_order <- rownames(mat)[hc$order]
  d$row_index <- match(d$pathway, row_order)
  d$comparison <- factor(d$comparison, levels = comparisons)
  d$coef_display <- pmax(pmin(d$coef, style$heatmap$limit), -style$heatmap$limit)
  d$star <- ifelse(d$tested, r3_q_stars(d$formal_q), "")
  labels <- d |>
    dplyr::group_by(pathway) |>
    dplyr::summarise(label = dplyr::first(label), .groups = "drop")
  labels$label <- make.unique(r3_html_symbols(labels$label), sep = " [duplicate ")
  label_vec <- stats::setNames(labels$label, labels$pathway)
  seg <- r3_hclust_segments(hc)
  p_dend <- ggplot2::ggplot(seg) +
    ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend), linewidth = .30, colour = "grey25", lineend = "square") +
    ggplot2::scale_y_continuous(limits = c(.5, length(row_order) + .5), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(.04, .02))) +
    ggplot2::theme_void()
  p_heat <- ggplot2::ggplot(d, ggplot2::aes(comparison, row_index, fill = coef_display)) +
    ggplot2::geom_tile(colour = "white", linewidth = .30) +
    ggplot2::geom_text(ggplot2::aes(label = star), family = style$font$family, size = 2.4, colour = "black") +
    ggplot2::scale_fill_gradient2(
      low = style$heatmap$low, mid = style$heatmap$mid, high = style$heatmap$high,
      midpoint = 0, limits = c(-2.5, 2.5), breaks = c(-2.5, 0, 2.5),
      labels = c("-2.5", "0", "2.5"), na.value = style$heatmap$na, name = NULL,
      guide = ggplot2::guide_colorbar(title = NULL, ticks = TRUE, frame.colour = NA,
        barheight = grid::unit(18, "mm"), barwidth = grid::unit(4, "mm"))
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(row_order), labels = unname(label_vec[row_order]),
      limits = c(.5, length(row_order) + .5), expand = c(0, 0), position = "right"
    ) +
    ggplot2::labs(title = panel$title, x = NULL, y = NULL) +
    r3_theme(style, 7.5) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      axis.line = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7.2),
      axis.text.y = ggplot2::element_text(size = 6.5), axis.ticks = ggplot2::element_blank(),
      legend.position = "right", legend.justification = c(0, 1), legend.box.just = "top",
      plot.margin = ggplot2::margin(4, 4, 4, 0)
    )
  combined <- patchwork::wrap_plots(p_dend, p_heat, nrow = 1, widths = c(1.0, 6.2), guides = "collect")
  audit <- data.frame(order = seq_along(row_order), pathway = row_order, label = unname(label_vec[row_order]), stringsAsFactors = FALSE)
  fill_scale <- p_heat$scales$get_scales("fill")
  list(plot = combined, display_data = d, audits = list(row_cluster_order = audit,
       colour_legend = data.frame(break_value = fill_scale$breaks, label = fill_scale$labels),
       layout_policy = r3_layout_policy(p_heat, list(legend_position = "right", legend_justification = "0|1", colour_title_absent = is.null(fill_scale$name),
         colour_scale_continuous = inherits(fill_scale, "ScaleContinuous"),
         heatmap_background_border_absent = is.na(p_heat$theme$panel.background$colour),
         heatmap_panel_border_blank = inherits(p_heat$theme$panel.border, "element_blank"),
         heatmap_axis_line_blank = inherits(p_heat$theme$axis.line, "element_blank"))),
       column_order = data.frame(order = seq_along(comparisons), comparison = comparisons),
       dendrogram_segments = seg, row_cluster_merges = data.frame(merge_index = seq_len(nrow(hc$merge)),
         left = hc$merge[, 1], right = hc$merge[, 2], height = hc$height)),
       caption = "The frozen 29-feature by five-comparison matrix is shown. Rows are clustered using saved pooled coefficients with missing cells set to zero only for distance calculation; comparison order is fixed. Colour clipping at ±2.5 is display-only. Stars use the saved within-contrast REML padj.")
}

r3_plot_pathway_forest <- function(d, panel, style) {
  id <- panel$figure_id
  r3_columns(d, c("study", "coef", "ci_lower", "ci_upper", "pathway", "I2", "padj", "pooled"), id)
  d$coef <- as.numeric(d$coef); d$ci_lower <- as.numeric(d$ci_lower); d$ci_upper <- as.numeric(d$ci_upper)
  d$I2 <- as.numeric(d$I2); d$padj <- as.numeric(d$padj); d$pooled <- r3_as_logical(d$pooled)
  expected_rows <- 6L
  r3_need(nrow(d) == expected_rows && sum(d$pooled) == 1L && length(unique(d$pathway)) == 1L, paste(id, "frozen forest rows differ"))
  other <- d$study[!d$pooled]
  d$study <- factor(d$study, levels = c(d$study[d$pooled], rev(other)))
  d$type <- ifelse(d$pooled, "REML pooled", "Cohort")
  pooled <- d[d$pooled, , drop = FALSE]
  subtitle <- sprintf("%s; I² = %.1f%%; adjusted P = %.3g", panel$title, pooled$I2[1], pooled$padj[1])
  pathway_title <- stringr::str_wrap(r3_clean_pathway_taxon(unique(d$pathway)), width = 72)
  p <- ggplot2::ggplot(d, ggplot2::aes(coef, study, xmin = ci_lower, xmax = ci_upper)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = .35) +
    ggplot2::geom_errorbarh(ggplot2::aes(colour = type), height = .16, linewidth = .65) +
    ggplot2::geom_point(ggplot2::aes(shape = type, fill = type, colour = type), size = 3.1, stroke = .55) +
    ggplot2::scale_shape_manual(values = c("Cohort" = 22, "REML pooled" = 23)) +
    ggplot2::scale_fill_manual(values = c("Cohort" = "#2A7AB9", "REML pooled" = "#CE4459")) +
    ggplot2::scale_colour_manual(values = c("Cohort" = "#2A7AB9", "REML pooled" = "#CE4459")) +
    ggplot2::labs(title = pathway_title, subtitle = subtitle, x = "Log2 fold change (95% CI)", y = NULL) +
    r3_theme(style, 9) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none", plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = 9.5),
      plot.subtitle = ggplot2::element_text(hjust = 0, size = 8),
      axis.text.y = ggplot2::element_text(face = "plain")
    )
  list(plot = p, display_data = d, audits = list(forest_values = d),
       caption = "The representative feature is selected from formal REML padj<0.05 results by cohort-direction concordance, cohort count, maximum cohort standard error, padj and pathway ID. Every valid contributing cohort is retained. Pooled effects/CIs/padj/I2 are copied from saved RC3 results; cohort CIs use the original export formula effect +/-1.96SE. No meta-analysis is refitted and no cohort is removed after feature selection.")
}

r3_plot_pathway_categories <- function(d, panel, style) {
  r3_columns(d, c("comparison", "functional_category", "direction", "Freq"), panel$figure_id)
  # as.data.frame(table(...)) in the integrated producer supplies factor keys.
  # Compare/order their labels, never factor attributes or level codes.
  # Counts and all source rows are preserved; no categories are renamed/merged.
  for (field in c("comparison", "functional_category", "direction"))
    d[[field]] <- as.character(d[[field]])
  comparison_order <- c("T2D vs Control", "NT2D vs NC", "OT2D vs NT2D")
  category_totals <- d |>
    dplyr::group_by(functional_category) |>
    dplyr::summarise(total_features = sum(Freq), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(total_features), functional_category)
  category_order <- category_totals$functional_category
  expected_order <- c("Other/Unclassified", "Nucleotide metabolism", "Amino acid metabolism", "Vitamin and cofactor", "TCA and energy", "Cell wall and LPS", "Carbohydrate metabolism", "Fatty acid metabolism", "One carbon and SAM")
  r3_need(identical(category_order, expected_order) && identical(as.integer(category_totals$total_features), c(45L, 38L, 29L, 16L, 12L, 10L, 8L, 4L, 4L)), "Frozen category totals/order changed")
  direction_order <- c("Positive pooled effect", "Negative pooled effect")
  r3_need(setequal(unique(d$comparison), comparison_order) && all(d$functional_category %in% category_order) &&
    all(d$direction %in% c("Positive", "Negative")) && all(is.finite(d$Freq) & d$Freq >= 0 & d$Freq == as.integer(d$Freq)), "Invalid saved category/count")
  totals <- tapply(d$Freq, d$comparison, sum)
  r3_need(identical(as.integer(totals[comparison_order]), c(13L, 133L, 20L)) && sum(d$Freq) == 166L, "Fig6D totals must remain13/133/20")
  d$direction_display <- ifelse(d$direction == "Positive", direction_order[1], direction_order[2])
  key <- function(x) paste(x$comparison, x$functional_category, x$direction_display, sep = "|")
  r3_need(!anyDuplicated(key(d)), "Saved category grid contains duplicate cells; do not aggregate silently")
  plot_data <- dplyr::transmute(d,
    comparison = factor(comparison, levels = comparison_order),
    functional_category = factor(functional_category, levels = rev(category_order)),
    direction_display = factor(direction_display, levels = direction_order), Freq = Freq) |>
    tidyr::complete(comparison, functional_category, direction_display, fill = list(Freq = 0L))
  plot_data$display_zero_filled <- !key(plot_data) %in% key(d)
  r3_need(nrow(plot_data) == 54L && sum(plot_data$Freq) == 166L, "Display completion changed category counts")
  fills <- c("Positive pooled effect" = "#CE4459", "Negative pooled effect" = "#2A7AB9")
  borders <- colorspace::darken(fills, amount = .30)
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(functional_category, Freq, fill = direction_display, colour = direction_display)) +
    ggplot2::geom_col(width = .55, linewidth = .50) + ggplot2::coord_flip() +
    ggplot2::facet_grid(. ~ comparison, scales = "fixed", space = "fixed") +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_fill_manual(values = fills, limits = direction_order, drop = FALSE, name = "Direction") +
    ggplot2::scale_colour_manual(values = borders, limits = direction_order, drop = FALSE, name = "Direction") +
    ggplot2::guides(colour = "none", fill = ggplot2::guide_legend(order = 1,
      override.aes = list(colour = unname(borders[direction_order])))) +
    ggplot2::labs(title = panel$title, x = "Functional category", y = "Number of significant pathway features") +
    r3_theme(style, 8.5) + ggplot2::theme(
      strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold", size = 9),
      panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom", panel.spacing = grid::unit(1.5, "lines"), axis.text.y = ggplot2::element_text(size = 7.2))
  p <- r3_no_legend_frame(p)
  policy <- r3_layout_policy(p, list(facet_type = class(p$facet)[1],
    fixed_x_scale = !isTRUE(p$facet$params$free$x), fixed_y_scale = !isTRUE(p$facet$params$free$y),
    category_drop_false = identical(p$scales$get_scales("x")$drop, FALSE),
    legend_position = p$theme$legend.position,
    one_common_legend = identical(p$guides$guides$colour, "none") && identical(p$theme$legend.position, "bottom")))
  list(plot = p, display_data = d, audits = list(category_counts = d, completed_plot_data = plot_data,
    category_levels = data.frame(position = seq_along(category_order), category = category_order, total_features = category_totals$total_features),
    comparison_levels = data.frame(position = seq_along(comparison_order), comparison = comparison_order),
    direction_levels = data.frame(position = seq_along(direction_order), direction = direction_order), layout_policy = policy),
    caption = "Saved totals remain13 overall,133 NT2D-vs-NC and20 OT2D-vs-NT2D significant REML pathway features (total166). Missing category/direction display cells are zero-filled only in the separate completed plot-data audit. All three fixed-scale panels share the same nine functional categories and one common bottom legend; original counts are unchanged.")
}

r3_sankey_nodes <- function(values, order, axis, colours, gap, halfwidth) {
  count <- table(factor(values, levels = order))
  r3_need(gap * (length(count) - 1L) < 1, "Configured Sankey node gap is too large")
  heights <- as.numeric(count) / sum(count) * (1 - gap * (length(count) - 1L))
  lower <- c(0, head(cumsum(heights + gap), -1L))
  data.frame(axis = axis, node = order, assignments = as.integer(count), ymin = lower, ymax = lower + heights,
             xmin = axis - halfwidth, xmax = axis + halfwidth, colour = unname(colours[order]), stringsAsFactors = FALSE)
}

r3_plot_pathway_sankey <- function(d, panel, style) {
  full <- d; d <- r3_taxon_pathway_candidates(full)
  genus <- sort(unique(d$genus), method = "radix")
  lookup <- unique(d[, c("taxonomic_contributor", "genus")])
  r3_need(length(genus) == 4L && nrow(lookup) == 6L && !anyDuplicated(lookup$taxonomic_contributor), "SFig3 frozen genus/species mapping differs")
  s <- style$sankey
  hues <- stats::setNames(seq(s$genus_hue_start, s$genus_hue_end, length.out = length(genus) + 1L)[seq_along(genus)], genus)
  genus_colours <- stats::setNames(grDevices::hcl(h = hues, c = s$genus_chroma, l = s$genus_luminance), genus)
  species <- unlist(lapply(genus, function(g) sort(lookup$taxonomic_contributor[lookup$genus == g], method = "radix")), use.names = FALSE)
  species_colours <- character()
  for (g in genus) {
    names_g <- species[species %in% lookup$taxonomic_contributor[lookup$genus == g]]
    h <- unname(hues[g])
    shades <- if (length(names_g) == 1L) grDevices::hcl(h = h, c = 46, l = 70) else
      grDevices::colorRampPalette(c(grDevices::hcl(h = h, c = 26, l = 90), grDevices::hcl(h = h, c = 62, l = 54)))(length(names_g))
    species_colours <- c(species_colours, stats::setNames(shades, names_g))
  }
  pathway <- unique(d$base_pathway_ID)
  r3_need(length(pathway) == 42L, "SFig3 restored42 pathway IDs differ")
  dominant_genus <- vapply(pathway, function(p) which.max(table(factor(d$genus[d$base_pathway_ID == p], levels = genus))), integer(1))
  frequency <- vapply(pathway, function(p) sum(d$base_pathway_ID == p), integer(1))
  pathway <- pathway[order(dominant_genus, -frequency, pathway, method = "radix")]
  pathway_colours <- stats::setNames(
    grDevices::hcl(h = seq(s$pathway_hue_start, s$pathway_hue_end, length.out = length(pathway) + 1L)[seq_along(pathway)],
                   c = s$pathway_chroma, l = s$pathway_luminance), pathway
  )
  axes <- list(d$taxonomic_contributor, d$genus, d$base_pathway_ID)
  orders <- list(species, genus, pathway)
  palettes <- list(species_colours, genus_colours, pathway_colours)
  nodes <- lapply(1:3, function(k) r3_sankey_nodes(axes[[k]], orders[[k]], k, palettes[[k]], s$node_gap, s$node_halfwidth))
  ribbons <- list(); edges <- list(); index <- 0L
  for (k in 1:2) {
    left_used <- stats::setNames(rep(0, nrow(nodes[[k]])), nodes[[k]]$node)
    right_used <- stats::setNames(rep(0, nrow(nodes[[k + 1L]])), nodes[[k + 1L]]$node)
    for (i in seq_len(nrow(d))) {
      a <- axes[[k]][i]; b <- axes[[k + 1L]][i]
      left <- nodes[[k]][nodes[[k]]$node == a, , drop = FALSE]
      right <- nodes[[k + 1L]][nodes[[k + 1L]]$node == b, , drop = FALSE]
      lh <- (left$ymax - left$ymin) / left$assignments
      rh <- (right$ymax - right$ymin) / right$assignments
      y1 <- left$ymin + left_used[a] * lh; y2 <- right$ymin + right_used[b] * rh
      left_used[a] <- left_used[a] + 1; right_used[b] <- right_used[b] + 1
      t <- seq(0, 1, length.out = s$polygon_points); smooth <- 3 * t^2 - 2 * t^3
      x <- k + s$node_halfwidth + t * (1 - 2 * s$node_halfwidth)
      index <- index + 1L
      ribbons[[index]] <- data.frame(
        x = c(x, rev(x)),
        y = c(y1 + (y2 - y1) * smooth, rev(y1 + lh + (y2 + rh - y1 - lh) * smooth)),
        ribbon = index, colour = unname(genus_colours[d$genus[i]])
      )
      edges[[index]] <- data.frame(ribbon = index, feature = d$feature[i], source_axis = k,
                                   source = a, target = b, genus = d$genus[i], colour = unname(genus_colours[d$genus[i]]))
    }
  }
  node_df <- do.call(rbind, nodes)
  node_df$label_x <- ifelse(node_df$axis == 1, node_df$xmin - .045, ifelse(node_df$axis == 3, node_df$xmax + .045, node_df$axis))
  node_df$hjust <- ifelse(node_df$axis == 1, 1, ifelse(node_df$axis == 3, 0, .5))
  node_df$face <- ifelse(node_df$axis < 3, s$label_fontface, "plain")
  node_df$label_size <- ifelse(node_df$axis < 3, 3.2, 2.25)
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = do.call(rbind, ribbons), ggplot2::aes(x, y, group = ribbon, fill = colour), alpha = s$ribbon_alpha, linewidth = 0) +
    ggplot2::geom_rect(data = node_df, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = colour),
                       colour = "white", linewidth = .40, alpha = s$node_alpha) +
    ggplot2::geom_text(data = node_df, ggplot2::aes(label_x, (ymin + ymax) / 2, label = node, hjust = hjust, fontface = face, size = label_size),
                       family = style$font$family, colour = "#263238", show.legend = FALSE) +
    ggplot2::scale_fill_identity() + ggplot2::scale_size_identity() +
    # Coordinate-only limits keep every full label in the layer data; clipping is
    # disabled so the complete species labels can occupy reserved physical margins.
    ggplot2::scale_x_continuous(breaks = 1:3, labels = c("Species", "Genus", "MetaCyc pathway ID"), expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(s$display_xmin, s$display_xmax), clip = "off") +
    ggplot2::scale_y_continuous(breaks = NULL, expand = c(.01, .01)) +
    ggplot2::labs(title = panel$title, x = NULL, y = NULL) +
    ggplot2::theme_void(base_family = style$font$family, base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = .5, size = 11),
      axis.text.x = ggplot2::element_text(face = "bold", colour = "black", size = 9),
      plot.margin = ggplot2::margin(5, s$label_margin_right_mm, 5, s$label_margin_left_mm, unit = "mm")
    )
  list(plot = p, display_data = d, source_features = d$feature,
       audits = list(nodes = node_df, edges = do.call(rbind, edges),
         label_layout = data.frame(taxon_fontface = s$label_fontface, left_margin_mm = s$label_margin_left_mm,
           right_margin_mm = s$label_margin_right_mm, coordinate_clip = "off", taxonomy_axes = "Species;Genus;MetaCyc pathway ID"),
         candidate_selection_audit = r3_candidate_filter_audit(full)),
       caption = "All 64 eligible taxon-stratified pathway features were included. The selection contains six species, four genera and42 MetaCyc pathway IDs. The selection is identical to Figures5B and5C, with no top-N truncation. Every retained feature contributes one assignment to each transition. Ribbon width represents the number of retained taxon-stratified feature assignments and does not represent microbial abundance, metabolic flux, or causal contribution.")
}

r3_build_panel <- function(d, panel, style) {
  switch(
    panel$plot_family,
    pathway_volcano = r3_plot_pathway_volcano(d, panel, style),
    pathway_taxonomy_bar = r3_plot_pathway_taxonomy_bar(d, panel, style),
    pathway_bubble = r3_plot_pathway_bubble(d, panel, style),
    pathway_overlap = r3_plot_pathway_overlap(d, panel, style),
    pathway_heatmap = r3_plot_pathway_heatmap(d, panel, style),
    pathway_forest = r3_plot_pathway_forest(d, panel, style),
    pathway_categories = r3_plot_pathway_categories(d, panel, style),
    pathway_sankey = r3_plot_pathway_sankey(d, panel, style),
    stop("Unregistered Round3 plot family: ", panel$plot_family, call. = FALSE)
  )
}
