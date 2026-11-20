r3_theme <- function(style, base_size = NULL) {
  if (is.null(base_size)) base_size <- style$font$base_size
  ggplot2::theme_bw(base_size = base_size, base_family = style$font$family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      axis.title = ggplot2::element_text(face = "bold", colour = "black"),
      axis.text = ggplot2::element_text(colour = "black"),
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.35)
    )
}

r3_html_symbols <- function(x) {
  x <- gsub("&beta;", "β", x, fixed = TRUE)
  x <- gsub("&alpha;", "α", x, fixed = TRUE)
  x <- gsub("&gamma;", "γ", x, fixed = TRUE)
  x
}


r3_as_logical <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}


r3_candidate_filter_audit <- function(d) {
  required <- c("feature", "coef", "qval", "visualization_candidate", "taxonomic_level", "genus", "taxonomic_contributor", "base_pathway_ID", "pathway_description")
  r3_need(all(required %in% names(d)), paste("Complete candidate input missing:", paste(setdiff(required, names(d)), collapse = ", ")))
  r3_need(nrow(d) == 5945L && !anyDuplicated(d$feature), "Complete pooled candidate input must retain5945 unique features")
  coef <- as.numeric(d$coef); qval <- as.numeric(d$qval)
  r3_need(all(is.finite(coef)) && all(is.finite(qval) & qval >= 0 & qval <= 1), "Invalid frozen candidate coefficient/qval")
  named <- function(x) !is.na(x) & nzchar(as.character(x))
  a <- data.frame(feature = d$feature, coef = coef, qval = qval,
    saved_visualization_candidate = r3_as_logical(d$visualization_candidate),
    is_species = !is.na(d$taxonomic_level) & d$taxonomic_level == "species",
    qval_below_001 = qval < .01, absolute_coef_above_1 = abs(coef) > 1,
    excludes_unintegrated_unmapped = !grepl("UNINTEGRATED|UNMAPPED", d$feature, ignore.case = TRUE),
    named_genus = named(d$genus), named_species = named(d$taxonomic_contributor),
    named_pathway = named(d$base_pathway_ID), stringsAsFactors = FALSE)
  conditions <- a[setdiff(names(a), c("feature", "coef", "qval"))]
  r3_need(!anyNA(conditions), "Missing candidate eligibility flag")
  a$included <- rowSums(conditions) == ncol(conditions)
  a
}

r3_assert_candidate_counts <- function(d) {
  r3_need(nrow(d) == 64L && !anyDuplicated(d$feature), "Expected64 unique eligible taxon-stratified features")
  r3_need(length(unique(d$taxonomic_contributor)) == 6L && length(unique(d$genus)) == 4L && length(unique(d$base_pathway_ID)) == 42L, "Eligible feature family must have6 species/4 genera/42 pathways")
  r3_need(sum(d$coef > 0) == 37L && sum(d$coef < 0) == 27L, "Eligible effects must have37 positive/27 negative coefficients")
  species <- c("Bacteroides fragilis", "Escherichia coli", "Eubacterium eligens", "Eubacterium eligens CAG 72", "Faecalibacterium prausnitzii", "Bacteroides fragilis CAG 47")
  r3_need(setequal(unique(d$taxonomic_contributor), species) && identical(as.integer(table(factor(d$taxonomic_contributor, levels = species))), c(19L, 12L, 12L, 8L, 7L, 6L)), "Eligible species counts must be19/12/12/8/7/6")
  invisible(TRUE)
}

# Independent base-R reconstruction for validation; renderer uses dplyr order.
r3_expected_candidate_data <- function(d) {
  a <- r3_candidate_filter_audit(d); selected <- d[a$included, , drop = FALSE]
  selected <- selected[order(selected$qval, -abs(selected$coef), selected$feature, method = "radix"), , drop = FALSE]
  r3_assert_candidate_counts(selected); selected
}

r3_candidate_selection_rule <- function() c(
  "Candidates were restricted to species-stratified pathway features with",
  "MaAsLin2 full-model qval < 0.01 and |coefficient| > 1.",
  "UNINTEGRATED and UNMAPPED entries were excluded.",
  "All 64 eligible features were retained; no top-N truncation was applied.",
  "The saved visualization_candidate flag must be TRUE; genus, species contributor and base pathway ID must be nonempty.",
  "Order: qval ascending, absolute coefficient descending, then original feature ID ascending.",
  "All original statistical fields are retained; no P value or q value was recalculated.")

