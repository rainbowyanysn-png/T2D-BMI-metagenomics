# Descriptive display preparation only. No MaAsLin2/LMM/DESeq2/REML fitting.
publication_phylum_display <- function(full, saved, metadata, selected, tolerance) {
  assert_true(is.matrix(full) && is.numeric(full) && all(is.finite(full)) && all(full >= 0),
    "Full phylum abundance must be a finite nonnegative numeric matrix")
  ids <- as.character(metadata$sample_id)
  assert_true(identical(colnames(full), ids) && !anyDuplicated(ids) && !anyDuplicated(rownames(full)),
    "Full phylum sample/feature order differs from current-run metadata")
  taxa <- as.character(selected$feature)
  assert_true(length(taxa) == 4L && !anyDuplicated(taxa) && all(taxa %in% rownames(full)),
    "The saved top-four phylum identities must be preserved")
  assert_columns(saved, c("taxonomic_level", "sample_id", "feature", "raw_abundance", "analysis_group"), "Saved top-four abundance")
  saved <- saved[saved$taxonomic_level == "phylum", ]
  assert_true(nrow(saved) == length(ids) * 4L && !anyDuplicated(saved[, c("sample_id", "feature")]),
    "Saved phylum display key is incomplete or duplicated")
  assert_true(setequal(unique(saved$feature), taxa), "Saved top-four identities differ from the selected-taxa manifest")
  original <- unname(full[cbind(match(saved$feature, rownames(full)), match(saved$sample_id, ids))])
  assert_true(!anyNA(original) && max(abs(original - saved$raw_abundance)) <= tolerance,
    "Raw phylum reconstruction differs from saved abundant-taxa input; do not alter statistical inputs")
  j <- match(saved$sample_id, ids)
  assert_true(identical(as.character(saved$analysis_group), as.character(metadata$BMI25_group[j])),
    "Saved phylum groups differ from the frozen metadata")
  rest <- setdiff(rownames(full), taxa)
  # Crucially: sum the actual omitted phyla, never 100 - selected taxa.
  other <- if (length(rest)) colSums(full[rest, , drop = FALSE]) else rep(0, length(ids))
  display <- rbind(saved[, c("sample_id", "analysis_group", "feature", "raw_abundance")],
    data.frame(sample_id = ids, analysis_group = metadata$BMI25_group,
      feature = "Other phyla", raw_abundance = unname(other)))
  display$abundance_unit <- "Provider percentage relative abundance; not reclosed"
  list(display = display, full = data.frame(phylum = rownames(full), full, check.names = FALSE),
    totals = data.frame(sample_id = ids, full_phylum_total = unname(colSums(full)),
      selected_four_total = unname(colSums(full[taxa, , drop = FALSE])),
      actual_other_phyla_total = unname(other)),
    alignment = data.frame(check = "all_saved_top4_sample_values", n_values = nrow(saved),
      max_absolute_difference = max(abs(original - saved$raw_abundance)), tolerance = tolerance))
}

prepare_publication_display_data <- function(data, run, config) {
  require_project_packages(c("curatedMetagenomicData", "SummarizedExperiment"))
  check_frozen_versions(config)
  md <- as.data.frame(data$metadata)
  assert_true(nrow(md) == 1186L && identical(colnames(data$species), md$sample_id),
    "Publication display data must use the frozen 1186-sample panel")
  saved <- read_current_result(run, "tables/04b_abundant_taxa/abundant_taxa_sample_level_data.csv")
  selected <- read_current_result(run, "tables/04b_abundant_taxa/selected_top4_phyla.csv")
  retrieve <- prepare_frozen_resources(run, config)
  long <- retrieve(md, dataType = "relative_abundance", rownames = "long")
  a <- SummarizedExperiment::assay(long)
  assert_true(setequal(colnames(a), md$sample_id), "Raw composition sample set changed")
  a <- a[, md$sample_id, drop = FALSE]
  assert_true(all(grepl("(^|[|])p__", rownames(a))), "Missing explicit provider p__ phylum field")
  phylum <- sub("[|;].*$", "", sub("^.*p__", "", rownames(a)))
  assert_true(!anyNA(phylum) && all(nzchar(phylum)), "Invalid provider phylum annotation")
  full <- rowsum(a, phylum, reorder = TRUE)
  prepared <- publication_phylum_display(full, saved, md, selected, config$project$numeric_tolerance)
  directory <- "tables/01c_publication_display"
  safe_write_csv(prepared$display, path_in_run(run, directory, "BMI25_phylum_composition_sample_level.csv"))
  safe_write_csv(prepared$full, path_in_run(run, directory, "BMI25_all_phylum_abundance.csv"))
  safe_write_csv(prepared$totals, path_in_run(run, directory, "BMI25_phylum_total_audit.csv"))
  safe_write_csv(prepared$alignment, path_in_run(run, directory, "BMI25_top4_alignment.csv"))
  # Species composition is exactly the saved retained 200-feature panel.
  # 'saved' holds both ranks from the original abundant-taxa table.
  ids <- md$sample_id
  top8 <- read_current_result(run, "tables/04b_abundant_taxa/selected_top8_species.csv")$feature
  assert_true(length(top8) == 8L && all(top8 %in% rownames(data$species)), "Saved top-eight species panel changed")
  species <- saved[saved$taxonomic_level == "species", c("sample_id", "analysis_group", "feature", "raw_abundance")]
  assert_true(setequal(unique(species$feature), top8), "Saved top-eight identities differ from the selected-taxa manifest")
  ij <- cbind(match(species$feature, rownames(data$species)), match(species$sample_id, ids))
  assert_true(!anyNA(ij) && nrow(species) == 1186L * 8L && !anyDuplicated(species[, c("sample_id", "feature")]),
    "Saved top-eight species keys changed")
  assert_true(max(abs(data$species[ij] - species$raw_abundance)) <= config$project$numeric_tolerance,
    "Saved top-eight abundance differs from the retained panel")
  assert_true(identical(as.character(species$analysis_group), as.character(md$BMI25_group[ij[, 2]])),
    "Saved top-eight groups differ from the frozen metadata")
  other <- colSums(data$species[setdiff(rownames(data$species), top8), , drop = FALSE])
  species <- rbind(species, data.frame(sample_id = ids, analysis_group = md$BMI25_group,
    feature = "Other retained features", raw_abundance = unname(other)))
  species$abundance_unit <- "Provider percentage within retained 200-feature panel; not reclosed"
  safe_write_csv(species, path_in_run(run, directory, "BMI25_species_composition_sample_level.csv"))
}
