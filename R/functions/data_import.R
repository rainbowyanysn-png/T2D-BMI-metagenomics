import_harmonized_data <- function(run, config, expected) {
  require_project_packages(c("curatedMetagenomicData", "ExperimentHub", "SummarizedExperiment", "dplyr", "tibble"))
  check_frozen_versions(config)
  retrieve_samples <- prepare_frozen_resources(run, config)
  raw_metadata <- curatedMetagenomicData::sampleMetadata
  candidate <- raw_metadata |>
    dplyr::filter(disease %in% unlist(config$source_data$disease_values), !is.na(BMI), body_site == config$source_data$body_site) |>
    dplyr::mutate(disease = dplyr::case_when(disease %in% c("healthy", "control") ~ "control", disease == "T2D" ~ "T2D", TRUE ~ NA_character_), analysis_group = factor(dplyr::case_when(disease == "T2D" & BMI < 25 ~ "Non_Overweight_T2D", disease == "T2D" & BMI >= 25 ~ "Overweight_T2D", disease == "control" & BMI < 25 ~ "Non_Overweight_Control", disease == "control" & BMI >= 25 ~ "Overweight_Control", TRUE ~ NA_character_), levels = c("Non_Overweight_Control", "Overweight_Control", "Non_Overweight_T2D", "Overweight_T2D"))) |>
    dplyr::filter(!is.na(disease), !is.na(analysis_group))
  species_se <- retrieve_samples(candidate, dataType = "relative_abundance", rownames = "short")
  species_raw <- SummarizedExperiment::assay(species_se)
  raw_returned_metadata <- as.data.frame(SummarizedExperiment::colData(species_se)) |>
    tibble::rownames_to_column("sample_id")
  metadata <- raw_returned_metadata |>
    dplyr::mutate(disease = factor(ifelse(disease %in% c("healthy", "control"), "control", "T2D"), levels = c("control", "T2D")), BMI = suppressWarnings(as.numeric(BMI)), sex = factor(tolower(as.character(gender))), age = suppressWarnings(as.numeric(age))) |>
    dplyr::group_by(study_name) |>
    dplyr::mutate(age = dplyr::if_else(is.na(age), median(age, na.rm = TRUE), age)) |>
    dplyr::ungroup()
  cohort_counts <- metadata |>
    dplyr::count(study_name, analysis_group, name = "n") |>
    dplyr::filter(n >= config$cohorts$min_samples_per_group) |>
    dplyr::mutate(group_type = ifelse(grepl("Control", analysis_group), "Control", "T2D")) |>
    dplyr::group_by(study_name) |>
    dplyr::filter(dplyr::n_distinct(group_type) >= 2L) |>
    dplyr::ungroup()
  keep_cohorts <- sort(unique(as.character(cohort_counts$study_name)))
  metadata <- metadata |>
    dplyr::filter(as.character(study_name) %in% keep_cohorts) |>
    dplyr::mutate(BMI25_group = factor(as.character(analysis_group), levels = analysis_group_levels), BMI23_threshold = ifelse(as.character(study_name) %in% unlist(config$filters$bmi23_cohorts), 23, 25), BMI23_group = factor(dplyr::case_when(disease == "control" & BMI < BMI23_threshold ~ "Non_Overweight_Control", disease == "control" ~ "Overweight_Control", disease == "T2D" & BMI < BMI23_threshold ~ "Non_Overweight_T2D", TRUE ~ "Overweight_T2D"), levels = analysis_group_levels))
  species_raw <- species_raw[, metadata$sample_id, drop = FALSE]
  species_keep <- rowMeans(species_raw > 0) >= config$filters$global_prevalence
  species <- species_raw[species_keep, , drop = FALSE]; rownames(species) <- make.names(rownames(species))
  assert_true(!anyDuplicated(rownames(species)), "Species IDs collide after historical make.names; audit before continuing.")
  pathway_se <- retrieve_samples(metadata, dataType = "pathway_abundance")
  pathway_raw <- SummarizedExperiment::assay(pathway_se)[, metadata$sample_id, drop = FALSE]
  pathway_keep <- rowMeans(pathway_raw > 0) >= config$filters$global_prevalence
  pathways <- pathway_raw[pathway_keep, , drop = FALSE]; rownames(pathways) <- make.names(rownames(pathways), unique = TRUE)
  pathway_cohorts <- list()
  extraction_audit <- list()
  for (study in unique(as.character(metadata$study_name))) {
    md <- as.data.frame(metadata[metadata$study_name == study, , drop = FALSE])
    rownames(md) <- NULL
    se <- retrieve_samples(md, dataType = "pathway_abundance", counts = TRUE)
    raw <- SummarizedExperiment::assay(se)
    extracted_metadata <- as.data.frame(SummarizedExperiment::colData(se))
    assert_true(identical(colnames(raw), rownames(extracted_metadata)), paste(study, "original pathway extraction order"))
    pathway_cohorts[[study]] <- list(matrix = raw, metadata = extracted_metadata)
    common <- intersect(rownames(raw), rownames(pathway_raw))
    extraction_audit[[study]] <- data.frame(study_name = study, n_features_individual = nrow(raw), n_features_global = nrow(pathway_raw), n_samples = ncol(raw), identical_feature_order = identical(rownames(raw), rownames(pathway_raw)), max_value_difference_common = max(abs(raw[common, , drop = FALSE] - pathway_raw[common, colnames(raw), drop = FALSE]), na.rm = TRUE))
  }
  safe_write_csv(do.call(rbind, extraction_audit), path_in_run(run, "logs", "pathway_cohort_extraction_audit.csv"))
  assert_true(identical(colnames(species), metadata$sample_id), "Species sample order mismatch.")
  assert_true(identical(colnames(pathways), metadata$sample_id), "Pathway sample order mismatch.")
  assert_true(!anyDuplicated(metadata$sample_id), "Duplicate sample IDs.")
  assert_equal(nrow(metadata), expected$inputs$n_samples, "analytic samples")
  assert_equal(sum(metadata$disease == "T2D"), expected$inputs$n_t2d, "T2D samples")
  assert_equal(sum(metadata$disease == "control"), expected$inputs$n_control, "Control samples")
  assert_true(identical(sort(unique(as.character(metadata$study_name))), sort(unlist(config$cohorts$names))), "Frozen cohort list mismatch.")
  assert_equal(nrow(species), expected$inputs$n_species, "species after prevalence filtering")
  assert_equal(nrow(pathways), expected$inputs$n_pathways, "pathways after prevalence filtering")
  assert_equal(sum(metadata$disease == "T2D" & metadata$BMI < 25), expected$inputs$bmi25_t2d_non_overweight, "BMI25 NT2D")
  assert_equal(sum(metadata$disease == "T2D" & metadata$BMI >= 25), expected$inputs$bmi25_t2d_overweight, "BMI25 OT2D")
  assert_equal(sum(metadata$disease == "T2D" & metadata$BMI < metadata$BMI23_threshold), expected$inputs$bmi23_t2d_non_overweight, "BMI23 NT2D")
  assert_equal(sum(metadata$disease == "T2D" & metadata$BMI >= metadata$BMI23_threshold), expected$inputs$bmi23_t2d_overweight, "BMI23 OT2D")
  list(metadata = metadata, species = species, pathways = pathways, pathways_raw = pathway_raw, pathway_cohorts = pathway_cohorts, raw_metadata = raw_returned_metadata[match(metadata$sample_id, raw_returned_metadata$sample_id), , drop = FALSE], source = "curatedMetagenomicData")
}
