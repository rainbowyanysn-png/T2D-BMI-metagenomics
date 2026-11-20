run_import_harmonize_qc <- function(run, config, expected) {
  data <- import_harmonized_data(run, config, expected)
  dir.create(run$derived_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(data, file.path(run$derived_dir, "harmonized_inputs.rds"))
  safe_write_csv(data$metadata, path_in_run(run, "tables", "01_sample_manifest.csv"))
  cells <- as.data.frame(table(study_name = data$metadata$study_name, analysis_group = data$metadata$BMI25_group))
  names(cells)[3] <- "n_samples"
  safe_write_csv(cells, path_in_run(run, "tables", "01_BMI25_cohort_four_group_counts.csv"))
  safe_write_csv(data.frame(feature = rownames(data$species), prevalence = rowMeans(data$species > 0)), path_in_run(run, "tables", "01_species_prevalence.csv"))
  safe_write_csv(data.frame(feature = rownames(data$pathways), prevalence = rowMeans(data$pathways > 0)), path_in_run(run, "tables", "01_pathway_prevalence.csv"))
  data
}
