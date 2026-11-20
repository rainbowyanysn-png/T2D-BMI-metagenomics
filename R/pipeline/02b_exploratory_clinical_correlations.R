run_exploratory_clinical_correlations <- function(data, run, config) {
  directory <- "tables/02b_exploratory_clinical_correlations"
  raw <- as.character(data$metadata$sex)
  counts <- as.data.frame(table(raw_value = raw, useNA = "ifany"))
  names(counts)[2] <- "n_samples"
  counts$normalized_value <- tolower(trimws(as.character(counts$raw_value)))
  counts$numeric_code <- ifelse(counts$normalized_value == "female", 0,
                                ifelse(counts$normalized_value == "male", 1, NA))
  safe_write_csv(counts, path_in_run(run, directory, "Clinical_correlation_sex_coding_audit.csv"))
  message("Clinical correlations sex audit: ", paste(counts$raw_value, counts$n_samples, sep = "=", collapse = ", "))
  main_path <- current_result_path(run, "tables/02_primary_BMI25_species/Total_T2D_vs_Control_all.csv")
  before <- sha256_file(main_path)
  main <- utils::read.csv(main_path, check.names = FALSE)
  result <- clinical_correlation_data(data, main, config)
  safe_write_csv(result$selected, path_in_run(run, directory, "Clinical_correlation_selected_top20_features.csv"))
  safe_write_csv(result$all, path_in_run(run, directory, "Clinical_correlations_all_60_tests.csv"))
  safe_write_csv(result$benchmark_checks, path_in_run(run, directory, "Clinical_correlation_benchmark_checks.csv"))
  stopifnot(all(result$benchmark_checks$pass), before == sha256_file(main_path))
  base <- list(root = run$root, id = run$id, dir = file.path(run$root, "output", run$id))
  state <- latest_pipeline_state(base$dir)
  input <- checkpoint_record(base, state$completed$inputs)
  summary_path <- current_result_path(run, "tables/02_primary_BMI25_species_summary.csv")
  summary <- utils::read.csv(summary_path, check.names = FALSE)
  expected <- c(Total_T2D_vs_Control = 69L, NT2D_vs_NC_BMI25 = 32L, OT2D_vs_OC_BMI25 = 10L,
                OC_vs_NC_BMI25 = 4L, OT2D_vs_NT2D_BMI25 = 5L)
  stopifnot(!anyDuplicated(summary$comparison_id),
            identical(as.integer(summary$n_formal_q_lt_005[match(names(expected), summary$comparison_id)]),
                      unname(expected)))
  provenance <- data.frame(
    role = c("Current-run harmonized inputs", "Current total-T2D MaAsLin2 result", "Current five-contrast species summary"),
    project_file = c(input$input_rds, substring(main_path, nchar(run$root) + 2L),
                     substring(summary_path, nchar(run$root) + 2L)),
    sha256 = c(input$input_sha256, before, sha256_file(summary_path)),
    statistical_columns = c("Saved untransformed relative abundance; BMI, age, sex",
                            "Original full-model qval for selection", "Original full-model significant counts"))
  safe_write_csv(provenance, path_in_run(run, directory, "Clinical_correlation_input_manifest.csv"))
  validation <- data.frame(check = c("n_samples_1186", "species_panel_200", "top_features_20",
    "tests_60", "complete_BMI_age", "unique_feature_trait", "one_BH60_family",
    "six_numeric_benchmarks", "main_result_hash_unchanged", "five_main_counts_unchanged"), pass = TRUE)
  safe_write_csv(validation, path_in_run(run, directory, "Clinical_correlation_validation_summary.csv"))
  invisible(result)
}
