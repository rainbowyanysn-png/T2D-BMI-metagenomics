run_metadata_medication_audit <- function(data, run, config) {
  require_project_packages(c("dplyr", "stringr"))
  build_metadata_availability_audit(data, run, config, medication_mapper = medication_model_covariate)
  dictionary <- utils::read.csv(file.path(run$root, "config/medication_dictionary.csv"), check.names = FALSE)
  samples <- build_medication_sample_audit(as.data.frame(data$metadata), dictionary)
  assert_equal(nrow(samples), 1186L, "Medication audit samples")
  assert_true(setequal(samples$sample_id, data$metadata$sample_id) && !anyDuplicated(samples$sample_id), "Medication audit sample identity")
  write <- function(x, name) safe_write_csv(x, path_in_run(run, "tables/11_metadata_medication_audit/medication", name))
  write(samples, "BMI25_Medication_Audit_SampleManifest.csv")
  crosswalk <- dplyr::count(samples, treatment_raw, treatment_normalized, medication_record_status, medication_provisional_class, provisional_metformin_status, name = "n_samples", sort = TRUE)
  write(crosswalk, "BMI25_Medication_Audit_RawTreatment_Crosswalk.csv")
  for (band in c("BMI25_group", "BMI3_group")) {
    fields <- c("study_name", "disease", band, "medication_record_status", "medication_provisional_class", "provisional_metformin_status")
    counts <- dplyr::count(samples, dplyr::across(dplyr::all_of(fields)), name = "n_samples", sort = TRUE)
    write(counts, paste0("BMI25_Medication_Audit_Counts_ByStudy_Disease_", sub("_group", "", band), ".csv"))
  }
  summarize <- function(fields) {
    grouped <- dplyr::group_by(samples, dplyr::across(dplyr::all_of(fields)))
    out <- dplyr::summarise(grouped,
      n_samples = dplyr::n(),
      n_unknown_or_not_recorded = sum(medication_record_status == "Unknown_or_not_recorded"),
      n_explicit_no_antidiabetic_medication = sum(medication_record_status == "Explicit_no_antidiabetic_medication"),
      n_literal_no_treatment_code = sum(literal_no_treatment_code),
      n_recorded_text = sum(medication_record_status == "Recorded_text"),
      n_metformin_detected = sum(provisional_metformin_status == "Metformin_detected"),
      n_multiple_recognised_classes = sum(medication_provisional_class == "Multiple_recognised_antidiabetic_classes"),
      n_recorded_text_manual_check_required = sum(provisional_metformin_status == "No_metformin_keyword_detected_manual_check_required"),
      n_generic_antidiabetic_treatment_labels = sum(contains_generic_antidiabetic_label), .groups = "drop")
    out <- dplyr::mutate(out,
      pct_unknown_or_not_recorded = round(100 * n_unknown_or_not_recorded / n_samples, 1),
      pct_metformin_detected = round(100 * n_metformin_detected / n_samples, 1),
      pct_literal_no_treatment_code = round(100 * n_literal_no_treatment_code / n_samples, 1),
      pct_multiple_recognised_classes = round(100 * n_multiple_recognised_classes / n_samples, 1))
    dplyr::arrange(out, dplyr::across(dplyr::all_of(fields)))
  }
  for (band in c("", "BMI25_group", "BMI3_group")) {
    fields <- c("study_name", "disease", if (nzchar(band)) band)
    write(summarize(fields), paste0("BMI25_Medication_Audit_Completeness_ByStudy_Disease", if (nzchar(band)) paste0("_", sub("_group", "", band)), ".csv"))
  }
  invisible(samples)
}
