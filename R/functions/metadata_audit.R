# Rebuilt raw/model metadata audit; preserves original field and search definitions.
build_metadata_availability_audit <- function(data, run, config, medication_mapper) {
  assert_true(is.function(medication_mapper), "Pass the original medication mapping function explicitly")
  write_csv_utf8 <- function(x, filename) safe_write_csv(x, path_in_run(run, "tables/11_metadata_medication_audit/availability", filename))
is_blank <- function(x) {
  is.na(x) | trimws(as.character(x)) == ""
}

# "Unknown" remains distinct from missing. It is not counted as documented use
# or documented absence. The exact raw value is retained in the sample-level file.
is_unknown_text <- function(x) {
  value <- tolower(trimws(as.character(x)))
  !is_blank(x) & value %in% c(
    "unknown", "not known", "not available", "not reported", "missing",
    "na", "n/a", "unspecified", "not specified", "nd", "not assessed",
    "not measured", "not recorded", "not collected"
  )
}

field_status <- function(x) {
  out <- rep("documented", length(x))
  out[is_blank(x)] <- "missing_or_blank"
  out[is_unknown_text(x)] <- "unknown_or_unspecified"
  out
}

as_audit_text <- function(x) {
  if (is.list(x)) {
    return(vapply(x, function(value) {
      if (is.null(value) || !length(value) || all(is.na(value))) return(NA_character_)
      if (is.atomic(value) && length(value) == 1L) return(as.character(value))
      paste(capture.output(dput(value)), collapse = " ")
    }, character(1)))
  }
  as.character(x)
}

safe_col <- function(data, column_name) {
  if (is.na(column_name) || !nzchar(column_name) || !column_name %in% names(data)) {
    return(rep(NA_character_, nrow(data)))
  }
  as_audit_text(data[[column_name]])
}

bind_rows_or_empty <- function(rows, template) {
  if (!length(rows)) return(template)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Only the identifying groups are used as split keys. Descriptive columns can
# contain NA (an intentionally unmapped concept) without losing those rows.
split_audit <- function(data, columns) {
  if (anyNA(data[, columns, drop = FALSE])) stop("Missing audit grouping key.")
  key <- do.call(interaction, c(unname(data[, columns, drop = FALSE]),
                               list(drop = TRUE, lex.order = TRUE)))
  split(data, key, drop = TRUE)
}

validation_rows <- list()
assert_audit <- function(check_name, condition, detail) {
  validation_rows[[length(validation_rows) + 1L]] <<- data.frame(
    check = check_name, passed = isTRUE(condition), detail = detail
  )
  if (!isTRUE(condition)) {
    write_csv_utf8(do.call(rbind, validation_rows), "21_ValidationChecks.csv")
    stop("Audit validation failed: ", check_name, ". ", detail)
  }
}

percent <- function(numerator, denominator) {
  ifelse(denominator == 0L, NA_real_, round(100 * numerator / denominator, 1))
}

normalize_disease <- function(x) {
  value <- tolower(trimws(as.character(x)))
  out <- rep("other_or_unresolved", length(value))
  out[value == "t2d"] <- "T2D"
  out[value %in% c("control", "healthy")] <- "control"
  out
}
  expected_studies <- unlist(config$cohorts$names)
  expected_counts <- data.frame(study_name = rep(expected_studies, each = 2L), disease = rep(c("T2D", "control"), length(expected_studies)), n_samples = c(53L,43L,502L,153L,170L,174L,19L,18L,16L,38L))
  meta <- medication_mapper(as.data.frame(data$metadata))
  meta$age_numeric <- as.numeric(as.character(meta$age)); meta$BMI_numeric <- as.numeric(as.character(meta$BMI))
  meta$BMI_centered <- meta$BMI_numeric - mean(meta$BMI_numeric, na.rm = TRUE)
  meta$T2D_x_BMI_continuous <- as.integer(meta$disease == "T2D") * meta$BMI_centered
  meta$treatment_lower <- tolower(meta$treatment)
  raw_meta <- as.data.frame(data$raw_metadata)
  assert_equal(meta$sample_id, raw_meta$sample_id, "Raw/model audit sample order")
  assert_equal(nrow(meta), 1186L, "Metadata audit sample count")
  input_columns <- names(meta); input_n_columns <- ncol(meta)
  meta$disease_audit <- factor(normalize_disease(meta$disease), levels = c("control", "T2D"))
  # These provenance labels match the historical table schema; both objects are
  # freshly reconstructed above, never restored from an old workspace.
  raw_source_object <- "meta_df_raw"
  database_columns <- names(curatedMetagenomicData::sampleMetadata)
  all_columns <- union(input_columns, names(raw_meta))
  column_source <- function(column) {
    if (column %in% names(raw_meta)) return(raw_source_object)
    if (column %in% input_columns) return("tax_meta")
    "no_predefined_mapping"
  }
  audit_col <- function(column) {
    if (!is.na(column) && column %in% names(raw_meta)) return(safe_col(raw_meta, column))
    safe_col(meta, column)
  }
  field_dictionary <- utils::read.csv(file.path(run$root, "config/metadata_field_dictionary.csv"), check.names = FALSE, na.strings = c("", "NA"))
mapped_fields <- field_dictionary$metadata_column[!is.na(field_dictionary$metadata_column)]
assert_audit("all_predefined_fields_exist", all(mapped_fields %in% all_columns),
  paste("Absent mapped fields:", paste(setdiff(mapped_fields, all_columns), collapse = ", ")))
field_dictionary$field_match_status <- ifelse(is.na(field_dictionary$metadata_column),
  "no_predefined_exact_mapping", "exact_column_name_verified_in_saved_metadata")
field_dictionary$value_source_object <- vapply(field_dictionary$metadata_column,
  function(column) if (is.na(column)) "no_predefined_mapping" else column_source(column), character(1))

write_csv_utf8(field_dictionary, "01_FieldDictionary.csv")

# ---- 5. Long sample-level audit ------------------------------------------------
long_list <- lapply(seq_len(nrow(field_dictionary)), function(i) {
  field <- field_dictionary[i, , drop = FALSE]
  value <- audit_col(field$metadata_column)
  status <- if (is.na(field$metadata_column)) {
    rep("no_predefined_exact_mapping", nrow(meta))
  } else {
    field_status(value)
  }
  data.frame(
    sample_id = meta$sample_id,
    study_name = meta$study_name,
    disease = meta$disease_audit,
    analysis_group = if ("analysis_group" %in% names(meta)) as.character(meta$analysis_group) else NA_character_,
    variable_id = field$variable_id,
    metadata_column = field$metadata_column,
    domain = field$domain,
    main_model_status = field$main_model_status,
    value_source_object = field$value_source_object,
    record_status = status,
    raw_value = value,
    stringsAsFactors = FALSE
  )
})
sample_level_long <- do.call(rbind, long_list)
write_csv_utf8(sample_level_long, "04_SampleLevel_FieldStatus_Long.csv")

# ---- 6. Cohort-level and cohort-by-disease summaries -------------------------
summarize_availability <- function(data, group_columns) {
  grouping_columns <- c(
    group_columns,
    "variable_id"
  )
  split_data <- split_audit(data, grouping_columns)
  rows <- lapply(split_data, function(z) {
    base <- z[1, c(group_columns, "variable_id", "metadata_column", "domain", "main_model_status"), drop = FALSE]
    n_total <- nrow(z)
    n_documented <- sum(z$record_status == "documented")
    n_unknown <- sum(z$record_status == "unknown_or_unspecified")
    n_missing <- sum(z$record_status == "missing_or_blank")
    n_unmapped <- sum(z$record_status == "no_predefined_exact_mapping")
    cbind(
      base,
      n_total = n_total,
      n_documented = n_documented,
      pct_documented = percent(n_documented, n_total),
      n_unknown_or_unspecified = n_unknown,
      pct_unknown_or_unspecified = percent(n_unknown, n_total),
      n_missing_or_blank = n_missing,
      pct_missing_or_blank = percent(n_missing, n_total),
      value_source_object = z$value_source_object[1],
      n_no_predefined_mapping = n_unmapped,
      automated_availability = if (n_unmapped == n_total) {
        "no_predefined_exact_mapping_check_alias_search"
      } else if (n_documented == 0L) {
        "no_documented_values_in_this_stratum"
      } else {
        "candidate_data_present_requires_source_verification"
      },
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

by_cohort_disease <- summarize_availability(
  sample_level_long,
  c("study_name", "disease")
)
by_cohort <- summarize_availability(sample_level_long, c("study_name"))
assert_audit("cohort_summary_rows", nrow(by_cohort) == nrow(field_dictionary) * 5L,
  "All 30 concepts must have five cohort rows, including unmapped concepts.")
assert_audit("cohort_disease_summary_rows", nrow(by_cohort_disease) == nrow(field_dictionary) * 10L,
  "All 30 concepts must have ten cohort-by-disease rows.")
assert_audit("sample_long_rows", nrow(sample_level_long) == nrow(meta) * nrow(field_dictionary),
  "Exactly one row per included sample and predefined concept.")
assert_audit("cohort_summary_unique_keys", !anyDuplicated(by_cohort[c("study_name", "variable_id")]),
  "No duplicated cohort-concept rows.")
assert_audit("availability_categories_partition", all(with(by_cohort_disease,
  n_documented + n_unknown_or_unspecified + n_missing_or_blank + n_no_predefined_mapping == n_total)),
  "Record statuses must partition every stratum denominator.")
assert_audit("availability_denominators", all(by_cohort_disease$n_total ==
  expected_counts$n_samples[match(paste(by_cohort_disease$study_name, by_cohort_disease$disease),
    paste(expected_counts$study_name, expected_counts$disease))]),
  "Every variable uses the same frozen cohort-disease denominators.")
write_csv_utf8(by_cohort_disease, "02_FieldAvailability_ByCohortDisease.csv")
write_csv_utf8(by_cohort, "03_FieldAvailability_ByCohort.csv")

# ---- 7. Main-model covariate completeness ------------------------------------
model_covariates <- data.frame(
  variable_id = c("age_raw", "age_main", "age_interaction", "sex_source", "sex_main",
    "bmi_raw", "bmi_main", "bmi_interaction", "study_name"),
  metadata_column = c("age", "age", "age_numeric", "gender", "sex", "BMI", "BMI", "BMI_numeric", "study_name"),
  requested_object = c("raw", "model", "model", "raw", "model", "raw", "model", "model", "model"),
  expected_role = c("pre-cleaning saved source", "main-script age covariate", "interaction age covariate",
    "pre-cleaning saved source", "main-script sex covariate", "pre-cleaning saved source",
    "main-script BMI variable", "interaction BMI variable", "fixed effect"),
  stringsAsFactors = FALSE
)
model_covariate_rows <- lapply(seq_len(nrow(model_covariates)), function(i) {
  spec <- model_covariates[i, , drop = FALSE]
  source_data <- if (spec$requested_object == "raw") raw_meta else meta
  source_available <- !is.null(source_data) && spec$metadata_column %in% names(source_data)
  value <- if (source_available) safe_col(source_data, spec$metadata_column) else rep(NA_character_, nrow(meta))
  status <- if (source_available) field_status(value) else rep("source_not_available", nrow(meta))
  aggregate_data <- data.frame(
    study_name = meta$study_name,
    disease = meta$disease_audit,
    record_status = status,
    stringsAsFactors = FALSE
  )
  split_list <- split_audit(aggregate_data, c("study_name", "disease"))
  do.call(rbind, lapply(split_list, function(z) {
    data.frame(
      variable_id = spec$variable_id,
      metadata_column = spec$metadata_column,
      expected_role = spec$expected_role,
      source_object = if (spec$requested_object == "raw") raw_source_object else "tax_meta",
      source_available = source_available,
      study_name = z$study_name[1],
      disease = z$disease[1],
      n_total = nrow(z),
      n_documented = sum(z$record_status == "documented"),
      n_missing_or_unknown = sum(z$record_status != "documented"),
      pct_documented = percent(sum(z$record_status == "documented"), nrow(z)),
      stringsAsFactors = FALSE
    )
  }))
})
model_covariate_completeness <- do.call(rbind, model_covariate_rows)
rownames(model_covariate_completeness) <- NULL
write_csv_utf8(model_covariate_completeness, "05_ModelCovariate_Completeness.csv")

# ---- 8. Medication-specific exports ------------------------------------------
medication_fields <- c("antibiotics_current_use", "treatment", "antidiabetes_group")
medication_columns <- intersect(medication_fields, all_columns)
medication_sample <- meta[, c("sample_id", "study_name", "disease_audit"), drop = FALSE]
for (column in medication_columns) medication_sample[[column]] <- audit_col(column)
names(medication_sample)[names(medication_sample) == "disease_audit"] <- "disease"
medication_sample$treatment_record_status <- if ("treatment" %in% names(medication_sample)) {
  field_status(medication_sample$treatment)
} else {
  "field_absent_from_saved_metadata"
}
medication_sample$interpretation_guardrail <- "Raw text retained; missing or unknown must not be interpreted as documented non-use."
write_csv_utf8(medication_sample, "06_MedicationText_SampleLevel.csv")

if ("treatment" %in% names(medication_sample)) {
  medication_group_split <- split(
    medication_sample,
    interaction(medication_sample$study_name, medication_sample$disease, drop = TRUE)
  )
  medication_group_summary <- do.call(rbind, lapply(medication_group_split, function(z) {
    status <- z$treatment_record_status
    data.frame(
      study_name = z$study_name[1],
      disease = z$disease[1],
      n_total = nrow(z),
      n_documented_text = sum(status == "documented"),
      n_unknown_or_unspecified_text = sum(status == "unknown_or_unspecified"),
      n_missing_or_blank_text = sum(status == "missing_or_blank"),
      pct_documented_text = percent(sum(status == "documented"), nrow(z)),
      conclusion = "Text availability only; not a validated medication exposure or non-exposure classification.",
      stringsAsFactors = FALSE
    )
  }))
} else {
  medication_group_summary <- data.frame(
    study_name = expected_studies,
    disease = NA_character_,
    n_total = NA_integer_,
    n_documented_text = 0L,
    n_unknown_or_unspecified_text = 0L,
    n_missing_or_blank_text = 0L,
    pct_documented_text = NA_real_,
    conclusion = "treatment field absent from saved metadata.",
    stringsAsFactors = FALSE
  )
}
rownames(medication_group_summary) <- NULL
write_csv_utf8(medication_group_summary, "07_MedicationText_ByCohortDisease.csv")

# ---- 9a. Every saved column, including fields outside the initial dictionary --
normalize_search_text <- function(x) {
  value <- gsub("([a-z])([A-Z])", "\\1 \\2", as.character(x), perl = TRUE)
  value <- tolower(value)
  value[is.na(value)] <- ""
  trimws(gsub("[[:space:]]+", " ", gsub("[^[:alnum:].%]+", " ", value)))
}

field_inventory <- do.call(rbind, lapply(all_columns, function(column) {
  value <- audit_col(column)
  present <- field_status(value) == "documented"
  both_sources <- !is.null(raw_meta) && column %in% names(raw_meta) && column %in% input_columns
  raw_value <- if (both_sources) safe_col(raw_meta, column) else rep(NA_character_, nrow(meta))
  model_value <- if (both_sources) safe_col(meta, column) else rep(NA_character_, nrow(meta))
  same_value <- (is.na(raw_value) & is.na(model_value)) |
    (!is.na(raw_value) & !is.na(model_value) & raw_value == model_value)
  same_value[is.na(same_value)] <- FALSE
  source_data <- if (column_source(column) == "meta_df_raw") raw_meta else meta
  data.frame(
    metadata_column = column,
    storage_class = paste(class(source_data[[column]]), collapse = ";"),
    value_source_object = column_source(column),
    present_in_tax_meta = column %in% input_columns,
    present_in_meta_df_raw = !is.null(raw_meta) && column %in% names(raw_meta),
    database_name_evidence = if (!length(database_columns)) {
      "saved_database_schema_unavailable"
    } else if (column %in% database_columns) {
      "name_in_saved_combined_metadata"
    } else {
      "name_not_in_saved_combined_metadata_may_be_derived"
    },
    predefined_concepts = paste(field_dictionary$variable_id[
      !is.na(field_dictionary$metadata_column) & field_dictionary$metadata_column == column], collapse = ";"),
    n_total = nrow(meta), n_documented = sum(present),
    n_distinct_documented_values = length(unique(value[present])),
    n_raw_vs_model_differences = if (both_sources) sum(!same_value) else NA_integer_,
    n_raw_unrecorded_but_model_recorded = if (both_sources) sum(
      field_status(raw_value) != "documented" & field_status(model_value) == "documented") else NA_integer_,
    stringsAsFactors = FALSE
  )
}))
rownames(field_inventory) <- NULL
write_csv_utf8(field_inventory, "12_AllMetadataColumns_Inventory.csv")
write_csv_utf8(field_inventory[, c("metadata_column", "value_source_object", "present_in_tax_meta",
  "present_in_meta_df_raw", "database_name_evidence", "n_raw_vs_model_differences",
  "n_raw_unrecorded_but_model_recorded")], "20_FieldProvenanceAndSourceAlignment.csv")

all_field_availability <- function(group_columns) {
  groups <- split_audit(data.frame(study_name = as.character(meta$study_name),
    disease = meta$disease_audit, row_index = seq_len(nrow(meta))), group_columns)
  do.call(rbind, lapply(all_columns, function(column) {
    value <- audit_col(column)
    status <- field_status(value)
    do.call(rbind, lapply(groups, function(group) {
      indices <- group$row_index
      documented <- status[indices] == "documented"
      data.frame(group[1, group_columns, drop = FALSE], metadata_column = column,
        value_source_object = column_source(column), n_total = length(indices),
        n_documented = sum(documented), pct_documented = percent(sum(documented), length(indices)),
        n_unknown_or_unspecified = sum(status[indices] == "unknown_or_unspecified"),
        n_missing_or_blank = sum(status[indices] == "missing_or_blank"),
        n_distinct_documented_values = length(unique(value[indices][documented])),
        stringsAsFactors = FALSE)
    }))
  }))
}
all_fields_by_cohort <- all_field_availability("study_name")
all_fields_by_disease <- all_field_availability(c("study_name", "disease"))
rownames(all_fields_by_cohort) <- rownames(all_fields_by_disease) <- NULL
assert_audit("full_column_inventory_rows", nrow(field_inventory) == length(all_columns),
  "Every column in the aligned saved metadata has an inventory record.")
assert_audit("full_column_cohort_rows", nrow(all_fields_by_cohort) == length(all_columns) * 5L,
  "All saved columns have all five cohort denominators.")
assert_audit("full_column_disease_rows", nrow(all_fields_by_disease) == length(all_columns) * 10L,
  "All saved columns have all ten cohort-disease denominators.")
write_csv_utf8(all_fields_by_cohort, "13_AllFields_ByCohort.csv")
write_csv_utf8(all_fields_by_disease, "14_AllFields_ByCohortDisease.csv")

# ---- 9b. Broad concept/alias dictionary, applied to ALL column names ----------
# Candidate patterns deliberately include related fields. A match establishes a
# review lead, not equivalence of measurements. Original exact mappings remain.
patterns <- utils::read.csv(file.path(run$root, "config/metadata_concept_patterns.csv"), check.names = FALSE)
concept_patterns <- setNames(patterns$pattern, patterns$concept_id)

concept_dictionary <- data.frame(
  concept_id = names(concept_patterns), pattern = unname(concept_patterns),
  in_predefined_dictionary = names(concept_patterns) %in% field_dictionary$variable_id,
  evidence_rule = "Candidate retrieval only; verify phenotype, definition, units, timing, and individual-level availability.",
  stringsAsFactors = FALSE
)
alias_regression_examples <- c(hba1c = "HbA1c", fasting_glucose = "fastingGlucose",
  alcohol = "alcohol_numeric", insulin_resistance = "HOMA_IR",
  diabetes_duration = "age_T1D_diagnosis")
assert_audit("alias_separator_and_case_examples", all(vapply(names(alias_regression_examples),
  function(concept) grepl(concept_patterns[[concept]],
    normalize_search_text(alias_regression_examples[[concept]]), perl = TRUE), logical(1))),
  "Retrieve underscore/camel-case aliases while preserving disease acronyms such as T1D as review candidates.")
assert_audit("concept_search_covers_dictionary", all(field_dictionary$variable_id %in% names(concept_patterns)),
  "Every predefined concept is covered by the expanded search.")
write_csv_utf8(concept_dictionary, "15_ConceptSearch_Dictionary.csv")

candidate_template <- data.frame(concept_id = character(), metadata_column = character(),
  match_type = character(), value_source_object = character(), n_documented = integer(),
  n_cohorts_with_documented_values = integer(), context_note = character())
candidate_rows <- list()
for (concept in names(concept_patterns)) {
  exact_column <- field_dictionary$metadata_column[match(concept, field_dictionary$variable_id)]
  candidates <- all_columns[grepl(concept_patterns[[concept]], normalize_search_text(all_columns), perl = TRUE)]
  if (length(exact_column) && !is.na(exact_column)) candidates <- union(exact_column, candidates)
  for (column in candidates) {
    inventory_row <- match(column, field_inventory$metadata_column)
    context_note <- if (grepl("t1d|infant|gestational|birth|seroconversion|breastfeed", column, ignore.case = TRUE)) {
      "Potential other-disease or early-life field; cannot establish a T2D measurement."
    } else if (grepl("tnm|ajcc|recist|performance_status|pfs|tumor", column, ignore.case = TRUE)) {
      "Potential oncology or generic severity field; not a validated T2D severity measure."
    } else if (column %in% c("insulin_cat", "treatment", "history_of_therapy") && concept == "insulin_resistance") {
      "Medication/categorical insulin information is not an insulin-resistance measurement."
    } else {
      "Name match only; inspect original definitions and observed records."
    }
    candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
      concept_id = concept, metadata_column = column,
      match_type = if (length(exact_column) && !is.na(exact_column) && column == exact_column) {
        "predefined_exact_column"
      } else "alias_or_related_column_candidate_unverified",
      value_source_object = column_source(column),
      n_documented = field_inventory$n_documented[inventory_row],
      n_cohorts_with_documented_values = sum(all_fields_by_cohort$n_documented[
        all_fields_by_cohort$metadata_column == column] > 0L),
      context_note = context_note, stringsAsFactors = FALSE
    )
  }
}
column_candidates <- bind_rows_or_empty(candidate_rows, candidate_template)
write_csv_utf8(column_candidates, "16_ColumnName_Candidates.csv")

# ---- 9c. Preserve uncurated payloads and retrieve free-text evidence -----------
uncurated_value <- audit_col("uncurated_metadata")
uncurated_export <- data.frame(sample_id = meta$sample_id, study_name = meta$study_name,
  disease = meta$disease_audit, source_column = "uncurated_metadata",
  value_source_object = column_source("uncurated_metadata"),
  record_status = if ("uncurated_metadata" %in% all_columns) field_status(uncurated_value) else
    rep("column_absent_from_saved_metadata", nrow(meta)),
  raw_payload = uncurated_value, stringsAsFactors = FALSE)
write_csv_utf8(uncurated_export, "17_UncuratedMetadata_SampleLevel.csv")

free_text_columns <- intersect(c("uncurated_metadata", "treatment", "history_of_therapy",
  "study_condition", "disease_subtype", "disease_stage", "diet", "lifestyle",
  "population", "smoker", "alcohol"), all_columns)
evidence_template <- data.frame(sample_id = character(), study_name = character(),
  disease = character(), concept_id = character(), source_column = character(),
  value_source_object = character(), matched_term_normalized = character(),
  context_normalized = character(), interpretation = character())
evidence_rows <- list()
for (column in free_text_columns) {
  value <- audit_col(column)
  text <- normalize_search_text(value)
  eligible <- field_status(value) == "documented"
  for (concept in names(concept_patterns)) {
    matches <- regexpr(concept_patterns[[concept]], text, perl = TRUE)
    hits <- which(eligible & matches > 0L)
    if (!length(hits)) next
    match_length <- attr(matches, "match.length")[hits]
    start <- as.integer(matches[hits])
    evidence_rows[[length(evidence_rows) + 1L]] <- data.frame(
      sample_id = meta$sample_id[hits], study_name = meta$study_name[hits],
      disease = meta$disease_audit[hits], concept_id = concept, source_column = column,
      value_source_object = column_source(column),
      matched_term_normalized = substring(text[hits], start, start + match_length - 1L),
      context_normalized = substring(text[hits], pmax(1L, start - 100L), start + match_length + 120L),
      interpretation = "Keyword evidence only; may describe exposure, absence, eligibility, or unrelated context. Not a validated patient-level measurement.",
      stringsAsFactors = FALSE
    )
  }
}
text_evidence <- bind_rows_or_empty(evidence_rows, evidence_template)
write_csv_utf8(text_evidence, "18_FreeText_KeywordEvidence.csv")

# ---- 9d. Cohort search summary, preserving candidate-versus-verified status ---
search_summary <- do.call(rbind, lapply(names(concept_patterns), function(concept) {
  exact_column <- field_dictionary$metadata_column[match(concept, field_dictionary$variable_id)]
  has_exact <- length(exact_column) && !is.na(exact_column)
  alias_columns <- unique(column_candidates$metadata_column[
    column_candidates$concept_id == concept & column_candidates$match_type != "predefined_exact_column"])
  do.call(rbind, lapply(expected_studies, function(study) {
    sample_indices <- which(meta$study_name == study)
    exact_n <- if (has_exact) sum(field_status(audit_col(exact_column)[sample_indices]) == "documented") else NA_integer_
    alias_has_value <- rep(FALSE, length(sample_indices))
    for (column in alias_columns) alias_has_value <- alias_has_value |
      field_status(audit_col(column)[sample_indices]) == "documented"
    evidence <- text_evidence[text_evidence$concept_id == concept & text_evidence$study_name == study, , drop = FALSE]
    uncurated_hits <- length(unique(evidence$sample_id[evidence$source_column == "uncurated_metadata"]))
    text_hits <- length(unique(evidence$sample_id))
    category <- if (has_exact && exact_n > 0L) {
      "exact_field_has_records_definition_still_requires_verification"
    } else if (any(alias_has_value)) {
      "alias_or_related_field_has_records_manual_review_required"
    } else if (text_hits > 0L) {
      "free_text_keyword_only_manual_review_required"
    } else if (has_exact) {
      "exact_field_no_documented_values_no_additional_record_evidence"
    } else {
      "no_predefined_mapping_no_record_evidence_in_saved_search_scope"
    }
    data.frame(concept_id = concept, study_name = study, n_total = length(sample_indices),
      predefined_column = if (has_exact) exact_column else NA_character_,
      n_exact_documented = exact_n, alias_columns_to_review = paste(alias_columns, collapse = ";"),
      n_samples_with_alias_or_related_field_record = sum(alias_has_value),
      n_samples_with_any_free_text_keyword = text_hits,
      n_samples_with_uncurated_keyword = uncurated_hits, search_conclusion = category,
      search_scope = "Aligned saved tax_meta and meta_df_raw; not the full current database or original publications.",
      stringsAsFactors = FALSE)
  }))
}))
rownames(search_summary) <- NULL
assert_audit("concept_search_summary_rows", nrow(search_summary) == length(concept_patterns) * 5L,
  "All 36 concepts have five cohort-level search results.")
assert_audit("keyword_evidence_uses_included_samples", all(text_evidence$sample_id %in% meta$sample_id),
  "All text evidence is restricted to the frozen included sample set.")
write_csv_utf8(search_summary, "19_ConceptSearch_ByCohort.csv")

# ---- 9. Source-verification template -----------------------------------------
# This table intentionally has empty evidence fields. Completion requires the
# original cohort papers, supplements, or data dictionaries; metadata presence
# alone cannot establish variable comparability.
source_template <- do.call(rbind, lapply(expected_studies, function(study) {
  data.frame(
    study_name = study,
    variable_id = field_dictionary$variable_id,
    metadata_column = field_dictionary$metadata_column,
    domain = field_dictionary$domain,
    automated_availability = by_cohort$automated_availability[
      match(paste(study, field_dictionary$variable_id), paste(by_cohort$study_name, by_cohort$variable_id))
    ],
    alias_columns_to_review = search_summary$alias_columns_to_review[
      match(paste(field_dictionary$variable_id, study), paste(search_summary$concept_id, search_summary$study_name))],
    expanded_search_conclusion = search_summary$search_conclusion[
      match(paste(field_dictionary$variable_id, study), paste(search_summary$concept_id, search_summary$study_name))],
    source_article_or_data_dictionary = "",
    PMID_or_DOI = "",
    original_variable_definition = "",
    measurement_unit_or_coding = "",
    measurement_timing = "",
    population_or_disease_group_measured = "",
    inclusion_exclusion_or_comorbidity_detail = "",
    comparable_across_five_cohorts = "undetermined_pending_source_verification",
    proposed_use = "do_not_adjust_until_verified",
    evidence_grade = "unverified",
    reviewer_response_wording = "",
    stringsAsFactors = FALSE
  )
}))
write_csv_utf8(source_template, "08_SourceVerification_Template.csv")

# ---- 10. Draft manuscript-facing summary -------------------------------------
summary_rows <- lapply(seq_len(nrow(field_dictionary)), function(i) {
  field <- field_dictionary[i, , drop = FALSE]
  z <- by_cohort[by_cohort$variable_id == field$variable_id, , drop = FALSE]
  n_candidate_cohorts <- sum(z$n_documented > 0L)
  search_rows <- search_summary[search_summary$concept_id == field$variable_id, , drop = FALSE]
  data.frame(
    variable_id = field$variable_id,
    domain = field$domain,
    metadata_column = field$metadata_column,
    main_model_status = field$main_model_status,
    n_cohorts_with_any_documented_value = if (is.na(field$metadata_column)) NA_integer_ else n_candidate_cohorts,
    n_cohorts_total = length(expected_studies),
    n_cohorts_with_alias_record_candidates = sum(search_rows$n_samples_with_alias_or_related_field_record > 0L),
    n_cohorts_with_free_text_keyword_evidence = sum(search_rows$n_samples_with_any_free_text_keyword > 0L),
    automated_conclusion = if (is.na(field$metadata_column)) {
      "No predefined exact mapping; review expanded alias and free-text search before making any absence claim."
    } else if (n_candidate_cohorts == 0L) {
      "A metadata column exists but contains no documented values in the included five-cohort sample."
    } else if (n_candidate_cohorts < length(expected_studies)) {
      "Available in only a subset of cohorts or strata; do not claim uniform adjustment."
    } else {
      "Values occur in all cohorts, but definition, timing, units, and comparability remain unverified."
    },
    manuscript_safe_wording = if (field$main_model_status == "exposure_or_effect_modifier") {
      "BMI is included as a covariate, stratifier, or effect modifier according to the model; see the model-specific Methods."
    } else if (field$main_model_status == "species_T2D_internal_primary") {
      "Treatment-text medication grouping is included in T2D-internal species primary models only. It does not establish comprehensive medication ascertainment or metformin-naive status."
    } else if (field$main_model_status == "retired_not_valid") {
      "Historical derived grouping is retained for traceability only and is not a reliable multimedication adjustment."
    } else if (field$main_model_status %in% c("adjusted", "fixed_effect")) {
      "Included in the primary model as specified."
    } else {
      "Not included as a uniform pooled covariate; verify exact and candidate field evidence before reporting availability or residual confounding."
    },
    stringsAsFactors = FALSE
  )
})
manuscript_summary <- do.call(rbind, summary_rows)
write_csv_utf8(manuscript_summary, "09_ManuscriptAvailability_Summary_Draft.csv")

# ---- 11. Input audit and reproducibility manifest ----------------------------
study_audit <- as.data.frame(table(meta$study_name, meta$disease_audit), stringsAsFactors = FALSE)
names(study_audit) <- c("study_name", "disease", "n_samples")
write_csv_utf8(study_audit, "10_InclusionAndStudyAudit.csv")

assert_audit("unmapped_concepts_retained", all(vapply(
  field_dictionary$variable_id[is.na(field_dictionary$metadata_column)],
  function(concept) sum(by_cohort$variable_id == concept) == 5L, logical(1))),
  "Duration, SES, and comorbidity concepts remain visible in all five cohorts.")
assert_audit("source_template_rows", nrow(source_template) == nrow(field_dictionary) * 5L &&
  !anyNA(source_template$automated_availability), "Source template covers every cohort-concept combination.")
assert_audit("all_saved_input_columns_scanned", all(input_columns %in% field_inventory$metadata_column),
  paste("All", input_n_columns, "tax_meta columns inventoried before adding audit fields."))

  write_csv_utf8(do.call(rbind, validation_rows), "21_ValidationChecks.csv")
  safe_write_lines(c("All input objects were rebuilt from frozen raw resources in this run.",
    "Historical provenance labels tax_meta/meta_df_raw identify model/raw roles only.",
    "Raw age is audited before imputation; model completeness is reported separately.",
    "Alias hits and recorded text are candidates requiring source verification.",
    "Missing/unknown medication text does not establish documented non-use.",
    "No sample exclusion, model refitting, alias merging or new covariate adjustment."),
    path_in_run(run, "tables/11_metadata_medication_audit/availability/22_README.txt"))
  invisible(list(by_cohort = by_cohort, by_cohort_disease = by_cohort_disease))
}

