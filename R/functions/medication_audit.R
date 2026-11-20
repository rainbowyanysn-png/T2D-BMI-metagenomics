# Central prespecified drug grouping for species models and sensitivity comparisons.
medication_model_covariate <- function(metadata) {
  assert_columns(metadata, c("treatment", "disease"), "Medication model mapping")
  treatment_lower <- tolower(as.character(metadata$treatment))
  metadata$antidiabetes_group <- factor(dplyr::case_when(
    grepl("metformin", treatment_lower) ~ "Metformin",
    grepl("insulin", treatment_lower) ~ "Insulin",
    grepl("dppiv", treatment_lower) ~ "DPP-4_Inhibitors",
    as.character(metadata$disease) == "control" ~ "No_Medication",
    TRUE ~ "Other_Anti-diabetic_Drugs"
  ))
  metadata
}

medication_model_reference <- function(metadata) {
  values <- metadata$antidiabetes_group
  assert_true(length(unique(na.omit(values))) > 1L, "Medication categories missing")
  if ("No_Medication" %in% as.character(values)) "No_Medication" else names(which.max(table(values)))[1]
}

# Compatibility for internal archived callers; the mapping has one definition.
historical_medication_covariate <- medication_model_covariate

# Selected original sample-level transformation; no modelling or exclusion is performed.
build_medication_sample_audit <- function(metadata_final, dictionary) {
  `%>%` <- dplyr::`%>%`
  unknown_treatment_tokens <- dictionary$token[dictionary$token_type == "unknown"]
  explicit_no_antidiabetic_tokens <- dictionary$token[dictionary$token_type == "explicit_no"]
  metadata_final %>% dplyr::transmute(sample_id = as.character(sample_id), study_name = as.character(study_name), 
      disease = as.character(disease), BMI = suppressWarnings(as.numeric(BMI)), age = suppressWarnings(as.numeric(age)), 
      sex = as.character(sex), analysis_group = as.character(analysis_group), treatment_raw = as.character(treatment), 
      treatment_normalized = dplyr::if_else(is.na(treatment_raw), NA_character_, stringr::str_squish(stringr::str_to_lower(treatment_raw))), 
      BMI25_group = dplyr::case_when(BMI < 25 ~ "BMI_lt_25", BMI >= 25 ~ "BMI_ge_25", TRUE ~ NA_character_), 
      BMI3_group = dplyr::case_when(BMI < 25 ~ "BMI_lt_25", BMI >= 25 & BMI < 30 ~ "BMI_25_to_lt_30", BMI >= 
          30 ~ "BMI_ge_30", TRUE ~ NA_character_)) %>% dplyr::mutate(BMI25_group = factor(BMI25_group, 
      levels = c("BMI_lt_25", "BMI_ge_25")), BMI3_group = factor(BMI3_group, levels = c("BMI_lt_25", "BMI_25_to_lt_30", 
      "BMI_ge_30")), medication_text_unknown = is.na(treatment_normalized) | treatment_normalized %in% 
      unknown_treatment_tokens, explicit_no_antidiabetic_medication = !medication_text_unknown & treatment_normalized %in% 
      explicit_no_antidiabetic_tokens, literal_no_treatment_code = !medication_text_unknown & treatment_normalized == 
      "no", medication_record_status = dplyr::case_when(medication_text_unknown ~ "Unknown_or_not_recorded", 
      explicit_no_antidiabetic_medication ~ "Explicit_no_antidiabetic_medication", TRUE ~ "Recorded_text"), 
      contains_metformin = !medication_text_unknown & stringr::str_detect(treatment_normalized, "metformin"), 
      contains_insulin = !medication_text_unknown & stringr::str_detect(treatment_normalized, "insulin|glargine|detemir|degludec|lispro|aspart|glulisine|\\bnph\\b"), 
      contains_dpp4_inhibitor = !medication_text_unknown & stringr::str_detect(treatment_normalized, "dpp[_ -]?(4|iv)|sitagliptin|vildagliptin|saxagliptin|linagliptin|alogliptin"), 
      contains_sglt2_inhibitor = !medication_text_unknown & stringr::str_detect(treatment_normalized, "sglt\\s*[- ]?\\s*2|dapagliflozin|empagliflozin|canagliflozin|ertugliflozin"), 
      contains_glp1_receptor_agonist = !medication_text_unknown & stringr::str_detect(treatment_normalized, 
          "glp[_ -]?1|liraglutide|semaglutide|dulaglutide|exenatide|lixisenatide"), contains_sulfonylurea = !medication_text_unknown & 
          stringr::str_detect(treatment_normalized, "(^|;)su(;|$)|sulfonylurea|glimepiride|gliclazide|glipizide|glyburide|glibenclamide|tolbutamide"), 
      contains_tzd = !medication_text_unknown & stringr::str_detect(treatment_normalized, "thiazolidinedione|pioglitazone|rosiglitazone"), 
      contains_other_common_antidiabetic = !medication_text_unknown & stringr::str_detect(treatment_normalized, 
          "acarbose|voglibose|miglitol|repaglinide|nateglinide"), contains_generic_antidiabetic_label = !medication_text_unknown & 
          stringr::str_detect(treatment_normalized, "(^|;)antidiab(;|$)"), contains_antihypertensive = !medication_text_unknown & 
          stringr::str_detect(treatment_normalized, "antihta|ace_inhib|at2_inhib|beta_blockers|ca2_cbl|thiazidique|kspdiur|ctrblock"), 
      contains_lipid_lowering = !medication_text_unknown & stringr::str_detect(treatment_normalized, "antilipid|statin|statins|fibrate|ezetimibe"), 
      contains_antiplatelet_or_anticoagulant = !medication_text_unknown & stringr::str_detect(treatment_normalized, 
          "(^|;)asa(;|$)|clopidogrel|cumarine"), contains_ppi = !medication_text_unknown & stringr::str_detect(treatment_normalized, 
          "(^|;)ppi(;|$)")) %>% dplyr::mutate(n_recognised_antidiabetic_classes = as.integer(contains_metformin) + 
      as.integer(contains_insulin) + as.integer(contains_dpp4_inhibitor) + as.integer(contains_sglt2_inhibitor) + 
      as.integer(contains_glp1_receptor_agonist) + as.integer(contains_sulfonylurea) + as.integer(contains_tzd) + 
      as.integer(contains_other_common_antidiabetic), medication_provisional_class = dplyr::case_when(medication_text_unknown ~ 
      "Unknown_or_not_recorded", explicit_no_antidiabetic_medication ~ "Explicit_no_antidiabetic_medication", 
      literal_no_treatment_code ~ "Literal_no_treatment_code_manual_source_validation_required", n_recognised_antidiabetic_classes >= 
          2 ~ "Multiple_recognised_antidiabetic_classes", contains_metformin ~ "Metformin_detected_no_other_recognised_class", 
      contains_insulin ~ "Insulin_detected_no_other_recognised_class", contains_dpp4_inhibitor ~ "DPP4_inhibitor_detected_no_other_recognised_class", 
      contains_sglt2_inhibitor ~ "SGLT2_inhibitor_detected_no_other_recognised_class", contains_glp1_receptor_agonist ~ 
          "GLP1_receptor_agonist_detected_no_other_recognised_class", contains_sulfonylurea ~ "Sulfonylurea_detected_no_other_recognised_class", 
      contains_tzd ~ "TZD_detected_no_other_recognised_class", contains_other_common_antidiabetic ~ "Other_recognised_antidiabetic_detected", 
      contains_generic_antidiabetic_label ~ "Generic_antidiabetic_treatment_label_without_specific_class", 
      TRUE ~ "Recorded_text_without_recognised_antidiabetic_keyword"), provisional_metformin_status = dplyr::case_when(medication_text_unknown ~ 
      "Cannot_determine_unknown_or_not_recorded", contains_metformin ~ "Metformin_detected", explicit_no_antidiabetic_medication ~ 
      "Explicit_no_antidiabetic_medication", literal_no_treatment_code ~ "Literal_no_treatment_code_manual_source_validation_required", 
      TRUE ~ "No_metformin_keyword_detected_manual_check_required"), metformin_free_screening_status = dplyr::case_when(medication_text_unknown ~ 
      "Not_eligible_unknown_or_not_recorded", contains_metformin ~ "Not_eligible_metformin_detected", explicit_no_antidiabetic_medication ~ 
      "Candidate_explicit_no_antidiabetic_medication", literal_no_treatment_code ~ "Not_eligible_until_literal_no_code_is_source_validated", 
      TRUE ~ "Candidate_recorded_text_without_metformin_manual_check_required")) %>% dplyr::arrange(study_name, 
      disease, BMI, sample_id)
}

