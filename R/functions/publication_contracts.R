# Export contracts: these functions inspect saved values and never fit a model.
pub_caption_text <- function(x) {
  x<-pub_publication_text(x)
  x<-gsub("\\bqval\\b","q value",x,perl=TRUE)
  x<-gsub("\\bpval\\b","P value",x,perl=TRUE)
  x
}

pub_plot_language_check <- function(plot,id) {
  labels<-unlist(lapply(plot$labels,function(z)if(is.character(z))z else character()),use.names=FALSE)
  built<-ggplot2::ggplot_build(plot)
  drawn<-unlist(lapply(built$data,function(z)if("label"%in%names(z))as.character(z$label) else character()),use.names=FALSE)
  pub_language_check(c(labels,drawn),id)
  if(inherits(plot,"patchwork"))for(child in plot$patches$plots)pub_plot_language_check(child,id)
  invisible(TRUE)
}
pub_named_mapping <- function(values, observed, context, colours = FALSE) {
  values <- unlist(values, use.names = TRUE)
  keys <- names(values)
  assert_true(length(values) > 0L && length(keys) == length(values) &&
    !anyNA(keys) && all(nzchar(keys)) && !anyDuplicated(keys),
    paste("Mapping needs unique nonempty names:", context))
  levels <- unique(as.character(observed))
  assert_true(!anyNA(levels) && all(nzchar(levels)), paste("Missing category before mapping:", context))
  absent <- setdiff(levels, keys)
  assert_true(!length(absent), paste("Unmapped categories:", context, paste(absent, collapse = ", ")))
  if (colours) {
    valid <- tryCatch({ grDevices::col2rgb(unname(values)); TRUE }, error = function(e) FALSE)
    assert_true(valid, paste("Invalid colour specification:", context))
  }
  values
}

pub_palette <- function(values, observed, context) pub_named_mapping(values, observed, context, colours = TRUE)

pub_warning_is_fatal <- function(message) {
  grepl(paste(c("No shared levels", "Insufficient values in manual scale", "Removed .*rows?",
    "non-finite", "Computation failed", "Failed to fit", "Unknown colour", "Unknown color",
    "Ignoring unknown aesthetics", "Too few points"), collapse = "|"), message, ignore.case = TRUE)
}

pub_probability_contract <- function(x, fields, context, allow_missing = FALSE) {
  assert_columns(x, fields, context)
  for (field in fields) {
    v <- x[[field]]
    assert_true(is.numeric(v), paste("Probability column must be numeric:", context, field))
    valid <- is.finite(v) & v >= 0 & v <= 1
    if (allow_missing) valid <- valid | is.na(v)
    assert_true(all(valid), paste("Invalid probability value:", context, field))
  }
  invisible(TRUE)
}

pub_key_contract <- function(x, fields, context) {
  assert_columns(x, fields, context)
  assert_true(!anyNA(x[, fields, drop = FALSE]) && !anyDuplicated(x[, fields, drop = FALSE]),
    paste("Missing or duplicate key:", context, paste(fields, collapse = ", ")))
}

pub_export_data_contract <- function(b, config, export, run) {
  pub_clr_alignment_contract(run)
  pub_key_contract(b$md, "sample_id", "sample metadata")
  cfg<-config$pooled_pathway_maaslin
  assert_true(isTRUE(cfg$standardize) && identical(unlist(cfg$fixed_effects),c("disease","BMI","age","sex","study_name")) &&
    !length(cfg$random_effects) && cfg$normalization=="TSS" && cfg$transform=="LOG",
    "Pooled pathway configuration must retain the accepted fixed-study model")
  pub_named_mapping(group_display, b$md$BMI25_group, "BMI25 display groups")
  style <- read_presentation_config(config$project_root)
  pub_palette(style$groups$fill, unname(group_display[b$md$BMI25_group]), "four-group palette")
  pub_palette(config$publication$effect_colours, c("Positive", "Negative", "NS"), "effect palette")
  for (domain in c("species", "pooled")) {
    x <- b[[domain]]
    pub_key_contract(x, if (domain == "species") c("comparison", "feature") else "feature", domain)
    pub_probability_contract(x, c("pval", "qval", "formal_q"), domain)
    assert_true(identical(x$qval, x$formal_q), paste("Full-model formal q must equal saved qval:", domain))
    assert_true(all(is.finite(x$coef)) && all(is.finite(x$stderr)), paste("Nonfinite estimates:", domain))
  }
  pub_key_contract(b$reml, c("comparison", "pathway"), "REML")
  pub_probability_contract(b$reml, c("p_value", "padj"), "REML")
  pub_key_contract(b$interactions, c("BMI_definition", "domain", "feature"), "interactions")
  pub_probability_contract(b$interactions, c("pval", "formal_q"), "interactions")
  pub_key_contract(b$additive, "feature", "additive")
  pub_probability_contract(b$additive, c("pval_T2D", "pval_BMI_per5",
    "q_within_exposure_BH_T2D", "q_within_exposure_BH_BMI_per5"), "additive")
  pub_key_contract(b$clinical, c("feature", "clinical_trait"), "clinical correlations")
  pub_probability_contract(b$clinical, c("p_value", "q_bh_60"), "clinical correlations")
  assert_true(setequal(b$clinical$clinical_trait, c("BMI", "age", "sex")), "Unexpected clinical display traits")
  safe_write_csv(data.frame(check = c("sample_keys", "category_mappings", "species_and_pooled_full_q",
    "REML_family_fields", "interaction_keys", "additive_exposure_fields", "clinical_keys_and_probabilities"),
    status = "MATCHED_SAVED_DATA_CONTRACT"), file.path(export, "validation/export_data_contracts.csv"))
  invisible(TRUE)
}

pub_validate_panel <- function(plot, data, id, width, height) {
  assert_true(is.data.frame(data) && nrow(data) > 0L, paste("Empty plot-data:", id))
  assert_true(is.finite(width) && width > 0 && is.finite(height) && height > 0,
    paste("Invalid panel dimensions:", id))
  # Build a plot before opening any output device. Mapping contracts also run
  # inside each subplot constructor, including multi-plot abundant-taxa panels.
  ggplot2::ggplot_build(plot)
  pub_plot_language_check(plot,id)
  invisible(TRUE)
}


# Reject pre-alignment CLR checkpoints before any current publication export.
pub_clr_alignment_contract <- function(run) {
  read<-function(n)pub_read(run,paste0("tables/08_CLR_sensitivity/",n))
  contracts<-read("CLR_model_contracts.csv");md<-read("CLR_sample_metadata.csv")
  original<-pub_read(run,"tables/01_sample_manifest.csv")
  assert_true(nrow(contracts)==16L&&!anyDuplicated(contracts$result_file),"CLR requires 16 aligned model contracts")
  assert_true(nrow(md)==1186L&&!anyDuplicated(md$sample_id)&&identical(md$sample_id,original$sample_id),"CLR aligned sample identities differ")
  assert_true(identical(as.character(md$BMI23_group),as.character(original$BMI23_group))&&identical(md$BMI23_threshold,original$BMI23_threshold),"CLR BMI23 differs from harmonized cohort-adapted groups")
  expected_drug<-as.character(medication_model_covariate(original)$antidiabetes_group)
  assert_true(identical(as.character(md$antidiabetes_group),expected_drug),"CLR medication grouping differs from primary mapping")
  assert_true(all(md$BMI23_binary==as.integer(md$BMI_numeric>=md$BMI23_threshold)),"CLR BMI23 binary is not cohort-adapted")
  assert_true(all(contracts$clr_contract_version=="RC4_CLR_ALIGNED_20261119")&&all(contracts$design_rank==contracts$design_columns),"Old or rank-deficient CLR model contract")
  expected_ids<-c("Total_T2D_vs_control","Lean_T2D_vs_Lean_Control","Overweight_T2D_vs_Overweight_Control","Overweight_vs_Lean_Controls","Overweight_vs_Lean_T2D","BMI25_Formal_Interaction","ContinuousBMI_Formal_Interaction","BMI23_Formal_Interaction")
  assert_true(setequal(paste(contracts$analysis_domain,contracts$analysis_id),as.vector(outer(c("Taxonomy","Pathway"),expected_ids,paste))),"CLR retained model set differs")
  for(i in seq_len(nrow(contracts))) {
    c<-contracts[i,,drop=FALSE];d<-read(c$result_file)
    adjusted<-c$analysis_domain=="Taxonomy"&&c$analysis_id=="Overweight_vs_Lean_T2D"
    assert_true(identical(c$medication_adjusted,adjusted)&&identical(grepl("antidiabetes_group",c$model_formula,fixed=TRUE),adjusted),"CLR treatment adjustment assigned to wrong model")
    if(adjusted)assert_true(c$n_samples==760L&&c$medication_reference==medication_model_reference(medication_model_covariate(original[original$disease=="T2D",,drop=FALSE])),"CLR T2D medication reference/sample count differs")
    if(c$analysis_id=="BMI23_Formal_Interaction")assert_true(startsWith(c$result_file,"BMI23_adapted_")&&grepl("Cohort-adapted",c$BMI_definition,fixed=TRUE),"Uniform BMI23 result cannot enter adapted CLR export")
    assert_columns(d,c("clr_contract_version","medication_adjusted","medication_reference","model_formula","target_term"),"Aligned CLR saved results")
    assert_true(nrow(d)==c$n_features_tested&&all(d$n_samples==c$n_samples)&&all(d$clr_contract_version==c$clr_contract_version)&&all(d$medication_adjusted==adjusted)&&all(d$model_formula==c$model_formula)&&all(d$target_term==c$target_term),"CLR result/contract mismatch")
  }
  invisible(contracts)
}
