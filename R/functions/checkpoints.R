# Immutable, project-generated checkpoints. Legacy statistical results are never inputs.
invocation_stamp <- function() paste0(format(Sys.time(), "%Y%m%d_%H%M%OS3"), "_", Sys.getpid())

relative_project_path <- function(path) {
  assert_true(is.character(path) && length(path) == 1L && nzchar(path), "Missing checkpoint path")
  assert_true(!grepl("^([A-Za-z]:|/|\\\\)|(^|[/\\\\])[.][.]([/\\\\]|$)", path), "Checkpoint path must be relative without parent traversal")
  path
}

latest_pipeline_state <- function(run_dir) {
  paths <- list.files(file.path(run_dir, "checkpoints"), pattern = "^State_[0-9]+[.]yml$", full.names = TRUE)
  if (!length(paths)) return(NULL)
  yaml::read_yaml(tail(sort(paths), 1L))
}

append_pipeline_state <- function(run, state) {
  folder <- file.path(run$dir, "checkpoints"); dir.create(folder, showWarnings = FALSE)
  previous <- list.files(folder, pattern = "^State_[0-9]+[.]yml$")
  number <- if (length(previous)) max(as.integer(sub("^State_([0-9]+)[.]yml$", "\\1", previous))) + 1L else 1L
  target <- file.path(folder, sprintf("State_%06d.yml", number))
  assert_true(!file.exists(target), "Refusing to overwrite a state snapshot")
  yaml::write_yaml(state, target)
  invisible(state)
}

module_function_dependencies <- function(module, main) {
  extras <- list(
    species_loco = c("species_loco_specs","species_loco_metadata","species_loco_primary_gate","species_loco_primary_significant","species_loco_summarise_all","species_loco_feature_summary","species_loco_acceptance","medication_model_covariate","medication_model_reference"),
    publication_display = c("prepare_publication_display_data", "publication_phylum_display",
      "prepare_frozen_resources", "cached_resource_path", "hub_cache_info", "read_current_result", "current_result_path"),
    pooled_pathways = c("run_pooled_pathway_maaslin", "audited_maaslin2", "categorical_references"),
    clinical_correlations = c("clinical_correlation_data", "clinical_feature_display", "current_result_path", "read_current_result"),
    raw_beta = "historical_beta_clr",
    manuscript_summary = "manuscript_summary_tables",
    sample_flow = c("prepare_frozen_resources", "cached_resource_path", "hub_cache_info"),
    alpha_details = "diversity_data",
    medication_sensitivity = c("read_current_result", "current_result_path"),
    primary = c("medication_model_covariate", "medication_model_reference"),
    clr = c("clr_prepared_metadata", "clr_pair_model_spec", "medication_model_covariate", "medication_model_reference"),
    abundant_taxa = c("abundant_taxa_prepare", "abundant_taxa_fit", "prepare_frozen_resources", "cached_resource_path", "hub_cache_info"),
    inputs = c("import_harmonized_data", "prepare_frozen_resources", "cached_resource_path", "hub_cache_info"),
    interactions = c("run_one_formal_interaction", "interaction_metadata"),
    primary_diversity = c("diversity_data", "historical_beta_clr", "replay_primary_jitter_rng"),
    diversity_interactions = c("diversity_data", "historical_beta_clr"),
    bmi23_species = c("medication_model_covariate", "medication_model_reference"),
    bmi23_diversity = c("diversity_data", "historical_beta_clr"),
    metadata_audit = c("build_metadata_availability_audit", "build_medication_sample_audit", "medication_model_covariate"),
    figures = c("clinical_export_figure", "raw_aitchison_figures", "abundant_taxa_figures", "publication_theme", "species_plot_label", "save_publication_figure", "read_current_result", "current_result_path", "read_presentation_config", "check_export_package_versions", "saved_p_label", "alpha_reference_panel", "pcoa_reference_panel", "permanova_reference_panel", "species_heatmap_data"),
    tables = c("clinical_export_workbook", "write_publication_workbooks", "workbook_group", "publication_column_headers", "verify_workbook_values", "read_presentation_config", "check_export_package_versions"))
  fit_modules <- c("species_loco","primary", "additive", "interactions", "bmi23_species", "bmi23_interactions", "no_medication")
  unique(c(main, extras[[module]], if (module %in% fit_modules) c("audited_maaslin2", "categorical_references")))
}

module_signature <- function(module, main, scope, run, config, expected, input_hash, upstream) {
  functions <- module_function_dependencies(module, main)
  code <- setNames(lapply(functions, function(name) paste(deparse(get(name, envir = scope, inherits = TRUE), width.cutoff = 500L, control = c("keepInteger", "niceNames")), collapse = "\n")), functions)
  keys <- list(
    species_loco = c("maaslin_normalization","maaslin_transform","primary_maaslin_standardize","primary_seed"),
    primary = c("maaslin_normalization", "maaslin_transform", "primary_maaslin_standardize", "primary_seed"),
    additive = c("maaslin_normalization", "maaslin_transform", "additive_maaslin_standardize", "additive_seed", "additive_bmi_center", "additive_bmi_scale"),
    interactions = c("maaslin_normalization", "maaslin_transform", "interaction_maaslin_standardize", "bmi25_interaction_seed", "continuous_interaction_seed"),
    pathways = c("pathway_cohort_prevalence", "pathway_min_features", "pathway_size_factors", "pathway_dispersion_fit", "pathway_meta_method"),
    primary_diversity = c("primary_seed", "diversity_permutations", "diversity_czm_warning", "diversity_czm_delete"),
    diversity_interactions = c("continuous_interaction_seed", "diversity_permutations", "diversity_czm_warning", "diversity_czm_delete"),
    bmi23_species = "bmi23_maaslin_standardize", bmi23_interactions = "interaction_maaslin_standardize",
    bmi23_pathways = c("bmi23_pathway_levels", "pathway_cohort_prevalence", "pathway_min_features", "pathway_size_factors", "pathway_dispersion_fit", "pathway_meta_method"),
    bmi23_diversity = c("bmi23_pathway_levels", "bmi23_diversity_seed", "diversity_permutations"),
    no_medication = c("primary_seed", "no_medication_maaslin_standardize"))
  parameters <- if (module == "inputs") list(source = config$source_data[c("provider", "curatedMetagenomicData_version", "body_site", "disease_values", "abundance_types")], cohorts = config$cohorts, filters = config$filters, expected = expected$inputs, levels = analysis_group_levels) else list(models = config$models[keys[[module]]], alpha = config$project$alpha)
  if (module == "medication_sensitivity") parameters$upstream <- upstream
  if (module == "clinical_correlations") {parameters$clinical_correlations <- config$clinical_correlations; parameters$upstream <- upstream}
  if (module == "pooled_pathways") parameters$pooled_pathway_maaslin <- config$pooled_pathway_maaslin
  if (module == "raw_beta") parameters$raw_aitchison <- config$raw_aitchison
  if (module == "abundant_taxa") parameters$abundant_taxa <- config$abundant_taxa
  if (module == "publication_display") {
    parameters$source <- config$source_data
    parameters$alignment_tolerance <- config$project$numeric_tolerance
    parameters$upstream <- upstream[intersect(names(upstream), c("inputs", "abundant_taxa"))]
  }
  if (module == "species_loco") {parameters$prevalence <- config$filters$global_prevalence;parameters$upstream <- upstream;parameters$cohorts <- config$cohorts;parameters$acceptance <- expected$species_loco}
  if (module == "primary") parameters$prevalence <- config$filters$global_prevalence
  expectation_keys <- list(primary = c("primary_bmi25_species", "primary_bmi25_tested"), additive = "additive", interactions = "interactions", pathways = "primary_bmi25_pathways")
  parameters$expected <- if (module == "inputs") expected$inputs else expected[expectation_keys[[module]]]
  if (module == "bmi23_pathways") parameters$collation <- config$output$historical_collation
  if (module == "manuscript_summary") parameters$upstream <- upstream
  if (module %in% c("tables", "figures")) {if (module == "figures") parameters$output <- config$output; parameters$display <- group_display; parameters$upstream <- upstream}
  if (module %in% c("tables", "figures")) files_presentation <- c("config/presentation.yml", "config/export_package_versions.csv") else files_presentation <- character()
  files <- c("config/package_versions.csv", if (module == "inputs") c("config/source_resources.csv", "config/package_input_objects.csv"), if (module == "metadata_audit") c("config/metadata_field_dictionary.csv", "config/metadata_concept_patterns.csv", "config/medication_dictionary.csv"))
  if (module %in% c("abundant_taxa", "sample_flow", "publication_display")) {
    files <- c(files, "config/source_resources.csv", "config/package_input_objects.csv")
  }
  files <- c(files, files_presentation)
  file_hashes <- setNames(vapply(file.path(run$root, files), sha256_file, character(1)), files)
  digest::digest(list(code = code, parameters = parameters, files = file_hashes, input = input_hash, R = as.character(getRversion())), algo = "sha256")
}

checkpoint_record <- function(run, relative) {
  file <- file.path(run$dir, relative_project_path(relative))
  assert_file(file, "Module checkpoint")
  yaml::read_yaml(file)
}

verify_module_checkpoint <- function(run, record) {
  assert_true(identical(record$status, "COMPUTATION_PASS"), "Checkpoint is not a completed module")
  index_file <- file.path(run$dir, relative_project_path(record$artifact_index))
  assert_file(index_file)
  assert_equal(sha256_file(index_file), record$artifact_index_sha256, "Checkpoint index SHA-256")
  index <- utils::read.csv(index_file, check.names = FALSE)
  for (i in seq_len(nrow(index))) {
    file <- file.path(run$root, relative_project_path(index$project_file[i]))
    assert_file(file, "Checkpoint artifact")
    assert_equal(sha256_file(file), index$sha256[i], paste("Checkpoint artifact SHA-256", index$project_file[i]))
  }
  invisible(TRUE)
}

read_verified_input_checkpoint <- function(run, record) {
  verify_module_checkpoint(run, record)
  relative <- relative_project_path(record$input_rds)
  assert_true(startsWith(relative, paste0("data_derived/", run$id, "/")) && basename(relative) == "harmonized_inputs.rds", "Only this run's rebuilt input checkpoint may be restored")
  file <- file.path(run$root, relative)
  assert_equal(sha256_file(file), record$input_sha256, "Rebuilt input checkpoint SHA-256")
  readRDS(file)
}

save_module_checkpoint <- function(base_run, module_run, module, signature) {
  files <- list.files(file.path(module_run$dir, "tables"), recursive = TRUE, full.names = TRUE)
  if (module == "figures") files <- c(files, list.files(file.path(module_run$dir, "figures"), recursive = TRUE, full.names = TRUE), list.files(file.path(module_run$dir, "plot_data"), recursive = TRUE, full.names = TRUE))
  rds <- if (module == "inputs") file.path(module_run$derived_dir, "harmonized_inputs.rds") else NULL
  files <- c(files, rds)
  assert_true(length(files) > 0L && all(file.exists(files)), "Completed module must have saved artifacts")
  index <- data.frame(project_file = substring(files, nchar(base_run$root) + 2L), sha256 = vapply(files, sha256_file, character(1)))
  index_path <- file.path(module_run$dir, "checkpoint_artifacts.csv")
  safe_write_csv(index, index_path)
  record <- list(schema_version = 1L, status = "COMPUTATION_PASS", independent_frozen_validation = "NOT_IMPLIED_BY_CHECKPOINT", module = module,
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), signature = signature,
    module_directory = substring(module_run$dir, nchar(base_run$dir) + 2L),
    artifact_index = substring(index_path, nchar(base_run$dir) + 2L), artifact_index_sha256 = sha256_file(index_path))
  if (module == "inputs") {record$input_rds <- substring(rds, nchar(base_run$root) + 2L); record$input_sha256 <- sha256_file(rds)}
  target <- file.path(module_run$dir, "COMPLETED.yml")
  assert_true(!file.exists(target), "Refusing to overwrite completion record")
  yaml::write_yaml(record, target)
  substring(target, nchar(base_run$dir) + 2L)
}

completed_result_directories <- function(run) {
  assert_rc4_run_identity(run)
  state <- latest_pipeline_state(run$dir)
  assert_true(!is.null(state), "RC4 result inputs require a valid checkpoint state; flat historical results are not accepted")
  setNames(vapply(state$completed, function(relative) file.path(run$dir, relative_project_path(checkpoint_record(run, relative)$module_directory)), character(1)), names(state$completed))
}

parse_pipeline_arguments <- function(arguments) {
  out <- list(selected = "all_candidate", resume = NULL, from = NULL, until = NULL)
  i <- 1L
  while (i <= length(arguments)) {
    token <- arguments[i]
    if (token %in% c("--resume", "--from", "--until")) {
      assert_true(i < length(arguments) && !startsWith(arguments[i + 1L], "--"), paste("Missing value for", token))
      out[[sub("^--", "", token)]] <- arguments[i + 1L]; i <- i + 2L
    } else {
      assert_true(i == 1L && !startsWith(token, "--"), paste("Unknown argument", token))
      out$selected <- token; i <- i + 1L
    }
  }
  assert_true(is.null(out$from) || !is.null(out$resume), "--from requires --resume RUN_ID")
  out
}
