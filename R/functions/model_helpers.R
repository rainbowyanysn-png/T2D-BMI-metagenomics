required_packages <- c("curatedMetagenomicData", "ExperimentHub", "SummarizedExperiment", "TreeSummarizedExperiment", "dplyr", "tibble", "tidyr", "Maaslin2", "vegan", "permute", "lme4", "lmerTest", "emmeans", "zCompositions", "compositions", "DESeq2", "metafor", "ggplot2", "ggrepel", "patchwork", "yaml", "digest", "sva", "performance", "limma", "stringi", "stringr")

# Transparent audit wrapper: arguments are forwarded unchanged to frozen MaAsLin2.
audited_maaslin2 <- function(..., standardize) {
  arguments <- list(..., standardize = standardize)
  assert_true(is.logical(standardize) && length(standardize) == 1L && !is.na(standardize), "Explicit standardize is mandatory")
  metadata <- as.data.frame(arguments$input_metadata)
  metadata$sample_id <- rownames(metadata)
  fixed <- arguments$fixed_effects
  design_data <- metadata
  for (reference in arguments$reference) {
    pair <- strsplit(reference, ",", fixed = TRUE)[[1]]
    assert_true(length(pair) == 2L && pair[2] %in% as.character(metadata[[pair[1]]]), "Invalid explicit model reference")
    design_data[[pair[1]]] <- stats::relevel(factor(design_data[[pair[1]]]), ref = pair[2])
  }
  assert_true(!anyNA(design_data[, fixed, drop = FALSE]), "Missing fixed effects before fitting")
  design <- stats::model.matrix(stats::reformulate(fixed), data = design_data)
  assert_true(qr(design)$rank == ncol(design), "Rank-deficient model design")
  assert_equal(nrow(design), nrow(metadata), "Model design samples")
  input <- arguments$input_data
  if (identical(colnames(input), metadata$sample_id)) features <- rownames(input) else {
    assert_equal(rownames(input), metadata$sample_id, "Model input sample order")
    features <- colnames(input)
  }
  directory <- file.path(arguments$output, "model_contract")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  safe_write_csv(metadata[, c("sample_id", fixed), drop = FALSE], file.path(directory, "sample_covariates.csv"))
  safe_write_csv(data.frame(feature = features), file.path(directory, "input_feature_ids.csv"))
  safe_write_csv(data.frame(sample_id = metadata$sample_id, design, check.names = FALSE), file.path(directory, "unstandardized_design.csv"))
  contract <- data.frame(parameter = c("fixed_effects", "random_effects", "references", "standardize", "normalization", "transform", "min_prevalence", "Maaslin2_version"), value = c(paste(fixed, collapse = " + "), paste(arguments$random_effects, collapse = ";"), paste(arguments$reference, collapse = ";"), as.character(standardize), arguments$normalization, arguments$transform, as.character(arguments$min_prevalence), as.character(utils::packageVersion("Maaslin2"))))
  safe_write_csv(contract, file.path(directory, "parameters.csv"))
  do.call(Maaslin2::Maaslin2, arguments)
}

check_frozen_versions <- function(config) {
  versions <- c(curatedMetagenomicData = config$source_data$curatedMetagenomicData_version, Maaslin2 = config$models$Maaslin2_version)
  for (package in names(versions)) assert_true(as.character(utils::packageVersion(package)) == versions[[package]], paste("Frozen package version required:", package, versions[[package]]))
}

require_project_packages <- function(packages = required_packages) {
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  assert_true(!length(absent), paste0("Missing required packages; do not install at runtime: ", paste(absent, collapse = ", ")))
}

read_config <- function(root) yaml::read_yaml(file.path(root, "config", "analysis_config.yml"))
read_expected <- function(root) yaml::read_yaml(file.path(root, "config", "expected_results.yml"))

categorical_references <- function(metadata, fixed_effects, largest_only = FALSE) {
  terms <- character()
  for (variable in intersect(c("sex", "study_name"), fixed_effects)) {
    values <- as.character(metadata[[variable]]); values <- values[!is.na(values)]
    reference <- if (!largest_only && variable == "sex" && "female" %in% values) "female" else names(sort(table(values), decreasing = TRUE))[1]
    terms <- c(terms, paste(variable, reference, sep = ","))
  }
  terms
}

