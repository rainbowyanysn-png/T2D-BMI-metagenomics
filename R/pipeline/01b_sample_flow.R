# Recomputes supplementary evidence from the passed current-run original inputs.
run_sample_flow <- function(data, run, config) {
  root <- run$root
  read <- function(x) read_current_result(run, file.path("tables", x))
  saved <- read("01_sample_manifest.csv")
  expected <- read_expected(root)$inputs
  tolerance <- config$project$numeric_tolerance
  compare <- function(x, y, label) {
    assert_true(isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE)), label)
  }
  calculate <- function() {
    require_project_packages(c("curatedMetagenomicData", "SummarizedExperiment", "dplyr", "tibble"))
    check_frozen_versions(config)
    retrieve <- prepare_frozen_resources(run, config)
    md <- curatedMetagenomicData::sampleMetadata
    stool <- md[md$disease %in% unlist(config$source_data$disease_values) &
                  md$body_site %in% config$source_data$body_site, , drop = FALSE]
    candidate <- stool[!is.na(stool$BMI), , drop = FALSE]
    candidate$disease <- ifelse(candidate$disease %in% c("healthy", "control"), "control", "T2D")
    candidate$analysis_group <- factor(dplyr::case_when(
      candidate$disease == "T2D" & candidate$BMI < 25 ~ "Non_Overweight_T2D",
      candidate$disease == "T2D" & candidate$BMI >= 25 ~ "Overweight_T2D",
      candidate$disease == "control" & candidate$BMI < 25 ~ "Non_Overweight_Control",
      candidate$disease == "control" & candidate$BMI >= 25 ~ "Overweight_Control"
    ), levels = analysis_group_levels)
    candidate <- candidate[!is.na(candidate$analysis_group), , drop = FALSE]
    se <- retrieve(candidate, dataType = "relative_abundance", rownames = "short")
    returned <- tibble::rownames_to_column(as.data.frame(SummarizedExperiment::colData(se)), "sample_id")
    counts <- dplyr::count(returned, study_name, analysis_group, name = "n_samples")
    qualified <- dplyr::filter(counts, n_samples >= config$cohorts$min_samples_per_group)
    qualified$group_type <- ifelse(grepl("Control", qualified$analysis_group), "Control", "T2D")
    keep <- qualified |>
      dplyr::group_by(study_name) |>
      dplyr::filter(dplyr::n_distinct(group_type) >= 2L) |>
      dplyr::pull(study_name) |>
      as.character() |>
      unique()
    final <- returned[returned$study_name %in% keep, , drop = FALSE]
    stopifnot(setequal(keep, unlist(config$cohorts$names)),
              identical(as.character(final$sample_id), as.character(saved$sample_id)),
              identical(as.character(final$analysis_group), as.character(saved$BMI25_group)),
              nrow(final) == expected$n_samples)
    counts$cell_meets_minimum <- counts$n_samples >= config$cohorts$min_samples_per_group
    counts$cohort_retained <- counts$study_name %in% keep
    steps <- data.frame(
      step = c("Eligible stool metadata", "BMI recorded and BMI25 group defined",
               "Relative-abundance sample returned", "Five cohorts retained, all their samples",
               "BMI25 classification", "Cohort-adapted BMI23 classification"),
      n_samples = c(nrow(stool), nrow(candidate), nrow(returned), nrow(final), nrow(final), nrow(final))
    )
    list(BMI25_BMI23_pre_QC_flow = steps, BMI25_pre_QC_cells = counts)
  }

  result <- calculate()
  for (file in names(result)) safe_write_csv(result[[file]],
    path_in_run(run, "tables/01b_sample_flow", paste0(file, ".csv")))
  invisible(result)
}
