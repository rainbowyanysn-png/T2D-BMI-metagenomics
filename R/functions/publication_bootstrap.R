# One dependency list for both publication export and its validation entry.
load_publication_support <- function(root, scope) {
  files <- c("assertions", "model_helpers", "medication_audit", "species_loco_helpers", "source_resources", "group_display", "checkpoints",
    "clinical_correlations", "figure_exports", "figure_panels", "publication_theme", "pathway_annotation",
    "publication_contracts", "submission_helpers", "publication_workbook_validation", "publication_table_helpers", "publication_figure_designs",
    "publication_figure_validation", "publication_panel_helpers", "publication_legends", "pathway_forest_selection", "final_integration_tables", "publication_refinement_adapter", "publication_export")
  for (file in files) {
    path <- file.path(root, "R/functions", paste0(file, ".R"))
    if (!file.exists(path)) stop("Missing publication dependency: ", path, call. = FALSE)
    source(path, local = scope)
  }
  for (file in c("12b_publication_tables.R", "13b_publication_panels.R"))
    source(file.path(root, "R/pipeline", file), local = scope)
  required <- c("sha256_file", "verify_module_checkpoint", "publication_data", "pub_palette",
    "pub_warning_is_fatal", "pub_read_csv_snapshot", "pub_validate_panel", "pub_save_panel", "publication_export",
    "pub_export_data_contract", "pub_verify_workbook_sheet", "pub_validate_saved_workbook",
    "pub_clear_empty_drawing_placeholders", "pub_select_volcano_labels", "pub_validate_saved_panel",
    "pub_write_combined_legends", "run_publication_tables", "run_publication_panels")
  missing <- required[!vapply(required, exists, logical(1), envir = scope, inherits = TRUE, mode = "function")]
  if (length(missing)) stop("Publication functions not loaded: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(files)
}
