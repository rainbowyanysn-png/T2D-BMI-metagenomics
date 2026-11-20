#!/usr/bin/env Rscript
local({
  entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (!grepl("^([A-Za-z]:[/\\\\]|/|\\\\\\\\)", entry)) entry <- file.path(getwd(), entry)
  root <- dirname(entry)
  source(file.path(root, "R/functions/project_paths.R"), local = environment())
  activate_existing_project_library(root)
  for (name in c("assertions", "checkpoints", "pipeline_runner")) source(file.path(root, "R/functions", paste0(name, ".R")), local = environment())
  args <- commandArgs(TRUE)
  assert_true(length(args) == 1L && grepl("^Run_[0-9]{8}_[0-9]{6}$", args[1]), "Usage: Rscript status.R RUN_ID")
  state <- latest_pipeline_state(file.path(root, "output", args[1]))
  assert_true(!is.null(state), "No checkpoint state found; older pre-checkpoint V3 runs are not resumable")
  modules <- c("inputs", pipeline_registry()$module)
  active <- if (is.null(state$active_module)) "" else state$active_module
  print(data.frame(module = modules, status = ifelse(modules %in% names(state$completed), "COMPLETED_CHECKPOINT", ifelse(modules == active, "FAILED_OR_ACTIVE", "PENDING"))), row.names = FALSE)
  cat("\nRun status:", state$status, "\n")
  if (!is.null(state$last_error)) cat("Last error:", state$last_error, "\n")
  cat("\nCOMPLETED_CHECKPOINT means computation and local assertions passed; it does not imply independent frozen-result validation.\n")
})
