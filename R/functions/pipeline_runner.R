pipeline_registry <- function(root) utils::read.csv(file.path(root,"config/pipeline_modules.csv"),stringsAsFactors=FALSE)

run_project_pipeline <- function(root, selected = "all_candidate", resume = NULL, from = NULL, until = NULL) {
  registry <- pipeline_registry(root)
  sequence <- c("inputs", registry$module)
  manual_only <- character()
  plan <- if (identical(selected, "all_candidate")) setdiff(sequence, manual_only) else if (identical(selected, "publication_candidate"))
    sequence else unique(c("inputs", selected))
  assert_true(all(plan %in% sequence), "Unknown module selection")
  if (!is.null(until)) {assert_true(until %in% plan, "--until is outside this plan"); plan <- plan[seq_len(match(until, plan))]}
  if (!is.null(from)) assert_true(!is.null(resume) && from %in% plan, "--from requires a resumed run and a module in this plan")
  config <- read_config(root); expected <- read_expected(root)
  if (is.null(resume)) {
    run <- run_context(root)
    state <- list(schema_version = 1L, run_id = run$id, completed = list(), status = "CREATED")
    append_pipeline_state(run, state)
  } else {
    assert_true(grepl("^Run_[0-9]{8}_[0-9]{6}$", resume), "Invalid resume Run ID")
    run <- list(root = root, id = resume, dir = file.path(root, "output", resume), derived_dir = file.path(root, "data_derived", resume))
    assert_rc4_run_identity(run)
    state <- latest_pipeline_state(run$dir)
    assert_true(!is.null(state), "This run predates the checkpoint engine or has no state. Start a new run; existing outputs are preserved.")
    assert_equal(state$run_id, run$id, "Checkpoint run identity")
  }
  lock <- file.path(run$dir, ".running.lock")
  assert_true(dir.create(lock, showWarnings = FALSE), "Run is locked. If an earlier R process was killed, first confirm it has stopped, then remove only this run's .running.lock directory.")
  safe_write_lines(c(paste("PID", Sys.getpid()), invocation_stamp()), file.path(lock, "owner.txt"))
  on.exit({if (startsWith(lock, paste0(run$dir, "/")) && basename(lock) == ".running.lock") unlink(lock, recursive = TRUE)}, add = TRUE)
  invocation <- file.path(run$dir, "logs", paste0("Invocation_", invocation_stamp()))
  dir.create(invocation, recursive = TRUE)
  if (!is.null(from)) {
    invalidate <- sequence[seq.int(match(from, sequence), length(sequence))]
    state$completed[intersect(names(state$completed), invalidate)] <- NULL
    state$status <- "RESTART_FROM"; state$restart_from <- from
    append_pipeline_state(run, state)
  }
  for (file in unique(c("01_import_harmonize_qc.R", registry$file[registry$module %in% plan]))) source(file.path(root, "R/pipeline", file), local = environment())
  timings <- list(); warning_rows <- list(); data <- NULL
  step <- function(label, action) {
    started <- Sys.time(); result_status <- "FAIL"; error_text <- ""
    message("START ", label, " ", format(started))
    on.exit({
      finished <- Sys.time()
      timings[[length(timings) + 1L]] <<- data.frame(step = label, started_utc = format(started, tz = "UTC", usetz = TRUE), finished_utc = format(finished, tz = "UTC", usetz = TRUE), elapsed_seconds = as.numeric(difftime(finished, started, units = "secs")), status = result_status, error = error_text)
      message(result_status, " ", label, " ", format(finished))
    })
    value <- tryCatch(withCallingHandlers(action(), warning = function(w) {
      warning_rows[[length(warning_rows) + 1L]] <<- data.frame(step = label, warning = conditionMessage(w))
    }), error = function(e) {error_text <<- conditionMessage(e); stop(e)})
    result_status <- "PASS"; value
  }
  status <- tryCatch({
    step("environment", function() {
      require_project_packages(); check_frozen_versions(config)
      versions <- utils::read.csv(file.path(root, "config/package_versions.csv"))
      installed <- vapply(versions$package, function(p) as.character(utils::packageDescription(p, fields = "Version")), character(1))
      assert_true(identical(unname(installed), versions$version), "Installed DESCRIPTION versions differ from frozen environment")
      safe_write_csv(versions, file.path(invocation, "package_versions.csv"))
    })
    source(file.path(root, "tests/validate_no_legacy_dependency.R"), local = environment())
    safe_write_csv(validate_no_legacy_dependency(root), file.path(invocation, "no_legacy_dependency.csv"))
    code_files <- c(list.files(file.path(root, "R"), pattern = "[.]R$", recursive = TRUE, full.names = TRUE), list.files(file.path(root, "config"), full.names = TRUE))
    code_files <- code_files[!grepl("^reference_", basename(code_files))]
    safe_write_csv(data.frame(file = substring(code_files, nchar(root) + 2L), sha256 = vapply(code_files, sha256_file, character(1))), file.path(invocation, "code_and_config_hashes.csv"))
    for (module in plan) {
      main <- if (module == "inputs") "run_import_harmonize_qc" else registry$function_name[match(module, registry$module)]
      input_record <- if ("inputs" %in% names(state$completed)) checkpoint_record(run, state$completed$inputs) else NULL
      input_hash <- if (module != "inputs" && !is.null(input_record)) input_record$input_sha256 else NULL
      earlier <- sequence[seq_len(match(module, sequence) - 1L)]
      upstream_names <- intersect(earlier, names(state$completed))
      upstream <- setNames(lapply(state$completed[upstream_names], function(path) checkpoint_record(run, path)$artifact_index_sha256), upstream_names)
      signature <- module_signature(module, main, environment(), run, config, expected, input_hash, upstream)
      if (module %in% names(state$completed)) {
        record <- checkpoint_record(run, state$completed[[module]])
        assert_equal(record$signature, signature, paste("Completed module", module, "has changed code/input/settings. Outputs preserved; explicitly restart from the earliest affected module using --resume RUN_ID --from MODULE, or start a new run"))
        verify_module_checkpoint(run, record)
        message("SKIP VERIFIED CHECKPOINT ", module)
        timings[[length(timings) + 1L]] <- data.frame(step = module, started_utc = NA_character_, finished_utc = NA_character_, elapsed_seconds = 0, status = "SKIPPED_VERIFIED_CHECKPOINT", error = "")
        next
      }
      if (module %in% c("tables", "figures", "manuscript_summary")) {
        required <- setdiff(earlier, c("tables", "figures", "manuscript_summary"))
        assert_true(all(required %in% names(state$completed)), paste("Export requires completed statistical modules:", paste(setdiff(required, names(state$completed)), collapse = ", ")))
      } else if (module != "inputs" && is.null(data)) {
        assert_true(!is.null(input_record), "The original-data preparation checkpoint is missing")
        data <- read_verified_input_checkpoint(run, input_record)
      }
      attempt_index <- 1L
      repeat {
        attempt <- sprintf("a%03d",attempt_index)
        if (!dir.exists(file.path(run$dir,"modules",module,attempt)) && !dir.exists(file.path(run$derived_dir,module,attempt))) break
        attempt_index <- attempt_index + 1L
      }
      module_run <- run
      module_run$base_dir <- run$dir
      module_run$dir <- file.path(run$dir, "modules", module, attempt)
      module_run$derived_dir <- file.path(run$derived_dir, module, attempt)
      module_run$result_dirs <- completed_result_directories(run)
      dir.create(module_run$dir, recursive = TRUE)
      state$status <- "RUNNING"; state$active_module <- module; state$last_error <- NULL
      append_pipeline_state(run, state)
      step(module, function() {
        fn <- get(main, inherits = TRUE)
        args <- if (module %in% c("tables", "figures", "manuscript_summary")) list(run = module_run, config = config) else if (module == "inputs") list(run = module_run, config = config, expected = expected) else list(data = data, run = module_run, config = config)
        if (module != "inputs" && "expected" %in% names(formals(fn))) args$expected <- expected
        value <- do.call(fn, args)
        if (module == "inputs") data <<- value
        invisible(NULL)
      })
      record_path <- save_module_checkpoint(run, module_run, module, signature)
      state$completed[[module]] <- record_path
      state$status <- "MODULE_COMPLETED"; state$active_module <- NULL
      append_pipeline_state(run, state)
    }
    safe_write_lines("COMPUTATION_PASS_PENDING_INDEPENDENT_VALIDATION", file.path(invocation, "execution_status.txt"))
    state$status <- "REQUESTED_PLAN_COMPLETED"; state$active_module <- NULL
    append_pipeline_state(run, state)
    0L
  }, error = function(e) {
    state$status <- "FAILED"; state$last_error <- conditionMessage(e)
    append_pipeline_state(run, state)
    safe_write_lines(c("FAIL", conditionMessage(e)), file.path(invocation, "execution_status.txt"))
    message(conditionMessage(e))
    1L
  })
  safe_write_lines(capture.output(sessionInfo()), file.path(invocation, "sessionInfo.txt"))
  if (length(timings)) safe_write_csv(do.call(rbind, timings), file.path(invocation, "step_timings.csv"))
  safe_write_csv(if (length(warning_rows)) do.call(rbind, warning_rows) else data.frame(step = character(), warning = character()), file.path(invocation, "warnings.csv"))
  if (exists("v3_file_access_report", envir = .GlobalEnv, inherits = FALSE)) safe_write_csv(get("v3_file_access_report", .GlobalEnv)(), file.path(invocation, "file_access_guard.csv"))
  message("Run ID: ", run$id, "\nLogs: ", invocation, "\nResume: Rscript --vanilla run_all_development.R --resume ", run$id)
  list(status = status, run = run)
}
