publication_export <- function(root, run_id, mode = "all") {
  assert_true(grepl("^Run_[0-9]{8}_[0-9]{6}$", run_id) && mode %in% c("all", "figures", "tables"),
    "Use a valid Run ID and all|figures|tables")
  readiness <- publication_preflight(root,run_id,mode)
  assert_true(readiness$pass,"Export blocked before creating an output directory. Fix all listed preflight blockers, then rerun preflight_publication.R.")
  export <- file.path(root, "output", paste0("PublicationExport_", invocation_stamp()))
  assert_true(!dir.exists(export), "Refusing to reuse an export directory")
  folders <- c("figures/main_panels", "figures/supplementary_panels", "figures/legends", "figures/previews",
    "figures/captions", "tables/main", "tables/supplementary", "tables/csv", "tables/archive",
    "plot_data", "manifests", "validation", "logs")
  for (folder in folders) dir.create(file.path(export, folder), recursive = TRUE, showWarnings = FALSE)
  message("Publication export directory: ", export)
  status <- warnings <- list()
  config <- state <- NULL
  run <- list(root = root, id = run_id, dir = file.path(root, "output", run_id))
  run$access <- new.env(parent = emptyenv())
  run$access$files <- character(); run$access$hashes <- character()
  attempt <- function(id, kind, action) {
    start <- Sys.time(); error <- NULL
    message("START ", kind, " ", id)
    tryCatch(withCallingHandlers(action(), warning = function(w) {
      text <- conditionMessage(w)
      fatal <- kind %in% c("figure_preparation", "figure") && pub_warning_is_fatal(text)
      warnings[[length(warnings) + 1L]] <<- data.frame(item = id, kind = kind, fatal = fatal, warning = text)
      if (fatal) stop(paste("Critical plotting warning:", text), call. = FALSE)
      invokeRestart("muffleWarning")
    }), error = function(e) { error <<- conditionMessage(e) })
    label <- if (!is.null(error)) "FAIL" else switch(kind,
      figure_preparation = "PREPARED_NOT_EXPORTED", preflight = "CHECKED",
      audit = "AUDIT_WRITTEN", "WRITTEN_REQUIRES_USER_REVIEW")
    status[[length(status) + 1L]] <<- data.frame(item = id, kind = kind, status = label,
      detail = if (is.null(error)) "" else error,
      elapsed_seconds = as.numeric(difftime(Sys.time(), start, units = "secs")))
    message(label, " ", id, if (!is.null(error)) paste0(": ", error) else "")
    invisible(is.null(error))
  }
  require_step <- function(ok) assert_true(ok, "Export preflight failed; inspect export_status.csv and export_gaps.csv")
  tryCatch({
    require_step(attempt("environment_and_config", "preflight", function() {
      require_project_packages(c("ggplot2", "ggrepel", "patchwork", "dplyr", "openxlsx", "readxl", "xml2", "yaml", "digest", "systemfonts", "colorspace", "scales", "tidyr", "stringr"))
      if (mode %in% c("all","tables")) pub_stream_reader_preflight(root)
      config <<- read_config(root)
      config$project_root <<- root
      config$publication <<- yaml::read_yaml(file.path(root, "config/publication_panel_config.yml"))
      config$output$font_family <<- config$publication$font_family
      config$output$font_size <<- config$publication$font_size
      check_export_package_versions(root, c("openxlsx", "readxl", "xml2", "systemfonts"))
      assert_true("Times New Roman" %in% systemfonts::system_fonts()$family,
        "Times New Roman is absent; silent font substitution is prohibited")
      assert_file(file.path(root, relative_project_path(config$publication$panel_registry)), "Panel registry")
    }))
    require_step(attempt("saved_module_checkpoints", "preflight", function() {
      state <<- latest_pipeline_state(run$dir)
      assert_true(!is.null(state), "Selected run has no checkpoint state")
      assert_true(!dir.exists(file.path(run$dir, ".running.lock")), "The selected Run is still running; finish it before exporting")
      required <- publication_required_modules(mode)
      assert_true(all(required %in% names(state$completed)), paste("Complete these modules first:",
        paste(setdiff(required, names(state$completed)), collapse = ", ")))
      for (m in required) verify_module_checkpoint(run, checkpoint_record(run, state$completed[[m]]))
      run$result_dirs <<- completed_result_directories(run)[required]
    }))
    if (mode %in% c("all", "figures")) require_step(attempt("composition_display_data", "preflight", function() {
      assert_true("publication_display" %in% names(state$completed), "Saved composition display checkpoint is missing")
      verify_module_checkpoint(run, checkpoint_record(run,state$completed$publication_display))
      run$result_dirs[["publication_display"]] <<- file.path(run$dir,checkpoint_record(run,state$completed$publication_display)$module_directory)
      message("Using verified saved display data. Export does not start the analysis runner.")
    }))
    b <- publication_data(run, config)
    require_step(attempt("statistical_counts_and_data_contracts", "preflight", function() {
      final_beta_sources(run)
      publication_numeric_checks(b, export)
      pub_export_data_contract(b, config, export, run)
    }))
    final_manuscript_beta_text(run,export)
    safe_write_csv(b$candidates, file.path(export, "plot_data/All_pooled_pathway_candidates_and_selection.csv"))
    if (mode %in% c("all", "tables")) run_publication_tables(b, run, config, export, attempt)
    if (mode %in% c("all", "figures")) run_publication_panels(b, run, config, export, attempt)
  }, error = function(e) {
    status[[length(status) + 1L]] <<- data.frame(item = "export_aborted", kind = "export", status = "FAIL",
      detail = conditionMessage(e), elapsed_seconds = NA_real_)
    safe_write_lines(conditionMessage(e), file.path(export, "logs/fatal_error.txt"))
    message(conditionMessage(e))
  })
  # Audit errors are recorded separately; no on.exit handler can mask the
  # original error or turn an incomplete export into a successful return code.
  attempt("input_immutability", "audit", function() {
    before <- run$access$hashes
    files <- names(before)
    after <- vapply(files, function(f) if (file.exists(f)) sha256_file(f) else NA_character_, character(1))
    checks <- data.frame(input_path = substring(files, nchar(root) + 2L), before_sha256 = unname(before),
      after_sha256 = unname(after), unchanged = !is.na(after) & unname(before) == unname(after))
    safe_write_csv(checks, file.path(export, "validation/input_immutability.csv"))
    assert_true(all(checks$unchanged), "A saved input changed during export; publication package is invalid")
  })
  attempt("provenance_manifests", "audit", function() {
    sources <- unique(run$access$files)
    if (!is.null(state) && "inputs" %in% names(state$completed)) {
      record <- checkpoint_record(run, state$completed$inputs)
      sources <- unique(c(sources, file.path(root, record$input_rds)))
    }
    if (!is.null(config)) {
      paths <- file.path(root, unlist(config$publication$source_review))
      sources <- unique(c(sources, paths[file.exists(paths)]))
      yaml::write_yaml(config, file.path(export, "manifests/run_configuration.yml"))
    }
    safe_write_csv(data.frame(input_path = substring(sources, nchar(root) + 2L),
      sha256 = vapply(sources, sha256_file, character(1))), file.path(export, "manifests/input_sha256.csv"))
    if (length(run$result_dirs) && "inputs" %in% names(run$result_dirs)) {
      resources <- file.path(run$result_dirs[["inputs"]], "logs/source_resource_verification.csv")
      assert_file(resources, "Current-run raw-resource verification")
      safe_write_csv(utils::read.csv(resources, check.names = FALSE), file.path(export, "manifests/raw_resource_sha256.csv"))
    }
    packages <- unique(c(required_packages, "openxlsx", "systemfonts", "colorspace", "scales", "stringr", "tidyr", "png"))
    safe_write_csv(data.frame(package = packages, version = vapply(packages, function(p)
      if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_, character(1))),
      file.path(export, "manifests/package_versions.csv"))
    code <- c(list.files(file.path(root, "R"), pattern = "[.]R$", recursive = TRUE, full.names = TRUE),
      file.path(root, c("export_publication.R", "run_publication_candidate.R", "config/analysis_config.yml",
        "config/publication_panel_config.yml", "config/publication_panels_manifest.csv", "config/pathway_function_dictionary.csv", "renv.lock")))
    safe_write_csv(data.frame(code_path = substring(code, nchar(root) + 2L), sha256 = vapply(code, sha256_file, character(1))),
      file.path(export, "manifests/code_sha256.csv"))
    safe_write_csv(data.frame(R_version = as.character(getRversion()), platform = R.version$platform,
      locale = Sys.getlocale(), run_id = run_id, code_version = "RC4 publication candidate; pending full-run acceptance",
      submission_output_revision = "RC4_submission16_20261121",
      OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS"), OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS")),
      file.path(export, "manifests/execution_context.csv"))
  })
  attempt("session_and_readme", "audit", function() {
    safe_write_lines(capture.output(sessionInfo()), file.path(export, "sessionInfo.txt"))
    safe_write_lines(c(paste("Run:", run_id), paste("Mode:", mode),
      "Saved statistical results are immutable inputs. Composition display data have a separate current-run checkpoint.",
      "All panel preparations are checked before image export. Critical plot warnings fail the relevant task.",
      "Written files require source-content and visual review; this is not production admission.",
      "See export_status.csv, export_gaps.csv, validation/panel_preparation.csv and input_immutability.csv.",
      "No figure composites or excluded exploratory analyses are exported."), file.path(export, "README.md"))
  })
  if (length(warnings)) {
    w <- pub_bind(warnings)
    key <- paste(w$item, w$kind, w$fatal, w$warning, sep = "\r")
    rows <- !duplicated(key); w <- w[rows, , drop = FALSE]
    w$occurrences <- as.integer(table(key)[key[rows]])
  } else w <- data.frame(item = character(), kind = character(), fatal = logical(), warning = character(), occurrences = integer())
  safe_write_csv(w, file.path(export, "logs/warnings.csv"))
  stat <- pub_bind(status)
  safe_write_csv(stat, file.path(export, "export_status.csv"))
  safe_write_csv(stat[stat$status == "FAIL", , drop = FALSE], file.path(export, "export_gaps.csv"))
  files <- list.files(export, recursive = TRUE, full.names = TRUE)
  safe_write_csv(data.frame(file = substring(files, nchar(export) + 2L), bytes = file.info(files)$size,
    sha256 = vapply(files, sha256_file, character(1))), file.path(export, "output_manifest.csv"))
  message("Export status and gaps: ", file.path(export, "export_status.csv"), " ; ", file.path(export, "export_gaps.csv"))
  invisible(list(status = if (any(stat$status == "FAIL")) 1L else 0L, directory = export))
}

# Read-only export readiness. Collect all independent failures before any Export
# directory is created. This never repairs state or invokes an analysis runner.
publication_required_modules <- function(mode="all") {
  c("inputs","primary","species_loco","clinical_correlations","pooled_pathways","additive",
    "interactions","pathways","primary_diversity","diversity_interactions","bmi23_species",
    "bmi23_interactions","bmi23_pathways","bmi23_diversity","no_medication","clr",
    "metadata_audit","sample_flow","alpha_details","raw_beta","medication_sensitivity",
    "abundant_taxa",if(mode%in%c("all","figures"))"publication_display")
}

publication_preflight <- function(root,run_id,mode="all",inspect_data=TRUE,check_environment=TRUE,verbose=TRUE) {
  assert_true(grepl("^Run_[0-9]{8}_[0-9]{6}$",run_id)&&mode%in%c("all","figures","tables"),"Invalid export Run/mode")
  checks<-list();config<-state<-NULL
  run<-list(root=root,id=run_id,dir=file.path(root,"output",run_id))
  run$access<-new.env(parent=emptyenv());run$access$files<-character();run$access$hashes<-character()
  check<-function(id,action){
    error<-"";ok<-tryCatch({action();TRUE},error=function(e){error<<-conditionMessage(e);FALSE})
    checks[[length(checks)+1L]]<<-data.frame(check=id,pass=ok,detail=error)
    if(verbose)message(if(ok)"READY " else "BLOCKED ",id,if(nzchar(error))paste0(": ",error))
    invisible(ok)
  }
  check("configuration",function(){
    config<<-read_config(root);config$project_root<<-root
    config$publication<<-yaml::read_yaml(file.path(root,"config/publication_panel_config.yml"))
    config$output$font_family<<-config$publication$font_family;config$output$font_size<<-config$publication$font_size
    assert_true(config$publication$expected_workbooks==16L,"Expected16 workbook plans")
    registry<-utils::read.csv(file.path(root,relative_project_path(config$publication$panel_registry)))
    assert_true(nrow(registry)==40L&&!anyDuplicated(registry$figure_id),"Expected40 unique panel registrations")
  })
  if(check_environment) {
    check("R_and_export_packages",function(){
      assert_true(as.character(getRversion())=="4.5.1","Frozen R version must be4.5.1")
      require_project_packages(c("ggplot2","ggrepel","patchwork","dplyr","openxlsx","readxl","xml2","yaml","digest","systemfonts","colorspace","scales","tidyr","stringr","jsonlite"))
      check_export_package_versions(root,c("openxlsx","readxl","xml2","systemfonts"))
    })
    check("font",function()assert_true("Times New Roman"%in%systemfonts::system_fonts()$family,"Times New Roman is absent"))
    if(mode%in%c("all","tables"))check("independent_workbook_reader",function()pub_stream_reader_preflight(root))
  }
  check("Run_identity",function()assert_rc4_run_identity(run))
  check("execution_lock",function(){
    lock<-file.path(run$dir,".running.lock")
    owner<-file.path(lock,"owner.txt")
    detail<-if(file.exists(owner))paste(readLines(owner,warn=FALSE),collapse="; ") else "no owner record"
    assert_true(!dir.exists(lock),paste("Execution lock exists:",detail,"; inspect process and checkpoint state together; never delete blindly"))
  })
  check("checkpoint_state",function(){
    state<<-latest_pipeline_state(run$dir)
    assert_true(!is.null(state)&&identical(state$run_id,run_id),"Missing or incorrect checkpoint state")
    assert_true(is.null(state$active_module),paste("Interrupted/active checkpoint state:",state$status,"module",state$active_module,
      "; a missing lock does not mean the interrupted plan completed"))
  })
  dirs<-character();required<-publication_required_modules(mode)
  for(m in required)check(paste0("checkpoint:",m),local({module<-m;function(){
    assert_true(!is.null(state)&&module%in%names(state$completed),paste("Not registered complete:",module))
    rel<-state$completed[[module]]
    assert_true(grepl(paste0("^modules/",module,"/a[0-9]{3}/COMPLETED[.]yml$"),rel),paste("Invalid registered module path:",module))
    record<-checkpoint_record(run,rel)
    assert_true(identical(record$module,module),paste("Checkpoint module identity differs:",module))
    verify_module_checkpoint(run,record)
    dirs[[module]]<<-file.path(run$dir,relative_project_path(record$module_directory))
  }}))
  if(!is.null(config))for(type in c("cohort","confounder"))check(paste0("source_review:",type),local({kind<-type;function()pub_source_review(root,config,kind,run)}))
  complete<-all(required%in%names(dirs))
  if(complete&&!is.null(config)&&inspect_data) {
    run$result_dirs<-dirs
    b<-NULL
    readable<-check("all_saved_result_sources",function(){b<<-publication_data(run,config);final_beta_sources(run)})
    if(readable) {
      check("saved_numeric_and_model_contracts",function(){
        tmp<-tempfile("publication_readonly_contract_");dir.create(tmp);on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
        publication_numeric_checks(b,tmp);pub_export_data_contract(b,config,tmp,run)
      })
      check("all16_table_sources_and_S2_annotations",function(){
        plan<-publication_table_plan(b,run,config)
        assert_true(length(plan)==nrow(submission_table_manifest()),"Expected16 table plans")
        submission_validate_table_identity(data.frame(table_id=names(plan)))
        submission_validate_s7(root,b$candidates,plan$Supplementary_Table_S7$sheets$Candidates_and_selection)
        for(id in names(plan)){
          item<-plan[[id]];if(is.function(item))item<-item()
          assert_true(is.list(item$sheets)&&length(item$sheets)>0L,paste("Missing sheet plan:",id))
          for(sheet in item$sheets)assert_true(ncol(as.data.frame(sheet))>0L,paste("Missing sheet columns:",id))
        }
      })
    }
  }
  result<-do.call(rbind,checks)
  if(verbose)message(if(all(result$pass))"EXPORT PREFLIGHT PASS: no statistical model or file export was run." else
    "EXPORT PREFLIGHT FAIL: all detected blockers are listed above. No PublicationExport directory was created.")
  list(pass=all(result$pass),checks=result,config=config,run=run,state=state)
}
