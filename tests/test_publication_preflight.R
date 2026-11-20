#!/usr/bin/env Rscript
# Read-only readiness tests. Never creates a real Export or resumes an analysis.
local({
  if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));root<-normalizePath(dirname(dirname(entry)),winslash="/",mustWork=TRUE)
  source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
  source(file.path(root,"R/functions/publication_bootstrap.R"),local=environment());load_publication_support(root,environment())
  stopifnot(length(publication_required_modules("all"))==23L,length(publication_required_modules("tables"))==22L)
  # Independent missing-module failures must all be reported, not hidden by the first one.
  fixture<-tempfile("preflight_fixture_");dir.create(file.path(fixture,"config"),recursive=TRUE)
  on.exit(unlink(fixture,recursive=TRUE),add=TRUE)
  file.copy(file.path(root,"config/publication_panel_config.yml"),file.path(fixture,"config/publication_panel_config.yml"))
  file.copy(file.path(root,"config/publication_panels_manifest.csv"),file.path(fixture,"config/publication_panels_manifest.csv"))
  real_read<-read_config
  # Read config from the actual code checkout; all Run paths stay inside the fixture.
  read_config<-function(unused)real_read(root)
  z<-publication_preflight(fixture,"Run_20000101_000000","all",inspect_data=FALSE,check_environment=FALSE,verbose=FALSE)
  stopifnot(!z$pass,sum(grepl("^checkpoint:",z$checks$check)&!z$checks$pass)==23L,!dir.exists(file.path(fixture,"output")))
  e<-suppressWarnings(suppressMessages(tryCatch({publication_export(fixture,"Run_20000101_000000","all");NULL},error=function(e)e)))
  stopifnot(inherits(e,"error"),!dir.exists(file.path(fixture,"output")))
  # Interrupted state remains a blocker even without a lock directory.
  run_dir<-file.path(fixture,"output/Run_20000101_000000");dir.create(file.path(run_dir,"checkpoints"),recursive=TRUE)
  writeLines("RC4_CURRENT_PROJECT_RUN",file.path(run_dir,"RC4_RUN_IDENTITY.txt"))
  yaml::write_yaml(list(run_id="Run_20000101_000000",completed=list(),status="RUNNING",active_module="raw_beta"),file.path(run_dir,"checkpoints/State_000001.yml"))
  z<-publication_preflight(fixture,"Run_20000101_000000","all",inspect_data=FALSE,check_environment=FALSE,verbose=FALSE)
  stopifnot(!z$checks$pass[z$checks$check=="checkpoint_state"],!dir.exists(file.path(run_dir,".running.lock")))
  # Export must no longer start the runner or create a timestamped output before preflight.
  implementation<-paste(deparse(publication_export),collapse="\n")
  stopifnot(!grepl("run_publication_updates|system2",implementation),
    regexpr("publication_preflight",implementation)[1]<regexpr("dir.create",implementation)[1])
  # Optional current-Run saved-summary test. All generated CSVs live in TEMP;
  # it never calls the analysis runner or changes a checkpoint.
  args<-commandArgs(TRUE)
  if(length(args)) {
    stopifnot(length(args)==2L,grepl("^Run_[0-9]{8}_[0-9]{6}$",args[1]),grepl("^State_[0-9]{6}[.]yml$",args[2]))
    selected<-list(root=root,id=args[1],dir=file.path(root,"output",args[1]))
    current<-latest_pipeline_state(selected$dir);old<-yaml::read_yaml(file.path(selected$dir,"checkpoints",args[2]));view<-current$completed
    for(m in setdiff(publication_required_modules("all"),names(view)))view[[m]]<-old$completed[[m]]
    selected$result_dirs<-setNames(vapply(view,function(path)file.path(selected$dir,checkpoint_record(selected,path)$module_directory),character(1)),names(view))
    for(path in view)verify_module_checkpoint(selected,checkpoint_record(selected,path))
    old_summary<-file.path(selected$result_dirs[["medication_sensitivity"]],"tables/09b_medication_sensitivity")
    target<-selected;target$dir<-file.path(fixture,"temporary_summary")
    source(file.path(root,"R/pipeline/09b_medication_sensitivity.R"),local=environment())
    run_medication_sensitivity(NULL,target,real_read(root))
    new_summary<-file.path(target$dir,"tables/09b_medication_sensitivity")
    names_old<-list.files(old_summary,pattern="[.]csv$");names_new<-list.files(new_summary,pattern="[.]csv$")
    stopifnot(length(names_old)==16L,identical(names_old,names_new),
      identical(unname(vapply(file.path(old_summary,names_old),sha256_file,character(1))),unname(vapply(file.path(new_summary,names_new),sha256_file,character(1)))))
    cat("PASS all16 TEMP medication-summary CSVs byte-identical to baseline; no models or P/q calculations\n")
  }
  cat("PASS complete blocker collection; no Export directory on failure\n")
  cat("PASS orphaned RUNNING state rejected even after lock removal\n")
  cat("PASS exporter is saved-result-only and does not call analysis runner\n")
})
