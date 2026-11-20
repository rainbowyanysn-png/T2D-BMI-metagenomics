#!/usr/bin/env Rscript
# Author-invoked engineering recovery; default is read-only PLAN.
# No MaAsLin2/DESeq2/REML/CLR/diversity/abundant-taxa function is called.
local({
  if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
  root<-normalizePath(dirname(entry),winslash="/",mustWork=TRUE)
  source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
  source(file.path(root,"R/functions/analysis_bootstrap.R"),local=environment());load_analysis_support(root,environment())
  source(file.path(root,"R/functions/publication_bootstrap.R"),local=environment());load_publication_support(root,environment())
  args<-commandArgs(TRUE)
  assert_true(length(args)>=2L,"Use RUN_ID State_NNNNNN.yml [--apply --refresh-medication-summary]")
  id<-args[1];baseline_name<-args[2];flags<-args[-(1:2)]
  assert_true(grepl("^Run_[0-9]{8}_[0-9]{6}$",id)&&grepl("^State_[0-9]{6}[.]yml$",baseline_name),"Explicit Run and baseline snapshot required")
  assert_true(!anyDuplicated(flags)&&all(flags%in%c("--apply","--refresh-medication-summary")),"Unknown/duplicate recovery option")
  apply<-"--apply"%in%flags;refresh<-"--refresh-medication-summary"%in%flags
  assert_true(!refresh||apply,"--refresh-medication-summary requires --apply")
  run<-list(root=root,id=id,dir=file.path(root,"output",id),derived_dir=file.path(root,"data_derived",id))
  assert_rc4_run_identity(run);assert_true(!dir.exists(file.path(run$dir,".running.lock")),"Run is locked; no recovery permitted")
  state<-latest_pipeline_state(run$dir);assert_true(!is.null(state),"Missing Run state")
  baseline_file<-file.path(run$dir,"checkpoints",baseline_name);assert_file(baseline_file,"Explicit recovery baseline")
  baseline<-yaml::read_yaml(baseline_file);assert_true(identical(baseline$run_id,id),"Baseline belongs to another Run")
  snapshots<-sort(list.files(file.path(run$dir,"checkpoints"),pattern="^State_[0-9]+[.]yml$",full.names=TRUE))
  head_file<-tail(snapshots,1);head_hash<-sha256_file(head_file)
  assert_true(baseline_name<basename(head_file),"Recovery baseline must precede current state")
  registry<-pipeline_registry(root);sequence<-c("inputs",registry$module)
  for(file in unique(c("01_import_harmonize_qc.R",registry$file)))source(file.path(root,"R/pipeline",file),local=environment())
  config<-read_config(root);expected<-read_expected(root)
  contract<-utils::read.csv(file.path(root,"config/statistical_implementation_sha256.csv"))
  assert_true(all(vapply(file.path(root,contract$file),sha256_file,character(1))==contract$sha256),"Frozen statistical implementation differs")
  required<-publication_required_modules("all");missing<-setdiff(required,names(state$completed))
  # This recovery is deliberately bounded to the interrupted publication tail.
  permitted<-c("raw_beta","medication_sensitivity","abundant_taxa","publication_display")
  assert_true(all(missing%in%permitted),paste("Recovery cannot replace missing core analyses:",paste(setdiff(missing,permitted),collapse=", ")))
  view<-state$completed
  for(m in names(view))verify_module_checkpoint(run,checkpoint_record(run,view[[m]]))
  input<-checkpoint_record(run,view$inputs)
  signature<-function(module,completed){
    earlier<-sequence[seq_len(match(module,sequence)-1L)]
    names_upstream<-intersect(earlier,names(completed))
    upstream<-setNames(lapply(completed[names_upstream],function(path)checkpoint_record(run,path)$artifact_index_sha256),names_upstream)
    module_signature(module,registry$function_name[match(module,registry$module)],environment(),run,config,expected,input$input_sha256,upstream)
  }
  # Verify exact primary/no-covariate sources of the old medication summary.
  same_medication_inputs<-function(completed,record){
    context<-run;context$result_dirs<-setNames(vapply(completed,function(path)file.path(run$dir,checkpoint_record(run,path)$module_directory),character(1)),names(completed))
    folder<-file.path(run$dir,record$module_directory,"tables/09b_medication_sensitivity")
    for(bmi in c("BMI25","BMI23"))for(kind in c("adjusted","unadjusted")){
      relative<-if(kind=="adjusted")if(bmi=="BMI25")"02_primary_BMI25_species/OT2D_vs_NT2D_BMI25_all.csv" else
        "07_BMI23_sensitivity/OT2D_vs_NT2D_BMI23_all.csv" else paste0("09_no_medication_covariate/",bmi,"_all.csv")
      actual<-read_current_result(context,file.path("tables",relative))
      if(kind=="unadjusted") {
        reference<-read_current_result(context,file.path("tables",if(bmi=="BMI25")"02_primary_BMI25_species/OT2D_vs_NT2D_BMI25_all.csv" else "07_BMI23_sensitivity/OT2D_vs_NT2D_BMI23_all.csv"))
        assert_true(!anyDuplicated(actual$feature)&&setequal(actual$feature,reference$feature),"Medication source feature set differs")
        actual<-actual[match(reference$feature,actual$feature),,drop=FALSE]
      }
      # Replay only the saved-summary's documented row alignment and CSV writing.
      # Its input/output P/q/coefficients are never recalculated or approximated.
      original<-file.path(folder,paste0(bmi,"_",kind,"_all.csv"))
      staged<-tempfile("medication_source_copy_",fileext=".csv")
      tryCatch({safe_write_csv(actual,staged);assert_true(identical(sha256_file(staged),sha256_file(original)),
        paste("Medication-summary input differs:",bmi,kind))},finally=unlink(staged))
    }
    invisible(TRUE)
  }
  plan<-list();pending<-NULL
  for(m in intersect(sequence,missing)){
    assert_true(m%in%names(baseline$completed),paste("Baseline does not explicitly register:",m))
    rel<-baseline$completed[[m]];record<-checkpoint_record(run,rel)
    assert_true(identical(record$module,m)&&grepl(paste0("^modules/",m,"/a[0-9]{3}/COMPLETED[.]yml$"),rel),"Baseline module identity mismatch")
    assert_true(identical(record$module_directory,dirname(rel))&&startsWith(record$artifact_index,paste0(dirname(rel),"/")),"Baseline checkpoint path binding differs")
    verify_module_checkpoint(run,record);now<-signature(m,view)
    same<-identical(now,record$signature)
    action<-if(same)"REUSE_EXACT_SIGNATURE_AND_ARTIFACTS" else "BLOCKED_SIGNATURE_MISMATCH"
    if(!same&&m=="medication_sensitivity"){
      same_medication_inputs(view,record);action<-"REFRESH_SAVED_MEDICATION_SUMMARY_ONLY";pending<-list(module=m,baseline=record)
    }
    plan[[m]]<-data.frame(module=m,baseline_checkpoint=rel,recorded_signature=record$signature,current_signature=now,
      exact_signature_match=same,artifacts_verified=TRUE,action=action)
    if(same)view[[m]]<-rel
  }
  report<-if(length(plan))do.call(rbind,plan) else data.frame(module=character(),action=character())
  print(report[,intersect(c("module","baseline_checkpoint","exact_signature_match","action"),names(report)),drop=FALSE],row.names=FALSE)
  assert_true(!any(report$action=="BLOCKED_SIGNATURE_MISMATCH"),"A model/input signature differs; no state recovery or model rerun permitted by this command")
  if(!apply){cat("PLAN ONLY. No state or statistical result changed. Exact-signature reuse is not yet registered.\n");quit(status=0L)}
  assert_true(is.null(pending)||refresh,"Medication summary needs an explicit --refresh-medication-summary flag; it will not be rebound with a mismatched signature")
  if(!length(missing)) {
    message("No missing publication checkpoints: no recovery state or summary written.")
    x<-publication_preflight(root,id,"all");quit(status=if(x$pass)0L else 1L)
  }
  # No existing file is overwritten. Lock plus head hash protect the state transition.
  lock<-file.path(run$dir,".running.lock");assert_true(dir.create(lock,showWarnings=FALSE),"Concurrent process acquired the Run lock")
  safe_write_lines(c(paste("PID",Sys.getpid()),"CHECKPOINT_RECONCILIATION_NO_MODEL_FITTING"),file.path(lock,"owner.txt"))
  release_own_lock<-function(){
    target<-normalizePath(lock,winslash="/",mustWork=FALSE)
    assert_true(identical(target,paste0(normalizePath(run$dir,winslash="/"),"/.running.lock")),"Refusing cleanup outside the selected Run lock")
    owner<-file.path(lock,"owner.txt")
    if(file.exists(owner)&&identical(readLines(owner,warn=FALSE)[1],paste("PID",Sys.getpid())))unlink(lock,recursive=TRUE)
  }
  on.exit(release_own_lock(),add=TRUE)
  assert_true(identical(sha256_file(tail(sort(list.files(file.path(run$dir,"checkpoints"),pattern="^State_[0-9]+[.]yml$",full.names=TRUE)),1)),head_hash),"Run state changed during recovery planning")
  if(!is.null(pending)){
    m<-pending$module;index<-1L
    while(dir.exists(file.path(run$dir,"modules",m,sprintf("a%03d",index))))index<-index+1L
    module_run<-run;module_run$base_dir<-run$dir;module_run$dir<-file.path(run$dir,"modules",m,sprintf("a%03d",index))
    module_run$result_dirs<-setNames(vapply(view,function(path)file.path(run$dir,checkpoint_record(run,path)$module_directory),character(1)),names(view))
    dir.create(module_run$dir,recursive=TRUE)
    message("Refreshing saved medication comparison only; no model or P/q computation.")
    run_medication_sensitivity(data=NULL,run=module_run,config=config)
    old<-file.path(run$dir,pending$baseline$module_directory,"tables/09b_medication_sensitivity")
    fresh<-file.path(module_run$dir,"tables/09b_medication_sensitivity")
    names_old<-list.files(old,pattern="[.]csv$");names_new<-list.files(fresh,pattern="[.]csv$")
    assert_true(identical(names_old,names_new)&&length(names_old)>0L,"Medication summary file family differs")
    assert_true(identical(unname(vapply(file.path(old,names_old),sha256_file,character(1))),unname(vapply(file.path(fresh,names_new),sha256_file,character(1)))),
      "Refreshed medication tables differ byte-for-byte from the verified baseline; no recovery state committed")
    view[[m]]<-save_module_checkpoint(run,module_run,m,signature(m,view))
  }
  # Keep all current references, including CLR a003; never roll back to baseline.
  assert_true(all(vapply(names(state$completed),function(m)identical(view[[m]],state$completed[[m]]),logical(1))),"Recovery attempted to replace a current completed checkpoint")
  event<-file.path(run$dir,"checkpoints",paste0("Reconciliation_",invocation_stamp(),".csv"))
  if(nrow(report))safe_write_csv(report,event)
  state$completed<-view;state$status<-"CHECKPOINTS_RECONCILED_PENDING_EXPORT_PREFLIGHT";state$active_module<-NULL
  state$recovery_baseline_snapshot<-baseline_name;state$recovery_event<-substring(event,nchar(run$dir)+2L)
  append_pipeline_state(run,state)
  release_own_lock()
  message("Checkpoint reconciliation committed. No statistical model or saved result was overwritten.")
  x<-publication_preflight(root,id,"all")
  quit(status=if(x$pass)0L else 1L)
})
