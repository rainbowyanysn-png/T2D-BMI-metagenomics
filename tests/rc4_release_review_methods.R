# Bind export provenance without mistaking identical rerun attempts for different statistics.
# Full checkpoint signatures AND every CLR table/metadata hash must match; no basename-only fallback.
review_export_clr_binding <- function(root,run,export,current) {
  inputs<-review_read(file.path(export,"manifests/input_sha256.csv"));assert_columns(inputs,c("input_path","sha256"),"Export input provenance")
  prefix<-paste0("output/",run$id,"/modules/clr/")
  selected<-inputs[startsWith(inputs$input_path,prefix),,drop=FALSE]
  assert_true(nrow(selected)>0L,"Selected export has no CLR input provenance")
  remainder<-substring(selected$input_path,nchar(prefix)+1L)
  attempts<-unique(sub("/.*$","",remainder))
  assert_true(length(attempts)==1L&&grepl("^a[0-9]{3}$",attempts),"Export contains mixed/invalid CLR attempts")
  exported<-file.path(run$dir,"modules/clr",attempts)
  current_record<-yaml::read_yaml(file.path(current,"COMPLETED.yml"))
  export_record<-yaml::read_yaml(file.path(exported,"COMPLETED.yml"))
  assert_true(identical(current_record$module,"clr")&&identical(export_record$module,"clr")&&
    is.character(current_record$signature)&&length(current_record$signature)==1L&&nzchar(current_record$signature)&&
    identical(current_record$signature,export_record$signature),"CLR attempts have different model/input signatures; re-export current Run")
  verify_module_checkpoint(run,current_record);verify_module_checkpoint(run,export_record)
  files<-list.files(file.path(current,"tables/08_CLR_sensitivity"),pattern="[.]csv$",full.names=TRUE)
  assert_true(length(files)>0L,"Current CLR checkpoint has no result tables")
  names_expected<-paste0("tables/08_CLR_sensitivity/",basename(files))
  suffix<-substring(remainder,nchar(attempts)+2L)
  assert_true(!anyDuplicated(suffix)&&setequal(suffix,names_expected),"Export CLR provenance is incomplete, duplicated or has extra files")
  evidence<-list()
  for(i in seq_along(files)){
    j<-match(names_expected[i],suffix);export_file<-file.path(root,relative_project_path(selected$input_path[j]))
    current_hash<-sha256_file(files[i]);export_hash<-sha256_file(export_file)
    assert_true(identical(current_hash,export_hash)&&identical(export_hash,selected$sha256[j]),paste("CLR current/export content hash differs:",basename(files[i]),"; re-export current Run"))
    evidence[[i]]<-data.frame(result_file=basename(files[i]),current_attempt=basename(current),export_attempt=attempts,current_sha256=current_hash,export_sha256=export_hash,export_manifest_sha256=selected$sha256[j],current_model_signature=current_record$signature,export_model_signature=export_record$signature,same_signature=TRUE,byte_identical=TRUE)
  }
  list(export_attempt=attempts,relation=if(identical(basename(current),attempts))"SAME_ATTEMPT"else"BYTE_IDENTICAL_COMPLETE_ATTEMPTS_WITH_IDENTICAL_SIGNATURE",evidence=do.call(rbind,evidence))
}

# Author-run review/freeze functions. Definitions only; never called by analysis/export.
review_read <- function(f) {assert_file(f,"Review input");utils::read.csv(f,check.names=FALSE,stringsAsFactors=FALSE,na.strings=c("NA",""))}
review_run <- function(root,id) {
 assert_true(grepl("^Run_[0-9]{8}_[0-9]{6}$",id),"Invalid RC4 Run identifier")
 run<-list(root=root,id=id,dir=file.path(root,"output",id));run$result_dirs<-completed_result_directories(run)
 run$access<-new.env(parent=emptyenv());run$access$files<-character();run$access$hashes<-character();run
}
review_new_external <- function(root,path) {
 path<-normalizePath(path,winslash="/",mustWork=FALSE);r<-tolower(normalizePath(root,winslash="/",mustWork=TRUE))
 assert_true(!file.exists(path)&&!dir.exists(path),"Review/archive destination already exists; refusing overwrite")
 assert_true(!identical(tolower(path),r)&&!startsWith(tolower(path),paste0(r,"/")),"Review/reference files must be outside the public project")
 assert_true(dir.create(path,recursive=TRUE),"Cannot create review destination");path
}
review_write <- function(x,f) safe_write_csv(x,f)
review_key <- function(d) paste(d$analysis_domain[1],sub("_UnadjustedPendingMedicationAudit$","",d$analysis_id[1]),sep="|")
review_changed_key <- function(k) k%in%c("Taxonomy|Overweight_vs_Lean_T2D","Taxonomy|BMI23_Formal_Interaction","Pathway|BMI23_Formal_Interaction")
review_validate_decisions <- function(d,summary) {
 assert_columns(d,c("result_file","source_sha256","decision","accepted_significant_count","reviewer","review_date","rationale","visual_review"),"CLR author decisions")
 assert_true(nrow(d)==16L&&!anyDuplicated(d$result_file)&&setequal(d$result_file,summary$result_file),"Author decisions must cover all16 unique CLR results")
 d<-d[match(summary$result_file,d$result_file),,drop=FALSE]
 assert_true(identical(d$source_sha256,summary$source_sha256),"Decision hashes differ from prepared CLR snapshot")
 assert_true(all(d$decision=="ACCEPT")&&all(d$visual_review=="CONFIRMED"),"Author acceptance/visual review is incomplete")
 for(n in c("reviewer","review_date","rationale"))assert_true(all(!is.na(d[[n]])&nzchar(trimws(d[[n]]))),paste("Missing author field:",n))
 assert_true(all(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$",d$review_date)),"Review dates require YYYY-MM-DD")
 counts<-suppressWarnings(as.numeric(d$accepted_significant_count))
 assert_true(all(is.finite(counts))&&all(counts==floor(counts))&&all(counts==summary$n_significant),"Explicit author significant counts differ from reviewed saved results")
 invisible(TRUE)
}
review_prepare <- function(root,run_id,export_id,diagnostics,baseline_attempt,out) {
 run<-review_run(root,run_id);assert_true(grepl("^PublicationExport_[0-9]{8}_[0-9]{6}[.][0-9]+_[0-9]+$",export_id),"Invalid export identifier")
 export<-file.path(root,"output",export_id);ctx<-review_read(file.path(export,"manifests/execution_context.csv"));assert_true(nrow(ctx)==1L&&identical(ctx$run_id,run_id),"Export/Run identity mismatch")
 contracts<-pub_clr_alignment_contract(run);state<-latest_pipeline_state(run$dir);verify_module_checkpoint(run,checkpoint_record(run,state$completed$clr))
 assert_true(grepl("^a[0-9]{3}$",baseline_attempt),"Baseline must be an explicit attempt ID")
 current<-run$result_dirs[["clr"]];baseline<-file.path(run$dir,"modules/clr",baseline_attempt)
 export_binding<-review_export_clr_binding(root,run,export,current)
 workbook_evidence<-review_workbook_evidence(export)
 assert_true(!identical(normalizePath(current,winslash="/"),normalizePath(baseline,winslash="/")),"Baseline cannot be current CLR attempt")
 oldfiles<-list.files(file.path(baseline,"tables/08_CLR_sensitivity"),pattern="_all[.]csv$",full.names=TRUE);assert_true(length(oldfiles)==16L,"Explicit historical attempt must have16 CLR families")
 old<-lapply(oldfiles,review_read);names(old)<-vapply(old,review_key,character(1));assert_true(!anyDuplicated(names(old)),"Ambiguous historical CLR family IDs")
 assert_true(grepl("^REML_Diagnostics_[0-9]{8}_[0-9]{6}[.][0-9]+_[0-9]+$",diagnostics),"Explicit diagnostic directory ID required")
 diag<-file.path(root,"output",diagnostics);all_diag<-review_read(file.path(diag,"all_diagnostic_refits.csv"));warnings<-review_read(file.path(diag,"feature_level_warnings.csv"));ih<-review_read(file.path(diag,"read_only_input_hashes.csv"))
 assert_true(nrow(ih)>0L&&!anyNA(ih$unchanged)&&all(ih$unchanged),"Diagnostic input immutability evidence is missing/failed")
 assert_true(any(grepl(paste0("/output/",run_id,"/"),gsub("\\\\","/",ih$file),fixed=TRUE)),"Diagnostic provenance does not identify requested Run")
 for(i in seq_len(nrow(ih)))assert_true(identical(sha256_file(ih$file[i]),ih$after_sha256[i]),paste("Diagnostic source changed:",ih$file[i]))
 out<-review_new_external(root,out);review_write(export_binding$evidence,file.path(out,"CLR_export_checkpoint_binding.csv"));dir.create(file.path(out,"clr_snapshot"));dir.create(file.path(out,"baseline_clr_evidence"));dir.create(file.path(out,"publication_review_evidence"))
 evidence<-c(file.path(export,c("tables/supplementary/Supplementary_Table_S2.xlsx","tables/supplementary/Supplementary_Table_S13.xlsx")),workbook_evidence)
 for(id in c("SFig4B_sensitivity_CLR_species","SFig4C_sensitivity_CLR_pathways"))evidence<-c(evidence,file.path(export,"figures/previews",paste0(id,".png")),file.path(export,"figures/supplementary_panels",paste0(id,".pdf")),file.path(export,"plot_data",paste0(id,".csv")),file.path(export,"figures/captions",paste0(id,".txt")))
 for(f in evidence){assert_file(f,"CLR publication review evidence");assert_true(file.copy(f,file.path(out,"publication_review_evidence",basename(f)),overwrite=FALSE),"Cannot archive publication review evidence")}
 for(f in oldfiles)assert_true(file.copy(f,file.path(out,"baseline_clr_evidence",basename(f)),overwrite=FALSE),"Cannot archive historical CLR evidence")
 sources<-unique(c(file.path(export,"manifests/input_sha256.csv"),evidence,file.path(current,"tables/08_CLR_sensitivity",contracts$result_file),file.path(current,"tables/08_CLR_sensitivity/CLR_model_contracts.csv"),file.path(current,"tables/08_CLR_sensitivity/CLR_sample_metadata.csv"),oldfiles,file.path(diag,c("all_diagnostic_refits.csv","feature_level_warnings.csv","read_only_input_hashes.csv")),file.path(export,"figure_index.csv")))
 before<-vapply(sources,sha256_file,character(1));summary<-comparison<-list()
 fields<-c("coef_CLR","stderr_moderated","t_statistic","pval_CLR","qval_CLR_BH")
 for(i in seq_len(nrow(contracts))){
  f<-contracts$result_file[i];source<-file.path(current,"tables/08_CLR_sensitivity",f);d<-review_read(source);k<-review_key(d);o<-old[[k]];assert_true(!is.null(o),paste("Historical CLR family absent:",k))
  changed<-review_changed_key(k);same_ids<-setequal(d$feature,o$feature)&&!anyDuplicated(d$feature)&&!anyDuplicated(o$feature)
  same_order<-identical(d$feature,o$feature);numeric_same<-same_ids;max_difference<-0
  if(same_ids){o<-o[match(d$feature,o$feature),,drop=FALSE];for(n in fields){assert_columns(o,n,"Historical scientific fields");delta<-abs(d[[n]]-o[[n]]);max_difference<-max(max_difference,max(delta,na.rm=TRUE));numeric_same<-numeric_same&&isTRUE(all.equal(d[[n]],o[[n]],tolerance=1e-8,check.attributes=FALSE))};numeric_same<-numeric_same&&identical(sign(d$coef_CLR),sign(o$coef_CLR))&&identical(d$qval_CLR_BH<.05,o$qval_CLR_BH<.05)}
  assert_true(isTRUE(all.equal(d$qval_CLR_BH,p.adjust(d$pval_CLR,"BH"),tolerance=1e-8,check.attributes=FALSE)),paste("CLR saved BH family differs:",f))
  assert_true(file.copy(source,file.path(out,"clr_snapshot",f),overwrite=FALSE),"Cannot snapshot CLR result")
  summary[[i]]<-data.frame(result_file=f,source_sha256=sha256_file(source),family_key=k,changed_model=changed,n_tested=nrow(d),n_significant=sum(d$qval_CLR_BH<.05),positive=sum(d$qval_CLR_BH<.05&d$coef_CLR>0),negative=sum(d$qval_CLR_BH<.05&d$coef_CLR<0),n_samples=unique(d$n_samples),model_formula=unique(d$model_formula),target_term=unique(d$target_term))
  comparison[[i]]<-data.frame(result_file=f,changed_model=changed,same_feature_set=same_ids,same_feature_order=same_order,scientific_fields_same=numeric_same,max_absolute_difference=max_difference,old_significant=sum(o$qval_CLR_BH<.05),new_significant=sum(d$qval_CLR_BH<.05))
 }
 summary<-do.call(rbind,summary);comparison<-do.call(rbind,comparison)
 assert_true(sum(summary$changed_model)==3L,"Expected exactly three changed CLR families")
 review_write(summary,file.path(out,"CLR_result_summary.csv"));review_write(comparison,file.path(out,"CLR_old_new_scientific_comparison.csv"))
 decisions<-summary[c("result_file","source_sha256")];decisions$decision<-"PENDING";decisions$accepted_significant_count<-NA_integer_;decisions$reviewer<-decisions$review_date<-decisions$rationale<-"";decisions$visual_review<-"PENDING";review_write(decisions,file.path(out,"CLR_AUTHOR_DECISIONS.csv"))
 for(n in c("CLR_model_contracts.csv","CLR_sample_metadata.csv"))assert_true(file.copy(file.path(current,"tables/08_CLR_sensitivity",n),file.path(out,"clr_snapshot",n),overwrite=FALSE),"Cannot copy CLR provenance")
 assert_columns(warnings,c("module","comparison","feature","warning","formal_significant","entered_representative_forest","finite_effect_and_se","valid_nonnegative_tau2","within_original_absolute_tolerance"),"Feature-level REML warnings")
 key<-paste(warnings$module,warnings$comparison,warnings$feature,sep="|");warnings$review_key<-vapply(key,function(x)digest::digest(x,algo="sha256",serialize=FALSE),character(1))
 # Never infer harmlessness from tau2=0, nonsignificance or agreement alone.
 warnings$risk_class<-ifelse(is.na(warnings$finite_effect_and_se)|is.na(warnings$valid_nonnegative_tau2)|is.na(warnings$within_original_absolute_tolerance)|!warnings$finite_effect_and_se|!warnings$valid_nonnegative_tau2|!warnings$within_original_absolute_tolerance,"NUMERICAL_DIFFERENCE_OR_INVALID_ESTIMATE",ifelse(warnings$formal_significant|warnings$entered_representative_forest,"FORMAL_OR_FOREST_FEATURE_REQUIRES_FOCUSED_REVIEW","DIAGNOSTIC_REPRODUCED_WARNING_REQUIRES_REVIEW"))
 warnings$historical_attribution<-"DIAGNOSTIC_REPLAY; ORIGINAL_AGGREGATE_LOG_HAS_NO_FEATURE_ID"
 review_write(warnings,file.path(out,"REML_warning_classification.csv"))
 aggregate_logs<-list.files(file.path(run$dir,"logs"),pattern="^warnings[.]csv$",recursive=TRUE,full.names=TRUE)
 logged<-unlist(lapply(aggregate_logs,function(f){d<-review_read(f);if(nrow(d))as.character(d$step[grepl("Fisher scoring|local maximum|Setting tau",d$warning,ignore.case=TRUE)])else character()}))
 counts<-data.frame(module=c("pathways","bmi23_pathways"),aggregate_log_count=as.integer(table(factor(logged,levels=c("pathways","bmi23_pathways")))),diagnostic_replay_count=as.integer(table(factor(warnings$module,levels=c("pathways","bmi23_pathways")))))
 counts$count_agreement<-counts$aggregate_log_count==counts$diagnostic_replay_count;counts$exact_historical_feature_attribution<-"NOT_PROVABLE_FROM_AGGREGATE_LOG"
 review_write(counts,file.path(out,"REML_warning_count_comparison.csv"))
 decisions<-unique(warnings[c("review_key","module","comparison","feature","risk_class")]);decisions$decision<-"PENDING";decisions$reviewer<-decisions$review_date<-decisions$rationale<-decisions$evidence_file<-"";review_write(decisions,file.path(out,"REML_AUTHOR_DECISIONS.csv"))
 review_write(all_diag,file.path(out,"REML_full_replay_stability.csv"))
 after<-vapply(sources,sha256_file,character(1));assert_true(identical(before,after),"Review input changed during preparation")
 review_write(data.frame(file=sources,sha256=unname(after)),file.path(out,"review_input_hashes.csv"))
 yaml::write_yaml(list(schema_version=1,source_run=run_id,source_export=export_id,current_clr_attempt=basename(current),export_clr_attempt=export_binding$export_attempt,export_attempt_relation=export_binding$relation,baseline_clr_attempt=baseline_attempt,diagnostics=diagnostics,statistical_contract_sha256=sha256_file(file.path(root,"config/statistical_implementation_sha256.csv")),snapshot_summary_sha256=sha256_file(file.path(out,"CLR_result_summary.csv")),scientific_comparison_sha256=sha256_file(file.path(out,"CLR_old_new_scientific_comparison.csv"))),file.path(out,"review_context.yml"))
 safe_write_lines(c("AUTHOR REVIEW PREPARATION ONLY; no result accepted or frozen.","Inspect all16 CLR families; three have changed models and13 require preserved scientific fields. Fill explicit counts, ACCEPT, CONFIRMED, reviewer/date/rationale only after review.","REML classes are triage, not harmless/resolved classifications. Original warning attribution remains limited; manuscript review/evidence must be supplied.","freeze-clr requires a new external destination and complete signed CLR decisions. archive-reml records review without clearing the existing warning release gate."),file.path(out,"README_REVIEW.txt"))
 message("Review prepared; fill author CSVs after inspection: ",out)
}
review_freeze_clr <- function(root,review,destination) {
 ctx<-yaml::read_yaml(file.path(review,"review_context.yml"));assert_true(identical(ctx$statistical_contract_sha256,sha256_file(file.path(root,"config/statistical_implementation_sha256.csv"))),"Statistical implementation changed after review preparation")
 summary<-review_read(file.path(review,"CLR_result_summary.csv"));assert_true(identical(ctx$snapshot_summary_sha256,sha256_file(file.path(review,"CLR_result_summary.csv"))),"Prepared CLR summary changed")
 decisions<-review_read(file.path(review,"CLR_AUTHOR_DECISIONS.csv"));review_validate_decisions(decisions,summary)
 assert_true(identical(ctx$scientific_comparison_sha256,sha256_file(file.path(review,"CLR_old_new_scientific_comparison.csv"))),"Prepared scientific comparison changed");comparison<-review_read(file.path(review,"CLR_old_new_scientific_comparison.csv"));assert_true(all(comparison$scientific_fields_same[!comparison$changed_model])&&all(comparison$same_feature_order[!comparison$changed_model]),"One of13 unchanged CLR scientific families differs; cannot freeze")
 source_run<-review_run(root,ctx$source_run);contracts<-pub_clr_alignment_contract(source_run)
 for(f in list.files(file.path(review,"clr_snapshot"),full.names=TRUE))assert_true(identical(sha256_file(f),sha256_file(file.path(source_run$result_dirs[["clr"]],"tables/08_CLR_sensitivity",basename(f)))),"CLR snapshot/provenance differs from current Run")
 for(i in seq_len(nrow(summary))){f<-summary$result_file[i];assert_true(identical(sha256_file(file.path(review,"clr_snapshot",f)),summary$source_sha256[i])&&identical(sha256_file(file.path(source_run$result_dirs[["clr"]],"tables/08_CLR_sensitivity",f)),summary$source_sha256[i]),paste("Reviewed result source/snapshot changed:",f))}
 destination<-review_new_external(root,destination);dir.create(file.path(destination,"clr"));index<-list()
 snapshot<-list.files(file.path(review,"clr_snapshot"),pattern="[.]csv$",full.names=TRUE)
 for(f in snapshot){name<-basename(f);to<-file.path(destination,"clr",name);assert_true(file.copy(f,to,overwrite=FALSE),"Reference copy failed");assert_true(identical(sha256_file(f),sha256_file(to)),"Reference copy hash mismatch");index[[name]]<-data.frame(module="clr",relative_file=paste0("tables/08_CLR_sensitivity/",name),reference_file=paste0("clr/",name),sha256=sha256_file(to))}
 review_write(do.call(rbind,index),file.path(destination,"reference_index.csv"))
 assert_true(file.copy(file.path(review,"CLR_AUTHOR_DECISIONS.csv"),file.path(destination,"author_decisions.csv"),overwrite=FALSE),"Cannot archive author review")
 yaml::write_yaml(ctx,file.path(destination,"reference_provenance.yml"))
 # Freeze only these acceptance expectations; never change a statistical implementation/hash.
 plan<-yaml::read_yaml(file.path(root,"config/release_validation_plan.yml"));j<-which(vapply(plan$modules,function(z)z$module_id=="clr",logical(1)));assert_true(length(j)==1L,"Ambiguous CLR validation module")
 changed<-summary[summary$changed_model,,drop=FALSE]
 for(i in seq_len(nrow(changed))){rel<-paste0("tables/08_CLR_sensitivity/",changed$result_file[i]);hit<-which(vapply(plan$modules[[j]]$expected_output_files,function(z)identical(z$file,rel),logical(1)));assert_true(length(hit)==1L,"Changed CLR validation file missing");plan$modules[[j]]$expected_output_files[[hit]]$expected_significant<-as.integer(changed$n_significant[i]);plan$modules[[j]]$expected_output_files[[hit]]$acceptance_status<-"ACCEPTED_CLR_AUTHOR_REVIEW";plan$modules[[j]]$expected_significant[[rel]]<-as.integer(changed$n_significant[i])}
 cfg<-list(schema_version=1,status="ACCEPTED_CLR_AUTHOR_REVIEW",required_for_computation=FALSE,reference_index="reference_index.csv",reference_index_sha256=sha256_file(file.path(destination,"reference_index.csv")),review_decisions_sha256=sha256_file(file.path(destination,"author_decisions.csv")),statistical_contract_sha256=ctx$statistical_contract_sha256,source_run=ctx$source_run,source_export=ctx$source_export)
 backups<-file.path(destination,"activation_before");dir.create(backups)
 activation_files<-c("config/release_validation_plan.yml","config/clr_review_reference.yml","public_file_manifest.csv","public_file_sha256.csv")
 for(f in activation_files){to<-file.path(backups,f);dir.create(dirname(to),recursive=TRUE,showWarnings=FALSE);assert_true(file.copy(file.path(root,f),to,overwrite=FALSE),"Cannot back up validation activation")}
 tryCatch({
   yaml::write_yaml(plan,file.path(root,"config/release_validation_plan.yml"));yaml::write_yaml(cfg,file.path(root,"config/clr_review_reference.yml"))
   source(file.path(root,"R/functions/public_release_files.R"),local=environment());public_release_manifests(root)
 },error=function(e){for(f in activation_files)file.copy(file.path(backups,f),file.path(root,f),overwrite=TRUE);stop("CLR reference activation failed; configuration/inventory restored where possible: ",conditionMessage(e),call.=FALSE)})
 safe_write_lines(c("Reviewed CLR baseline for future independent runs; freezing is NOT independent reproduction.","Only the three changed CLR expected counts were activated after signed author decisions;13 unchanged scientific families were checked. No production admission flag changed.","Reference files must never be read by analysis/export. Use a separate validation process on a NEW Run; source Run self-comparison is rejected."),file.path(destination,"README_REFERENCE.txt"))
 message("Reviewed CLR reference frozen externally: ",destination,"; all other release gates unchanged.")
}
review_archive_reml <- function(root,review,destination) {
 d<-review_read(file.path(review,"REML_AUTHOR_DECISIONS.csv"));w<-review_read(file.path(review,"REML_warning_classification.csv"))
 assert_true(!anyDuplicated(d$review_key)&&setequal(d$review_key,unique(w$review_key)),"REML review coverage differs")
 assert_true(all(d$decision%in%c("REVIEWED_WITH_LIMITATIONS","REQUIRES_FURTHER_DIAGNOSIS","REJECT_RELEASE")),"Each warning feature needs an explicit review decision; no automatic harmless class")
 for(n in c("reviewer","review_date","rationale"))assert_true(all(!is.na(d[[n]])&nzchar(trimws(d[[n]]))),paste("Missing REML review field:",n))
 destination<-review_new_external(root,destination)
 for(f in c("REML_AUTHOR_DECISIONS.csv","REML_warning_classification.csv","REML_full_replay_stability.csv","review_context.yml","review_input_hashes.csv"))assert_true(file.copy(file.path(review,f),file.path(destination,f),overwrite=FALSE),"REML archive copy failed")
 evidence<-list()
 for(i in seq_len(nrow(d)))if(!is.na(d$evidence_file[i])&&nzchar(d$evidence_file[i])){f<-normalizePath(d$evidence_file[i],winslash="/",mustWork=TRUE);to<-file.path(destination,paste0("evidence_",i,"_",basename(f)));assert_true(file.copy(f,to,overwrite=FALSE),"Evidence copy failed");evidence[[length(evidence)+1L]]<-data.frame(review_key=d$review_key[i],archived_file=basename(to),sha256=sha256_file(to))}
 if(length(evidence))review_write(do.call(rbind,evidence),file.path(destination,"evidence_index.csv"))
 safe_write_lines(c("REML REVIEW ARCHIVE; RELEASE WARNING GATE NOT CLEARED.","Classification/rationale/evidence are archived, not automatic proof that Fisher scoring reached a global optimum. No formal estimate replaced; no admission flag changed."),file.path(destination,"README_REML_REVIEW.txt"))
 message("REML classification review archived: ",destination,". Existing warning gate remains unchanged.")
}
review_check_clr_reference <- function(root,run,reference,out) {
 cfg<-yaml::read_yaml(file.path(root,"config/clr_review_reference.yml"));assert_true(identical(cfg$status,"ACCEPTED_CLR_AUTHOR_REVIEW"),"CLR reference has not been author-reviewed/frozen")
 assert_true(identical(cfg$statistical_contract_sha256,sha256_file(file.path(root,"config/statistical_implementation_sha256.csv"))),"CLR accepted reference was frozen under different statistical code")
 provenance<-yaml::read_yaml(file.path(reference,"reference_provenance.yml"));assert_true(identical(provenance$source_run,cfg$source_run)&&identical(provenance$statistical_contract_sha256,cfg$statistical_contract_sha256),"Reference provenance differs from accepted configuration");assert_true(!identical(provenance$source_run,run$id),"Source Run cannot certify its own independent CLR reproduction")
 assert_true(identical(sha256_file(file.path(reference,cfg$reference_index)),cfg$reference_index_sha256),"CLR reference index hash differs")
 assert_true(identical(sha256_file(file.path(reference,"author_decisions.csv")),cfg$review_decisions_sha256),"CLR author review archive hash differs")
 index<-review_read(file.path(reference,cfg$reference_index));assert_true(!anyDuplicated(index$relative_file),"Duplicate CLR reference entry")
 plan<-read_release_validation_plan(root);j<-Filter(function(x)x$module_id=="clr",plan$modules)[[1]]
 assert_true(all(vapply(j$expected_output_files,function(x)x$file,character(1))%in%index$relative_file),"CLR reference missing required result/provenance files")
 for(i in seq_len(nrow(index))){relative_project_path(index$relative_file[i]);relative_project_path(index$reference_file[i]);current<-file.path(run$result_dirs[["clr"]],index$relative_file[i]);frozen<-file.path(reference,index$reference_file[i]);assert_true(identical(sha256_file(frozen),index$sha256[i]),"CLR reference artifact hash differs");release_compare_frames(review_read(current),review_read(frozen),paste0("aligned_CLR_reference_",i),out,plan$absolute_tolerance,0)}
 invisible(TRUE)
}

review_workbook_evidence <- function(export) {
  file<-file.path(export,"validation/workbook_ooxml_integrity.csv")
  x<-review_read(file)
  required<-submission_table_manifest()$table_id
  assert_columns(x,c("workbook","final_pass","worksheet_range_pass","streaming_readback_pass"),"Workbook XML/streaming acceptance")
  assert_true(nrow(x)==length(required)&&!anyDuplicated(x$workbook)&&setequal(x$workbook,required)&&
    all(!is.na(x$final_pass)&x$final_pass)&all(!is.na(x$worksheet_range_pass)&x$worksheet_range_pass)&
    all(!is.na(x$streaming_readback_pass)&x$streaming_readback_pass),"New review requires all16 workbooks with XML ranges and independent streaming checks passed")
  evidence<-file
  for(id in required) {
    ranges<-file.path(export,"validation/workbooks",paste0(id,"__ranges.csv"))
    stream<-file.path(export,"validation/workbooks",paste0(id,"__streaming.json"))
    r<-review_read(ranges)
    assert_columns(r,c("range_pass","streaming_pass","expected_rows","expected_columns","read_rows","read_columns"),paste(id,"range receipt"))
    assert_true(nrow(r)>0L&&all(r$range_pass&r$streaming_pass)&&all(r$read_rows==r$expected_rows)&all(r$read_columns==r$expected_columns),paste("Incomplete sheet streaming evidence:",id))
    s<-jsonlite::read_json(stream,simplifyVector=FALSE)
    workbook<-file.path(export,"tables",if(id%in%c("Main_Table_1","Analysis_Summary_Table"))"main" else "supplementary",paste0(id,".xlsx"))
    assert_true(isTRUE(s$pass_status)&&isTRUE(s$read_only)&&identical(s$reader,"openpyxl")&&
      identical(s$workbook_sha256,sha256_file(workbook)),paste("Invalid or stale independent reader receipt:",id))
    assert_true(length(s$sheets)==nrow(r)&&identical(vapply(s$sheets,`[[`,character(1),"sheet"),as.character(r$sheet))&&
      all(vapply(s$sheets,function(z)isTRUE(z$pass_status)&&isTRUE(z$cell_type_and_value_pass),logical(1))),
      paste("Independent reader sheet coverage differs:",id))
    evidence<-c(evidence,ranges,stream)
  }
  evidence
}
