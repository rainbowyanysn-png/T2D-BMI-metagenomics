#!/usr/bin/env Rscript
if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
local({
 entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));stopifnot(length(entry)==1L)
 root<-normalizePath(dirname(entry),winslash="/",mustWork=TRUE)
 source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
 source(file.path(root,"R/functions/publication_bootstrap.R"),local=environment());load_publication_support(root,environment())
 source(file.path(root,"tests/publication_output_review_helpers.R"),local=environment())
 source(file.path(root,"tests/release_validation_config.R"),local=environment())
 source(file.path(root,"tests/release_validation_methods.R"),local=environment())
 source(file.path(root,"tests/validate_final_integration_saved.R"),local=environment())
 args<-commandArgs(TRUE);opts<-list();i<-1L
 while(i<=length(args)){
  assert_true(args[i]%in%c("--run","--export","--reference-root","--clr-reference-root")&&i<length(args),"Use --run RUN --export EXPORT [--reference-root DIRECTORY] [--clr-reference-root DIRECTORY]")
  key<-sub("^--","",args[i]);assert_true(is.null(opts[[key]]),"Duplicate argument");opts[[key]]<-args[i+1L];i<-i+2L
 }
 assert_true(!is.null(opts$run)&&grepl("^Run_[0-9]{8}_[0-9]{6}$",opts$run),"Invalid new Run ID")
 assert_true(!is.null(opts$export)&&grepl("^PublicationExport_[0-9]{8}_[0-9]{6}([.][0-9]+)?_[0-9]+$",opts$export),"Invalid export ID")
 run<-list(root=root,id=opts$run,dir=file.path(root,"output",opts$run));run$result_dirs<-completed_result_directories(run)
 run$access<-new.env();run$access$files<-character();run$access$hashes<-character()
 out<-file.path(root,"output",paste0("ReleaseValidation_",invocation_stamp()));dir.create(out)
 rows<-new.env();rows$data<-list();plan<-read_release_validation_plan(root)
 release_internal_checks(root,run,plan,out,rows,export=file.path(root,"output",opts$export))
 config<-release_publication_config(root)
 release_check("publication numeric and model-family contracts",function(){b<-publication_data(run,config);publication_numeric_checks(b,out)},rows)
 Sys.setenv(T2DBMI_VALIDATION_EXPORT=opts$export)
 release_check("16 workbooks and manifest-defined panels",function(){validate_publication_output_directories(root,opts$run,out,config)},rows)
 release_check("final S3/S4 and forest integration contracts",function(){validate_final_integration_saved(root,run,file.path(root,"output",opts$export))},rows)
 internal<-do.call(rbind,rows$data);internal_pass<-all(internal$pass)
 saved_pass<-all(internal$pass[internal$check!="statistical warning classification"])
 warning_pass<-all(internal$pass[internal$check=="statistical warning classification"])
 ref<-opts[["reference-root"]];if(is.null(ref))ref<-Sys.getenv("T2DBMI_REFERENCE_ROOT")
 ref_start<-length(rows$data)
 if(nzchar(ref))release_check("external reference process",function()release_reference_checks(root,run,ref,plan,out,rows),rows)
 clr_ref<-opts[["clr-reference-root"]];if(is.null(clr_ref))clr_ref<-Sys.getenv("T2DBMI_CLR_REFERENCE_ROOT")
 clr_start<-length(rows$data)
 if(nzchar(clr_ref))release_check("independent aligned CLR reference comparison",function(){source(file.path(root,"tests/rc4_release_review_methods.R"),local=environment());review_check_clr_reference(root,run,clr_ref,out)},rows)
 results<-do.call(rbind,rows$data);safe_write_csv(results,file.path(out,"release_checks.csv"))
 external<-if(!nzchar(ref))"NOT RUN"else if(all(results$pass[seq.int(ref_start+1L,clr_start)]))"PASS"else"FAIL"
 clr_external<-if(!nzchar(clr_ref))"NOT RUN"else if(all(results$pass[seq.int(clr_start+1L,nrow(results))]))"PASS (NEW-RUN CLR ONLY)"else"FAIL"
 safe_write_csv(data.frame(scope=c("saved-result internal validation","statistical warning release gate","independent reference comparison","independent aligned CLR reference comparison","manual Excel and visual review","new-computer no-cache test"),
  status=c(if(saved_pass)"PASS"else"FAIL",if(warning_pass)"PASS"else"UNRESOLVED; RELEASE BLOCKED",external,clr_external,"PENDING AUTHOR REVIEW","NOT EXECUTED BY VALIDATOR")),file.path(out,"validation_scope.csv"))
 safe_write_lines(c(paste("saved-result internal validation:",if(saved_pass)"PASS"else"FAIL"),paste("statistical warning gate:",if(warning_pass)"PASS"else"UNRESOLVED"),paste("internal validation:",if(internal_pass)"PASS"else"FAIL"),paste("independent reference comparison:",external),
  "No model was fitted by this validator. No statistic, threshold or submitted artifact was replaced."),file.path(out,"status.txt"))
 safe_write_lines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
 message("Release validation report: ",out)
 quit(status=if(all(results$pass))0L else 1L)
})
