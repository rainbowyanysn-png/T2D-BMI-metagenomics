#!/usr/bin/env Rscript
# Static dependency and clean-transfer inventory only. Never invokes analysis or export.
local({
 if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
 entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));root<-normalizePath(dirname(entry),winslash="/",mustWork=TRUE)
 source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
 source(file.path(root,"R/functions/assertions.R"),local=environment());source(file.path(root,"R/functions/public_release_files.R"),local=environment())
 args<-commandArgs(TRUE);assert_true(length(args)%in%c(0L,2L),"Use audit_rc4_dependencies.R [--out NEW_EXTERNAL_DIRECTORY]")
 if(length(args))assert_true(args[1]=="--out","Expected --out")
 out<-if(length(args))args[2]else file.path(dirname(root),paste0(basename(root),"_Audit"),paste0("DependencyAudit_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid()))
 out<-normalizePath(out,winslash="/",mustWork=FALSE);canonical<-tolower(normalizePath(root,winslash="/",mustWork=TRUE))
 assert_true(!dir.exists(out)&&!file.exists(out)&&!startsWith(tolower(out),paste0(canonical,"/"))&&!identical(tolower(out),canonical),"Audit destination must be new and outside public package")
 dir.create(out,recursive=TRUE);payload<-public_release_selection(root)
 utils::write.csv(payload,file.path(out,"clean_test_payload_files.csv"),row.names=FALSE)
 utils::write.csv(data.frame(path=c("public_file_manifest.csv","public_file_sha256.csv","output/.gitkeep","data_derived/.gitkeep"),role=c("payload inventory metadata","payload hash metadata","empty output directory marker","empty derived directory marker")),file.path(out,"clean_test_additional_metadata.csv"),row.names=FALSE)
 entries<-c("run_all.R","run_publication_candidate.R","run_all_development.R","export_publication.R","validate_release.R")
 registry<-utils::read.csv(file.path(root,"config/pipeline_modules.csv"),stringsAsFactors=FALSE)
 assert_true(all(entries%in%payload$file)&&nrow(registry)==22L,"Missing formal entry or registered modules")
 # Conservative source inventory: retain all internal R code, plus validation support.
 files<-payload$file[grepl("[.]R$",payload$file)]
 issues<-calls<-list();old<-"FinalIntegration|FigureRefinement|IntegrationPreview|CleanTest|LOCO_Trial|MedicationRestored|1015返修|1115返修|T2D_BMI_RC2"
 walk<-function(expr,file){
  if(missing(expr))return(invisible(NULL))
  if(is.call(expr)){
   op<-paste(deparse(expr[[1]]),collapse="")
   if(op%in%c("source","sys.source","readRDS","load","setwd","install.packages"))calls[[length(calls)+1L]]<<-data.frame(file=file,call=op,expression=paste(deparse(expr,width.cutoff=200L),collapse=" "),scope=if(startsWith(file,"tests/"))"validation_or_fixture"else"runtime_or_support")
   # Literal personal absolute paths are different from path-recognition regexes.
   for(i in seq_along(expr)[-1])walk(expr[[i]],file)
  }else if(is.expression(expr)||is.pairlist(expr)){lapply(as.list(expr),function(e){if(!missing(e))walk(e,file)})}
  else if(is.character(expr))for(value in expr){
   if(grepl("^[A-Za-z]:[/\\\\]|^/Users/|^/home/[^/]+/",value))issues[[length(issues)+1L]]<<-data.frame(file=file,issue="LITERAL_PERSONAL_ABSOLUTE_PATH",value=value)
   if(grepl(old,value,ignore.case=TRUE)&&!startsWith(file,"tests/")&&file!="audit_rc4_dependencies.R")issues[[length(issues)+1L]]<<-data.frame(file=file,issue="LEGACY_PROJECT_NAME_LITERAL_REVIEW",value=value)
  }
 }
 for(f in files)walk(parse(file.path(root,f),encoding="UTF-8",keep.source=FALSE),f)
 calls<-if(length(calls))do.call(rbind,calls)else data.frame(file=character(),call=character(),expression=character(),scope=character())
 issues<-if(length(issues))do.call(rbind,issues)else data.frame(file=character(),issue=character(),value=character())
 utils::write.csv(calls,file.path(out,"source_io_call_inventory.csv"),row.names=FALSE);utils::write.csv(issues,file.path(out,"literal_path_findings.csv"),row.names=FALSE)
 # Existing strict analysis-code scan complements the AST literal inventory.
 source(file.path(root,"tests/validate_no_legacy_dependency.R"),local=environment());legacy<-tryCatch(validate_no_legacy_dependency(root),error=function(e)e)
 legacy_ok<-!inherits(legacy,"error");if(legacy_ok)utils::write.csv(legacy,file.path(out,"analysis_legacy_dependency_checks.csv"),row.names=FALSE)
 model_files<-unique(c("R/pipeline/01_import_harmonize_qc.R",paste0("R/pipeline/",registry$file)))
 assert_true(all(model_files%in%payload$file),"A registered pipeline dependency was excluded")
 for(boot in c("R/functions/analysis_bootstrap.R","R/functions/publication_bootstrap.R"))assert_file(file.path(root,boot))
 # Names explicitly loaded by current fixed bootstraps are checked without evaluating them.
 bootdeps<-list()
 for(boot in c("R/functions/analysis_bootstrap.R","R/functions/publication_bootstrap.R")){
  expressions<-parse(file.path(root,boot),encoding="UTF-8");strings<-character()
  collect<-function(x){if(missing(x))return(invisible(NULL));if(is.character(x))strings<<-c(strings,x)else if(is.call(x)||is.expression(x)||is.pairlist(x))lapply(as.list(x),function(y){if(!missing(y)&&!is.symbol(y))collect(y)})};collect(expressions)
  lists<-list()
  find_lists<-function(x){if(missing(x))return();if(is.call(x)){if(identical(x[[1]],as.name("c"))&&all(vapply(as.list(x)[-1],is.character,logical(1)))){v<-unlist(as.list(x)[-1]);if("assertions"%in%v)lists[[length(lists)+1L]]<<-v};lapply(as.list(x),function(y){if(!missing(y)&&!is.symbol(y))find_lists(y)})}else if(is.expression(x)||is.pairlist(x))lapply(as.list(x),function(y){if(!missing(y)&&!is.symbol(y))find_lists(y)})};find_lists(expressions)
  assert_true(length(lists)==1L,paste("Unresolved explicit bootstrap dependency vector:",boot))
  names<-unique(lists[[1]])
  missing_files<-setdiff(paste0("R/functions/",names,".R"),payload$file)
  assert_true(!length(missing_files),paste("Bootstrap dependency missing:",paste(missing_files,collapse=", ")))
  selected<-intersect(paste0("R/functions/",names,".R"),payload$file)
  bootdeps[[boot]]<-data.frame(loader=boot,dependency=selected,present=TRUE)
 }
 utils::write.csv(do.call(rbind,bootdeps),file.path(out,"bootstrap_dependency_inventory.csv"),row.names=FALSE)
 dependency<-data.frame(module=c("inputs",registry$module),script=c("R/pipeline/01_import_harmonize_qc.R",paste0("R/pipeline/",registry$file)),included=TRUE)
 utils::write.csv(dependency,file.path(out,"registered_module_dependencies.csv"),row.names=FALSE)
 source(file.path(root,"tests/release_access_guard.R"),local=environment());stopifnot(identical(rc3_guard_filesystem_paths(c("","https://experimenthub.bioconductor.org")),character()))
 saved<-utils::read.csv(file.path(root,"public_file_manifest.csv"),stringsAsFactors=FALSE)
 manifest_ok<-identical(payload$file,saved$file)&&identical(payload$sha256,saved$sha256)
 forbidden<-grepl("(^|/)(audit|output|data_derived|renv|reference_results)/|(^|/)(preview|replace_|resize_)",payload$file)
 statistical<-utils::read.csv(file.path(root,"config/statistical_implementation_sha256.csv"));stat_ok<-all(vapply(file.path(root,statistical$file),function(f)digest::digest(file=f,algo="sha256"),character(1))==statistical$sha256)
 checks<-data.frame(check=c("formal entries/modules included","all R syntax parsed","no literal personal paths or legacy names outside validation fixtures","analysis legacy dependency scan","manifest matches current code","development/results excluded","statistical hashes unchanged","network/anonymous connection classifier"),pass=c(TRUE,TRUE,nrow(issues)==0L,legacy_ok,manifest_ok,!any(forbidden),stat_ok,TRUE))
 utils::write.csv(checks,file.path(out,"dependency_checks.csv"),row.names=FALSE)
 writeLines(c("STATIC DEPENDENCY AUDIT; NOT AN INDEPENDENT CLEAN-RUN CERTIFICATE.","Source/I/O inventory includes dynamic bootstrap/module source calls. Fixed bootstraps and the complete registered module set are inventoried; static analysis is conservative and cannot prove arbitrary dynamic execution.","No models, raw-data download, image export or workbook export ran. Network and anonymous connections were classified, not contacted.","Clean-test payload list describes code to transfer; it does not copy old results. Explicit validation-only reference access and diagnostic/review operations are separate from formal calculation.","Next steps: author review/freeze, prepare_transfer.R --zip, fresh-computer empty output/cache setup, renv restore, static_rc4_candidate, full run_all.R --candidate, export validation and independent reference comparison. Admission flags unchanged."),file.path(out,"README_DEPENDENCY_AUDIT.txt"))
 writeLines(c("# R commands after setting project_dir on the new computer:","renv::restore(project=project_dir, prompt=FALSE)","rc4(\"tests/static_rc4_candidate.R\")","rc4(\"run_all.R\", \"--candidate\")","# On interruption: rc4(\"run_all_development.R\", c(\"--resume\", NEW_RUN_ID))","# Re-export only: rc4(\"export_publication.R\", c(NEW_RUN_ID, \"all\"))","# Saved validation: rc4(\"validate_release.R\", c(\"--run\",NEW_RUN_ID,\"--export\",NEW_EXPORT_ID))","# New-run CLR reference only: append c(\"--clr-reference-root\",EXTERNAL_ACCEPTED_CLR_REFERENCE)."),file.path(out,"CLEAN_TEST_COMMANDS.R"))
 message("Dependency/clean-test inventory: ",out)
 if(!all(checks$pass))quit(status=1L)
})
