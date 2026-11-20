#!/usr/bin/env Rscript
# Saved submission output only. No model fits, result rewriting or release admission.
local({
 if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
 entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));root<-normalizePath(dirname(entry),winslash="/",mustWork=TRUE)
 source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
 source(file.path(root,"R/functions/publication_bootstrap.R"),local=environment());load_publication_support(root,environment())
 source(file.path(root,"tests/publication_output_review_helpers.R"),local=environment())
 args<-commandArgs(TRUE);assert_true(length(args)%in%c(2L,3L),"Use RUN_ID TABLE_EXPORT_ID [FIGURE_EXPORT_ID]")
 assert_true(grepl("^Run_[0-9]{8}_[0-9]{6}$",args[1])&&all(grepl("^PublicationExport_[0-9]{8}_[0-9]{6}[.][0-9]+_[0-9]+$",args[-1])),"Invalid Run/export identity")
 run<-list(root=root,id=args[1],dir=file.path(root,"output",args[1]));run$result_dirs<-completed_result_directories(run)
 run$access<-new.env(parent=emptyenv());run$access$files<-character();run$access$hashes<-character()
 export<-file.path(root,"output",args[2]);figure_export<-if(length(args)==3L)file.path(root,"output",args[3]) else export
 for(folder in unique(c(export,figure_export))){c<-read.csv(file.path(folder,"manifests/execution_context.csv"));assert_true(nrow(c)==1L&&c$run_id==run$id,"Export is bound to another Run")}
 out<-file.path(root,"output",paste0("SubmissionValidation_",invocation_stamp()));dir.create(out)
 message("Submission validation: ",out)
 error<-NULL
 tryCatch({
  result<-submission_validate_saved(root,run,export,figure_export);safe_write_csv(result,file.path(out,"S7_panel_source_sets.csv"))
  tables<-read.csv(file.path(export,"table_index.csv"))
  for(i in seq_len(nrow(tables))){
   id<-tables$table_id[i];dictionary<-read.csv(file.path(export,"tables/csv",paste0(id,"__data_dictionary.csv")))
   sheets<-setNames(lapply(unique(dictionary$sheet),function(s){
    z<-dictionary[dictionary$sheet==s,];classes<-ifelse(z$storage_class%in%c("numeric","integer","logical","character"),z$storage_class,"character")
    publication_read_export_csv(file.path(export,"tables/csv",paste0(id,"__",s,".csv")),setNames(classes,z$column))
   }),unique(dictionary$sheet))
   target<-file.path(export,tables$workbook_path[i]);assert_true(sha256_file(target)==tables$workbook_sha256[i],"Workbook index hash differs")
   pub_validate_saved_workbook(target,sheets,id,out,root)
  }
  pub_finalize_workbook_integrity(out)
  figures<-read.csv(file.path(figure_export,"figure_index.csv"));assert_true(nrow(figures)==40L&&!anyDuplicated(figures$figure_id),"Expected40 unchanged panels")
  manifest<-read.csv(file.path(figure_export,"output_manifest.csv"),check.names=FALSE)
  for(field in c("PDF_path","TIFF_path","PNG_path","plot_data_path","caption_path"))for(path in figures[[field]]) {
    k<-match(path,manifest$file);assert_true(!is.na(k)&&identical(sha256_file(file.path(figure_export,path)),manifest$sha256[k]),paste("Figure artifact identity differs:",path))
  }
  safe_write_csv(data.frame(role=c("tables","figures"),export=c(args[2],basename(figure_export)),run_id=run$id),file.path(out,"export_source_mapping.csv"))
 },error=function(e)error<<-conditionMessage(e))
 safe_write_lines(c(if(is.null(error))"SUBMISSION_OUTPUT_PASS" else "SUBMISSION_OUTPUT_FAIL",if(!is.null(error))error,
   "Only submission output was checked. REML warning gate, reference validation and author release admission remain separate."),file.path(out,"status.txt"))
 if(!is.null(error))stop(error,call.=FALSE)
 message("PASS16 workbooks,40 panel registry, S7 source identities and streaming read-back. Report: ",out)
})
