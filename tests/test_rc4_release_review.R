#!/usr/bin/env Rscript
# Pure review-decision tests; no actual result is accepted/frozen and no model is run.
local({
 if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
 entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE));root<-normalizePath(dirname(dirname(entry)),winslash="/",mustWork=TRUE)
 source(file.path(root,"R/functions/assertions.R"),local=environment());source(file.path(root,"tests/rc4_release_review_methods.R"),local=environment())
 summary<-data.frame(result_file=paste0("fixture",1:16,".csv"),source_sha256=rep(paste(rep("a",64),collapse=""),16),n_significant=0:15)
 d<-summary[c("result_file","source_sha256")];d$decision<-"ACCEPT";d$accepted_significant_count<-0:15;d$reviewer<-"TEST_FIXTURE";d$review_date<-"2026-11-19";d$rationale<-"TEST_ONLY_NOT_REAL_APPROVAL";d$visual_review<-"CONFIRMED"
 review_validate_decisions(d,summary)
 reject<-function(z){e<-tryCatch({review_validate_decisions(z,summary);NULL},error=function(e)e);stopifnot(inherits(e,"error"))}
 z<-d;z$decision[1]<-"PENDING";reject(z);z<-d;z$accepted_significant_count[1]<-99;reject(z);z<-d;z$reviewer[1]<-"";reject(z);z<-d;z$source_sha256[1]<-"wrong";reject(z);reject(d[-1,]);z<-d;z$result_file[2]<-z$result_file[1];reject(z)
 stopifnot(sum(review_changed_key(c("Taxonomy|Overweight_vs_Lean_T2D","Taxonomy|BMI23_Formal_Interaction","Pathway|BMI23_Formal_Interaction","Pathway|Overweight_vs_Lean_T2D")))==3L)
 # New review receipts: complete rectangular streaming evidence is mandatory.
 source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
 source(file.path(root,"R/functions/source_resources.R"),local=environment())
 tmp<-tempfile("review_workbook_receipts_");dir.create(file.path(tmp,"validation/workbooks"),recursive=TRUE);on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
 source(file.path(root,"R/functions/submission_helpers.R"),local=environment())
 ids<-submission_table_manifest()$table_id
 integrity<-data.frame(workbook=ids,final_pass=TRUE,worksheet_range_pass=TRUE,streaming_readback_pass=TRUE)
 utils::write.csv(integrity,file.path(tmp,"validation/workbook_ooxml_integrity.csv"),row.names=FALSE)
 for(id in ids){
   folder<-file.path(tmp,"tables",if(id%in%ids[1:2])"main" else "supplementary");dir.create(folder,recursive=TRUE,showWarnings=FALSE)
   book<-file.path(folder,paste0(id,".xlsx"));writeLines("RECEIPT_FIXTURE_ONLY",book)
   utils::write.csv(data.frame(sheet="Fixture",range_pass=TRUE,streaming_pass=TRUE,expected_rows=2L,expected_columns=2L,read_rows=2L,read_columns=2L),file.path(tmp,"validation/workbooks",paste0(id,"__ranges.csv")),row.names=FALSE)
   jsonlite::write_json(list(pass_status=TRUE,read_only=TRUE,reader="openpyxl",workbook_sha256=sha256_file(book),sheets=list(list(sheet="Fixture",pass_status=TRUE,cell_type_and_value_pass=TRUE))),file.path(tmp,"validation/workbooks",paste0(id,"__streaming.json")),auto_unbox=TRUE)
 }
 stopifnot(length(review_workbook_evidence(tmp))==33L)
 bad_receipt<-file.path(tmp,"validation/workbooks/Main_Table_1__ranges.csv")
 utils::write.csv(data.frame(sheet="Fixture",range_pass=TRUE,streaming_pass=TRUE,expected_rows=2L,expected_columns=2L,read_rows=1L,read_columns=2L),bad_receipt,row.names=FALSE)
 e<-tryCatch({review_workbook_evidence(tmp);NULL},error=function(e)e);stopifnot(inherits(e,"error"))
 cat("PASS complete streaming receipts required; incomplete read rejected\n")
 files<-c("review_rc4_release.R","audit_rc4_dependencies.R","tests/rc4_release_review_methods.R","validate_release.R","tests/release_validation_methods.R")
 for(f in files)parse(file.path(root,f),encoding="UTF-8")
 cat("PASS fixture-only author acceptance checks; pending/missing/hash/count/coverage defects rejected\nPASS exact three changed-family classification\nPASS review/validation/audit R syntax\nNo real acceptance, reference freeze, diagnostic fit or dependency workflow executed.\n")
})
