#!/usr/bin/env Rscript
# Temporary fixtures only. Does not export a Run, fit a model or accept results.
local({
  if(.Platform$OS.type=="windows"&&!isTRUE(l10n_info()[["UTF-8"]]))Sys.setlocale("LC_CTYPE","English_United States.utf8")
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
  root<-normalizePath(dirname(dirname(entry)),winslash="/",mustWork=TRUE)
  source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
  source(file.path(root,"R/functions/publication_bootstrap.R"),local=environment());load_publication_support(root,environment())
  base_tmp<-tempfile("workbook_contract_");dir.create(base_tmp);on.exit(unlink(base_tmp,recursive=TRUE),add=TRUE)
  tmp<-file.path(base_tmp,"工作簿 fixture");dir.create(tmp)
  reject<-function(action,pattern=NULL){e<-tryCatch({action();NULL},error=function(e)e);stopifnot(inherits(e,"error"));if(!is.null(pattern))stopifnot(grepl(pattern,conditionMessage(e),fixed=TRUE))}
  values<-data.frame(feature=c("00123","species: Escherichia coli",NA,"NA"),
    coef=c(-.123456789123456,0,NA,2.5),pval=c(1e-120,.049999,.050001,NA),
    logical_flag=c(TRUE,FALSE,NA,TRUE),note=c("first\nsecond","",NA,"literal"),check.names=FALSE)
  header<-values[FALSE,,drop=FALSE]
  sheets<-list(Normal=values,Header_only=header)
  row<-pub_workbook("Temporary_fixture","Workbook fixture",sheets,"Temporary test data",tmp,list(project_root=root))
  file<-file.path(tmp,row$workbook_path)
  range<-utils::read.csv(file.path(tmp,"validation/workbooks/Temporary_fixture__ranges.csv"))
  stopifnot(all(range$range_pass),range$declared_range[range$sheet=="Normal"]=="A1:E5",
    range$declared_range[range$sheet=="Header_only"]=="A1:E1")
  stream<-jsonlite::read_json(file.path(tmp,"validation/workbooks/Temporary_fixture__streaming.json"))
  stopifnot(isTRUE(stream$pass_status),isTRUE(stream$read_only))
  # Bad declaration is created in a TEMPORARY fixture before save. No old XLSX is patched.
  w<-openxlsx::createWorkbook();openxlsx::addWorksheet(w,"Normal");openxlsx::writeData(w,"Normal",values,keepNA=TRUE,na.string="NA")
  w$worksheets[[1]]$dimension <- '<dimension ref="A1"/>'
  pub_clear_empty_drawing_placeholders(w)
  bad<-file.path(tmp,"negative_A1.xlsx");openxlsx::saveWorkbook(w,bad)
  bad_range<-pub_sheet_dimension_report(bad,list(Normal=values));stopifnot(!bad_range$range_pass)
  reject(function()pub_stream_read_workbook(bad,list(Normal=values),"negative_A1",tmp,root),"streaming/type/value validation failed")
  reject(function()pub_validate_saved_workbook(bad,list(Normal=values),"negative_A1_saved",tmp,root),"Worksheet declared/actual/expected range differs")
  # Display formats must not change the numeric Excel values or their cell types.
  numeric<-data.frame(direction_agreement=c(.968,NA),descriptive_coefficient_r=c(.968123456789,NA),pval=c(1e-90,NA))
  out<-pub_workbook("Supplementary_Table_S13","Temporary S13",list(Concordance=numeric),"Temporary",tmp,list(project_root=root))
  z<-openxlsx::read.xlsx(file.path(tmp,out$workbook_path),sheet="Concordance",skipEmptyRows=FALSE,na.strings=NULL)
  stopifnot(isTRUE(all.equal(suppressWarnings(as.numeric(z$direction_agreement)),numeric$direction_agreement,tolerance=1e-12)),
    isTRUE(all.equal(suppressWarnings(as.numeric(z$descriptive_coefficient_r)),numeric$descriptive_coefficient_r,tolerance=1e-12)))
  # The independent streaming checker above additionally requires numeric cell types;
  # openxlsx itself may infer a character column when a missing marker is present.
  styles<-pub_xml_part(file.path(tmp,out$workbook_path),"xl/styles.xml")
  formats<-xml2::xml_attr(xml2::xml_find_all(styles,"//*[local-name()='numFmt']"),"formatCode")
  stopifnot(all(c("0.0%","0.000","0.000000E+00")%in%formats))
  # Source identity annotations preserve original values, including after row shuffling.
  lr<-data.frame(npar=c(8L,5L,5L,8L),Chisq=c(2,NA,NA,3),metric=c("Shannon","Shannon","Richness","Richness"),check.names=FALSE)
  lr[["Pr(>Chisq)"]]<-c(.2,NA,NA,.1)
  labeled<-pub_bmi23_alpha_identity(lr,"model_comparison","BMI23_Alpha_Group_Overall_LRT.csv")
  stopifnot(identical(labeled[,names(lr)],lr),all(labeled$model_id[lr$npar==5]=="BMI23_alpha_null"))
  wrong<-lr;wrong$npar[1]<-7L;reject(function()pub_bmi23_alpha_identity(wrong,"model_comparison","fixture"))
  pair<-data.frame(metric=rep(c("Shannon","Richness"),each=6),contrast=rep(paste0("group",1:6," - reference"),2),p.value=seq_len(12)/20)
  annotated<-pub_bmi23_alpha_identity(pair,"pairwise_comparison","BMI23_Alpha_Pairwise_BH.csv")
  stopifnot(identical(annotated[,names(pair)],pair),all(grepl("BH-adjusted",annotated$p_adjustment,fixed=TRUE)),!"p.adj"%in%names(annotated),
    "model_id"%in%names(pub_reader_table(annotated,preserve_model_identity=TRUE)),grepl("Tukey",pub_column_definition("p.value","Supplementary_Table_S2","Pairwise"),fixed=TRUE))
  # Actually empty worksheets have no populated cells, unlike header-only tables.
  empty<-openxlsx::createWorkbook();openxlsx::addWorksheet(empty,"Empty")
  pub_set_written_dimensions(empty,list(Empty=data.frame()));pub_clear_empty_drawing_placeholders(empty)
  empty_file<-file.path(tmp,"true_empty.xlsx");openxlsx::saveWorkbook(empty,empty_file)
  stopifnot(pub_sheet_dimension_report(empty_file,list(Empty=data.frame()))$range_pass)
  large_shape<-as.data.frame(matrix(0,nrow=47701L,ncol=18L))
  wide_header<-as.data.frame(matrix(numeric(),nrow=0L,ncol=28L))
  stopifnot(pub_written_range(large_shape)=="A1:R47702",pub_written_range(wide_header)=="A1:AB1")
  # The formal p.adj threshold is exact, even inside the numeric read-back tolerance.
  crossed<-data.frame(p.adj=.05-1e-14);observed<-data.frame(p.adj=.05+1e-14)
  reject(function()pub_verify_workbook_sheet(crossed,observed,"threshold","Pairwise",tmp),"Excel read-back mismatch")
  cat("PASS dynamic S13-sized and beyond-Z/header-only ranges; exact p.adj state checks\n")
  cat("PASS rectangular/header-only range, streaming headers/cell types/identifiers/logicals/NA/scientific values\n")
  cat("PASS incorrect A1 declaration rejected by XML and independent streaming reader\n")
  cat("PASS S13 numeric display only; S2 source identities/adjustments without numeric changes\n")
  cat("No formal Run/Export/author decision/reference changed. Temporary fixtures removed.\n")
})
