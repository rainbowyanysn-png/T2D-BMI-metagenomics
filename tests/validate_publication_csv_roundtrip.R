# Tiny CSV/Excel fixtures for independent verification; no study export or model.
local({
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
  if(!grepl("^([A-Za-z]:|/)",entry))entry<-file.path(getwd(),entry)
  root<-dirname(dirname(entry))
  source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
  source(file.path(root,"R/functions/publication_bootstrap.R"),local=environment());load_publication_support(root,environment())
  source(file.path(root,"tests/publication_output_review_helpers.R"),local=environment())
  folder<-tempfile("publication_csv_fixture_");dir.create(folder)
  x<-data.frame(sample_id=c("Ctrl_sample_1","sample_2","sample_3","sample_4","sample_5","sample_6"),
    text=c("",NA_character_,"NA","left,,right\n,next",'He said "NA", then left',"__T2DBMI_UNQUOTED_MISSING__"),
    qval=c(1e-30,NA_real_,.049,.051,0,1),flag=c(TRUE,NA,FALSE,TRUE,FALSE,TRUE),
    leading_missing=c(NA_character_,"",NA_character_,"NA","last",""),check.names=FALSE)
  classes<-vapply(x,function(v)class(v)[1],character(1))
  csv<-file.path(folder,"fixture.csv");safe_write_csv(x,csv)
  broken<-utils::read.csv(csv,check.names=FALSE,colClasses=classes,na.strings="")
  assert_true(is.na(broken$text[1]),"Fixture must reproduce the quoted-empty reader defect")
  expected<-publication_read_export_csv(csv,classes)
  for(name in names(x))assert_true(identical(x[[name]],expected[[name]]),paste("CSV round-trip changed",name))
  # Missing fields at the first/last column and CRLF record boundaries.
  boundary<-data.frame(first=c(NA_character_,"","NA"),last=c("",NA_character_,"NA"))
  boundary_csv<-file.path(folder,"boundary.csv");safe_write_csv(boundary,boundary_csv)
  back<-publication_read_export_csv(boundary_csv,c(first="character",last="character"))
  assert_true(identical(boundary,back),"Leading/trailing missing fields changed")
  crlf_csv<-file.path(folder,"crlf.csv")
  writeBin(charToRaw(paste0(paste(c('"first","last"',',""','"",','"NA","NA"'),collapse="\r\n"),"\r\n")),crlf_csv)
  back<-publication_read_export_csv(crlf_csv,c(first="character",last="character"))
  assert_true(identical(boundary,back),"CRLF missing fields changed")
  unicode<-data.frame(text=c("\u03b1\u2013\u03b2\n\u4e2d\u6587","",NA_character_,"NA"))
  unicode_csv<-file.path(folder,"unicode.csv")
  # Create exact UTF-8 bytes even when the test process inherits an ASCII locale.
  writeBin(charToRaw(enc2utf8(paste0('"text"\r\n"',unicode$text[1],'"\r\n""\r\n\r\n"NA"\r\n'))),unicode_csv)
  back<-publication_read_export_csv(unicode_csv,c(text="character"))
  assert_true(identical(unicode,back),"Unicode or single-column fields changed")
  # Verify against the independently read workbook, without broadening the
  # comparator to equate empty strings and missing values in numeric columns.
  wb<-openxlsx::createWorkbook();openxlsx::addWorksheet(wb,"Results")
  openxlsx::writeData(wb,"Results",x,keepNA=TRUE,na.string="NA")
  pub_clear_empty_drawing_placeholders(wb)
  path<-file.path(folder,"fixture.xlsx");openxlsx::saveWorkbook(wb,path)
  observed<-openxlsx::read.xlsx(path,sheet="Results",check.names=FALSE,skipEmptyRows=FALSE,
    skipEmptyCols=FALSE,na.strings=NULL)
  pub_verify_workbook_sheet(expected,observed,"fixture","Results",folder)
  changed<-observed;changed$qval[3]<-.051
  error<-tryCatch({pub_verify_workbook_sheet(expected,changed,"changed","Results",folder);NULL},error=function(e)e)
  assert_true(inherits(error,"error"),"Changed numeric/significance state was not rejected")
  changed<-observed;changed$qval[1]<-0
  error<-tryCatch({pub_verify_workbook_sheet(expected,changed,"zeroed","Results",folder);NULL},error=function(e)e)
  assert_true(inherits(error,"error"),"A tiny nonzero probability became zero without rejection")
  message("PASS: quoted empty strings, true missing values, literal NA, escaped/multiline CSV, CRLF, typed Excel round trip and strict numeric rejection")
})
