# Small synthetic workbooks only; never exports real study tables or fits models.
local({
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
  if(!grepl("^([A-Za-z]:|/)",entry))entry<-file.path(getwd(),entry)
  root<-dirname(dirname(entry))
  source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
  source(file.path(root,"R/functions/publication_bootstrap.R"),local=environment());load_publication_support(root,environment())
  require_project_packages(c("openxlsx","readxl","xml2","dplyr","digest"))
  folder<-tempfile("publication_workbook_fixture_");dir.create(folder)
  raw<-openxlsx::createWorkbook();openxlsx::addWorksheet(raw,"fixture")
  openxlsx::writeData(raw,"fixture",data.frame(value=1:2))
  raw_file<-file.path(folder,"openxlsx_default_probe.xlsx");openxlsx::saveWorkbook(raw,raw_file)
  before<-pub_ooxml_integrity(raw_file)
  assert_true(before$dangling>0L,"The pinned openxlsx empty-drawing regression was not reproduced; inspect writer version")
  x<-data.frame(feature=c("NA",NA_character_,"taxon_A","taxon_B"),
    qval=c(1e-30,NA_real_,.049,.051),flag=c(TRUE,NA,FALSE,TRUE))
  config<-read_config(root)
  index<-pub_workbook("Fixture","Export contract fixture",list(Results=x,
    Model_and_q_definition=data.frame(model="Prespecified model",family="Full-model q values")),
    "Synthetic data only.",folder,config)
  report<-utils::read.csv(file.path(folder,"validation/workbooks/Fixture__integrity.csv"))
  assert_true(all(report$final_pass) && report$dangling_relationship_count==0L,"Fixture OOXML checks failed")
  table1<-data.frame(Characteristic=c("N","Study cohort, n"),Overall=c("1186",""),
    NC=c("307",""),OC=c("119",""),NT2D=c("138",""),OT2D=c("622",""))
  pub_workbook("Main_Table_1","Participant characteristics",list(Characteristics=table1),
    "Values are descriptive.",folder,config,main=TRUE)
  main_file<-file.path(folder,"tables/main/Main_Table_1.xlsx")
  doc<-pub_xml_part(main_file,"xl/worksheets/sheet1.xml")
  setup<-xml2::xml_find_first(doc,"//*[local-name()='pageSetup']")
  assert_true(xml2::xml_attr(setup,"orientation")=="landscape" && xml2::xml_attr(setup,"fitToWidth")=="1",
    "Main Table 1 landscape/width contract failed")
  book<-pub_xml_part(main_file,"xl/workbook.xml")
  assert_true(length(xml2::xml_find_all(book,"//*[local-name()='definedName' and @name='_xlnm.Print_Area']"))==1L,
    "Main Table 1 print area missing")
  assert_true(identical(readxl::excel_sheets(main_file),c("Characteristics","Notes")),"Notes must be last")
  pub_finalize_workbook_integrity(folder)
  message("PASS: default writer defect reproduced; in-memory adapter, OOXML, R/readxl round trips and print setup checked. Fixtures: ",folder)
})
