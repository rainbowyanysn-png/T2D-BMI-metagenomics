validate_no_legacy_dependency <- function(root) {
 files<-list.files(file.path(root,"R"),pattern="[.]R$",recursive=TRUE,full.names=TRUE)
 rows<-lapply(files,function(f){
  lines<-readLines(f,warn=FALSE)
  bad<-grep("reference_results_manifest|internal_audit|submission_materials|FinalIntegration|FigureRefinement|IntegrationPreview|CleanTest|LOCO_Trial|MedicationRestored|T2D_BMI_RC2|[A-Za-z]:[/\\\\]|Run_[0-9]{8}_[0-9]{6}|load[(]",lines)
  if(basename(f)!="checkpoints.R")bad<-union(bad,grep("readRDS[(]",lines))
  # Path recognition syntax is not a hardcoded personal absolute path.
  bad<-bad[!grepl("grepl[(]|normalizePath",lines[bad])]
  data.frame(file=substring(f,nchar(root)+2L),pass=!length(bad),lines=paste(bad,collapse=";"))
 })
 result<-do.call(rbind,rows);assert_true(all(result$pass),"Legacy/private input dependency in public analysis code");result
}
