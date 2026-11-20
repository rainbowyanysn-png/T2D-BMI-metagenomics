# Reader-facing language/column projection fixtures only.
local({
  entry<-sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))
  if(!grepl("^([A-Za-z]:|/)",entry))entry<-file.path(getwd(),entry)
  root<-dirname(dirname(entry))
  source(file.path(root,"R/functions/project_paths.R"),local=environment());activate_existing_project_library(root)
  source(file.path(root,"R/functions/publication_bootstrap.R"),local=environment());load_publication_support(root,environment())
  fails<-function(text)inherits(tryCatch({pub_language_check(text,"fixture");NULL},error=function(e)e),"error")
  forbidden<-c("frozen","HISTORICAL","old result","revised analysis","final model","audit only",
    "removed medication","retired_not_valid","candidate_data_present_requires_source_verification",
    "awaiting_user","pending medication","validation")
  assert_true(all(vapply(forbidden,fails,logical(1))),"Forbidden language was not rejected")
  pub_language_check("Prespecified threshold; medication-adjusted sensitivity model; not reported in the database","allowed")
  x<-data.frame(feature=c("taxon_A","taxon_B"),coef=c(-.21,.38),pval=c(1e-30,.049),qval=c(.002,.051),
    formal_q=c(.002,.051),sample_id=c("Ctrl_sample_1","Ctrl_sample_2"),model_id="internal",source_result_file="internal")
  y<-pub_reader_table(x)
  assert_true(identical(x$coef,y$coef)&&identical(x$pval,y$pval)&&identical(x$qval,y$qval),"Display projection altered statistics")
  assert_true(identical(x$sample_id,y$sample_id)&&identical(x$feature,y$feature),"Display terminology altered sample or feature IDs")
  assert_true(!any(c("model_id","source_result_file")%in%names(y)),"Internal columns escaped projection")
  pub_language_check(pub_caption_text("Frozen result; MaAsLin2 qval; continuous_BMI; Ctrl"),"display conversion")
  message("PASS: reader language and numeric-preserving projection fixtures")
})
