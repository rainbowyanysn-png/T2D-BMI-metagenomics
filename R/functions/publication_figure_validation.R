# File-format inspection only. This does not refit a model or recalculate a q value.
pub_tiff_properties <- function(path) {
  con<-file(path,"rb");on.exit(close(con),add=TRUE)
  bytes<-readBin(con,"raw",2);order<-rawToChar(bytes)
  assert_true(order%in%c("II","MM"),"Invalid TIFF byte order")
  endian<-if(order=="II")"little" else "big"
  u16<-function()readBin(con,"integer",1,size=2,signed=FALSE,endian=endian)
  u32<-function(){z<-as.double(readBin(con,"integer",1,size=4,endian=endian));if(z<0)z+2^32 else z}
  assert_true(u16()==42L,"Only classic TIFF is supported by this validator")
  offset<-u32();seek(con,offset,"start");n<-u16();values<-list()
  for(i in seq_len(n)){
    tag<-u16();type<-u16();count<-u32();raw<-readBin(con,"raw",4);back<-seek(con)
    numeric_bytes<-as.integer(raw);if(endian=="big")numeric_bytes<-rev(numeric_bytes)
    value<-sum(numeric_bytes*256^(0:3))
    if(type==3L && count==1L)value<-if(endian=="little")sum(as.integer(raw[1:2])*256^(0:1)) else sum(as.integer(raw[1:2])*256^(1:0))
    if(type==5L && count==1L){seek(con,value,"start");numerator<-u32();denominator<-u32();value<-numerator/denominator;seek(con,back,"start")}
    values[[as.character(tag)]]<-value
  }
  list(compression=values[["259"]],x_dpi=values[["282"]],y_dpi=values[["283"]],resolution_unit=values[["296"]])
}

pub_validate_saved_panel <- function(row,export) {
  fields<-c("PDF_path","TIFF_path","PNG_path","plot_data_path","caption_path","selection_rule_path")
  files<-file.path(export,unlist(row[fields],use.names=FALSE))
  present<-all(file.exists(files)) && all(file.info(files)$size>0)
  props<-pub_tiff_properties(file.path(export,row$TIFF_path))
  caption<-readLines(file.path(export,row$caption_path),warn=FALSE)
  pub_language_check(c(caption,row$title,row$selection_rule),row$figure_id)
  label_ok<-TRUE
  d<-utils::read.csv(file.path(export,row$plot_data_path),check.names=FALSE)
  if("label_selected"%in%names(d)){
    take<-which(d$label_selected)
    label_ok<-all(d$qval[take]<.05) && !anyNA(d$label_rank[take]) &&
      identical(sort(as.integer(d$label_rank[take])),seq_along(take)) &&
      all(is.na(d$label_rank[!d$label_selected]))
  }
  expected_dpi<-if(row$figure_id=="SFig3_taxon_pathway_assignment_Sankey")600 else 300
  assert_true(row$dpi==expected_dpi,"Figure index TIFF DPI does not match the panel contract")
  pass<-present && identical(as.integer(props$compression),5L) &&
    isTRUE(abs(props$x_dpi-expected_dpi)<.01) && isTRUE(abs(props$y_dpi-expected_dpi)<.01) &&
    isTRUE(props$resolution_unit==2) && row$background=="white" && label_ok
  check<-data.frame(figure_id=row$figure_id,all_files_nonempty=present,x_dpi=props$x_dpi,y_dpi=props$y_dpi,
    tiff_LZW=props$compression==5,white_background_contract=row$background=="white",
    label_selection_pass=label_ok,publication_language_pass=TRUE,final_pass=pass)
  target<-file.path(export,"validation/figures",paste0(row$figure_id,"__integrity.csv"))
  safe_write_csv(check,target)
  assert_true(pass,paste("Figure format/content contract failed:",row$figure_id))
  invisible(check)
}

pub_finalize_figure_integrity <- function(export) {
  files<-list.files(file.path(export,"validation/figures"),pattern="__integrity[.]csv$",full.names=TRUE)
  report<-if(length(files))pub_bind(lapply(files,utils::read.csv,check.names=FALSE)) else data.frame(figure_id=character(),final_pass=logical())
  safe_write_csv(report,file.path(export,"validation/figure_file_integrity.csv"))
  invisible(report)
}

pub_validate_panel_counts <- function(index,config) {
  registry<-utils::read.csv(file.path(config$project_root,config$publication$panel_registry),stringsAsFactors=FALSE)
  assert_true(!anyDuplicated(index$figure_id)&&nrow(index)==nrow(registry)&&setequal(index$figure_id,registry$figure_id),"Exported panels differ from the final manifest")
  assert_true(!any(grepl("^SFig6|^Fig2F|^SFig5|^SFig4[A-Z]_diversity_diagnostics",index$figure_id)),"Excluded or obsolete supplementary panel is present")
  required_sensitivity<-c("SFig4A_sensitivity_BMI23_species","SFig4B_sensitivity_CLR_species","SFig4C_sensitivity_CLR_pathways")
  assert_true(all(required_sensitivity%in%index$figure_id),"Final sensitivity panels are incomplete")
  invisible(TRUE)
}
