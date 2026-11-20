# Independent post-export verification only. Exporters never source this file.
publication_read_export_csv <- function(file, classes) {
  # safe_write_csv uses write.csv(na=""): a quoted "" is an actual empty
  # string; an UNQUOTED empty field is missing. na.strings="" erases that
  # distinction. Protect quoted fields, including escaped quotes/newlines,
  # and mark only unquoted empties before asking R to read the typed data.
  con<-base::file(file,"rb")
  text<-tryCatch(readChar(con,nchars=file.info(file)$size,useBytes=TRUE),finally=close(con))
  marker<-"__T2DBMI_UNQUOTED_MISSING__"
  while(grepl(marker,text,fixed=TRUE))marker<-paste0(marker,"_")
  pattern<-'"(?:[^"]|"")*"(*SKIP)(*F)|(?<=,)(?=,|\\r?\\n|\\z)|(?m)^(?=,|\\r?\\n)'
  text<-gsub(pattern,marker,text,perl=TRUE)
  # textConnection supplies its own last line terminator. Retaining the CSV's
  # terminator would create an extra blank record when blank.lines.skip=FALSE.
  text<-sub("\\r?\\n\\z","",text,perl=TRUE)
  utf8<-validUTF8(text)
  if(utf8)Encoding(text)<-"UTF-8"
  connection<-textConnection(text,encoding=if(utf8)"UTF-8" else "bytes")
  on.exit(close(connection),add=TRUE)
  utils::read.csv(connection,check.names=FALSE,colClasses=classes,na.strings=marker,
    blank.lines.skip=FALSE,comment.char="",encoding=if(utf8)"UTF-8" else "unknown")
}

publication_export_directory <- function(root,run_id,kind) {
  pinned<-Sys.getenv("T2DBMI_VALIDATION_EXPORT")
  if(nzchar(pinned)){
    assert_true(!grepl("[/\\\\]|[.][.]",pinned),"Pinned export must be a directory name")
    folder<-file.path(root,"output",pinned)
    meta<-utils::read.csv(file.path(folder,"manifests/execution_context.csv"))
    assert_true(nrow(meta)==1L && identical(meta$run_id,run_id),"Pinned export belongs to a different Run")
    return(folder)
  }
  folders<-sort(list.dirs(file.path(root,"output"),recursive=FALSE,full.names=TRUE),decreasing=TRUE)
  folders<-folders[grepl("^PublicationExport_",basename(folders))]
  for(folder in folders){
    context<-file.path(folder,"manifests/execution_context.csv");readme<-file.path(folder,"README.md")
    if(!file.exists(context)||!file.exists(readme))next
    meta<-utils::read.csv(context,check.names=FALSE)
    mode<-sub("^Mode: ","",grep("^Mode: ",readLines(readme,warn=FALSE),value=TRUE))
    if(nrow(meta)==1L && identical(meta$run_id,run_id) && length(mode)==1L && mode%in%c("all",kind))return(folder)
  }
  stop(paste("No matching publication export for",kind,"and",run_id),call.=FALSE)
}

validate_publication_output_directories <- function(root,run_id,out,config) {
  dirs<-setNames(vapply(c("tables","figures"),function(k)publication_export_directory(root,run_id,k),character(1)),c("tables","figures"))
  safe_write_csv(data.frame(kind=names(dirs),directory=unname(dirs)),file.path(out,"export_directories_checked.csv"))
  for(folder in unique(dirs)){
    gaps<-utils::read.csv(file.path(folder,"export_gaps.csv"),check.names=FALSE)
    warnings<-utils::read.csv(file.path(folder,"logs/warnings.csv"),check.names=FALSE)
    assert_true(!nrow(gaps),paste("Export has unresolved failures:",folder))
    assert_true(!nrow(warnings),paste("Export has warnings requiring explanation:",folder))
    immutable<-utils::read.csv(file.path(folder,"validation/input_immutability.csv"))
    assert_true(all(immutable$unchanged),"Export input-immutability check failed")
  }
  tables<-utils::read.csv(file.path(dirs[["tables"]],"table_index.csv"),check.names=FALSE)
  submission_validate_table_identity(tables)
  workbook_out<-file.path(out,"workbook_checks");dir.create(workbook_out,recursive=TRUE)
  for(i in seq_len(nrow(tables))){
    id<-tables$table_id[i];target<-file.path(dirs[["tables"]],tables$workbook_path[i])
    assert_true(identical(sha256_file(target),tables$workbook_sha256[i]),paste("Workbook hash changed:",id))
    dictionary<-utils::read.csv(file.path(dirs[["tables"]],"tables/csv",paste0(id,"__data_dictionary.csv")),check.names=FALSE)
    sheet_names<-unique(dictionary$sheet)
    assert_true(tail(sheet_names,1L)=="Notes" && !sheet_names[1L]%in%c("README","Notes"),"Results must precede Notes")
    sheets<-setNames(lapply(sheet_names,function(name){
      types<-dictionary[dictionary$sheet==name,]
      classes<-setNames(ifelse(types$storage_class%in%c("numeric","integer","logical","character"),types$storage_class,"character"),types$column)
      publication_read_export_csv(file.path(dirs[["tables"]],"tables/csv",paste0(id,"__",name,".csv")),classes)
    }),sheet_names)
    # Windows ZIP libraries may misdecode a valid Chinese path. Read a byte-identical
    # ASCII temporary copy; never rewrite or repair the submitted workbook.
    staged<-tempfile("mr_workbook_",fileext=".xlsx")
    assert_true(all(as.integer(charToRaw(staged))<128L),"Validation TEMP path must be ASCII")
    assert_true(file.copy(target,staged,overwrite=FALSE),"Cannot stage read-only workbook")
    assert_true(identical(sha256_file(staged),sha256_file(target)),"Staged workbook hash differs")
    tryCatch(pub_validate_saved_workbook(staged,sheets,id,workbook_out),finally=unlink(staged))
  }
  pub_finalize_workbook_integrity(workbook_out)
  figures<-utils::read.csv(file.path(dirs[["figures"]],"figure_index.csv"),check.names=FALSE,na.strings="")
  pub_validate_panel_counts(figures,config)
  file_manifest<-utils::read.csv(file.path(dirs[["figures"]],"output_manifest.csv"),check.names=FALSE)
  fields<-c("PDF_path","TIFF_path","PNG_path","plot_data_path","caption_path","selection_rule_path")
  checked<-list()
  for(i in seq_len(nrow(figures))){
    z<-figures[i,,drop=FALSE]
    paths<-unlist(z[,fields],use.names=FALSE)
    for(p in paths){
      index<-match(p,file_manifest$file)
      assert_true(!is.na(index) && identical(sha256_file(file.path(dirs[["figures"]],p)),file_manifest$sha256[index]),
        paste("Figure artifact changed:",p))
    }
    properties<-pub_tiff_properties(file.path(dirs[["figures"]],z$TIFF_path))
    pub_language_check(readLines(file.path(dirs[["figures"]],z$caption_path),warn=FALSE),z$figure_id)
    expected_dpi<-if(z$figure_id=="SFig3_taxon_pathway_assignment_Sankey")600 else 300
    assert_true(z$dpi==expected_dpi,"Figure index TIFF DPI contract differs")
    checked[[i]]<-data.frame(figure_id=z$figure_id,x_dpi=properties$x_dpi,y_dpi=properties$y_dpi,
      LZW=properties$compression==5,pass=properties$compression==5 &&
        abs(properties$x_dpi-expected_dpi)<.01 && abs(properties$y_dpi-expected_dpi)<.01)
  }
  result<-pub_bind(checked);safe_write_csv(result,file.path(out,"figure_format_checks.csv"))
  assert_true(all(result$pass),"Figure file-format checks failed")
  legend_names<-paste0(gsub(" ","_",unique(figures$manuscript_figure)),"_legend.txt")
  for(name in legend_names){
    file<-file.path(dirs[["figures"]],"figures/combined_legends",name)
    assert_file(file);pub_language_check(readLines(file,warn=FALSE),name)
  }
  invisible(dirs)
}
