# OOXML checks read ZIP entries directly. Saved workbooks are never unpacked or patched.
pub_xml_part <- function(file, member) {
  con <- unz(file, member, open = "rb")
  on.exit(close(con), add = TRUE)
  xml2::read_xml(con)
}

pub_ooxml_resolve <- function(source, target) {
  target <- utils::URLdecode(sub("#.*$", "", target))
  parts <- strsplit(if (startsWith(target, "/")) substring(target, 2L) else
    paste0(if (nzchar(source)) paste0(dirname(source), "/") else "", target), "/", fixed = TRUE)[[1]]
  stack <- character()
  for (p in parts) {
    if (p %in% c("", ".")) next
    if (p == "..") {
      if (!length(stack)) stop("OOXML target escapes package root", call. = FALSE)
      stack <- head(stack, -1L)
    } else stack <- c(stack, p)
  }
  paste(stack, collapse = "/")
}

pub_ooxml_integrity <- function(file) {
  require_project_packages("xml2")
  parts <- utils::unzip(file, list = TRUE)$Name
  relationships <- list()
  for (member in parts[grepl("[.]rels$", parts)]) {
    doc <- pub_xml_part(file, member)
    nodes <- xml2::xml_find_all(doc, "//*[local-name()='Relationship']")
    source <- if (member == "_rels/.rels") "" else sub("/_rels/([^/]+)[.]rels$", "/\\1", member)
    for (node in nodes) {
      if (identical(xml2::xml_attr(node, "TargetMode"), "External")) next
      target <- xml2::xml_attr(node, "Target")
      resolved <- pub_ooxml_resolve(source, target)
      relationships[[length(relationships) + 1L]] <- data.frame(relationship_part = member,
        relationship_id = xml2::xml_attr(node, "Id"), target = target, resolved_part = resolved,
        type = xml2::xml_attr(node, "Type"), exists = resolved %in% parts)
    }
  }
  links <- if (length(relationships)) pub_bind(relationships) else data.frame(
    relationship_part = character(), relationship_id = character(), target = character(),
    resolved_part = character(), type = character(), exists = logical())
  types <- pub_xml_part(file, "[Content_Types].xml")
  declared <- xml2::xml_attr(xml2::xml_find_all(types, "//*[local-name()='Override']"), "PartName")
  missing <- setdiff(sub("^/", "", declared), parts)
  book <- pub_xml_part(file, "xl/workbook.xml")
  sheet_order <- xml2::xml_attr(xml2::xml_find_all(book, "//*[local-name()='sheet']"), "name")
  errors <- sum(vapply(parts[grepl("^xl/worksheets/sheet[0-9]+[.]xml$", parts)], function(p)
    length(xml2::xml_find_all(pub_xml_part(file, p), "//*[local-name()='c' and @t='e']")), integer(1)))
  list(sheet_order = sheet_order, links = links, dangling = sum(!links$exists),
    missing_parts = length(unique(c(links$resolved_part[!links$exists], missing))), formula_errors = errors)
}

pub_clear_empty_drawing_placeholders <- function(wb) {
  # openxlsx 4.2.8.1 pre-populates drawing/VML links on an empty worksheet.
  # These exports contain data tables only. Reject actual drawing/comment objects,
  # then remove only their empty placeholders IN MEMORY, before saveWorkbook().
  assert_true(identical(as.character(utils::packageVersion("openxlsx")), "4.2.8.1"),
    "Recheck the table-only OOXML writer adapter before changing openxlsx version")
  for (i in seq_along(wb$worksheets)) {
    assert_true(!length(wb$drawings[[i]]) && !length(wb$comments[[i]]),
      "Table-only writer cannot discard actual drawings or comments")
    links <- wb$worksheets_rels[[i]]
    empty <- grepl('/relationships/(drawing|vmlDrawing)"', links)
    wb$worksheets_rels[[i]] <- links[!empty]
    wb$worksheets[[i]]$drawing <- character()
    wb$worksheets[[i]]$legacyDrawing <- character()
  }
  # The package also pre-declares absent drawing parts in Content_Types.
  # Remove those declarations before serialization, after the table-only guards.
  wb$Content_Types <- wb$Content_Types[!grepl('PartName="/xl/drawings/', wb$Content_Types, fixed = TRUE)]
  invisible(wb)
}

pub_record_workbook_integrity <- function(export, row) {
  path <- file.path(export, "validation/workbooks", paste0(row$workbook,"__integrity.csv"))
  safe_write_csv(row, path)
}

pub_finalize_workbook_integrity <- function(export) {
  files<-list.files(file.path(export,"validation/workbooks"),pattern="__integrity[.]csv$",full.names=TRUE)
  report<-if(length(files))pub_bind(lapply(files,utils::read.csv,check.names=FALSE)) else data.frame(workbook=character(),final_pass=logical())
  safe_write_csv(report,file.path(export,"validation/workbook_ooxml_integrity.csv"))
  invisible(report)
}

pub_validate_saved_workbook <- function(target, sheets, id, export, root=NULL) {
  report <- data.frame(workbook = id, sheet_count = length(sheets),
    expected_sheet_order = paste(names(sheets), collapse = " | "), actual_sheet_order = "",
    dangling_relationship_count = NA_integer_, missing_internal_part_count = NA_integer_,
    R_readback_pass = FALSE, third_party_readback_pass = FALSE, numeric_content_pass = FALSE,
    formula_error_count = NA_integer_, worksheet_range_pass = FALSE, streaming_readback_pass = FALSE,
    publication_language_pass = FALSE, final_pass = FALSE,
    error = "", stringsAsFactors = FALSE)
  tryCatch({
    integrity <- pub_ooxml_integrity(target)
    report$actual_sheet_order <- paste(integrity$sheet_order, collapse = " | ")
    report$dangling_relationship_count <- integrity$dangling
    report$missing_internal_part_count <- integrity$missing_parts
    report$formula_error_count <- integrity$formula_errors
    safe_write_csv(integrity$links, file.path(export, "validation/workbooks", paste0(id, "__relationships.csv")))
    assert_true(identical(names(sheets), integrity$sheet_order), paste("Sheet order differs:", id))
    assert_true(integrity$dangling == 0L && integrity$missing_parts == 0L && integrity$formula_errors == 0L,
      paste("Invalid workbook OOXML or formula cells:", id))
    ranges <- pub_sheet_dimension_report(target,sheets)
    safe_write_csv(ranges,file.path(export,"validation/workbooks",paste0(id,"__xml_ranges.csv")))
    assert_true(all(ranges$range_pass),paste("Worksheet declared/actual/expected range differs:",id))
    report$worksheet_range_pass <- TRUE
    require_project_packages("readxl")
    assert_true(identical(readxl::excel_sheets(target), names(sheets)), "Third-party sheet order differs")
    for (name in names(sheets)) {
      x <- as.data.frame(sheets[[name]], check.names = FALSE)
      back <- openxlsx::read.xlsx(target, sheet = name, check.names = FALSE, sep.names = " ",
        skipEmptyRows = FALSE, skipEmptyCols = FALSE, na.strings = NULL)
      pub_verify_workbook_sheet(x, back, id, name, export)
      other <- as.data.frame(readxl::read_xlsx(target, sheet = name, col_types = "text",
        .name_repair = "minimal", trim_ws = FALSE, na = character()), check.names = FALSE)
      assert_true(nrow(other) == nrow(x) && identical(names(other), names(x)),
        paste("Third-party workbook shape differs:", id, name))
      pub_verify_workbook_sheet(x, other, id, name, export)
      pub_language_check(c(names(x), unlist(x, use.names = FALSE)), paste(id, name))
    }
    report$R_readback_pass <- report$numeric_content_pass <- TRUE
    report$third_party_readback_pass <- report$publication_language_pass <- TRUE
    streaming <- pub_stream_read_workbook(target,sheets,id,export,root)
    ranges <- merge(ranges,streaming[,c("sheet","read_rows","read_columns","streaming_pass")],by="sheet",sort=FALSE)
    ranges <- ranges[match(names(sheets),ranges$sheet),,drop=FALSE]
    safe_write_csv(ranges,file.path(export,"validation/workbooks",paste0(id,"__ranges.csv")))
    report$streaming_readback_pass <- TRUE
    report$final_pass <- TRUE
  }, error = function(e) report$error <<- conditionMessage(e))
  pub_record_workbook_integrity(export, report)
  assert_true(report$final_pass, paste("Workbook validation failed:", id, report$error,
    "See validation/workbook_ooxml_integrity.csv"))
  invisible(report)
}

# All changes occur in the in-memory workbook, before saveWorkbook().
pub_written_range <- function(x) {
  if (!ncol(x)) {
    assert_true(!nrow(x), "A worksheet with rows but no columns is invalid")
    return("A1")
  }
  paste0("A1:", openxlsx::int2col(ncol(x)), nrow(x) + 1L)
}

pub_set_written_dimensions <- function(wb, sheets) {
  assert_true(identical(as.character(utils::packageVersion("openxlsx")), "4.2.8.1"),
    "Recheck the worksheet-dimension adapter before changing openxlsx version")
  assert_true(identical(names(wb), names(sheets)), "Workbook/sheet order differs before dimension assignment")
  for (i in seq_along(sheets))
    wb$worksheets[[i]]$dimension <- paste0('<dimension ref="', pub_written_range(sheets[[i]]), '"/>')
  invisible(wb)
}

pub_sheet_dimension_report <- function(target, sheets) {
  book <- pub_xml_part(target, "xl/workbook.xml")
  nodes <- xml2::xml_find_all(book, "//*[local-name()='sheet']")
  links <- pub_xml_part(target, "xl/_rels/workbook.xml.rels")
  relations <- xml2::xml_find_all(links, "//*[local-name()='Relationship']")
  rel_id <- xml2::xml_attr(relations, "Id")
  entries <- lapply(seq_along(sheets), function(i) {
    # Relationship IDs, rather than numeric worksheet filenames, identify sheets.
    attrs <- xml2::xml_attrs(nodes[[i]])
    rid <- unname(attrs[grepl("(^|:)id$", names(attrs))])
    assert_true(length(rid) == 1L && sum(rel_id == rid) == 1L, "Ambiguous worksheet relationship")
    member <- pub_ooxml_resolve("xl/workbook.xml", xml2::xml_attr(relations[[match(rid, rel_id)]], "Target"))
    doc <- pub_xml_part(target, member)
    declared <- xml2::xml_attr(xml2::xml_find_first(doc, "//*[local-name()='dimension']"), "ref")
    cells <- xml2::xml_find_all(doc, "//*[local-name()='sheetData']/*[local-name()='row']/*[local-name()='c'][*[local-name()='v' or local-name()='is' or local-name()='f']]")
    refs <- xml2::xml_attr(cells, "r")
    actual <- "A1"
    if (length(refs)) {
      cols <- openxlsx::col2int(sub("[0-9]+$", "", refs))
      rows <- as.integer(sub("^[A-Z]+", "", refs))
      actual <- paste0(openxlsx::int2col(min(cols)), min(rows), ":", openxlsx::int2col(max(cols)), max(rows))
    }
    expected <- pub_written_range(sheets[[i]])
    data.frame(sheet=names(sheets)[i], worksheet_part=member, declared_range=declared,
      actual_range=actual, expected_range=expected,
      expected_rows=if(ncol(sheets[[i]]))nrow(sheets[[i]])+1L else 0L,
      expected_columns=ncol(sheets[[i]]),
      range_pass=!is.na(declared) && identical(declared,expected) && identical(actual,expected))
  })
  pub_bind(entries)
}

pub_stream_reader_script <- function(root = NULL) {
  if (is.null(root)) {
    arg <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE))
    assert_true(length(arg)==1L,"Supply the project root for workbook streaming validation")
    candidates <- c(dirname(arg), dirname(dirname(arg)))
    found <- candidates[file.exists(file.path(candidates,"tests/workbook_stream_read.py"))]
    assert_true(length(found)>0L,"Missing tests/workbook_stream_read.py")
    root <- found[1]
  }
  script <- file.path(root,"tests/workbook_stream_read.py")
  assert_file(script,"Independent streaming workbook reader")
  script
}

pub_stream_reader_preflight <- function(root = NULL) {
  python <- Sys.getenv("T2DBMI_VALIDATION_PYTHON")
  if(!nzchar(python)) {
    available <- unname(Sys.which(c("python","python3")))
    available <- available[nzchar(available)]
    assert_true(length(available)>0L,paste("Workbook acceptance requires Python with openpyxl.",
      "Set T2DBMI_VALIDATION_PYTHON to its executable; no package is installed at runtime."))
    python <- available[1]
  }
  script <- pub_stream_reader_script(root)
  log <- tempfile("stream_probe_",fileext=".txt");on.exit(unlink(log),add=TRUE)
  code <- suppressWarnings(system2(python,c("-B",shQuote(script),"--probe"),stdout=log,stderr=log))
  assert_true(identical(as.integer(code),0L),paste("Independent workbook reader is unavailable.",
    "Set T2DBMI_VALIDATION_PYTHON to Python with openpyxl; probe output:",paste(readLines(log,warn=FALSE),collapse=" | ")))
  invisible(list(python=python,script=script))
}

pub_stream_read_workbook <- function(target,sheets,id,export,root=NULL) {
  require_project_packages("jsonlite")
  reader <- pub_stream_reader_preflight(root)
  spec <- list(sheets=lapply(names(sheets),function(name) {
    x<-sheets[[name]]
    list(name=name,expected_range=pub_written_range(x),rows=if(ncol(x))nrow(x)+1L else 0L,
      columns=lapply(names(x),function(n) {
        v<-x[[n]];kind<-if(is.numeric(v))"numeric" else if(is.logical(v))"logical" else "character"
        if(!kind%in%c("numeric","logical"))v<-as.character(v)
        list(name=n,kind=kind,values=as.list(v),missing=as.list(is.na(v)))
      }))
  }))
  path<-tempfile("stream_expected_",fileext=".json");on.exit(unlink(path),add=TRUE)
  jsonlite::write_json(spec,path,auto_unbox=TRUE,na="null",null="null",digits=NA)
  report<-file.path(export,"validation/workbooks",paste0(id,"__streaming.json"))
  log<-file.path(export,"validation/workbooks",paste0(id,"__streaming.log"))
  dir.create(dirname(report),recursive=TRUE,showWarnings=FALSE)
  code<-suppressWarnings(system2(reader$python,c("-B",shQuote(reader$script),"--workbook",shQuote(target),
    "--spec",shQuote(path),"--report",shQuote(report)),stdout=log,stderr=log))
  assert_file(report,paste("Independent reader did not produce a report; see",log))
  result<-jsonlite::read_json(report,simplifyVector=FALSE)
  rows<-lapply(result$sheets,function(z)data.frame(sheet=z$sheet,declared_range=z$declared_range,
    actual_range=z$actual_range,expected_range=z$expected_range,
    read_rows=z$read_rows,read_columns=z$read_columns,streaming_pass=isTRUE(z$pass_status),
    error=paste(unlist(z$errors),collapse=" | ")))
  if(length(rows))safe_write_csv(pub_bind(rows),file.path(export,"validation/workbooks",paste0(id,"__streaming_ranges.csv")))
  assert_true(identical(as.integer(code),0L)&&isTRUE(result$pass_status),
    paste("Independent streaming/type/value validation failed; see",report))
  invisible(pub_bind(rows))
}
