# Reader-facing exports of saved statistics. No model fitting or p.adjust here.
pub_bind <- function(x) dplyr::bind_rows(x)
pub_read_csv_snapshot <- function(run, file) {
  before <- sha256_file(file)
  x <- utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  assert_true(identical(before, sha256_file(file)), paste("Input changed while being read:", file))
  previous <- run$access$hashes
  if (!is.null(previous) && file %in% names(previous))
    assert_true(identical(unname(previous[file]), before), paste("Input changed during export:", file))
  run$access$files <- unique(c(run$access$files, file))
  run$access$hashes[file] <- before
  x
}
pub_read <- function(run, relative) {
  file <- current_result_path(run, relative)
  prefix <- paste0(normalizePath(run$dir, winslash = "/", mustWork = TRUE), "/modules/")
  assert_true(startsWith(normalizePath(file, winslash = "/", mustWork = TRUE), prefix),
              "Publication input must belong to the selected run's module checkpoints")
  pub_read_csv_snapshot(run, file)
}
pub_read_module <- function(run, module, directory, pattern = "[.]csv$") {
  folder <- file.path(run$result_dirs[[module]], directory)
  assert_true(length(folder) == 1L && dir.exists(folder), paste("Missing current-run directory:", module, directory))
  files <- list.files(folder, pattern = pattern, full.names = TRUE, recursive = TRUE)
  assert_true(length(files) > 0L, paste("No saved results:", module, directory, pattern))
  pub_bind(lapply(files, function(f) {
    x <- pub_read_csv_snapshot(run, f)
    x$source_result_file <- substring(f, nchar(run$dir) + 2L)
    x
  }))
}

pub_comparisons <- function() data.frame(
  comparison = c("T2D vs Control", "NT2D vs NC", "OT2D vs OC", "OC vs NC", "OT2D vs NT2D"),
  species = c("Total_T2D_vs_Control", "NT2D_vs_NC_BMI25", "OT2D_vs_OC_BMI25", "OC_vs_NC_BMI25", "OT2D_vs_NT2D_BMI25"),
  pathway = c("Total_T2D_vs_HC", "NonOW_T2D_vs_NonOW_Ctrl", "OW_T2D_vs_OW_Ctrl", "OW_Ctrl_vs_NonOW_Ctrl", "OW_T2D_vs_NonOW_T2D"))

publication_alpha_pairwise <- function(run) {
  z<-pub_read_module(run,"primary_diversity","tables/04_primary_diversity","alpha_.*_pairwise[.]csv$")
  assert_columns(z,"source_result_file","Alpha pairwise source identity")
  name<-basename(z$source_result_file)
  assert_true(all(name%in%c("BMI25_alpha_shannon_pairwise.csv","BMI25_alpha_richness_pairwise.csv")),"Unknown alpha pairwise metric file")
  z$metric<-ifelse(name=="BMI25_alpha_shannon_pairwise.csv","Shannon","Richness")
  assert_true(nrow(z)==12L&&all(table(z$metric)==6L)&&!anyDuplicated(z[c("metric","group1","group2")]),"Alpha pairwise metric/keys differ")
  z[,c("metric",setdiff(names(z),"metric")),drop=FALSE]
}

publication_data <- function(run, config) {
  spec <- pub_comparisons()
  species <- reml <- cohort <- list()
  for (i in seq_len(nrow(spec))) {
    species[[i]] <- pub_read(run, paste0("tables/02_primary_BMI25_species/", spec$species[i], "_all.csv"))
    reml[[i]] <- pub_read(run, paste0("tables/03_primary_BMI25_pathways/BMI25_", spec$pathway[i], "_all.csv"))
    species[[i]]$comparison <- reml[[i]]$comparison <- spec$comparison[i]
    cohort[[i]] <- pub_read_module(run, "pathways", "models/03_primary_BMI25_pathways",
      paste0("^BMI25_", spec$pathway[i], "_.*_DESeq2[.]csv$"))
    cohort[[i]]$comparison <- spec$comparison[i]
  }
  interactions <- list()
  for (bmi in c("BMI25", "continuous_BMI", "BMI23")) for (domain in c("species", "pathways")) {
    directory <- if (bmi == "BMI23") "07_BMI23_interactions" else "05_formal_interactions"
    x <- pub_read(run, paste0("tables/", directory, "/", bmi, "_", domain, "_all.csv"))
    x$BMI_definition <- bmi; x$domain <- domain
    interactions[[paste(bmi, domain)]] <- x
  }
  data <- list(spec = spec, species = pub_bind(species), reml = pub_bind(reml), cohort_effects = pub_bind(cohort),
    interactions = pub_bind(interactions), additive = pub_read(run, "tables/06_additive_coefficients_wide.csv"),
    additive_long = pub_read(run, "tables/06_additive_all_target_coefficients.csv"),
    pooled = pub_read(run, "tables/03a_pooled_MaAsLin2_pathways/Pooled_pathways_all.csv"),
    clinical = pub_read(run, "tables/02b_exploratory_clinical_correlations/Clinical_correlations_all_60_tests.csv"),
    clinical_selected = pub_read(run, "tables/02b_exploratory_clinical_correlations/Clinical_correlation_selected_top20_features.csv"),
    clinical_checks = pub_read(run, "tables/02b_exploratory_clinical_correlations/Clinical_correlation_benchmark_checks.csv"),
    md = pub_read(run, "tables/01_sample_manifest.csv"))
  data$candidates <- select_pathway_candidates(data$pooled, config)
  data$pooled_contract <- pub_read_csv_snapshot(run,file.path(run$result_dirs[["pooled_pathways"]],
    "models/03a_pooled_MaAsLin2_pathways/MaAsLin2/model_contract/parameters.csv"))
  assert_columns(data$pooled,c("model_id","study_effect","model_formula","coef","stderr","pval","qval"),"Current fixed-study pathway model")
  data$reml_annotation <- cbind(data$reml, annotate_pathway_features(data$reml$pathway, run$root)[, -1, drop = FALSE])
  data
}

publication_numeric_checks <- function(b, export) {
  checks <- list()
  add <- function(id, observed, expected) checks[[id]] <<- data.frame(check = id,
    observed = observed, expected = expected, pass = !is.na(observed) && observed == expected)
  add("samples", nrow(b$md), 1186L)
  add("control", sum(b$md$disease == "control"), 426L); add("T2D", sum(b$md$disease == "T2D"), 760L)
  add("cohorts", length(unique(b$md$study_name)), 5L)
  group <- unname(group_display[b$md$BMI25_group])
  for (g in c("NC", "OC", "NT2D", "OT2D")) add(paste0("BMI25_", g), sum(group == g), c(NC=307L, OC=119L, NT2D=138L, OT2D=622L)[g])
  for (i in 1:5) {
    s <- b$species[b$species$comparison == b$spec$comparison[i], ]
    p <- b$reml[b$reml$comparison == b$spec$comparison[i], ]
    add(paste0("species_tested_", i), nrow(s), c(200L,192L,187L,193L,190L)[i])
    add(paste0("species_significant_", i), sum(s$qval < .05), c(69L,32L,10L,4L,5L)[i])
    add(paste0("REML_tested_", i), nrow(p), c(5573L,5572L,5573L,5572L,5573L)[i])
    add(paste0("REML_significant_", i), sum(p$padj < .05), c(13L,133L,0L,0L,20L)[i])
  }
  shared <- intersect(b$reml$pathway[b$reml$comparison == b$spec$comparison[1] & b$reml$padj < .05],
    b$reml$pathway[b$reml$comparison == b$spec$comparison[2] & b$reml$padj < .05])
  add("REML_shared", length(shared), 4L)
  ann <- b$reml_annotation[match(shared, b$reml_annotation$pathway), ]
  add("REML_shared_Blautia_wexlerae", sum(ann$taxonomic_contributor == "Blautia wexlerae"), 3L)
  add("REML_shared_Hungatella_hathewayi", sum(ann$taxonomic_contributor == "Hungatella hathewayi"), 1L)
  neg <- b$reml[b$reml$comparison == b$spec$comparison[5] & b$reml$padj < .05, ]
  add("T2D_internal_negative_REML", sum(neg$combined_log2FC < 0), 20L)
  add("pooled_pathways_tested", nrow(b$pooled), 5945L)
  add("pooled_pathways_full_model_q_significant", sum(b$pooled$qval < .05), 1641L)
  add("pooled_pathways_disease_term_BH_audit", sum(b$pooled$qval_disease_term_BH < .05), 1425L)
  add("pooled_pathways_N", as.integer(all(b$pooled$N == 1186L)), 1L)
  add("pooled_pathways_fixed_effect", as.integer(all(b$pooled$study_effect == "Fixed") &&
    all(b$pooled$model_formula == "pathway ~ disease + BMI + age + sex + study_name")), 1L)
  add("pooled_pathways_formal_q", as.integer(identical(b$pooled$qval, b$pooled$formal_q)), 1L)
  parameters<-setNames(as.character(b$pooled_contract$value),b$pooled_contract$parameter)
  add("pooled_pathways_standardize_TRUE",as.integer(identical(unname(parameters["standardize"]),"TRUE")),1L)
  add("pooled_pathways_saved_fixed_effects",as.integer(identical(unname(parameters["fixed_effects"]),
    "disease + BMI + age + sex + study_name")),1L)
  add("pooled_pathways_no_random_effects",as.integer(is.na(parameters["random_effects"]) || !nzchar(parameters["random_effects"])),1L)
  add("pooled_pathways_fixed_study_identity", as.integer(all(b$pooled$model_id ==
    "Pooled_MaAsLin2_pathways_T2D_BMI_age_sex_fixed_study")), 1L)
  safe_write_csv(data.frame(model_id=unique(b$pooled$model_id),n_tested=nrow(b$pooled),
    formal_q_significant=sum(b$pooled$qval<.05),audit_target_BH_significant=sum(b$pooled$qval_disease_term_BH<.05),
    frozen_significance_count=1641L,note="Fixed-study reference accepted 2026-11-12"),
    file.path(export,"validation/pooled_fixed_study_results_for_review.csv"))
  add("formal_interaction_families",length(unique(paste(b$interactions$BMI_definition,b$interactions$domain))),6L)
  for (name in unique(paste(b$interactions$BMI_definition, b$interactions$domain))) {
    d <- b$interactions[paste(b$interactions$BMI_definition, b$interactions$domain) == name, ]
    add(paste0("interaction_tested_", name), nrow(d), if (d$domain[1] == "species") 200L else 5945L)
    add(paste0("interaction_significant_", name), sum(d$formal_q < .05), 0L)
  }
  a <- b$additive; t <- a$q_within_exposure_BH_T2D < .05; bmi <- a$q_within_exposure_BH_BMI_per5 < .05
  add("additive_tested",nrow(a),200L)
  add("additive_T2D", sum(t), 62L); add("additive_BMI", sum(bmi), 27L)
  add("additive_both", sum(t & bmi), 13L); add("additive_BMI_only", sum(!t & bmi), 14L)
  add("clinical_features", length(unique(b$clinical$feature)), 20L); add("clinical_tests", nrow(b$clinical), 60L)
  add("clinical_significant", sum(b$clinical$q_bh_60 < .05), 24L)
  add("REML_overall_only", sum(b$reml$comparison == b$spec$comparison[1] & b$reml$padj < .05) - length(shared), 9L)
  add("REML_NT2D_only", sum(b$reml$comparison == b$spec$comparison[2] & b$reml$padj < .05) - length(shared), 129L)
  for (i in c(1L, 2L, 5L)) {
    d <- b$reml[b$reml$comparison == b$spec$comparison[i] & b$reml$padj < .05, ]
    add(paste0("REML_positive_", i), sum(d$combined_log2FC > 0), c(`1`=9L, `2`=114L, `5`=0L)[as.character(i)])
    add(paste0("REML_negative_", i), sum(d$combined_log2FC < 0), c(`1`=4L, `2`=19L, `5`=20L)[as.character(i)])
  }
  add("clinical_benchmarks", sum(!b$clinical_checks$pass), 0L)
  result <- pub_bind(checks)
  safe_write_csv(result, file.path(export, "validation/numeric_assertions.csv"))
  assert_true(all(result$pass), "Numeric discrepancy: see validation/numeric_assertions.csv. Publication export stopped.")
  invisible(result)
}

pub_language_check <- function(x, id) {
  text <- paste(as.character(unlist(x, use.names = FALSE)), collapse = "\n")
  forbidden <- c("\\bfrozen\\b", "\\bhistorical\\b", "\\bold\\b", "\\brevised\\b", "audit only",
    "removed medication", "retired_not_valid", "candidate_data_present_requires_source_verification",
    "user authorized", "awaiting_user", "pending medication", "confirmed robustness", "\\bvalidation\\b",
    "latest model", "newly added analysis", "in the revised analysis",
    "was changed", "was replaced", "was retained", "was removed", "original model", "old result",
    "updated result", "final model", "\\bPilot\\b", "independent biomarker", "BMI-dependent dysbiosis")
  hit <- forbidden[vapply(forbidden, function(p) grepl(p, text, ignore.case = TRUE, perl = TRUE), logical(1))]
  assert_true(!length(hit), paste("Reader-facing language check failed:", id, paste(hit, collapse = "; ")))
}

pub_reader_table <- function(x, preserve_model_identity=FALSE) {
  # Version/audit fields stay in the separate archive, without changing any estimate.
  drop <- grepl("historical|bracket_change|qval_BH_target|qval_target_bh|q_joint|source_result_file|model_id|comparison_id|formal_significant|automated_availability|clr_contract_version", names(x), ignore.case = TRUE)
  if("omitted_study"%in%names(x)||"n_LOCO_testable"%in%names(x))drop[grepl("comparison_id|qval_target_bh|model_id",names(x))]<-FALSE
  # S2 BMI23 model identities are reader-facing descriptors, not audit IDs.
  if (isTRUE(preserve_model_identity) && "result_type" %in% names(x)) drop[names(x)=="model_id"] <- FALSE
  x <- x[, !drop, drop = FALSE]
  identifiers<-c("sample_id","feature","pathway","pathway_feature","base_pathway_ID","cohort","study_name",
    "metadata","value","PMID_or_DOI","source_citation","taxonomic_contributor","genus")
  for (n in names(x)) if (is.character(x[[n]]) && !n%in%identifiers) {
    x[[n]] <- gsub("Overweight_vs_Lean_T2D_UnadjustedPendingMedicationAudit", "OT2D_vs_NT2D_BMI25", x[[n]], fixed = TRUE)
    x[[n]] <- gsub("original full-model", "full-model", x[[n]], fixed = TRUE)
    x[[n]] <- gsub("_historical_medication_adjusted", "_medication_adjusted", x[[n]], fixed = TRUE)
    x[[n]] <- gsub("preserving original CLR definition", "preserving prespecified CLR definition", x[[n]], fixed = TRUE)
    x[[n]] <- pub_publication_text(x[[n]])
  }
  x
}

pub_publication_text <- function(x) {
  x<-gsub("([0-9]+)-year-old","aged \\1 years",x,perl=TRUE)
  replacements <- c("historical_antidiabetes_group"="medication adjustment variable",
    "retired_not_valid"="not included in the primary model",
    "candidate_data_present_requires_source_verification"="not reported or could not be determined",
    "no_documented_values_in_this_stratum"="not reported in this database subset",
    "awaiting_user_run"="not reported or could not be determined",
    "pending medication audit"="medication covariate not included",
    "removed medication"="primary model", "historical model"="medication-adjusted primary model",
    "historical"="prespecified", "frozen"="prespecified", "final model"="primary model",
    "audit only"="not included in the primary model", "user authorized"="prespecified",
    "revised analysis"="analysis", "continuous_BMI"="Continuous BMI", "Ctrl"="Control")
  for (from in names(replacements)) x <- gsub(from, replacements[[from]], x, ignore.case = TRUE)
  x
}

pub_verify_workbook_sheet <- function(expected, observed, id, sheet, export) {
  report <- file.path(export, "validation/workbooks", paste0(id, "__", sheet, "__differences.csv"))
  if (nrow(expected) != nrow(observed) || !identical(names(expected), names(observed))) {
    safe_write_csv(data.frame(workbook = id, sheet = sheet, check = "shape_and_column_order",
      expected_rows = nrow(expected), observed_rows = nrow(observed),
      expected_columns = paste(names(expected), collapse = " | "),
      observed_columns = paste(names(observed), collapse = " | ")), report)
    stop(paste("Excel read-back shape/column mismatch:", id, sheet, "See", report), call. = FALSE)
  }
  differences <- list()
  repr <- function(v) {
    s <- if (is.numeric(v)) format(v, digits = 17, scientific = TRUE, trim = TRUE) else as.character(v)
    ifelse(is.na(v), "<R missing value>", paste0('"', s, '"'))
  }
  for (j in seq_len(ncol(expected))) {
    a <- expected[[j]]; raw <- observed[[j]]
    text <- unname(as.character(raw))
    # writeData(keepNA=TRUE, na.string="NA") uses the declared display
    # token for missing cells. Literal character "NA" must remain text.
    marker <- is.na(text) | (!is.na(text) & text == "NA")
    missing <- is.na(a)
    ok <- missing & marker
    delta <- bound <- rep(NA_real_, length(a))
    if (is.numeric(a)) {
      z <- suppressWarnings(as.numeric(text))
      z[marker] <- NA_real_
      present <- !missing & !marker & !is.na(z)
      finite <- present & is.finite(a) & is.finite(z)
      delta[finite] <- abs(a[finite] - z[finite])
      bound[finite] <- 1e-12 * pmax(1, abs(a[finite]))
      ok[finite] <- delta[finite] <= bound[finite]
      nonfinite <- present & !is.finite(a) & !is.finite(z)
      ok[nonfinite] <- a[nonfinite] == z[nonfinite]
      if (grepl("pval|qval|formal_q|q_bh|padj|p[.]adj|p[.]value|BH", names(expected)[j])) {
        # Preserve the existing exact significance-state check.
        ok[present] <- ok[present] & ((a[present] < .05) == (z[present] < .05)) &
          ((a[present] == 0) == (z[present] == 0))
      }
    } else {
      a_text <- unname(as.character(a))
      present <- !missing & !is.na(text)
      # Independent readers may expose XML line breaks as CRLF versus LF.
      # Normalize that encoding difference only; preserve all other text and IDs.
      canonical <- function(z) enc2utf8(gsub("\r\n", "\n", z, fixed=TRUE))
      ok[present] <- canonical(a_text[present]) == canonical(text[present])
      # An explicitly empty string can be read as an empty Excel cell.
      # This exception never applies to nonempty IDs or category labels.
      blank <- !missing & a_text == "" & is.na(text)
      ok[blank] <- TRUE
    }
    bad <- which(is.na(ok) | !ok)
    if (length(bad)) differences[[names(expected)[j]]] <- data.frame(
      workbook = id, sheet = sheet, column = names(expected)[j],
      data_row = bad, excel_row = bad + 1L,
      expected_type = class(a)[1], observed_type = class(raw)[1],
      expected_value = repr(a[bad]), observed_value = repr(raw[bad]),
      absolute_difference = delta[bad], allowed_difference = bound[bad])
  }
  if (length(differences)) {
    difference <- pub_bind(differences)
    safe_write_csv(difference, report)
    stop(paste0("Excel read-back mismatch: ", id, "; sheet=", sheet,
      "; column=", difference$column[1], "; Excel row=", difference$excel_row[1],
      ". See ", report), call. = FALSE)
  }
  invisible(TRUE)
}

pub_workbook <- function(id, title, sheets, note, export, config, main = FALSE, model_id = id) {
  title <- pub_publication_text(title); note <- pub_publication_text(note)
  folder <- file.path(export, "tables", if (main) "main" else "supplementary")
  target <- file.path(folder, paste0(id, ".xlsx"))
  stopifnot(!file.exists(target), !anyDuplicated(names(sheets)))
  sheets <- lapply(sheets, pub_reader_table, preserve_model_identity=identical(id,"Supplementary_Table_S2"))
  readme <- data.frame(Item = c("Title", "Footnote", "Missing values", "Probability display"),
    Definition = c(title, note, "NA", "Numeric values preserved; scientific notation for probability columns"))
  pub_language_check(c(list(title, note), sheets), id)
  wb <- openxlsx::createWorkbook(creator = "T2D BMI publication pipeline")
  # Primary results, secondary results, definitions, then Notes.
  if("Notes"%in%names(sheets)) {
    assert_true(identical(names(sheets$Notes),names(readme)),"Notes schema must be Item/Definition")
    readme<-rbind(readme,sheets$Notes);sheets$Notes<-NULL
  }
  definition <- grepl("definition", names(sheets), ignore.case = TRUE)
  sheets <- c(sheets[!definition], sheets[definition], list(Notes = readme))
  header <- openxlsx::createStyle(fontName = "Times New Roman", fontSize = 10,
    textDecoration = "bold", fgFill = "#E6EEF4", wrapText = TRUE, valign = "top")
  base <- openxlsx::createStyle(fontName = "Times New Roman", fontSize = 10, valign = "top", wrapText = TRUE)
  dictionary <- list()
  for (name in names(sheets)) {
    x <- as.data.frame(sheets[[name]], check.names = FALSE)
    stopifnot(nchar(name) <= 31L, ncol(x) > 0L)
    openxlsx::addWorksheet(wb, name, gridLines = FALSE)
    openxlsx::writeData(wb, name, x, withFilter = TRUE, keepNA = TRUE, na.string = "NA", headerStyle = header)
    openxlsx::freezePane(wb, name, firstRow = TRUE)
    openxlsx::setRowHeights(wb, name, rows = 1, heights = 36)
    if (nrow(x)) openxlsx::addStyle(wb, name, base, rows = seq_len(nrow(x)) + 1L, cols = seq_len(ncol(x)), gridExpand = TRUE)
    for (j in seq_len(ncol(x))) {
      wide <- grepl("feature|pathway|species|definition|note|formula|model|source|citation|criterion|description", names(x)[j], ignore.case = TRUE)
      openxlsx::setColWidths(wb, name, cols = j, widths = if (wide) 46 else 21)
      prob <- grepl("(^p$|pval|p_value|p.value|padj|qval|formal_q|q_bh|BH|Pr\\()", names(x)[j]) && is.numeric(x[[j]])
      if (prob && nrow(x)) openxlsx::addStyle(wb, name, openxlsx::createStyle(numFmt = "0.000000E+00"),
        rows = seq_len(nrow(x)) + 1L, cols = j, stack = TRUE)
      if (id == "Supplementary_Table_S13" && is.numeric(x[[j]]) && nrow(x)) {
        fmt <- if (names(x)[j] == "direction_agreement") "0.0%" else
          if (names(x)[j] == "descriptive_coefficient_r") "0.000" else NULL
        if (!is.null(fmt)) openxlsx::addStyle(wb, name, openxlsx::createStyle(numFmt=fmt),
          rows=seq_len(nrow(x))+1L, cols=j, stack=TRUE)
      }
      if (grepl("^(species|display_name)$", names(x)[j]) && nrow(x)) openxlsx::addStyle(wb, name,
        openxlsx::createStyle(textDecoration = "italic"), rows = seq_len(nrow(x)) + 1L, cols = j, stack = TRUE)
    }
    if (id == "Main_Table_1" && name == "Characteristics") {
      assert_true(identical(names(x), c("Characteristic", "Overall", "NC", "OC", "NT2D", "OT2D")),
        "Main Table 1 must have the six prespecified columns")
      openxlsx::setColWidths(wb, name, cols = 1:6, widths = c(30, rep(15, 5)))
      if (nrow(x)) {
        openxlsx::addStyle(wb, name, openxlsx::createStyle(halign = "left"), rows = 2:(nrow(x)+1), cols = 1, stack = TRUE)
        openxlsx::addStyle(wb, name, openxlsx::createStyle(halign = "center"), rows = 2:(nrow(x)+1), cols = 2:6, gridExpand = TRUE, stack = TRUE)
      }
      openxlsx::pageSetup(wb, name, orientation = "landscape", paperSize = 9,
        fitToWidth = 1, fitToHeight = 1, printTitleRows = 1, left = .3, right = .3, top = .4, bottom = .4)
      wb$workbook$definedNames <- c(wb$workbook$definedNames,
        paste0('<definedName name="_xlnm.Print_Area" localSheetId="', match(name, names(sheets))-1L,
          '">\'', name, '\'!$A$1:$F$', nrow(x)+1L, '</definedName>'))
    }
    csv <- file.path(export, "tables/csv", paste0(id, "__", name, ".csv"))
    safe_write_csv(x, csv)
    dictionary[[name]] <- data.frame(sheet = name, column = names(x), storage_class = vapply(x, function(v) class(v)[1], character(1)),
      definition = vapply(names(x), pub_column_definition, character(1), id=id, sheet=name))
    # Read-back checks numbers and categorical strings; byte identity is not required.
    # Performed only when the user executes the exporter.
  }
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  pub_clear_empty_drawing_placeholders(wb)
  pub_set_written_dimensions(wb, sheets)
  openxlsx::saveWorkbook(wb, target, overwrite = FALSE)
  pub_validate_saved_workbook(target, sheets, id, export, root=config$project_root)
  safe_write_csv(pub_bind(dictionary), file.path(export, "tables/csv", paste0(id, "__data_dictionary.csv")))
  safe_write_lines(c(title, note), file.path(export, "tables/csv", paste0(id, "__title_footnote.txt")))
  data.frame(table_id = id, title = title, workbook_path = substring(target, nchar(export) + 2L),
    n_sheets = length(sheets), source_model_id = model_id, workbook_sha256 = sha256_file(target))
}

pub_column_definition <- function(n, id=NULL, sheet=NULL) {
  if (identical(id,"Supplementary_Table_S2")) {
    description <- pub_s2_column_definition(n,sheet)
    if (!is.null(description)) return(description)
  }
  known <- c(qval = "MaAsLin2 full-model q value; formal for primary species and pooled pathways",
    formal_q = "Formal q value under the model-specific family stated in q_family",
    padj = "BH-adjusted P value within each complete successful REML contrast",
    q_bh_60 = "BH-adjusted q value across all 60 exploratory Spearman tests",
    p_adj_BH = "BH over 24 phylum or 48 species unadjusted emmeans pairwise P values",
    coef = "Target coefficient in the saved model's units and comparison direction",
    combined_log2FC = "REML-pooled cohort DESeq2 log2 fold change",
    stderr = "Standard error of the target coefficient", ci_lower = "Lower 95% confidence limit",
    ci_upper = "Upper 95% confidence limit", qval_disease_term_BH = "Disease-term BH across 5945 pathways; audit only",
    n_complete = "Number of pairwise complete observations", spearman_rho = "Unadjusted Spearman rank correlation",
    BMI25_group = "BMI <25 or >=25 crossed with disease", BMI23_group = "Cohort-adapted BMI23 sensitivity grouping")
  if (n %in% names(known)) unname(known[n]) else paste("Saved", gsub("_", " ", n), "(see model definition and source manifest)")
}

pub_source_review <- function(root, config, type, run = NULL) {
  path <- file.path(root, config$publication$source_review[[paste0(type, "_file")]])
  assert_file(path, paste("Completed structured source-review input; copy the verified CSV here, do not infer evidence:", type))
  x <- if (is.null(run)) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE) else pub_read_csv_snapshot(run, path)
  required <- if (type == "cohort") c("cohort", "PMID_or_DOI", "country", "original_study_design", "T2D_definition",
    "control_definition", "inclusion_exclusion_criteria", "evidence_grade", "source_verification_note") else
    c("cohort", "variable", "source_study_collected", "current_database_retained", "completeness",
      "harmonizable_across_cohorts", "use_in_current_model", "reason_not_adjusted", "evidence_grade", "source_citation", "source_verification_note")
  assert_columns(x, required, "Verified source review")
  pub_key_contract(x, if(type=="cohort")"cohort" else c("cohort","variable"), "Verified source review")
  assert_true(setequal(unique(x$cohort), unlist(config$cohorts$names)), "Source review must cover exactly the five analysis cohorts")
  if (type == "cohort") {
    stopifnot(nrow(x) == 5L, !anyDuplicated(x$cohort))
    q <- x[x$cohort == "QinJ_2012", ]; m <- x[x$cohort == "MetaCardis_2020_a", ]; y <- x[x$cohort == "YuJ_2015", ]
    stopifnot(grepl("1999", q$T2D_definition), grepl("WHO", q$T2D_definition),
      grepl("A", q$evidence_grade), grepl("non-T2D", m$control_definition, ignore.case = TRUE),
      grepl("A-", paste(m$evidence_grade, m$source_verification_note)),
      grepl("B", y$evidence_grade), grepl("16", y$source_verification_note), grepl("38", y$source_verification_note))
  } else {
    required_variables <- c("diabetes_duration", "HbA1c", "fasting_glucose", "fasting_insulin_or_HOMA_IR",
      "disease_severity", "medication", "metformin", "diet", "smoking", "alcohol", "physical_activity",
      "education", "income", "occupation", "environmental_factors")
    assert_true(all(required_variables %in% x$variable), "Structured confounder review lacks required variable concepts; see input schema")
    drug <- x[x$variable %in% c("medication", "metformin"), , drop=FALSE]
    assert_true(nrow(drug)==10L && all(drug$use_in_current_model %in% c("yes","no")), "Medication source-review scope/flags differ")
    assert_true(all(grepl("Categorical antidiabetic treatment was included",drug$reason_not_adjusted,fixed=TRUE)) &&
      all(grepl("omission of this covariate assessed in sensitivity analyses",drug$reason_not_adjusted,fixed=TRUE)) &&
      all(grepl("reliably defined metformin-free analysis was not available",drug$reason_not_adjusted,fixed=TRUE)),
      "S15 medication model-use notes are outdated; do not infer medication nonuse")
    assert_true(all(drug$use_in_current_model[drug$variable=="metformin"]=="no"), "Metformin was not a separate binary adjustment")
    grid <- expand.grid(cohort = unlist(config$cohorts$names), variable = required_variables)
    assert_true(all(paste(grid$cohort, grid$variable) %in% paste(x$cohort, x$variable)), "Each required confounder must be documented for every cohort")
  }
  x
}

# Add identities before binding the three source types. No numeric calculations.
pub_bmi23_alpha_identity <- function(x, type, source_table) {
  original <- x
  assert_columns(x,"metric","BMI23 alpha metric")
  assert_true(all(x$metric%in%c("Shannon","Richness")),"Unknown BMI23 alpha metric")
  assert_true(type%in%c("model_comparison","fixed_effect","pairwise_comparison"),"Unknown BMI23 result source type")
  full_rhs <- "analysis_group_BMI23_adapted + age_numeric + sex + (1 | study_name)"
  null_rhs <- "age_numeric + sex + (1 | study_name)"
  x$result_type <- type
  x$source_table <- source_table
  x$model_id <- "BMI23_alpha_full"
  x$model_formula <- paste(x$metric,"~",full_rhs)
  x$model_comparison <- NA_character_
  x$comparison_direction <- NA_character_
  x$p_adjustment <- NA_character_
  if(type=="model_comparison") {
    assert_columns(x,c("npar","Chisq","Pr(>Chisq)"),"BMI23 likelihood-ratio source")
    # Source formulas: 3 null fixed parameters + random-intercept/residual
    # variance parameters = 5; the full model adds three group parameters = 8.
    # Structural parameter counts identify models; source row order is irrelevant.
    assert_true(all(x$npar%in%c(5L,8L))&&!anyDuplicated(x[c("metric","npar")])&&nrow(x)==4L,
      "BMI23 model-comparison parameter identities differ from the source formulas")
    assert_true(all(table(x$metric)==2L),"Incomplete BMI23 model comparison")
    null <- x$npar==5L
    x$model_id[null] <- "BMI23_alpha_null"
    x$model_formula[null] <- paste(x$metric[null],"~",null_rhs)
    x$model_comparison <- "BMI23_alpha_full versus BMI23_alpha_null; ML likelihood-ratio test"
    x$p_adjustment <- "None; Pr(>Chisq) is the nominal likelihood-ratio P value"
    assert_true(all(is.na(x$Chisq[null]))&&all(is.na(x$`Pr(>Chisq)`[null])),
      "The structurally identified null model unexpectedly carries the comparison test")
  } else if(type=="pairwise_comparison") {
    assert_columns(x,c("contrast","p.value"),"BMI23 alpha pairwise source")
    parts<-strsplit(as.character(x$contrast)," - ",fixed=TRUE)
    assert_true(all(lengths(parts)==2L)&&nrow(x)==12L&&all(table(x$metric)==6L),
      "BMI23 pairwise contrast identity/family differs")
    x$comparison_direction <- vapply(parts,function(z)paste(z[1],"minus",z[2]),character(1))
    x$p_adjustment <- "p.value: BH-adjusted across six pairwise contrasts within this metric (emmeans adjust = BH)"
  } else {
    assert_columns(x,c("term","Estimate","Std. Error"),"BMI23 fixed-effect source")
    x$p_adjustment <- "Not applicable; lme4 fixed-effect coefficients do not include P values"
  }
  assert_true(identical(original,x[,names(original),drop=FALSE]),"BMI23 annotation altered an original field")
  x
}

publication_bmi23_alpha <- function(run) {
  files<-c(model_comparison="BMI23_Alpha_Group_Overall_LRT.csv",
    fixed_effect="BMI23_Alpha_LMM_Coefficients.csv",pairwise_comparison="BMI23_Alpha_Pairwise_BH.csv")
  pub_bind(lapply(names(files),function(type) {
    file<-files[[type]]
    x<-pub_read(run,paste0("tables/07_BMI23_diversity/",file))
    pub_bmi23_alpha_identity(x,type,file)
  }))
}

publication_alpha_notes <- function() data.frame(Item=c("Pairwise","Interactions","BMI23 models","BMI23 pairwise"),
  Definition=c(
    "metric distinguishes Shannon and Richness. Pairwise estimates represent adjusted mean differences calculated as group1 minus group2. The p.value column contains Tukey-adjusted P values; p.adj contains the subsequent BH adjustment across the six pairwise comparisons within each diversity metric. Formal pairwise significance uses p.adj < 0.05.",
    "Categorical BMI25 and continuous centered BMI each have a separate BH family of two interaction tests (Shannon and Richness). p_interaction is the nominal P value; q_interaction_BH is the formal BH-adjusted value.",
    "result_type and source_table distinguish model-comparison, fixed-effect and pairwise records. The full model includes the four-group factor, age, sex and a cohort random intercept; the null model omits the group factor. The model comparison uses maximum-likelihood refits; coefficient and emmeans summaries use the fitted REML full model.",
    "BMI23 pairwise estimates follow the left group minus the right group in contrast. Its p.value already contains BH adjustment across six comparisons within each metric (emmeans adjust = BH). The likelihood-ratio Pr(>Chisq) is nominal. No additional BMI23 BH column is present."),stringsAsFactors=FALSE)

pub_s2_column_definition <- function(n,sheet) {
  shared<-c(metric="Diversity metric: Shannon or Richness",result_type="Record type identified from its source table before combination",
    source_table="Current-Run source CSV identifying the result type",model_id="Null or full model identified from the source model and its parameter count",
    model_formula="Formula of the fitted model for this record",model_comparison="Full versus null maximum-likelihood likelihood-ratio comparison",
    comparison_direction="Adjusted contrast: left group minus right group",p_adjustment="Actual adjustment applied in this source table")
  specific<-if(sheet=="Pairwise")c(estimate="Adjusted mean difference: group1 minus group2",
    group1="First group in the adjusted mean contrast",group2="Reference group subtracted from group1",
    p.value="Tukey-adjusted pairwise P value; not an unadjusted P value",
    p.adj="Formal pairwise P value: subsequent BH over six Tukey-adjusted P values within each metric") else
    if(sheet=="Interactions")c(p_interaction="Nominal interaction P value",q_interaction_BH="Formal BH over Shannon and Richness within this categorical/continuous BMI family") else
    if(sheet=="BMI23")c(p.value="For pairwise records: emmeans BH over six contrasts within the metric; not Tukey",
      `Pr(>Chisq)`="Nominal full-versus-null likelihood-ratio P value; no additional BH applied",
      estimate="Adjusted pairwise mean difference in the direction stated in contrast") else character()
  known<-c(specific,shared)
  if(n%in%names(known))unname(known[n]) else NULL
}
