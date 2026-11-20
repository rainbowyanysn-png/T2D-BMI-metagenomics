# User-executable export-only regression checks. No model fitting or real-data export.
local({
  entry <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  stopifnot(length(entry) == 1L)
  if (!grepl("^([A-Za-z]:|/)", entry)) entry <- file.path(getwd(), entry)
  root <- dirname(dirname(entry))
  source(file.path(root, "R/functions/project_paths.R"), local = environment())
  activate_existing_project_library(root)
  for (file in c("assertions", "source_resources", "model_helpers", "publication_contracts",
    "publication_display_data", "publication_table_helpers"))
    source(file.path(root, "R/functions", paste0(file, ".R")), local = environment())
  require_project_packages(c("dplyr", "openxlsx"))
  fails <- function(action) inherits(tryCatch({action(); NULL}, error = function(e) e), "error")
  expect <- function(value, label) assert_true(value, paste("Export regression failed:", label))

  # A named colour vector with the wrong semantic keys must fail before plotting.
  wrong <- c(Negative = "#2166AC", "#9985B5", Positive = "#B2182B")
  groups <- c("Overall only", "Shared", "NT2D vs NC only")
  expect(fails(function() pub_palette(wrong, groups, "fixture wrong palette")), "wrong palette rejected")
  right <- setNames(unname(wrong), groups)
  expect(identical(pub_palette(right, groups, "fixture correct palette"), right), "explicit palette preserved")
  expect(fails(function() pub_palette(right, c(groups, NA_character_), "missing category")), "missing category rejected")

  # Recorded rounded taxa can sum to 100.00006. Preserve them exactly and sum
  # omitted taxa directly; never loosen a 100-complement threshold or reclose.
  phyla <- c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria")
  values <- c(2.14546, 37.58639, 46.63219, 13.63602)
  full <- matrix(values, ncol = 1L, dimnames = list(phyla, "sample_A"))
  md <- data.frame(sample_id = "sample_A", BMI25_group = "Non_Overweight_Control")
  saved <- data.frame(taxonomic_level = "phylum", sample_id = "sample_A", feature = phyla,
    analysis_group = "Non_Overweight_Control", raw_abundance = values)
  selected <- data.frame(feature = phyla)
  before <- digest::digest(saved, algo = "sha256")
  prepared <- publication_phylum_display(full, saved, md, selected, 1e-8)
  expect(identical(prepared$display$raw_abundance[1:4], values), "rounded input values unchanged")
  expect(prepared$display$raw_abundance[5] == 0, "no omitted phyla means actual other abundance zero")
  expect(prepared$totals$full_phylum_total > 100, "total not forced to 100")
  expect(identical(before, digest::digest(saved, algo = "sha256")), "source data unchanged")
  full2 <- rbind(full, Unclassified = .003)
  prepared2 <- publication_phylum_display(full2, saved, md, selected, 1e-8)
  expect(prepared2$display$raw_abundance[5] == .003, "actual omitted abundance retained")
  changed <- saved; changed$raw_abundance[1] <- changed$raw_abundance[1] + .01
  expect(fails(function() publication_phylum_display(full, changed, md, selected, 1e-8)), "abundance mismatch rejected")

  # Excel fixtures distinguish literal NA text, declared missing cells, logical
  # values and tiny probabilities. No model results are calculated.
  expected <- data.frame(text = c("NA", NA_character_, "", "NC"),
    qval = c(1e-30, NA_real_, .049, .051), flag = c(TRUE, NA, FALSE, TRUE))
  folder <- tempfile("publication_export_contracts_"); dir.create(folder)
  wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, "fixture")
  openxlsx::writeData(wb, "fixture", expected, keepNA = TRUE, na.string = "NA")
  path <- file.path(folder, "fixture.xlsx"); openxlsx::saveWorkbook(wb, path, overwrite = FALSE)
  back <- openxlsx::read.xlsx(path, sheet = "fixture", check.names = FALSE,
    skipEmptyRows = FALSE, skipEmptyCols = FALSE, na.strings = NULL)
  pub_verify_workbook_sheet(expected, back, "fixture", "fixture", folder)
  corrupted <- back; corrupted$qval[3] <- .051
  expect(fails(function() pub_verify_workbook_sheet(expected, corrupted, "corrupted", "fixture", folder)),
    "changed numeric/significance content rejected")
  zeroed <- back; zeroed$qval[1] <- 0
  expect(fails(function() pub_verify_workbook_sheet(expected, zeroed, "zeroed_probability", "fixture", folder)),
    "tiny nonzero probability cannot become zero")
  expect(pub_warning_is_fatal("No shared levels found between names(values) and data"), "palette warning is fatal")
  expect(pub_warning_is_fatal("Removed 2 rows containing non-finite values"), "dropped rows are fatal")
  message("Export-only regression checks completed. Temporary fixture files: ", folder)
})
