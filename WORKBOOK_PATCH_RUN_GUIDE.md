# Workbook methods, worksheet ranges and independent read-back

This patch changes table annotations, Excel display styles and validation only.
It does not fit models, recalculate probabilities, update CLR references or clear
the REML release gate. Existing Run, Export and review directories stay unchanged.

## Validation dependency

The frozen R environment is unchanged. Workbook acceptance additionally uses a
separate Python interpreter with openpyxl (tested with Python 3.12.14 and openpyxl
3.1.5). It opens workbooks with `read_only=True`, respects the declared sheet
dimensions and checks every expected header and cell. It never repairs XLSX files
or resets dimensions. No Python/R package is installed during export or validation.

Set `T2DBMI_VALIDATION_PYTHON` to the interpreter on the computer running the
commands. If unset, the code checks `python`/`python3` on PATH. A missing dependency
fails in export preflight, before writing workbooks. On a new computer this
validation dependency must be provisioned separately from renv.

## R commands

Use a normal R console. If it shows `Browse[1]>`, enter `Q` first. Replace the
Python path below with the interpreter available on that computer.

```r
project_dir <- "Y:/研究生相关课题/分析类/胖瘦DM/T2D_BMI_Publication_RC4"
Sys.setenv(T2DBMI_VALIDATION_PYTHON = "C:/path/to/python.exe")
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
rc4 <- function(script, args = character()) {
  code <- system2(rscript, c("--vanilla", shQuote(file.path(project_dir, script)), shQuote(args)))
  if (code != 0L) stop("Step incomplete; keep the printed error and report path.", call. = FALSE)
  invisible(code)
}

# Temporary fixtures only; they are deleted after the test.
rc4("tests/test_publication_workbook_contract.R")
rc4("tests/test_review_export_binding.R")
rc4("tests/test_rc4_release_review.R")

# Read existing completed checkpoints; no statistical model is run.
# Use all so the new review has one complete, traceable set of 16 workbooks
# and 40 panels. The figure code and parameters are unchanged by this patch.
before_exports <- list.dirs(file.path(project_dir,"output"), recursive=FALSE, full.names=FALSE)
rc4("export_publication.R", c("Run_20261119_143714", "all"))
export_id <- setdiff(list.dirs(file.path(project_dir,"output"), recursive=FALSE, full.names=FALSE), before_exports)
export_id <- export_id[grepl("^PublicationExport_", export_id)]
stopifnot(length(export_id)==1L)
print(export_id)
```

The exporter prints `output/PublicationExport_<new timestamp>`. It never reuses
an existing Export directory. To inspect only tables first, substitute `tables`
for `all`; prepare the release review from a subsequent complete `all` export.

The commands above select only the directory created by that invocation, without
choosing the newest old directory. Continue in the same R session:

```r
stopifnot(grepl("^PublicationExport_[0-9]{8}_[0-9]{6}[.][0-9]+_[0-9]+$", export_id))
review_dir <- file.path(dirname(project_dir), paste0(
  "RC4_ReleaseReview_", format(Sys.time(), "%Y%m%d_%H%M%S")))
stopifnot(!dir.exists(review_dir))
rc4("review_rc4_release.R", c(
  "prepare", "--run", "Run_20261119_143714",
  "--export", export_id,
  "--diagnostics", "REML_Diagnostics_20261119_171725.766_32664",
  "--baseline-attempt", "a001", "--out", review_dir
))
```

This prepares review materials only. Both new author decision tables start at
`PENDING`. No previous ACCEPT/CONFIRMED decision is inherited. The previous
`RC4_ReleaseReview_20261119` is not changed. Do not run `freeze-clr` or `archive-reml`
until the corresponding author review is complete.

## Outputs to inspect

- `tables/supplementary/Supplementary_Table_S2.xlsx`: original sheet names and
  records; Notes explains BMI25 Tukey then within-metric BH, the two-metric
  interaction families and BMI23's distinct procedures. BMI23 descriptors are
  added separately to each source table before binding. LRT null/full identity
  uses the model formulas' parameter counts (5/8), never a row-number guess.
- `tables/csv/Supplementary_Table_S2__data_dictionary.csv`: sheet-specific meanings
  of nominal/adjusted probability columns and result/model identity columns.
- `tables/supplementary/Supplementary_Table_S13.xlsx`: direction agreement is a
  percentage and descriptive coefficient correlation has three displayed decimal
  places. Underlying Excel cells remain numeric and retain the saved values.
- `validation/workbooks/<workbook>__xml_ranges.csv`: original OOXML declaration,
  occupied cell extent and expected written range, including header-only tables.
- `validation/workbooks/<workbook>__ranges.csv`: adds independent streaming row and
  column counts. Every range and value check must pass.
- `validation/workbooks/<workbook>__streaming.json`: independent reader/version,
  workbook SHA-256, per-sheet cell-type/value/range checks. Text checks normalize
  only the XML reader's CRLF/LF representation; numeric tolerance is unchanged.
- `validation/workbook_ooxml_integrity.csv`: aggregated workbook acceptance.
- `export_status.csv`, `export_gaps.csv`: overall export status and unresolved gaps.
- The new review's `CLR_export_checkpoint_binding.csv`: current Attempt comes
  from completed checkpoint state. An older export retains its actual a002 source;
  a new export records whichever Attempt it really read. Across Attempts, complete
  file coverage, identical model signatures and byte-identical relevant files are
  all required. Matching significant counts alone is insufficient.

Inspect SFig4B/SFig4C in the new export as well. Their existing styles, including
SFig4C's annotation layout, are unchanged. This patch does not complete the CLR
reference freeze, REML author review, manuscript figure review or new-computer
empty-cache test. Those remain separate release conditions.
