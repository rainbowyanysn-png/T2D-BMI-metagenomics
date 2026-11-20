# RC4 submission-output revision: RC4_submission16_20261121

Only submission output and verification changed. Frozen statistical files and
results, CLR reference, REML archive, ethics statements and all40 panel functions,
selection criteria, styles and dimensions are preserved. No analysis is required
to apply this patch. Actual new export and author review remain pending.

## Table identity and manuscript references

| Previous content | Current submission |
|---|---|
| Two main tables | Unchanged |
| S1-S13 | Unchanged content identities and numbering; S7 display metadata corrected |
| S14 Medication-adjustment sensitivity | Internal archive only; remove manuscript/legend/response-letter citations |
| S15 Clinical and confounder availability | S14; update its manuscript, legend and response-letter citations |

The formal set contains16 workbooks: two main plus S1-S14. New S14 has
Source_and_database, Database_completeness and Notes. It cannot contain the
retired Models/Paired_coefficients/Concordance medication workbook. The shared
submission_table_manifest() is the numbering contract; exported indices carry
content_id, with separate formal and retired identity manifests.

Source evidence is read unchanged and scoped medication/metformin wording is
mapped only in the exported view. Collection/retention/completeness/counts,
evidence grades, citations, missing-treatment interpretation, non-binary
metformin limitation and unadjusted pathway-meta/interaction scope are preserved.

S7 retains all68 eligible stratified candidate rows in the current Run, with64
displayed species-stratified contributions (six species, four genera,42 pathway
IDs). Four additional genus-level candidates remain explicitly unshown. Tables
and the existing figure renderers use the same r3_candidate_filter_audit() and
r3_expected_candidate_data() implementation. Source matching is by complete
feature ID. S7 uses current flags and display_selection_rank; old Top-N ranks
remain only in the internal complete-source archive. Genus bars aggregate the64
contributions; they are not64 individual bars. No coefficient or P/q value is
recalculated. Inferred functional potential is not measured flux or causation.

## R: new complete export and submission-only validation

Run from a normal console (exit Browse with Q). No file chooser is used.

```r
local({
  project_dir <- "Y:/研究生相关课题/分析类/胖瘦DM/T2D_BMI_Publication_RC4"
  python <- file.path(Sys.getenv("USERPROFILE"),
    ".cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe")
  stopifnot(file.exists(python)) # Other computers must supply their own Python.
  Sys.setenv(T2DBMI_VALIDATION_PYTHON=python)
  rscript <- file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe" else "Rscript")
  rc4 <- function(script,args=character()) {
    code <- system2(rscript,c("--vanilla",shQuote(file.path(project_dir,script)),shQuote(args)))
    if(code!=0L)stop("Step incomplete; later steps stopped. Keep the printed report.",call.=FALSE)
  }
  before <- list.dirs(file.path(project_dir,"output"),recursive=FALSE,full.names=FALSE)
  rc4("export_publication.R",c("Run_20261119_143714","all"))
  export_id <- setdiff(list.dirs(file.path(project_dir,"output"),recursive=FALSE,full.names=FALSE),before)
  export_id <- export_id[grepl("^PublicationExport_",export_id)]
  stopifnot(length(export_id)==1L)
  message("New export: ",file.path(project_dir,"output",export_id))
  rc4("validate_submission_output.R",c("Run_20261119_143714",export_id))
})
```

This creates one new PublicationExport directory, then a SubmissionValidation
directory. The latter contains S7_panel_source_sets.csv, workbook XML/range and
streaming/type/value reports, export_source_mapping.csv and status.txt.
The source Run and historical Exports remain unchanged.

For tables only, replace all with tables. Inside the same local block above,
replace the validator call by the following; it is not a standalone command:

```r
# After a successful tables-only export and explicit selection of its new ID:
rc4("validate_submission_output.R",c("Run_20261119_143714",export_id,
    "PublicationExport_20261120_014747.282_2096"))
```

The mapping records the two different Export IDs. This is not a new full visual
acceptance. For a complete `all` Export only, the following can be added inside that same
local block after export; it is not a standalone command:

```r
rc4("validate_release.R",c("--run","Run_20261119_143714","--export",export_id))
```

Do not pass a tables-only Export to that single-Export full validator; use the
specialized validator with the explicit figure-source mapping instead.

The broader validator still reports unresolved REML warnings and unexecuted
independent comparisons honestly; it may return a nonzero status for those
conditions. Submission-only PASS does not clear those gates. Do not compare this
source Run to the CLR reference frozen from itself as independent validation.

## Candidate transfer and GitHub steps (author executes)

```r
local({
project_dir <- "Y:/研究生相关课题/分析类/胖瘦DM/T2D_BMI_Publication_RC4"
rscript <- file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe" else "Rscript")
rc4 <- function(script,args=character()) {
  code <- system2(rscript,c("--vanilla",shQuote(file.path(project_dir,script)),shQuote(args)))
  if(code!=0L)stop("Candidate transfer incomplete.",call.=FALSE)
}
rc4("build_public_package.R","--manifest")
destination <- file.path(dirname(project_dir),paste0(
  "T2D_BMI_RC4_Submission16_",format(Sys.time(),"%Y%m%d_%H%M%S")))
rc4("prepare_transfer.R",c("--zip","--dest",destination))
message("New candidate: ",destination,".zip")
})
```

Keep the old test package. The new timestamped folder/ZIP carries this revision
note and RC4_VERSION.yml; output/data_derived start empty. An old-package model
test may continue, but its17-workbook code and this16-workbook output revision
are distinct acceptance runs. Record code-manifest hashes and Run/Export IDs.

For GitHub, review changed_files.csv/diffs from this patch audit, update the
listed public code/configuration/documentation and regenerated inventory, then
make an author-reviewed commit/tag. Do not add output, derived caches, external
CLR references, author decision forms or private audits. No commit, push,
repository publication or production admission was performed by this patch.

Repository tests include retired/misnamed medication workbook, wrong content ID,
legacy Top20/30 flags, missing/duplicate/wrong display IDs and altered scientific
values. XML ranges, independent openpyxl read_only streaming and original numeric
read-back checks remain active. Extraction-only empty-output tests must not be
run against this used working project; they apply to a new transfer candidate.
