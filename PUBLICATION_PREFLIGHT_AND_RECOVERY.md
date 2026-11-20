# Complete preflight and recovery of an interrupted publication tail

An export is not evidence that the current Run is complete. A later
`--from clr` invocation withdrew downstream completion registrations and was
interrupted during `raw_beta`. Removing its orphaned lock could not restore the
missing completion registrations. The exporter now checks the entire state,
every required checkpoint, artifact hashes, environment, source-review inputs,
saved numerical/model contracts and all16 table source plans before creating an
Export directory. Failures are listed together in the console. No repeated empty
`PublicationExport_*` folders are created by a preflight failure.

Exports now read saved composition data directly. They never invoke
`run_publication_updates.R` or the pipeline runner. Preparing or recovering a Run
and exporting it are separate operations.

## This interrupted Run

Run: `Run_20261119_143714`. Explicit earlier complete snapshot:
`checkpoints/State_000068.yml`. The current snapshot remains authoritative; it is
not rolled back to the earlier one.

Read-only planning verified:

- `raw_beta`, `abundant_taxa`, `publication_display`: their explicitly registered
  a002 checkpoints have matching current code/input/parameter signatures and
  complete unchanged artifacts. They can be registered again without model fits.
- `medication_sensitivity`: its stored signature differs because it includes
  upstream checkpoint identities. It must not be rebound with that mismatched
  signature. Its four saved primary/sensitivity input copies reproduce exactly,
  including the source implementation's feature alignment. An explicitly
  requested refresh calls only the existing saved-result comparison function,
  not the four underlying models. P values, q values and coefficients are copied,
  not calculated. All16 refreshed CSVs must be byte-identical to the verified
  baseline before the new checkpoint/state is committed.

`reconcile_publication_checkpoints.R` defaults to PLAN ONLY. Apply requires the
explicit flags below. It checks the frozen implementation, preserves every
current completed reference (including CLR a003), selects candidates only from
the specified prior state, verifies signatures and every artifact, and appends a
new state and provenance event. It never chooses the largest Attempt directory.
Any other signature mismatch blocks recovery. Historical states, outputs, failed
Exports and author review tables are not overwritten.

## Author-run R commands

Exit `Browse[1]>` with `Q` before running commands. Keep the same frozen R
environment and `T2DBMI_VALIDATION_PYTHON` setting described in
`WORKBOOK_PATCH_RUN_GUIDE.md`.

```r
project_dir <- "Y:/研究生相关课题/分析类/胖瘦DM/T2D_BMI_Publication_RC4"
python <- file.path(Sys.getenv("USERPROFILE"),
  ".cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe")
stopifnot(file.exists(python)) # On another computer, supply its own interpreter.
Sys.setenv(T2DBMI_VALIDATION_PYTHON=python)
rscript <- file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe" else "Rscript")
rc4 <- function(script,args=character()) {
  code <- system2(rscript,c("--vanilla",shQuote(file.path(project_dir,script)),shQuote(args)))
  if(code!=0L)stop("Step incomplete. Read all printed blockers; do not start export.",call.=FALSE)
  invisible(code)
}

# One explicit recovery; no statistical model is fitted.
rc4("reconcile_publication_checkpoints.R",c(
  "Run_20261119_143714","State_000068.yml",
  "--apply","--refresh-medication-summary"
))
# The recovery already ends with a complete read-only export preflight.
# An independent repeat is available, without generating any Export directory:
# rc4("preflight_publication.R",c("Run_20261119_143714","all"))

# Only after recovery and preflight both succeed:
before <- list.dirs(file.path(project_dir,"output"),recursive=FALSE,full.names=FALSE)
rc4("export_publication.R",c("Run_20261119_143714","all"))
export_id <- setdiff(list.dirs(file.path(project_dir,"output"),recursive=FALSE,full.names=FALSE),before)
export_id <- export_id[grepl("^PublicationExport_",export_id)]
stopifnot(length(export_id)==1L)
message("New complete export: ",file.path(project_dir,"output",export_id))
```

Recovery creates no new top-level Run or Export folder. It adds only a new
medication-summary attempt, a reconciliation event and a checkpoint state inside
the selected Run. The expensive existing raw-Aitchison and abundant-taxa results
are reused only after exact signature and integrity verification. If a source or
signature changed, stop and inspect the evidence; never relax the gate.

After a successful export, prepare a new review using the actual `export_id` and
a new external review directory as described in `WORKBOOK_PATCH_RUN_GUIDE.md`.
All author decisions remain PENDING. CLR reference freeze, REML review, visual
review and new-computer testing remain separate release conditions.
