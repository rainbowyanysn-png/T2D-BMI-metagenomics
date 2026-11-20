# RC4 publication candidate: run guide

Status: candidate awaiting complete independent run and author review. Extraction/static PASS is not statistical reproduction, visual acceptance or production approval. All admission flags retain their actual pending/false values. No model or image was run during extraction.

## Environment preparation (separate from analysis)
Use R 4.5.1. Install the renv management tool if absent, then restore the unmodified renv.lock. Restore is an explicit, network-enabled preparation step; analysis never installs/updates packages. Frozen analysis packages include curatedMetagenomicData 3.16.1 and Maaslin2 1.22.0; other versions are fixed by renv.lock and config/package_versions.csv. Supply Times New Roman and required system libraries (e.g. Cairo/font support) according to your OS; do not silently substitute a font.

From R, set project_dir to this copied RC4 folder on your own computer (any path/drive; no old checkout required):

```r
# Edit this value once; do not choose an earlier project.
project_dir <- normalizePath("path/to/T2D_BMI_Publication_RC4", winslash="/", mustWork=TRUE)
if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv")
renv::restore(project=project_dir, prompt=FALSE)
exe <- file.path(R.home("bin"), if (.Platform$OS.type=="windows") "Rscript.exe" else "Rscript")
rc4 <- function(script, args=character()) {
  code <- system2(exe, c("--vanilla", shQuote(file.path(project_dir,script)), shQuote(args)))
  if (code!=0L) stop("Step incomplete; retain console errors and printed log path.",call.=FALSE)
  invisible(code)
}
rc4("tests/static_rc4_candidate.R")
# Author-run full computation and all publication exports, ONLY when ready:
rc4("run_all.R", "--candidate")
```

For command-line users, every entry supports its absolute script path from an unrelated working directory. The root comes from the script location; working directory objects are not used. No .RData workspace is loaded.

## Recovery, export and validation

```r
# Set identifiers to this RC4 run/export, printed by the pipeline.
run_id <- "Run_YYYYMMDD_HHMMSS"
export_id <- "PublicationExport_YYYYMMDD_HHMMSS.mmm_PID"
rc4("run_all_development.R", c("--resume",run_id))
# After resumed analysis completes (or to repeat exports without refitting):
rc4("export_publication.R", c(run_id,"all"))
# Or c(run_id,"tables") / c(run_id,"figures"). Each creates a NEW export folder.
rc4("validate_release.R", c("--run",run_id,"--export",export_id))
# Optional reference-root is supplied ONLY to that independent validator.
# Without external references, their absence must remain explicitly unverified.
```

A new Run writes RC4_RUN_IDENTITY.txt. Creation, resume, export and saved-result validation require RC4 identity plus checkpoint provenance; RC3/V3 Runs are rejected. Use only a Run created in this RC4 checkout. The legacy clean_rc3_run_pending admission key is retained as historical schema compatibility; clean_rc4_run_pending controls RC4's production gate. No pending flag was cleared.

## Raw resources and first use without caches

The package contains licensed source documentation/mappings, not previous Run results or fitted input matrices. ExperimentHub resource IDs/titles/versions/content SHA-256 are fixed by data_manifest.csv and config/source_resources.csv. Current-run acquisition uses those exact resources, verifies hashes and logs source/acquisition status. Missing cache is allowed to trigger exact-ID downloads; failed network/resource/version/hash checks stop without replacements. The unmodified source_resources.R also checks packaged provider input-object hashes.

A computer with no cache needs Internet both for environment restore and first raw-data acquisition. No old Run, derived RDS, plot-data, CleanTest, V3, FigureRefinement or IntegrationPreview is needed. For genuine no-cache acceptance, use a fresh OS user/VM and an RC4 copy with empty output/data_derived. Restore packages, run static checks and then --candidate; retain source_resource_verification.csv and download logs. Do not copy caches or intermediate data from another project. This no-cache statistical test has NOT been executed during extraction.

## Retained analyses and outputs

Input module plus 22 registered analysis/preparation modules are unchanged, including medication-adjusted T2D-internal species, full 25-model species LOCO, clinical Spearman, pooled fixed-study pathways, DESeq2+REML, interactions/additive, diversity, BMI23/CLR/no-medication sensitivity, exploratory abundant taxa, audits and display preparation. Frozen statistical definitions, BH families, references, cohort rules and random seeds are unchanged.

Outputs are 40 individual image panels (PDF/TIFF/PNG, plot-data/captions) and 16 workbooks: two main and S1-S14. S3 has eight sheets; S4 adds six full LOCO sheets. The retired medication sensitivity workbook remains internal and distinguishes adjusted main and covariate-omission sensitivity; S15 includes verified source notes AND database completeness from this same Run. No REVIEW worksheet/file or preview process is part of formal S15. Images retain latest integrated style and dimensions, including SFig2's saved BH brackets, full 64-contribution pathway display, two five-cohort forests, 600-dpi SFig3 and three SFig4 sensitivity scatters. No manual composite is made.

## Publication packaging

Public selection/copying share config/public_release_files.yml. `Rscript --vanilla build_public_package.R --manifest` regenerates actual payload hashes. Formal `--copy NEW_DESTINATION` still requires release admission and a nonexistent destination. Candidate extraction does not grant that admission. Inventory metadata accompanies the payload but cannot hash itself. Development preview/replacement scripts, source audit records and old outputs/caches are excluded. Extraction/static audit files remain outside this package.


Saved-result validation fixes and unresolved numerical-warning procedure: see [VALIDATION_REPAIR.md](VALIDATION_REPAIR.md). Post-run checks use tests/test_release_validation_repair.R and validate_release.R; do not remove results to satisfy extraction-time empty-directory tests. The independent REML diagnostic refits only when explicitly invoked by the author and does not change formal estimates.


## CLR alignment correction 2026-11-19
See CLR_ALIGNMENT_GUIDE.md before resuming/exporting. The author-approved CLR implementation supersedes its previous statistical hash; the other 30 statistical implementation files are unchanged. Old CLR checkpoints are not valid current publication inputs. Three changed result families await rerun and review; no new model was run during this edit.

Current submission S14 is clinical/confounder availability (formerly S15). See SUBMISSION_OUTPUT_PATCH.md.
