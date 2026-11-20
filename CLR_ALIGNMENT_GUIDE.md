# CLR alignment correction 2026-11-19

Author-approved changes replace the current CLR implementation: T2D-internal species models reuse the primary categorical medication mapping/reference; CLR BMI23 interaction reuses the harmonized cohort-adapted BMI23 grouping. Other CLR methods (CZM, limma, filtering, target-family BH) are unchanged. Pathway and interaction models remain without medication. This is an alternative statistical framework, not a claim that transformation is the only methodological difference.

There is no legacy uniform-23 or unadjusted-species branch in the current CLR module. Historical Run/Attempt/Export files remain immutable evidence and are never selected as the aligned result. New results include model contracts and sample-group/covariate provenance. Publication export rejects old CLR checkpoints. S2 gains metric labels from source-file identities without recalculating values.

## Run locally after reviewing this code change

```r
project_dir <- "Y:/研究生相关课题/分析类/胖瘦DM/T2D_BMI_Publication_RC4"
exe <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
rc4 <- function(script, args = character()) {
  code <- system2(exe, c("--vanilla", shQuote(file.path(project_dir, script)), shQuote(args)))
  if (code != 0L) stop("Step incomplete; retain the printed error and report path.")
}
rc4("tests/test_clr_alignment.R", c("--metadata-run", "Run_20261119_143714"))
rc4("run_all_development.R", c("--resume", "Run_20261119_143714", "--from", "clr"))
rc4("export_publication.R", c("Run_20261119_143714", "all"))
# Use the NEW PublicationExport identifier printed by the preceding command:
# rc4("validate_release.R", c("--run", "Run_20261119_143714", "--export", "NEW_EXPORT_ID"))
```

The existing checkpoint engine starts a new CLR Attempt and reruns the subsequently registered supporting modules. It does not refit earlier primary species/pathway/LOCO/additive/interaction models. Old attempts/exports are preserved; the current checkpoint selects the new results. This is a current-RC4 partial rerun, not an independent no-cache test. A new computer still starts with empty output/data_derived and executes the complete candidate pipeline.

## Required review before release

Three changed result families have no invented expected significant count: BMI25 internal species, CLR adapted-BMI23 taxonomy interaction, and CLR adapted-BMI23 pathway interaction. The validator still checks their feature counts, sample/model identities and target-family BH exactly. It deliberately reports pending author result/reference freeze until those new results are reviewed. The other 13 CLR result families retain their frozen expectations.

After an author reviews new coefficients, q values, comparison directions and figure changes, replace the three explicitly pending expectations in config/release_validation_plan.yml with the accepted observed counts and status ACCEPTED, archive that review, and rebuild the independently stored reference bundle/index and its hash. Do not merely remove the pending gate or copy unreviewed observed counts into the configuration. The CLR reference bundle must also accommodate new provenance columns/file names for otherwise unchanged results; compare the original scientific numeric fields in those other 13 families before replacing their reference snapshots. The old external reference is not the approved reference for the three altered models. No run/count was accepted during this code edit.

Re-export 16 workbooks and 40 panels into a new folder. Review S13, S2, SFig4B internal-comparison facet, their captions and numeric summaries; review the two revised CLR BMI23 interaction tables. Other frozen scientific results should match. Update manuscript/Response Letter and the handoff only from accepted new results.

Previously repaired validation defects remain repaired. REML numerical warnings, their diagnostic/author review, external reference comparison, fresh-computer no-cache testing, LICENSE/CITATION and author manuscript/ethics review remain separate release requirements. This correction does not resolve or waive them. No publication-admission flag is enabled here.
