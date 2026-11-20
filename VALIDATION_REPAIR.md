# RC4 release validation repair and unresolved REML warnings

Validation-only changes: typed finite/nonnegative scalar tolerance reading (absolute1e-8, pooled absolute1e-10, pooled relative1e-7 unchanged), named three-row pooled numeric-check contract, per-feature validation-only BH family audits, complete publication YAML loading before publication_data, separated saved-result/internal/external/warning/manual/no-cache scopes.

No saved result, statistical model/input/filter/q definition or admission flag is changed. No full clean run is needed to rerun the validator. Example in R after defining rc4 as in RC4_RUN_GUIDE.md:

```r
rc4("tests/test_release_validation_repair.R")
rc4("validate_release.R", c("--run", "Run_20261119_143714",
  "--export", "PublicationExport_20261119_152725.242_13804"))
```

Replace example IDs on another computer with that RC4's completed run/export. Each validator execution creates a fresh ReleaseValidation directory. Internal saved-result PASS does not certify independent reference reproduction, author visual review or a fresh-computer no-cache run. Unresolved numerical warnings remain a failing release gate and keep a nonzero validator exit status. Do not change production_ready to bypass this.

## REML warnings remain unresolved

The current aggregate warning logs record13 BMI25 and18 BMI23 Fisher-scoring warnings, without contrast/feature IDs. BMI25 pooled CSV omits tau2; neither module saved REML model objects. Zero tau2/I2 cannot identify a warning feature. The validation report includes an unlocalized31-warning inventory and the full saved-pathway universe with actual coefficients/SE/P/padj/k/I2 and tau2 only where available. Forest membership is listed for the universe, not attributed to the warnings. No manuscript was supplied to identify discussion features.

The following is an AUTHOR-INITIATED diagnostic refit, NOT a saved-result check:

```r
rc4("diagnose_reml_warnings.R", c("--run", "Run_20261119_143714",
  "--export", "PublicationExport_20261119_152725.242_13804"))
```

This diagnostic was NOT executed during the repair. It reads only that RC4 Run's saved per-cohort DESeq2 effects/SE, preserving the original REML formula/default settings and full contrast scope, and captures every warning with feature/contrast and k/effect/SE/P/tau2/I2. It records original formal padj/significance, representative-forest membership, finite/variance checks and differences from saved estimates, with a diagnostic BH column only when the successful family matches the formal family. It saves warned/differing fitted objects for later inspection. No optimizer substitution or profile likelihood is run automatically and no formal estimate is replaced. The diagnostic output is a fresh output/REML_Diagnostics_timestamp directory.

Saved CSV serialization can differ slightly from the original in-memory precision. Diagnostic warning counts/identities may therefore differ from the original aggregate31 warnings; do not automatically claim the historical warnings are resolved. Boundary/finiteness checks do not prove a global likelihood optimum. A separately justified author-run profile/optimizer comparison may still be needed, retained diagnostically and never substituted into formal results without approval.

The extraction-time static_rc4_candidate check requires empty output/data_derived. Do not delete completed Runs to satisfy it. Use test_release_validation_repair plus validate_release for this post-run task. Packaging manifests are updated for changed code only. Prior failed reports and submitted workbooks/images remain unchanged.


## Diagnostic version string check

The diagnostic compares the loaded package raw DESCRIPTION Version against config/package_versions.csv, exactly as recorded in the original Run. R packageVersion() normalizes metafor 4.8-0 to 4.8.0 for display, so its formatted string is not compared directly to the raw frozen version. Both strings and the loaded namespace path are recorded in diagnostic_package_version_check.csv. Actual DESCRIPTION version mismatches still stop; no package is installed or changed. This preflight correction does not run or resolve the REML diagnostic itself.


## Shared complete input preflight

Diagnostic module directories are named character vectors, accessed through a checked [[module]] accessor (lists also supported). Full diagnosis and --preflight-only share complete preparation of all nine contrast families: sample/export identity, cohort file availability and schemas, study/feature identities, original validity rules and formal-feature coverage. The original REML calls/filters/settings and result-record logic are unchanged. Preflight performs no model fit and writes its report outside the public package. Run the same diagnostic command with --preflight-only to check inputs first; without that option it performs author-requested REML refits. Existing failed diagnostic folders are retained.

Diagnostic warned/differing fitted objects now use short fits/family-hash/feature-hash.rds paths to limit Windows path length. The full pathway/comparison identity and diagnostic_fit_file mapping remain in the diagnostic CSV. Hash-path collisions or existing files stop instead of overwriting.
