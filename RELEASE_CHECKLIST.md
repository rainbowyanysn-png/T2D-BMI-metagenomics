# Release checklist: not yet admitted

- Run tests/static_rc4_candidate.R; inspect exact static scope.
- Start the candidate with empty output/data_derived (the delivered folders contain only .gitkeep). Do not delete results from any other project.
- Rebuild all original inputs; verify sample/feature/cohort identities and frozen resource hashes.
- Complete all baseline modules plus 25 species LOCO fits.
- Reproduce all formal q-family/model contracts and LOCO acceptance values.
- Export 16 workbooks; S3 has eight sheets and S4 retains original sheets plus six LOCO sheets.
- Verify Figure6C automatic selection and the manifest-defined panel inventory; no SFig6 or Figure2F.
- Run validate_release.R against the generated Run/Export; external reference comparison is separate and must be truthfully marked when unavailable.
- Review all images, S3/S4 Excel values, source-review evidence, and generated manuscript/Response Letter text.
- Confirm repository URL/DOI, data permissions and institutional/journal ethics wording.
- Only after these checks may an author explicitly update admission flags and extract an upload package. Static PASS is not clean-run or publication PASS.


## Public file selection

The shared config/public_release_files.yml policy is used by both inventory generation and publication copying. Regenerate the inventory with `Rscript --vanilla build_public_package.R --manifest`. After clean-run validation and explicit author admission only, `Rscript --vanilla build_public_package.R --copy NEW_DESTINATION` extracts a new package; it refuses existing destinations and the closed release gate. Inventory files describe the public payload and are copied alongside it (self hashes are intentionally not recursive).

`preview_publication.R` is a local development-only entry that reads accepted saved checkpoints. It and preview replacements/resizing scripts, review backups, local outputs and audit directories are excluded from public payloads. The public independent workflow uses run_all.R, run_publication_candidate.R, run_all_development.R, export_publication.R and validate_release.R, with their included R/config/tests/data-source dependencies. It does not require CleanTest, FigureRefinement or IntegrationPreview. Do not copy the whole development folder as a release.

- RC4 remains a candidate. Require RC4 identity for all resume/export inputs. Keep extraction/static audits outside the code package. See RC4_RUN_GUIDE.md.
