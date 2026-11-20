# Methods and outputs

All baseline statistical implementation files retain their accepted-source hashes. Central analysis definitions and sample/feature assertions remain in config/analysis_config.yml, expected_results.yml and release_validation_plan.yml.

New species LOCO: five original BMI25 contrasts by five omitted cohorts; all 200 species as input, original strict prevalence/variance filtering, TSS/LOG, standardize TRUE, full-model MaAsLin2 qval, medication-adjusted T2D-internal comparison only. Complete 5000-row outputs, primary 120-record/79-feature summaries, design/sample/medication/exclusion audits and package/session/code hashes are saved in the current Run.

Tables: 16 workbooks, S3 beta-diversity eight-sheet normalized format, S4 with original main species results and six LOCO sheets, S5-S13 numbering unchanged; former medication S14 archived; confounder availability renumbered from S15 to S14.

Figures: inventory is defined only by config/publication_panels_manifest.csv. Current build has 40 panels. SFig1 cohort metrics are combined in an unnumbered shared-legend panel. No Figure2F, SFig6 or supplementary beta diagnostic figures are exported. Full beta results stay in S3 and full REML results in S8. SFig4 contains only three sensitivity coefficient scatters; old alpha-interaction images and all SFig5 images are omitted. Statistical tables remain unchanged.

Formal uncorrected Aitchison PERMDISP and new forest selection are validated directly from current-run saved estimates. Body/Response Letter update paragraphs are generated under each export/manuscript_updates for author insertion. See RC4_RUN_GUIDE.md and RC4_中文使用说明.md for current operation and output descriptions. Extraction/engineering audit records are intentionally outside this public package.
