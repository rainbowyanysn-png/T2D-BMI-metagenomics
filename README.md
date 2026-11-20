# T2D-BMI Publication RC4 — candidate

This is the complete code candidate extracted from the integrated analysis implementation. It awaits a full clean statistical run, all table/image checks and author review. No statistical or visual acceptance is claimed by code extraction.

Read [RC4_RUN_GUIDE.md](RC4_RUN_GUIDE.md) for R commands, environment restoration, first-use data downloads, resume, separate export and independent validation. Public execution requires only this package and obtainable frozen raw resources, not any old project/output.

R 4.5.1; curatedMetagenomicData 3.16.1; Maaslin2 1.22.0; other packages pinned by renv.lock. Main author-run command: `Rscript --vanilla run_all.R --candidate`. Static-only command: `Rscript --vanilla tests/static_rc4_candidate.R`.

22 registered modules plus inputs, 40 image panels, 16 workbooks. Statistical implementations and parameters match the source; Run identity and path protection are RC4-specific engineering. output/ and data_derived/ start empty. All release admission flags remain unapproved/pending. Independent reference comparison happens only after calculation, never as an analysis/plot input.

中文使用、打包、新电脑运行及分析/导出内容介绍：见 [RC4_中文使用说明.md](RC4_中文使用说明.md)。候选转交可运行 `Rscript --vanilla prepare_transfer.R --zip`；不会复制旧结果或改变正式发布准入。


Saved-result validation fixes and unresolved numerical-warning procedure: see [VALIDATION_REPAIR.md](VALIDATION_REPAIR.md). Post-run checks use tests/test_release_validation_repair.R and validate_release.R; do not remove results to satisfy extraction-time empty-directory tests. The independent REML diagnostic refits only when explicitly invoked by the author and does not change formal estimates.

修正CLR审阅/外部参考冻结、REML分类归档、正式入口依赖审查与新电脑洁净测试：见 [RELEASE_REVIEW_AND_CLEAN_TEST.md](RELEASE_REVIEW_AND_CLEAN_TEST.md)。所有实际审阅、冻结及完整运行均由作者执行；本轮仅编写代码。


Workbook export and independent streaming acceptance: see [WORKBOOK_PATCH_RUN_GUIDE.md](WORKBOOK_PATCH_RUN_GUIDE.md). This validation-only update preserves all statistical calculations and saved results; author review, CLR reference acceptance and REML release conditions remain separate.

Complete export preflight and explicit recovery of interrupted same-Run checkpoints: [PUBLICATION_PREFLIGHT_AND_RECOVERY.md](PUBLICATION_PREFLIGHT_AND_RECOVERY.md). Export no longer invokes the analysis runner; failed readiness checks do not create an Export directory.

Submission-output revision RC4_submission16_20261121: [SUBMISSION_OUTPUT_PATCH.md](SUBMISSION_OUTPUT_PATCH.md). Stable statistics/40 panels;16 formal workbooks; S7 current64-source display flags; former medication workbook internal only.
