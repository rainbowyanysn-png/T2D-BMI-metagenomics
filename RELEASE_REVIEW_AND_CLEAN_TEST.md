# RC4：修正CLR审阅、REML分类与新电脑洁净测试

本阶段只提供代码。程序尚未替作者审阅、接受或冻结任何结果；统计实现、模型、过滤、q值家族、历史输出及production_ready保持不变。正式统计和导出进程不会读取这些审阅包/外部参考。

## 1. 在R中定义命令函数

```r
project_dir <- "Y:/研究生相关课题/分析类/胖瘦DM/T2D_BMI_Publication_RC4"
exe <- file.path(R.home("bin"), if (.Platform$OS.type=="windows") "Rscript.exe" else "Rscript")
rc4 <- function(script, args=character()) {
  status <- system2(exe, c("--vanilla", shQuote(file.path(project_dir,script)), shQuote(args)))
  if (status!=0L) stop("本步未完成，请保留控制台错误和打印的报告目录。",call.=FALSE)
  invisible(status)
}
```

以下路径均在RC4之外；目标目录必须尚未存在。如已有同名目录，请换一个新名称，不删除旧目录。全部Run/Export/诊断ID显式指定，不自动选最新结果。

## 2. 生成审阅材料（只读结果，不拟合）

```r
review_dir <- "Y:/研究生相关课题/分析类/胖瘦DM/RC4_ReleaseReview_20261119"
rc4("tests/test_rc4_release_review.R") # 使用虚构小数据测试审阅拒绝条件，不接受真实结果
rc4("review_rc4_release.R", c("prepare",
  "--run", "Run_20261119_143714",
  "--export", "PublicationExport_20261119_202555.628_3588",
  "--diagnostics", "REML_Diagnostics_20261119_171725.766_32664",
  "--baseline-attempt", "a001",
  "--out", review_dir))
```

程序根据当前checkpoint选修正后的CLR（应为a002），核对模型、样本、药物分类和队列适配BMI23，再与显式指定的旧a001进行**审阅验证**。旧结果只作为差异证据，不作为拟合、绘图或新统计的输入。

产生的主要文件：

| 文件/目录 | 用途 |
|---|---|
| CLR_result_summary.csv | 16组当前结果的样本/特征/正负及显著数、公式、目标项、SHA-256 |
| CLR_old_new_scientific_comparison.csv | 三组改模型结果与其余13组的区分；科学字段/ID/顺序/显著状态对照 |
| clr_snapshot/ | 当前16组结果、模型契约和样本元数据副本 |
| baseline_clr_evidence/ | 历史CLR证据，不能成为修正结果的正式参考 |
| publication_review_evidence/ | 新S2、S13、SFig4B/C的PDF/PNG/plot-data/图注，便于核对图表变化 |
| CLR_AUTHOR_DECISIONS.csv | 所有行初始PENDING；作者填写审阅记录 |
| REML_warning_classification.csv | 按数值有效性/差异、正式显著与森林图涉及情况分层，不等于无害/已解决 |
| REML_warning_count_comparison.csv | 原聚合警告与诊断重放数量对照；数量相同也不证明历史逐特征一一对应 |
| REML_full_replay_stability.csv | 全范围诊断的错误、差异和检验家族/显著状态证据 |
| REML_AUTHOR_DECISIONS.csv | 逐警告特征的作者意见/理由/补充证据路径 |
| review_input_hashes.csv、review_context.yml | 来源绑定、代码/结果哈希、Run/Export/Attempt/诊断身份 |

需要重点审阅：BMI25 T2D内部物种增加分类药物协变量；两组CLR BMI23交互恢复队列适配规则。其余13组必须保持既有科学字段和ID顺序。药物删除敏感性不是排除用药者；CLR结论不能替代主模型或确证交互。

## 3. 作者确认后才冻结CLR参考

先打开审阅目录中的CLR_AUTHOR_DECISIONS.csv，逐行填写：

- decision：审阅接受后写ACCEPT，否则保留PENDING；不能批量预设接受。
- accepted_significant_count：明确填写作者核对的显著数，与结果汇总一致。不是为了显著数调整模型。
- reviewer、review_date（YYYY-MM-DD）、rationale：姓名/日期/理由。
- visual_review：核对新表和相关图片后写CONFIRMED。

保留CSV字段名称，不修改source_sha256或准备阶段的汇总、科学对照、快照文件。13组未改模型结果若存在实质差异，冻结必须停止，不改容差。

```r
clr_reference <- "Y:/研究生相关课题/分析类/胖瘦DM/RC4_CLR_Reference_20261119"
rc4("review_rc4_release.R", c("freeze-clr", "--review", review_dir, "--out", clr_reference))
```

此命令在新的外部目录复制**当前已审阅**CLR结果、元数据、作者记录，生成精确reference_index.csv和SHA-256。然后仅更新验证计划中三组结果的验收显著数及ACCEPTED_CLR_AUTHOR_REVIEW状态，绑定专用config/clr_review_reference.yml，并重建公开清单。激活前的配置备份保存在外部参考目录activation_before/，失败不静默激活。统计实现及其哈希不改，其他验收数字不改，发布准入不改。

这个包是今后独立新Run比较的基线，**冻结本身不是独立复现**。来源Run不能拿它自比并报告独立CLR复现通过。旧的全模块参考索引保持原样，不能将它当作三组修正CLR的权威结果；新的CLR参考只覆盖CLR，不能宣称所有分析都完成外部逐特征比较。

## 4. REML分类审阅归档（不自动放行）

检查REML_warning_classification.csv及完整诊断结果，必要时回查诊断fit对象和另行授权的profile/优化器比较。不能因为tau2=0、结果不显著或系数近似一致就认定警告无害。原警告日志没有逐特征ID，诊断重放使用序列化输入，历史精确归因仍有局限。作者需核对稿件讨论，程序不推断未提供的正文内容。

填写REML_AUTHOR_DECISIONS.csv：decision从REVIEWED_WITH_LIMITATIONS、REQUIRES_FURTHER_DIAGNOSIS、REJECT_RELEASE中选择；同时填写reviewer/date/rationale，可提供evidence_file实际路径。不要写“全部无害”。

```r
reml_archive <- "Y:/研究生相关课题/分析类/胖瘦DM/RC4_REML_Review_20261119"
rc4("review_rc4_release.R", c("archive-reml", "--review", review_dir, "--out", reml_archive))
```

归档复制分类、意见、哈希及作者提供的证据。**现有REML发布警告门禁不会因此清除。** 审阅归档、统计警告解决、作者正式发布批准是不同事项；仍需数值研究时保留待完成状态。

## 5. 正式入口依赖审查与洁净转交文件清单

```r
rc4("audit_rc4_dependencies.R")
```

结果写入包外同级RC4_Audit/DependencyAudit_时间戳/。包括：

- clean_test_payload_files.csv：实际应转交的代码/配置/说明/字典及大小/SHA-256。
- clean_test_additional_metadata.csv：清单元数据及两个空目录标记。
- registered_module_dependencies.csv：输入+22模块文件完整性。
- bootstrap_dependency_inventory.csv：分析/导出固定加载器的完整依赖。
- source_io_call_inventory.csv、literal_path_findings.csv：source及I/O调用、个人绝对路径/旧项目名字检查。
- dependency_checks.csv：静态检查结果；不会运行模型、下载资源或导出图表。
- CLEAN_TEST_COMMANDS.R：后续命令提示；需先定义rc4命令函数及本机实际路径。

依赖审查覆盖当前代码及正式入口，并保留验证专用的参考读取和路径防护测试夹具。动态source的固定bootstrap/模块注册表受检查；静态扫描不能证明任意运行路径，完整无缓存测试仍由新电脑执行。不会为了满足空目录检查而删除当前结果。

## 6. 打包并交给新电脑

```r
rc4("prepare_transfer.R", "--zip")
```

按同一公开选择规则创建新的候选ZIP，只含代码/配置/说明/原始资源与来源字典清单；output和data_derived为空，不包含审阅包、外部参考、诊断结果、旧Run或任何缓存。审阅资料/CLR参考如用于独立验证，另行交付至代码包之外，不作为分析输入。production_ready等保持真实状态。

## 7. 新电脑无缓存完整测试（由作者执行）

使用新OS用户/虚拟机或确认无数据缓存的环境；安装R4.5.1、Times New Roman及必要系统依赖，将ZIP解压到短路径，例如D:/T2D_BMI_RC4。不得复制旧Run、派生RDS或ExperimentHub缓存。

```r
project_dir <- normalizePath("D:/T2D_BMI_RC4", winslash="/", mustWork=TRUE)
stopifnot(as.character(getRversion())=="4.5.1")
if (!requireNamespace("renv",quietly=TRUE)) install.packages("renv")
renv::restore(project=project_dir,prompt=FALSE)  # 明确的环境准备，不是分析时安装
exe <- file.path(R.home("bin"),if (.Platform$OS.type=="windows") "Rscript.exe" else "Rscript")
rc4 <- function(script,args=character()) {
  status <- system2(exe,c("--vanilla",shQuote(file.path(project_dir,script)),shQuote(args)))
  if (status!=0L) stop("本步未完成，请保留错误和日志目录。",call.=FALSE)
  invisible(status)
}
rc4("tests/static_rc4_candidate.R") # 此时空目录才适用
rc4("audit_rc4_dependencies.R")
rc4("run_all.R","--candidate")    # 这里才完整重建统计与图表
```

首次环境和原始资源获取需要联网，严格按冻结ID/版本/SHA-256；不替换资源或升级模型包。来源下载不可用时停止并记录。

```r
# 以下编号换成本机RC4实际新Run/Export：
run_id <- "Run_YYYYMMDD_HHMMSS"
export_id <- "PublicationExport_YYYYMMDD_HHMMSS.mmm_PID"
# 仅当分析中断时：
# rc4("run_all_development.R",c("--resume",run_id))
# 完成后需要单独重导出时，不重拟合：
# rc4("export_publication.R",c(run_id,"all"))
rc4("validate_release.R",c("--run",run_id,"--export",export_id))
# 外部已接受的修正CLR参考另行放置，只有验证进程读取：
rc4("validate_release.R",c("--run",run_id,"--export",export_id,
  "--clr-reference-root","D:/RC4_CLR_Reference_20261119"))
```

--clr-reference-root只验证CLR，新Run不能与参考来源Run相同。旧全模块参考如另有已批准版本，仍通过--reference-root指定；其未更新的CLR条目不会自动获得豁免或被替换，缺少/不一致要如实报告。全模块参考比较与专用修正CLR比较分开记录。

17工作簿、40面板、当前Run数据来源、冻结值及每种q家族继续验收。REML数值门禁、人工图表审阅、来源许可/伦理说明和仓库URL/DOI保持各自状态；不因审阅归档、冻结或静态PASS就将全部发布条件改成通过。


## 等价CLR重跑Attempt与导出来源绑定

如果Run因重跑指向新的CLR Attempt，而既有导出来自上一Attempt，不按路径差异直接判统计不一致。审阅程序会核验两次COMPLETED检查点、相同模型/输入签名、完整CLR文件集合，以及每个结果和元数据在当前Attempt、导出Attempt、导出input_sha256清单中的三方SHA-256。全部一致才允许继续，并生成CLR_export_checkpoint_binding.csv记录两个Attempt和等价关系；否则仍要求从当前Run重新导出。不重新计算/替换统计值，不跳过来源验证。
