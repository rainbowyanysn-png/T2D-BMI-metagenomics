# RC4新电脑使用说明、分析与导出目录

适用对象：首次接收本项目、没有旧项目或原始数据缓存的使用者。版本：RC4发表候选代码。当前电脑已有运行与导出文件，不等于完成了新电脑无缓存复现、独立核验和作者最终审图；发布准入状态以config/module_admission.yml为准，不自动改为通过。

## 1. 给对方什么文件

交付一份按公开规则生成的代码候选包即可。不要直接压缩当前工作文件夹，因为它已经包含本机运行结果。

| 必须包含 | 用途 |
|---|---|
| R/、config/、tests/ | 全部分析/绘图/导表函数、配置、分类与注释字典和验证程序 |
| data_raw/ | 已核实的来源说明、必要映射和schema；不是预先计算的分析矩阵 |
| renv.lock、data_manifest.csv | 包版本和原始资源ID、版本、SHA-256 |
| 根目录正式R入口 | 新运行、恢复、单独导出与独立验证 |
| README、运行指南、版本/许可/引用文件 | 操作、状态和来源说明 |
| public_file_manifest.csv、public_file_sha256.csv | 实际公开载荷文件及校验值 |
| 空output/、空data_derived/ | 新电脑运行后再生成内容；仅保留.gitkeep |

不交付：本机Run、PublicationExport、ReleaseValidation、派生RDS/矩阵、ExperimentHub/renv缓存、审计目录、旧V3/CleanTest/FigureRefinement/IntegrationPreview、预览/替换/缩放脚本、.Rhistory、.RData和历史备份。分析所需原始公共资源在首次运行时按冻结ID取得；来源表和字典已随包包含。许可及数据引用仍按原始资源要求执行。

### 在当前电脑生成干净转交包

在R 4.5.1中设定当前RC4路径，定义下方命令函数；执行：

```r
rc4("prepare_transfer.R", "--zip")
```

这一步只复制公开代码和说明、验证文件哈希，不运行统计或导图。默认在RC4工作目录旁创建一个T2D_BMI_RC4_CodeCandidate_时间戳文件夹，以及同名ZIP；只需把ZIP发给对方，解压后即有完整代码。复制记录放在包外的RC4_Audit中。源工作目录里的统计结果保留原样。

如只需要文件夹：`rc4("prepare_transfer.R")`。指定新的目的地：`rc4("prepare_transfer.R", c("--zip","--dest","D:/Transfer/RC4_Candidate"))`；目的地不得已存在，也不能位于当前RC4项目里面。

这是候选代码转交，不代表正式发布通过。build_public_package.R --copy的正式准入门禁保留。不要手工改production_ready来打包。

## 2. 空电脑先准备什么

1. 安装R **4.5.1**；RStudio可选，代码不依赖RStudio。
2. 将ZIP完整解压到一个较短路径，例如D:/T2D_BMI_RC4；可以换盘符、改文件夹名称，但不改变内部结构。避免Windows路径过长。
3. 提供Times New Roman字体。图形导出会检查字体，不静默替换；按所用系统安装所需图形/编译依赖。如果恢复提示需要源码编译，按R安装提示补齐匹配的编译工具。
4. 首次恢复包和取得原始资源需要联网。其他电脑不需要本人的Y/W盘，也不需要任何旧Run。互联网不通或指定资源不可用时应停止并保留错误，不改用新版本资源。

### 首次环境恢复：只做一次，独立于分析

在R中执行，第一行改成解压后的实际路径。此处允许主动安装环境管理工具renv并恢复锁定环境；模型运行过程不会安装或更新包。

```r
project_dir <- normalizePath("D:/T2D_BMI_RC4", winslash="/", mustWork=TRUE)
stopifnot(as.character(getRversion()) == "4.5.1")
if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv")
renv::restore(project=project_dir, prompt=FALSE)
```

冻结核心版本：curatedMetagenomicData 3.16.1、Maaslin2 1.22.0；其余以renv.lock及版本清单为准。不要用更新包版本解决报错。如果设置过T2DBMI_R_LIBRARY，应确认它没有指向旧项目或另一套版本库；没有特殊需求时不设置它。

### 定义运行命令函数：每次重启R后重新执行

```r
# 改为这台电脑上解压后的路径；不需要setwd或文件选择弹窗。
project_dir <- normalizePath("D:/T2D_BMI_RC4", winslash="/", mustWork=TRUE)
rscript <- file.path(R.home("bin"),
  if (.Platform$OS.type=="windows") "Rscript.exe" else "Rscript")
rc4 <- function(script, args=character()) {
  code <- system2(rscript, c("--vanilla",
    shQuote(file.path(project_dir, script)), shQuote(args)))
  if (code != 0L) stop("本步未完成，请保留控制台报错和打印的日志路径。", call.=FALSE)
  invisible(code)
}
```

入口由脚本自身位置确定根目录。R环境可以为空，不需要预先加载任何对象；每次子进程使用--vanilla。可从其他工作目录运行。

## 3. 新电脑执行顺序

### 3.1 空输出目录的静态检查

```r
rc4("tests/static_rc4_candidate.R")
```

只核查代码、配置、依赖和空目录，不下载资源、不拟合模型或绘图。该检查含“候选包output/data_derived为空”的准入要求，适用于刚解压的代码包；已有正式结果后再次运行会在空目录项失败，不要为了让它通过删除已有结果。静态报告默认写到项目外的同级审计目录。

### 3.2 完整新运行并统一导出

```r
rc4("run_all.R", "--candidate")
```

此命令会从冻结原始资源重建数据，执行全部分析，之后导出16个工作簿和40个图片面板。正式运行可能耗时，不预估未经实测的耗时。首次缺缓存时按固定ExperimentHub ID/数据版本下载，逐一核对SHA-256，记录获取状态；失败不替换数据。

不要在运行时清理目录、修改参数或开启针对同一Run的第二个进程。失败时后续步骤会停止；保留Run编号和日志路径。

### 3.3 中断恢复与单独重导出

以下Run编号必须换成**本电脑、本RC4项目**控制台打印的编号。不要填写旧项目编号，也不要输入PublicationExport编号。

```r
run_id <- "Run_YYYYMMDD_HHMMSS"  # 替换后才运行
rc4("run_all_development.R", c("--resume", run_id))

# 恢复分析结束后，或只有导出报错时：只导出，不重拟合模型
rc4("export_publication.R", c(run_id, "all"))
# 二选一的局部重导出命令：
# rc4("export_publication.R", c(run_id, "tables"))
# rc4("export_publication.R", c(run_id, "figures"))
```

已验证且签名一致的模块跳过；未完成模块继续。检查点内容、代码/参数不符时停止，不能自动复用别的Run。RC4创建/恢复/导出均核查自身Run身份。每次重导出新建PublicationExport，不覆盖旧图表。

### 3.4 完成后独立验证

```r
run_id <- "Run_YYYYMMDD_HHMMSS"                   # 换为实际RC4 Run
export_id <- "PublicationExport_YYYYMMDD_HHMMSS.mmm_PID" # 换为实际导出
rc4("validate_release.R", c("--run",run_id,"--export",export_id))
```

此阶段核对保存结果、工作簿和图片，不重新拟合模型。可选参考目录仅提供给独立验证进程：在参数后添加c("--reference-root","实际参考目录")。缺少外部参考不能表述为逐特征已复现；必须如实保留未核验状态。分析和绘图进程不读取外部旧结果。

## 4. 保留的分析内容

输入模块＋22个分析/准备模块（确切顺序见config/pipeline_modules.csv）。

| 内容 | 模型/用途与正式校正口径 |
|---|---|
| 输入与QC | 原始资源、五队列、1186人；物种200、基础通路5945；原队列纳入/参考组/ID断言 |
| BMI25五个物种比较 | MaAsLin2，standardize TRUE，全模型qval；T2D内部主模型含分类降糖治疗协变量 |
| 跨队列LOCO | 五比较×删除五个队列，共25模型；完整200物种输入，保留排除/缺失标记和稳定性审计 |
| 物种—临床探索性相关 | 20物种×BMI/年龄/性别=60项Spearman，统一BH；不是因果或独立验证 |
| 总体个体通路模型 | disease+BMI+age+sex+study_name均固定效应；TSS/LOG、standardize TRUE；全模型qval |
| 队列内通路＋荟萃 | DESeq2＋REML；每对比所有成功合并特征内BH padj；保留队列效应和异质性 |
| 正式交互 | BMI25、总体中心化连续BMI、队列适配BMI23；standardize FALSE；目标项按定义/物种或通路分家族BH |
| 加性模型 | disease+BMI_per5+age+sex+study_name；standardize FALSE；T2D/BMI分别在200物种中BH |
| α多样性 | 主模型、交互、ICC和两两比较明细，沿用原模型/比较规则 |
| β多样性 | 原始CZM/CLR Aitchison、边际PERMANOVA、PERMDISP和逐队列结果；ComBat排序单列为描述性 |
| BMI23敏感性 | QinJ/YuJ用23，其他队列保留25；保留已批准的四组levels例外 |
| CLR敏感性 | 保留CZM零处理、CLR、原模型和目标项检验家族 |
| 删除药物协变量 | BMI25/BMI23 T2D内部敏感性；全模型qval，不是排除用药者或metformin-naive分析 |
| 高丰度菌门/物种 | Top4/Top8探索性LMM；未调整emmeans P按层级24/48项BH；作图读取保存BH，不重算 |
| 其余准备/审计 | 样本流程、元数据和药物可得性、来源核查、当前Run作图数据准备 |

不包含已归档的F/B、分层通路Pilot、极端BMI、功效/MDE、核心得分或旧预览拼接分析。完整统计参数见config/analysis_config.yml和model_specifications.csv；说明文字不能替代配置。

## 5. 16个工作簿：当前实际编号

注意：以本RC4代码为准，不使用早期设计草案的S编号。

| 文件 | 主要内容 |
|---|---|
| Main_Table_1.xlsx | 总体及BMI25四组人口学/队列构成描述 |
| Analysis_Summary_Table.xlsx | 分析模型、检验家族、测试/显著数量汇总 |
| Supplementary_Table_S1.xlsx | 五队列来源、筛选流程、四组人数、最终样本ID |
| S2 | α多样性模型、ICC、Tukey等比较与明细 |
| S3 | β多样性8个规范sheet；正式raw Aitchison、PERMDISP、交互及描述性ComBat分开 |
| S4 | 五个BMI25物种比较完整结果，加6个完整LOCO sheet |
| S5 | 60项探索性临床Spearman及BH结果 |
| S6 | 总体固定队列效应MaAsLin2全部5945通路特征 |
| S7 | taxon-stratified通路贡献候选、注释及展示选择 |
| S8 | 队列DESeq2、REML合并效应、置信区间及异质性 |
| S9 | 通路注释、功能分类及固定字典 |
| S10 | BMI25/连续BMI/BMI23正式物种与通路交互 |
| S11 | T2D＋BMI加性模型完整200物种及目标系数 |
| S12 | 队列适配BMI23敏感性、样本构成、交互和一致性 |
| S13 | CZM＋CLR完整敏感性及一致性结果 |
| S14 | 来源核查表＋本Run的Database_completeness；包含修正药物说明 |

S2–S14完整文件名均为Supplementary_Table_S编号.xlsx。各sheet另存CSV、数据字典、标题/脚注，并执行工作簿完整性检查。正式S14不含REVIEW文件名或审阅过程说明。

## 6. 40个单独图片面板

| 图组 | 面板数量 | 内容 |
|---|---:|---|
| Figure 1 | 1 | 研究设计和分析流程 |
| Figure 2 | 5 | Shannon、正式原始Aitchison总体/四组排序、逐队列PERMANOVA |
| Figure 3 | 5 | 门/物种组成、UpSet、lollipop、五比较物种效应热图 |
| Figure 4 | 6 | 五个物种火山图＋探索性临床Spearman热图 |
| Figure 5 | 3 | 总体通路火山图、菌属条形图、物种—通路气泡图 |
| Figure 6 | 5 | REML通路重叠、效应热图、两张五队列森林图、功能类别 |
| Figure 7 | 4 | 正式交互家族、BMI加性森林、双效应散点、结论汇总 |
| Supplementary Figure 1 | 5 | Richness、上下组合逐队列α多样性、描述性ComBat |
| Supplementary Figure 2 | 2 | Top4门/Top8物种紧凑箱线组合及保存的BH括号 |
| Supplementary Figure 3 | 1 | 全64条贡献的物种→菌属→通路桑基图 |
| Supplementary Figure 4 | 3 | BMI23物种、CLR物种、CLR通路一致性散点 |

合计29正文＋11补充=40面板。保留SFig1C上下组合和Top4/Top8分面组合，其余按当前manifest逐面板独立导出；不自动拼正文整版图，由作者组图。

每面板提供矢量PDF、LZW TIFF、300 dpi PNG、plot-data CSV及图注。普通TIFF为300 dpi，SFig3为600 dpi。SFig2尺寸180×47.5/180×90 mm，Fig6B为140×180，SFig3为240×180，SFig4为234×115。准确清单与尺寸以配置和figure_index.csv为准。正式阴性交互不能解释为确证BMI效应修饰；通路贡献不是实验验证机制。

## 7. 运行后到哪里找

```text
output/
  Run_日期_时间/
    RC4_RUN_IDENTITY.txt
    checkpoints/             # 模块完成状态与签名
    modules/                 # 本次完整统计结果及模块审计
    logs/                    # 每次调用、错误、警告、版本和耗时
  PublicationExport_新时间戳/
    tables/main/             # 两个主表工作簿
    tables/supplementary/    # S1–S14
    tables/csv/              # sheet CSV、字典、脚注
    figures/main_panels/     # 正文PDF/TIFF
    figures/supplementary_panels/ # 补充PDF/TIFF
    figures/previews/        # PNG
    figures/captions/        # 单面板图注
    figures/combined_legends/ # 图组图注
    plot_data/               # 实际作图数据、选择和括号审计
    validation/              # 导出文件检查
    manifests/               # 输入、代码、资源/版本哈希
    manuscript_updates/      # 正文/回复信β多样性段落草稿
    export_status.csv
    export_gaps.csv
  ReleaseValidation_新时间戳/ # 独立验收报告

data_derived/Run_日期_时间/   # 本次从原始数据重建的输入检查点
```

Run和PublicationExport是不同职责：前者保存统计计算，后者保存某一次导出；一次Run可以重导出多次。这不是重复拟合。看到export_gaps为空也不等于作者已完成科学内容和视觉审阅。

## 8. 报错后发什么

发送完整错误、Run/Export编号、控制台打印的日志目录、export_status/export_gaps（若已生成）。下载失败发送resource_acquisition_failure.txt或download_attempts日志；不要删除旧结果后盲目全跑。R包版本错误优先核对环境恢复，路径错误优先核对project_dir和文件完整解压。

不要把README中的示例Run_YYYYMMDD_HHMMSS原样运行；必须改成实际编号。本说明不使用readline、file.choose或choose.dir，避免文件选择窗口文字异常。

## 9. 新电脑验收及正式发布

新电脑验收须实际完成：环境恢复、空缓存下载记录、输入ID/样本/特征断言、所有模型、图表导出、独立验收和作者审图。没有运行过就不能标为PASS。不同操作系统允许既定容差内浮点差异，但ID/方向/显著状态等必须满足原有规则，不放宽阈值。

目前的候选转交包不包含本机旧结果；给对方一个ZIP和本说明即可。若另外给审稿人本机图表，应放在独立“结果材料”包，不混入可独立运行代码包。最终正式公开仍需确认来源许可、机构/期刊伦理说明和仓库URL/DOI，并按既有release admission门禁办理。

投稿输出修订 RC4_submission16_20261121：原药物敏感性S14仅内部归档；临床与混杂因素表由原S15改为S14，正式合计16个工作簿。详见SUBMISSION_OUTPUT_PATCH.md。
