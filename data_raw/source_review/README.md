# 已完成来源核查的结构化输入

这里接收已经完成的来源核查表，不重新检索、不自动生成来源证据、不改变样本。

把已核查内容转换为 UTF-8 CSV 后放入：

- `cohort_sources_verified.csv`：五队列各一行；字段见同目录 `.schema.csv`。
- `confounders_verified.csv`：五队列 × 变量概念；字段见同目录 `.schema.csv`。

`.schema.csv` 只有表头，是输入格式说明，不能改名冒充完成的核查数据。
两个真实输入尚未放入时，分析模块可以运行；S1、S15 导出会清楚报告缺少文件，整个投稿包返回非零状态。
这表示结构化输入未放置，不表示来源核查未完成。

队列代码固定为 KarlssonFH_2013、MetaCardis_2020_a、QinJ_2012、SankaranarayananK_2015、YuJ_2015。
保留 PMID/DOI、原研究诊断/纳排标准、证据等级和限制说明。QinJ_2012 使用1999 WHO标准；
MetaCardis的T2D映射A、non-T2D controls定义A-；YuJ_2015总体16例T2D和38例对照、个体映射无法独立重建的B级限定应保留。
边界注释异常仅记录，不触发样本删除。

混杂因素 `variable` 使用这些稳定概念ID：
`diabetes_duration`, `HbA1c`, `fasting_glucose`, `fasting_insulin_or_HOMA_IR`,
`disease_severity`, `medication`, `metformin`, `diet`, `smoking`, `alcohol`,
`physical_activity`, `education`, `income`, `occupation`, `environmental_factors`。
每个队列每个概念都需一行，可记录来源未采集、数据库未保留或不适合统一；不能把缺失当作没有暴露。
所有未知或限定状态照录已核查结论，代码不提升证据等级。
