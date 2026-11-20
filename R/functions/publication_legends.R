pub_write_combined_legends <- function(index,export) {
  summaries<-c(
    `Figure 1`="Study design. Analyses branch in parallel from the prepared study data. Species and pooled MaAsLin2 pathway features use full-model q values; cohort DESeq2 + REML uses BH within each contrast; interactions use target-term BH and the additive model uses separate exposure-specific BH families.",
    `Figure 2`="Community diversity at BMI25. Overall alpha models include disease, BMI, age, sex and a cohort random intercept; four-group models include group, age, sex and a cohort random intercept. Saved pairwise comparisons use Tukey-adjusted P values followed by BH across six comparisons per metric. Formal beta inference uses uncorrected CZM/CLR Aitchison distances, marginal PERMANOVA and PERMDISP; dispersion affects interpretation of ordination.",
    `Figure 3`="Descriptive composition and species associations. Formal species significance uses MaAsLin2 full-model q values across feature-by-covariate coefficients within each contrast (69/32/10/4/5 discoveries). The T2D-internal primary comparison adjusts for medication; its 190 tested species include five q<0.05 associations. Warm colours indicate positive effects, blue negative effects, and grey untested heatmap cells. Stars denote q<0.05, q<0.01 and q<0.001.",
    `Figure 4`="Species associations and exploratory clinical correlations. Volcanoes use MaAsLin2 full-model q values; warm/blue points denote positive/negative coefficients and pale grey nonsignificant features. Correlation colours encode rho; stars use one BH family containing all 60 tests. Correlations are pairwise-complete and unadjusted for covariates.",
    `Figure 5`="Pooled pathway-feature associations from MaAsLin2 with disease, BMI, age, sex and study as fixed effects, TSS normalization, LOG transformation and standardize TRUE. Significance is based on the MaAsLin2 full-model q value. Taxonomic assignments represent inferred functional potential.",
    `Figure 6`="Cohort DESeq2 and REML pathway-feature associations. BH correction includes all successfully pooled features separately within each comparison. Heatmap stars denote padj<0.05, <0.01 and <0.001; grey cells are untested. Forest intervals are 95% confidence intervals. Stratum-specific differences in significance do not establish effect modification.",
    `Figure 7`="BMI-stratified associations, formal interactions and additive models. Interaction BH families are separate for each BMI definition and feature domain. The additive model is disease + BMI_per5 + age + sex + study_name with standardize FALSE, with separate 200-species BH families for T2D and BMI. Warm/blue denote positive/negative effects; association classes distinguish T2D-only, BMI-only, Both and Neither. Statistical effect modification was not established.",
    `Supplementary Figure 1`="Additional diversity summaries. Alpha comparisons use their saved model-specific adjustments. Uncorrected Aitchison distances provide formal beta inference. ComBat ordination is descriptive only. BMI25 groups use the 25 kg/m² threshold.",
    `Supplementary Figure 2`="Secondary/exploratory abundant-taxa distributions and LMM contrasts, adjusted for age, sex and random cohort intercepts. Unadjusted emmeans pairwise P values form one 24-test phylum BH family and one 48-test species BH family. Brackets use p_adj_BH<0.05. Formal species inference remains MaAsLin2.",
    `Supplementary Figure 3`="HUMAnN taxonomic assignments of selected pooled MaAsLin2 pathway features. Assignment counts do not measure metabolic flux or establish an experimentally verified mechanism. Selection uses full-model q values.",
    `Supplementary Figure 4`="Sensitivity analyses using alternative BMI thresholds and CLR transformations. Each point represents a feature matched between the BMI25 main analysis and the corresponding sensitivity analysis. Blue and orange indicate concordant and discordant coefficient directions. Dark points were FDR significant in both analyses, whereas pale points were not jointly FDR significant. Facet annotations report matched feature count, direction concordance, descriptive coefficient correlation, and shared FDR-significant feature count. No correlation P value is calculated. Coefficients from CLR and main-model analyses may be on different scales.")
  folder<-file.path(export,"figures/combined_legends");dir.create(folder,recursive=TRUE,showWarnings=FALSE)
  for(figure in intersect(names(summaries),unique(index$manuscript_figure))){
    z<-index[index$manuscript_figure==figure,,drop=FALSE]
    assert_true(nrow(z)>0,paste("No panels for combined legend:",figure))
    z<-z[order(z$panel),,drop=FALSE]
    paragraphs<-vapply(seq_len(nrow(z)),function(i){
      caption<-readLines(file.path(export,z$caption_path[i]),warn=FALSE)
      paste0(if(nzchar(z$panel[i]))paste0("(",z$panel[i],") ") else "",paste(caption,collapse=" "))
    },character(1))
    text<-pub_caption_text(c(figure,summaries[[figure]],paragraphs))
    pub_language_check(text,paste(figure,"combined legend"))
    safe_write_lines(text,file.path(folder,paste0(gsub(" ","_",figure),"_legend.txt")))
  }
  invisible(TRUE)
}
