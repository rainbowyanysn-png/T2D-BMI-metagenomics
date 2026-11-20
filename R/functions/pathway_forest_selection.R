# Display ranking of current-run estimates only; no model/CI/q-value refit.
final_forest_rule <- function() paste(
  "Representative pathway features were selected for visualization from the formal REML-significant set using a predefined direction-concordance and cohort-coverage rule.",
  "Eligibility: padj<0.05, at least four contributing cohorts, finite pooled effect/CI/padj/I2 and valid cohort effects with positive finite SE; exclude UNINTEGRATED/UNMAPPED.",
  "Order: descending direction concordance, descending cohort count, ascending I2, ascending maximum cohort SE, ascending padj, then pathway ID.",
  "Every valid contributing cohort is displayed; this does not imply independent significance in every cohort.")

final_pathway_forest <- function(b,panel) {
  comparison<-if(panel=="C1")"T2D vs Control"else"NT2D vs NC"
  p<-b$reml[b$reml$comparison==comparison & !is.na(b$reml$padj) & b$reml$padj<.05,,drop=FALSE]
  c<-b$cohort_effects[b$cohort_effects$comparison==comparison,,drop=FALSE]
  assert_columns(p,c("pathway","combined_log2FC","ci_lower","ci_upper","padj","Heterogeneity_I2","n_studies","stderr"),"REML forest")
  assert_columns(c,c("pathway","study","log2FoldChange","lfcSE"),"Forest cohort effects")
  assert_true(nrow(p)==if(panel=="C1")13L else 133L,"Forest significant family changed")
  assert_true(!anyDuplicated(p$pathway)&&!anyDuplicated(c[c("pathway","study")]),"Duplicate forest keys")
  ranking<-do.call(rbind,lapply(seq_len(nrow(p)),function(i) {
    one<-c[c$pathway==p$pathway[i],,drop=FALSE]
    valid<-is.finite(one$log2FoldChange)&is.finite(one$lfcSE)&one$lfcSE>0
    one<-one[valid,,drop=FALSE];k<-nrow(one)
    assert_true(k==p$n_studies[i],paste("Saved contributing cohort membership differs:",p$pathway[i]))
    finite<-all(is.finite(unlist(p[i,c("combined_log2FC","ci_lower","ci_upper","padj","Heterogeneity_I2","stderr")],use.names=FALSE)))
    data.frame(comparison=comparison,pathway=p$pathway[i],n_studies=k,
      direction_concordance=if(k)mean(sign(one$log2FoldChange)==sign(p$combined_log2FC[i]))else NA_real_,
      I2=p$Heterogeneity_I2[i],maximum_lfcSE=if(k)max(one$lfcSE)else Inf,padj=p$padj[i],
      eligible=finite&&k>=4L&&!grepl("UNINTEGRATED|UNMAPPED",p$pathway[i],ignore.case=TRUE),display_rank=NA_integer_,selected=FALSE)
  }))
  eligible<-which(ranking$eligible);assert_true(length(eligible)>0L,"No eligible forest feature")
  eligible<-eligible[order(-ranking$direction_concordance[eligible],-ranking$n_studies[eligible],ranking$I2[eligible],
    ranking$maximum_lfcSE[eligible],ranking$padj[eligible],enc2utf8(ranking$pathway[eligible]),method="radix")]
  ranking$display_rank[eligible]<-seq_along(eligible);ranking$selected[eligible[1]]<-TRUE
  selected<-ranking$pathway[ranking$selected]
  # Acceptance assertions only, after ranking every candidate.
  if(panel=="C1")assert_true(identical(selected,"P162-PWY: L-glutamate degradation V (via hydroxyglutarate)"),paste("Automatic overall forest selection differs:",selected))else
    assert_true(identical(selected,"VALSYN-PWY: L-valine biosynthesis|g__Blautia.s__Blautia_wexlerae"),paste("Automatic NT2D forest selection differs:",selected))
  pooled<-p[p$pathway==selected,,drop=FALSE];one<-c[c$pathway==selected & is.finite(c$log2FoldChange)&is.finite(c$lfcSE)&c$lfcSE>0,,drop=FALSE]
  assert_true(nrow(one)==5L&&all(one$log2FoldChange>0)&&pooled$combined_log2FC>0,"Expected five positive contributing effects")
  values<-data.frame(study=c(one$study,"REML pooled"),coef=c(one$log2FoldChange,pooled$combined_log2FC),
    stderr=c(one$lfcSE,pooled$stderr),ci_lower=c(one$log2FoldChange-1.96*one$lfcSE,pooled$ci_lower),
    ci_upper=c(one$log2FoldChange+1.96*one$lfcSE,pooled$ci_upper),pathway=selected,I2=pooled$Heterogeneity_I2,padj=pooled$padj,pooled=c(rep(FALSE,5),TRUE))
  list(values=values,ranking=ranking,comparison=comparison,pooled=pooled,cohorts=one)
}
