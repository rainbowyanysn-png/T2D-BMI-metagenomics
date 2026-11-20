plot_cohort_permanova <- function(d,panel,root,paths,style) {
 required<-c("study_name","R2","F","p_value","n_samples","model_formula","n_permutations")
 r1_need(all(required%in%names(d)),"Frozen cohort PERMANOVA fields missing")
 order<-names(unlist(style$ordination$cohort_shapes))
 r1_need(nrow(d)==5L&&!anyDuplicated(d$study_name)&&setequal(d$study_name,order),"Expected exactly the five fixed cohorts")
 r1_need(all(d$n_permutations==999L)&&sum(d$n_samples)==1182L&&all(grepl("analysis_group",d$model_formula,fixed=TRUE)),"Frozen cohort model/permutation/sample definition differs")
 r1_need(all(is.finite(d$R2)&d$R2>=0&d$R2<=1)&&all(is.finite(d$p_value)&d$p_value>0&d$p_value<=1),"Invalid saved cohort R2/P")
 z<-d;z$display_study<-factor(d$study_name,levels=rev(order))
 z$display_label<-paste0(vapply(d$p_value,format_probability,character(1)),"; N = ",d$n_samples)
 offset<-max(d$R2)*.06
 p<-ggplot2::ggplot(z,ggplot2::aes(R2,display_study))+
  ggplot2::geom_segment(ggplot2::aes(x=0,xend=R2,yend=display_study),colour="#D9D9D9",linewidth=.35)+
  ggplot2::geom_point(colour="#2166AC",shape=16,size=2.1)+
  ggplot2::geom_text(ggplot2::aes(x=R2+offset,label=display_label),hjust=0,size=2.5,family=style$font$family)+
  ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(.02,.85)))+
  ggplot2::labs(x="Marginal PERMANOVA R²",y=NULL)+theme_round1(style)+
  ggplot2::theme(legend.position="none",axis.text.y=ggplot2::element_text(size=7))
 p<-r1_titles(p,panel,"",style)
 list(plot=p,annotation="Saved cohort-specific four-group marginal PERMANOVA; exact P values and sample sizes",
  statistics=d[required],n_samples=sum(d$n_samples),display_data=z,
  caption="Cohort-specific marginal PERMANOVA based on uncorrected CZM/CLR Aitchison distances with 999 permutations. Points show the saved analysis-group marginal R²; reference stems are not confidence intervals. Saved P values and cohort sample counts are shown. No ordination or permutation test was recomputed.")
}
