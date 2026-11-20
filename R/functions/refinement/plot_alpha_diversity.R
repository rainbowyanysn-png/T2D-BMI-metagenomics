plot_alpha_diversity <- function(d,panel,root,paths,style) {
 r1_need(nrow(d)==1186L&&!anyDuplicated(d$sample_id)&&all(is.finite(d$value)),"Alpha data must retain all1186 participants")
 value<-if(panel$metric=="Shannon")"shannon"else"richness"
 r1_need(identical(d$value,d[[value]]),"Saved alpha metric/value mismatch")
 total<-panel$plot_family=="alpha_total";cohort<-panel$plot_family=="alpha_cohort"
 palette<-r1_group_vectors(style,total);z<-d;z$group<-r1_factor(z$group,palette$order)
 ann<-r1_alpha_annotation(root,paths,panel$metric,panel$plot_family)
 box<-style$boxplot
 if(cohort){
  st<-style$cohort_boxplot;z$study_name<-factor(z$study_name,levels=unique(z$study_name))
  p<-ggplot2::ggplot(z,ggplot2::aes(study_name,value,fill=group,colour=group))+
   ggplot2::geom_boxplot(width=st$box_width,position=ggplot2::position_dodge(width=st$dodge_width),alpha=box$alpha,linewidth=box$linewidth,coef=style$whisker_coef,outlier.shape=NA)+
   ggplot2::geom_point(position=ggplot2::position_jitterdodge(jitter.width=st$jitter_width,jitter.height=0,dodge.width=st$dodge_width,seed=style$jitter_seed),size=st$point_size,alpha=st$point_alpha)
 }else{
  p<-ggplot2::ggplot(z,ggplot2::aes(group,value,fill=group,colour=group))+
   ggplot2::geom_boxplot(width=box$width,alpha=box$alpha,linewidth=box$linewidth,coef=style$whisker_coef,outlier.shape=NA)+
   ggplot2::geom_point(position=ggplot2::position_jitter(width=box$jitter_width,height=0,seed=style$jitter_seed),size=box$point_size,alpha=box$point_alpha)
 }
 p<-p+ggplot2::scale_fill_manual(values=palette$fill,limits=palette$order,drop=FALSE)+
  ggplot2::scale_colour_manual(values=palette$dark,limits=palette$order,drop=FALSE)+
  ggplot2::labs(title=panel$title,subtitle=ann$text,x=NULL,y=if(panel$metric=="Shannon")"Shannon diversity"else"Species richness",fill="Group",colour="Group")+
  theme_round1(style)+ggplot2::theme(legend.position=if(cohort)"right"else"none")
 if(cohort)p<-p+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=style$cohort_boxplot$x_text_angle,hjust=1))
 p<-r1_titles(p,panel,ann$text,style)
 list(plot=p,annotation=ann$text,statistics=ann$statistics,n_samples=nrow(d),display_data=z,
  caption="All participants are shown. Boxes indicate median/IQR; whiskers use1.5 x IQR. Saved pooled LMM P values and ICC; no pairwise brackets or new inference.")
}
