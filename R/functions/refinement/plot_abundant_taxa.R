r1_abundant_annotations <- function(d,panel,root,paths,style) {
 level<-if(panel$plot_family=="abundant_species")"species"else"phylum"
 nfeatures<-if(level=="species")8L else 4L;family_size<-6L*nfeatures
 saved<-r1_read(root,paths,"abundant_taxa_emmeans_pairwise_all.csv")
 required<-c("feature","taxonomic_level","model_formula","n_samples","n_features","family_size","group1","group2","p_raw","p_adj_BH","formal_significant","formal_significance_column")
 r1_need(all(required%in%names(saved)),"Saved abundant-taxa pairwise columns missing")
 a<-saved[saved$taxonomic_level==level,,drop=FALSE]
 features<-unique(as.character(d$feature))
 map<-c(Non_Overweight_Control="NC",Overweight_Control="OC",Non_Overweight_T2D="NT2D",Overweight_T2D="OT2D")
 r1_need(nrow(a)==family_size&&setequal(a$feature,features)&&all(table(a$feature)==6L),"Saved pairwise family does not match the plotted Top4/Top8")
 r1_need(all(a$n_samples==1186L)&&all(a$n_features==nfeatures)&&all(a$family_size==family_size),"Abundant-taxa saved family/sample counts differ")
 r1_need(all(a$model_formula=="Abundance ~ analysis_group + age + sex + (1 | study_name)"),"Unexpected saved abundant-taxa model")
 r1_need(all(a$formal_significance_column=="p_adj_BH")&&all(d$formal_significance_column=="p_adj_BH"),"Do not replace the frozen abundant-taxa significance definition")
 r1_need(all(is.finite(a$p_raw)&a$p_raw>=0&a$p_raw<=1)&&all(is.finite(a$p_adj_BH)&a$p_adj_BH>=0&a$p_adj_BH<=1),"Invalid saved raw/adjusted P")
 r1_need(!anyNA(a$formal_significant)&&all(a$formal_significant==(a$p_adj_BH<.05)),"Saved significance flags disagree with p_adj_BH")
 a$display_group1<-unname(map[a$group1]);a$display_group2<-unname(map[a$group2])
 group_order<-unlist(style$groups$order);a$display_x1<-match(a$display_group1,group_order);a$display_x2<-match(a$display_group2,group_order)
 r1_need(!anyNA(a$display_x1)&&!anyNA(a$display_x2)&&all(a$display_x1!=a$display_x2),"Unmapped saved pairwise groups")
 pair_key<-paste(a$feature,pmin(a$display_x1,a$display_x2),pmax(a$display_x1,a$display_x2),sep="|")
 r1_need(!anyDuplicated(pair_key),"Duplicate saved unordered pairwise comparison")
 s<-style$abundant_taxa
 r1_need(s$annotation_p_column=="p_adj_BH"&&s$annotation_alpha==.05,"Frozen bracket column/threshold changed")
 a$display_selected<-a$p_adj_BH<s$annotation_alpha
 a$display_label<-vapply(a$p_adj_BH,function(p)format_probability(p,"BH P"),character(1))
 a$display_y<-a$display_tip_y<-a$display_label_y<-NA_real_
 for(feature in features){
  values<-d$raw_abundance[as.character(d$feature)==feature]
  span<-diff(range(values));if(span==0)span<-max(abs(values),1)
  i<-which(a$feature==feature&a$display_selected)
  if(length(i)){
   i<-i[order(abs(a$display_x2[i]-a$display_x1[i]),a$display_x1[i],a$display_x2[i])]
   y<-max(values)+span*(s$bracket_base_fraction+s$bracket_step_fraction*(seq_along(i)-1L))
   a$display_y[i]<-y;a$display_tip_y[i]<-y-span*s$bracket_tip_fraction
   a$display_label_y[i]<-y+span*s$bracket_label_fraction
  }
 }
 a
}

plot_abundant_taxa <- function(d,panel,root,paths,style) {
 species<-panel$plot_family=="abundant_species";nfeatures<-if(species)8L else 4L
 r1_need(nrow(d)==1186L*nfeatures&&length(unique(d$feature))==nfeatures&&!anyDuplicated(d[,c("sample_id","feature")]),"Abundant-taxa input panel differs")
 r1_need(all(table(d$feature)==1186L)&&all(is.finite(d$raw_abundance)),"Abundant-taxa sample values incomplete")
 r1_need(all(grepl("percentage",d$transformation)),"Abundance units not established by saved input")
 map<-c(Non_Overweight_Control="NC",Overweight_Control="OC",Non_Overweight_T2D="NT2D",Overweight_T2D="OT2D")
 z<-d;palette<-r1_group_vectors(style,FALSE);z$group<-r1_factor(unname(map[as.character(d$analysis_group)]),palette$order)
 levels<-unique(z$feature);labels<-if(species)gsub("[._]+"," ",sub("^species[.]","",levels))else levels
 z$feature<-factor(z$feature,levels=levels)
 st<-style$abundant_taxa;box<-style$boxplot
 statistics<-r1_abundant_annotations(d,panel,root,paths,style)
 brackets<-statistics[statistics$display_selected,,drop=FALSE]
 brackets$feature<-factor(brackets$feature,levels=levels)
 p<-ggplot2::ggplot(z,ggplot2::aes(group,raw_abundance,fill=group,colour=group))+
  ggplot2::geom_boxplot(width=st$box_width,alpha=box$alpha,linewidth=box$linewidth,coef=style$whisker_coef,outlier.shape=NA)+
  ggplot2::geom_point(position=ggplot2::position_jitter(width=box$jitter_width,height=0,seed=style$jitter_seed),size=st$point_size,alpha=st$point_alpha)+
  ggplot2::facet_wrap(~feature,ncol=if(species)st$facet_ncol_species else st$facet_ncol_phylum,scales="free_y",labeller=ggplot2::as_labeller(setNames(labels,levels)))+
  ggplot2::scale_fill_manual(values=palette$fill,limits=palette$order,drop=FALSE)+ggplot2::scale_colour_manual(values=palette$dark,limits=palette$order,drop=FALSE)+
  ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(.03,.12)))+
  ggplot2::labs(title=NULL,subtitle=NULL,x=NULL,y="Relative abundance (%)")+
  theme_round1(style)+ggplot2::theme(strip.text=ggplot2::element_text(face=if(species)"italic"else"plain"),legend.position="none")
 if(nrow(brackets))p<-p+
  ggplot2::geom_segment(data=brackets,ggplot2::aes(x=display_x1,xend=display_x2,y=display_y,yend=display_y),inherit.aes=FALSE,colour=st$bracket_colour,linewidth=st$bracket_linewidth)+
  ggplot2::geom_segment(data=brackets,ggplot2::aes(x=display_x1,xend=display_x1,y=display_y,yend=display_tip_y),inherit.aes=FALSE,colour=st$bracket_colour,linewidth=st$bracket_linewidth)+
  ggplot2::geom_segment(data=brackets,ggplot2::aes(x=display_x2,xend=display_x2,y=display_y,yend=display_tip_y),inherit.aes=FALSE,colour=st$bracket_colour,linewidth=st$bracket_linewidth)+
  ggplot2::geom_text(data=brackets,ggplot2::aes(x=(display_x1+display_x2)/2,y=display_label_y,label=display_label),inherit.aes=FALSE,colour=st$bracket_colour,family=style$font$family,size=st$bracket_text_size,vjust=0)
 annotation<-paste0("Secondary/exploratory; saved BH-adjusted P < 0.05 brackets (",if(species)48 else 24," comparisons per taxonomic level)")
 list(plot=p,annotation=annotation,statistics=statistics,n_samples=1186L,display_data=z,
  caption="Saved high-abundance taxa and all sample-level percentage abundances are unchanged. Secondary/exploratory LMM pairwise brackets use saved p_adj_BH < 0.05 only: BH P denotes the BH-adjusted P value across all24 phylum or48 species comparisons. Full raw P, adjusted P and historical audit columns are saved in the annotation table. No LMM/emmeans refit or P/BH recalculation; formal species inference remains MaAsLin2.")
}
