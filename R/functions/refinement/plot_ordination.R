plot_ordination <- function(d,panel,root,paths,style) {
 r1_need(nrow(d)==1182L&&!anyDuplicated(d$sample_id)&&all(is.finite(d$PC1))&&all(is.finite(d$PC2)),"Saved PCoA must retain all1182 coordinates")
 total<-panel$plot_family%in%c("ordination_total","combat_total");combat<-startsWith(panel$plot_family,"combat")
 palette<-r1_group_vectors(style,total);z<-d;z$phenotype<-r1_factor(z$phenotype,palette$order)
 if(!combat){
  palette$fill<-unlist(if(total)style$ordination$raw_total_fill else style$ordination$raw_four_fill)
  r1_need(identical(names(palette$fill),palette$order),"Raw beta point colours do not match group levels")
  palette$dark<-setNames(colorspace::darken(unname(palette$fill),amount=style$colour_adjustment$darken_amount),palette$order)
 }
 shapes<-unlist(style$ordination$cohort_shapes);r1_need(all(z$study_name%in%names(shapes)),"Unmapped cohort shape")
 z$study_name<-factor(z$study_name,levels=names(shapes))
 if(combat){ann<-list(text="Descriptive visualization only",statistics=data.frame(inference="Not displayed for ComBat"))}else ann<-r1_beta_annotation(root,paths,total)
 fractions<-all(c("PC1_fraction","PC2_fraction")%in%names(d))
 if(fractions){
  r1_need(length(unique(d$PC1_fraction))==1L&&length(unique(d$PC2_fraction))==1L&&all(d$PC1_fraction>=0&d$PC1_fraction<=1)&all(d$PC2_fraction>=0&d$PC2_fraction<=1),"Invalid saved coordinate fractions")
  labels<-sprintf("PCoA %d (%.1f%%)",1:2,100*c(d$PC1_fraction[1],d$PC2_fraction[1]))
 }else{r1_need(combat,"Formal raw ordination requires saved axis fractions");labels<-c("PCoA 1","PCoA 2")}
 st<-style$ordination
 p<-ggplot2::ggplot(z,ggplot2::aes(PC1,PC2,colour=phenotype))+
  ggplot2::geom_point(ggplot2::aes(shape=study_name),size=st$point_size,alpha=st$point_alpha,key_glyph=ggplot2::draw_key_point)
 # Constant darkened group colours avoid adding an ellipse/line legend.
 for(g in palette$order)p<-p+ggplot2::stat_ellipse(data=z[z$phenotype==g,,drop=FALSE],
  ggplot2::aes(group=phenotype),colour=unname(palette$dark[g]),type=st$ellipse_type,level=st$ellipse_level,
  linetype=st$ellipse_linetype,linewidth=st$ellipse_linewidth,show.legend=FALSE)
 p<-p+ggplot2::scale_colour_manual(values=if(combat)palette$dark else palette$fill,limits=palette$order,drop=FALSE)+
  ggplot2::scale_shape_manual(values=shapes,drop=FALSE)+
  ggplot2::guides(colour=ggplot2::guide_legend(order=1,ncol=1,override.aes=list(shape=16,alpha=1,size=1.3)),shape=ggplot2::guide_legend(order=2,ncol=1,override.aes=list(colour="black",alpha=1,size=1.3)))+
  ggplot2::labs(title=panel$title,subtitle=ann$text,x=labels[1],y=labels[2],colour="Group",shape="Cohort")+
  theme_round1(style)+ggplot2::theme(legend.position=st$legend_position,
   plot.subtitle=ggplot2::element_text(size=st$subtitle_size,lineheight=st$subtitle_lineheight,hjust=.5,margin=ggplot2::margin(b=3)),
   legend.text=ggplot2::element_text(size=st$legend_text_size,lineheight=st$legend_lineheight),
   legend.title=ggplot2::element_text(size=st$legend_title_size,margin=ggplot2::margin(b=1)),
   legend.key.height=grid::unit(st$legend_key_height_mm,"mm"),legend.key.width=grid::unit(st$legend_key_width_mm,"mm"),
   legend.key.spacing.y=grid::unit(st$legend_row_spacing_mm,"mm"),legend.spacing.y=grid::unit(st$legend_guide_spacing_mm,"mm"),
   legend.box.spacing=grid::unit(st$legend_box_spacing_mm,"mm"),legend.margin=ggplot2::margin(0,0,0,0))
 caption<-if(combat)"ComBat-corrected coordinates are presented for descriptive visualization only. Formal beta-diversity inference was based on uncorrected CZM/CLR Aitchison distances with study-aware analyses."else
  paste("Saved uncorrected CZM/CLR Aitchison coordinates. Formal labels use the saved target-term marginal PERMANOVA, not full-model or study R2. There are1182 saved samples and999 within-study permutations. Ellipses are95% t-type descriptive envelopes drawn from saved coordinates; no ordination or permutation model was recomputed.",r1_dispersion_caption(root,paths))
 if(!fractions)caption<-paste(caption,"Axis fractions are absent from the saved ComBat plot-data; percentages are omitted rather than recalculated or borrowed.")
 p<-r1_titles(p,panel,ann$text,style)
 list(plot=p,annotation=ann$text,statistics=ann$statistics,n_samples=nrow(d),display_data=z,axis_fraction_status=if(fractions)"SAVED"else"NOT_SAVED",caption=caption)
}
