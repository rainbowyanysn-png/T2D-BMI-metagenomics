plot_lollipop <- function(d,panel,root,paths,style) {
 r2_columns(d,c("feature","coef","qval","phylum","label"))
 r1_need(nrow(d)==20L&&!anyDuplicated(d$feature)&&all(d$qval>0&d$qval<.05),"Saved Top 20 differs")
 count<-table(factor(d$phylum,levels=c("Firmicutes","Bacteroidetes","Proteobacteria")))
 r1_need(!anyNA(match(d$phylum,names(count)))&&identical(as.integer(count),c(17L,2L,1L)),"Saved phylum annotation differs")
 # Display order only: all saved Top20 rows and statistics remain unchanged.
 display_order<-order(-d$coef,seq_len(nrow(d)))
 z<-d;z$display_feature<-factor(z$feature,levels=rev(d$feature[display_order]));z$display_minus_log10_q<--log10(z$qval)
 s<-style$lollipop;palette<-unlist(s$phylum_palette)
 r1_need(all(z$phylum%in%names(palette)),"Lollipop colour not mapped")
 p<-ggplot2::ggplot(z,ggplot2::aes(y=display_feature,x=coef,colour=phylum))+
  ggplot2::geom_vline(xintercept=0,colour=s$zero_line_colour,linetype=s$zero_line_type,linewidth=s$zero_line_width)+
  ggplot2::geom_segment(ggplot2::aes(x=0,xend=coef,yend=display_feature),linewidth=s$segment_width,alpha=s$segment_alpha)+
  ggplot2::geom_point(ggplot2::aes(size=display_minus_log10_q),alpha=s$point_alpha)+
  ggplot2::scale_colour_manual(values=palette)+ggplot2::scale_size_continuous(range=unlist(s$point_size_range))+
  ggplot2::scale_y_discrete(labels=setNames(d$label,d$feature))+
  ggplot2::labs(x="MaAsLin2 coefficient",y=NULL,colour="Phylum",size="−log10(q value)")+theme_round1(style)+
  ggplot2::theme(axis.text.y=ggplot2::element_text(face="italic"),legend.position="right")
 panel$title<-"Top 20 T2D-associated Species by Absolute Effect Size";p<-r1_titles(p,panel,"",style)
 list(plot=p,display_data=d,audits=list(Fig3D_saved_order.csv=data.frame(source_position=seq_len(nrow(d)),feature=d$feature,phylum=d$phylum),
  Fig3D_display_order.csv=data.frame(position=seq_along(display_order),source_position=display_order,feature=d$feature[display_order],coef=d$coef[display_order])),
  caption="Saved RC3 Top 20 selection and phylum annotations are unchanged. Displayed from top to bottom in descending coefficient order (positive to negative); saved source order is retained in the source CSV and audited separately. No feature was reselected.")
}
