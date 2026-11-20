theme_round1 <- function(style) {
 f<-style$font
 ggplot2::theme_classic(base_family=f$family,base_size=f$base_size)+
 ggplot2::theme(panel.grid.major=ggplot2::element_blank(),panel.grid.minor=ggplot2::element_blank(),
  panel.border=ggplot2::element_blank(),axis.line=ggplot2::element_line(colour="black",linewidth=style$lines$axis_linewidth),
  axis.text=ggplot2::element_text(colour="black",size=f$axis_text_size),axis.title=ggplot2::element_text(size=f$axis_title_size),
  plot.title=ggplot2::element_text(size=f$title_size,face="bold",hjust=.5),plot.subtitle=ggplot2::element_text(size=f$subtitle_size,hjust=.5),
  legend.title=ggplot2::element_text(size=f$legend_title_size),legend.text=ggplot2::element_text(size=f$legend_text_size),
  legend.background=ggplot2::element_blank(),legend.key=ggplot2::element_blank(),strip.background=ggplot2::element_rect(fill=style$export$background,colour=NA),
  strip.text=ggplot2::element_text(size=f$facet_text_size),plot.background=ggplot2::element_rect(fill=style$export$background,colour=NA))
}

r1_titles <- function(p,panel,subtitle,style) {
 title<-paste(strwrap(panel$title,width=if(panel$width_mm<=style$layout$narrow_max_mm)style$layout$title_wrap_narrow else style$layout$title_wrap_wide),collapse="\n")
 subtitle<-gsub(" | ","\n",subtitle,fixed=TRUE)
 p+ggplot2::labs(title=title,subtitle=subtitle)
}
