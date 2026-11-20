# Display-only cleaning. Never evaluate a saved plotmath string as R code.
clean_species_label <- function(x) {
 x <- as.character(x)
 wrapped <- !is.na(x) & grepl("^italic\\('.*'\\)$", x)
 x[wrapped] <- substring(x[wrapped], 9L, nchar(x[wrapped])-2L)
 x <- sub("^species\\s*[:|]\\s*", "", x, ignore.case=TRUE, perl=TRUE)
 x <- sub("^s__", "", x, ignore.case=TRUE)
 trimws(gsub("_", " ", x, fixed=TRUE))
}

r2_volcano_data <- function(d, panel) {
 r2_columns(d,c("feature","coef","pval","qval","label_selected","label_rank","comparison"))
 r1_need(!any(c("label_display","display_class","display_minus_log10_q")%in%names(d)),"Saved source collides with display-only columns")
 r1_need(!anyDuplicated(d$feature),"Duplicate volcano feature IDs")
 r1_need(is.logical(d$label_selected)&&!anyNA(d$label_selected),"Saved label_selected must contain explicit TRUE/FALSE")
 expected <- switch(panel$panel,A=c(200,69,10),B=c(192,32,10),C=c(187,10,8),D=c(190,5,5),E=c(193,4,4))
 comparison <- switch(panel$panel,A="T2D vs Control",B="NT2D vs NC",C="OT2D vs OC",D="OT2D vs NT2D",E="OC vs NC")
 r1_need(length(expected)==3L&&all(d$comparison==comparison),"Volcano comparison or direction differs")
 r1_need(all(is.finite(d$coef))&&all(is.finite(d$qval))&&all(d$qval>0&d$qval<=1),"Invalid coefficient/qval; zero q needs explicit display policy")
 r1_need(identical(as.integer(c(nrow(d),sum(d$qval<.05),sum(d$label_selected))),as.integer(expected)),"Saved volcano counts differ")
 # Selection, row order and all source fields are kept; no top/rank selection occurs.
 labels <- if("label_text"%in%names(d)) d$label_text else if("label"%in%names(d))d$label else d$feature
 cleaned <- clean_species_label(labels)
 r1_need(all(!is.na(cleaned[d$label_selected])&nzchar(cleaned[d$label_selected])),"A selected saved label is empty")
 d$label_display <- ifelse(d$label_selected,cleaned,NA_character_)
 d$display_class <- ifelse(d$qval<.05&d$coef>0,"Positive",ifelse(d$qval<.05&d$coef<0,"Negative","NS"))
 d$display_minus_log10_q <- -log10(d$qval)
 valid <- d$label_display[!is.na(d$label_display)]
 r1_need(sum(!is.na(d$label_display))==sum(d$label_selected),"Label selection changed")
 r1_need(!any(grepl("^species\\s*[:|]|^s__",valid,ignore.case=TRUE,perl=TRUE)),"Uncleaned taxonomy prefix in display label")
 d
}

plot_species_volcano <- function(d,panel,root,paths,style) {
 z <- r2_volcano_data(d,panel);s<-style$volcano
 fill<-c(Positive=s$positive_fill,Negative=s$negative_fill,NS=s$nonsignificant_fill)
 border<-c(Positive=s$positive_border,Negative=s$negative_border,NS=s$nonsignificant_border)
 selected<-z[z$label_selected,,drop=FALSE]
 r1_need(all(is.finite(selected$label_rank))&&!anyDuplicated(selected$label_rank),"Selected labels need unique saved finite ranks")
 selected<-selected[order(selected$label_rank),,drop=FALSE] # saved priority only; no re-selection
 duplicate<-!is.na(z$label_display)&(duplicated(z$label_display)|duplicated(z$label_display,fromLast=TRUE))
 source_column<-if("label_text"%in%names(z))"label_text"else if("label"%in%names(z))"label"else"feature"
 audit<-data.frame(feature=z$feature,label_selected=z$label_selected,label_rank=z$label_rank,
  source_label_column=source_column,source_label=z[[source_column]],label_display=z$label_display,duplicate_display_name=duplicate,stringsAsFactors=FALSE)
 p<-ggplot2::ggplot(z,ggplot2::aes(x=coef,y=display_minus_log10_q))+
  ggplot2::geom_vline(xintercept=0,colour=s$zero_line_colour,linewidth=s$zero_line_width)+
  ggplot2::geom_hline(yintercept=-log10(.05),colour=s$threshold_line_colour,linewidth=s$threshold_line_width,linetype="dashed")+
  ggplot2::geom_point(ggplot2::aes(fill=display_class,colour=display_class),shape=s$point_shape,size=s$point_size,alpha=s$point_alpha,stroke=s$point_stroke)+
  ggplot2::scale_fill_manual(values=fill,drop=FALSE)+ggplot2::scale_colour_manual(values=border,drop=FALSE)+
  ggrepel::geom_text_repel(data=selected,ggplot2::aes(label=label_display),family=style$font$family,fontface="italic",
   size=s$label_size,box.padding=s$box_padding,point.padding=s$point_padding,min.segment.length=s$min_segment_length,
   segment.color=s$label_segment_colour,segment.size=s$label_segment_width,max.overlaps=Inf,seed=s$seed,na.rm=TRUE)+
  ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=.22))+ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(.03,.18)))+
  ggplot2::labs(x="MaAsLin2 coefficient",y="−log10(q value)")+theme_round1(style)+
  ggplot2::theme(legend.position=s$legend_position,plot.margin=ggplot2::margin(8,12,8,12))
 p<-r1_titles(p,panel,"",style)
 list(plot=p,display_data=z,audits=setNames(list(audit),paste0(panel$output_filename,"_label_audit.csv")),
  caption="Formal significance: saved MaAsLin2 full-model qval < 0.05. Saved label_selected and label_rank are unchanged. Labels alone have taxonomy prefixes removed; duplicate display names are audited and points are never merged. ggrepel segment.size is the supported line-width argument (0.25).")
}
