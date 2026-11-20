r2_stars <- function(q) ifelse(is.na(q),"",ifelse(q<.001,"***",ifelse(q<.01,"**",ifelse(q<.05,"*",""))))

r2_heatmap_layout <- function(d,features,columns,value,cluster_rows,cluster_cols,style) {
 r1_need(!anyDuplicated(paste(d$feature,d$display_column,sep="\r")),"Duplicate heatmap cell")
 r1_need(nrow(d)==length(features)*length(columns),"Heatmap grid incomplete")
 m<-matrix(NA_real_,length(features),length(columns),dimnames=list(features,columns))
 ij<-cbind(match(d$feature,features),match(d$display_column,columns));r1_need(!anyNA(ij),"Unknown heatmap row/column")
 m[ij]<-d[[value]]
 # ONLY the temporary clustering copy replaces NA; m and saved coefficients stay NA.
 clustering_copy<-m;clustering_copy[is.na(clustering_copy)]<-0
 cluster<-function(x)stats::hclust(stats::dist(x,method=style$clustering$distance),method=style$clustering$linkage)
 row_tree<-if(cluster_rows)cluster(clustering_copy)else NULL
 column_tree<-if(cluster_cols)cluster(t(clustering_copy))else NULL
 ri<-if(cluster_rows)row_tree$order else seq_along(features)
 ci<-if(cluster_cols)column_tree$order else seq_along(columns)
 list(matrix=m,row_order=features[ri],column_order=columns[ci],row_tree=row_tree,column_tree=column_tree)
}

r2_dendrogram_segments <- function(layout) {
 tree<-layout$row_tree;r1_need(!is.null(tree),"A true saved-display clustering tree is required")
 features<-rownames(layout$matrix);n<-length(features)
 # Leaves and tiles share the exact same feature-to-y mapping. No invented tree.
 leaf_y<-setNames(n-seq_len(n)+1L,layout$row_order)
 node_x<-node_y<-numeric(n-1L);parts<-vector("list",n-1L)
 for(i in seq_len(n-1L)){
  children<-tree$merge[i,]
  xx<-vapply(children,function(k)if(k<0)0 else node_x[k],numeric(1))
  yy<-vapply(children,function(k)if(k<0)unname(leaf_y[features[-k]])else node_y[k],numeric(1))
  x<--tree$height[i];node_x[i]<-x;node_y[i]<-mean(yy)
  parts[[i]]<-data.frame(x=c(xx[1],xx[2],x),xend=c(x,x,x),y=c(yy[1],yy[2],yy[1]),yend=c(yy[1],yy[2],yy[2]))
 }
 do.call(rbind,parts)
}

r2_heatmap_plot <- function(z,layout,labels,limits,s,style,panel,legend=NULL) {
 r1_need(length(limits)==2L&&limits[1]<0&&limits[2]>0,"Heatmap limits must bracket zero")
 n<-length(layout$row_order)
 z$display_y<-n-match(z$feature,layout$row_order)+1L
 z$display_x<-match(z$display_column,layout$column_order)
 segments<-r2_dendrogram_segments(layout)
 tree<-ggplot2::ggplot(segments)+
  ggplot2::geom_segment(ggplot2::aes(x=x,y=y,xend=xend,yend=yend),linewidth=.25,colour="#444444")+
  ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(.04,0)))+
  ggplot2::scale_y_continuous(limits=c(.5,n+.5),expand=c(0,0))+
  ggplot2::theme_void(base_family=style$font$family)+
  ggplot2::theme(plot.background=ggplot2::element_rect(fill="white",colour=NA),plot.margin=ggplot2::margin(5.5,0,5.5,0))
 p<-ggplot2::ggplot(z,ggplot2::aes(x=display_x,y=display_y,fill=display_value))+
  ggplot2::geom_tile(colour=s$border_colour,linewidth=s$border_width)+
  ggplot2::geom_text(ggplot2::aes(label=display_stars),size=s$star_size,family=style$font$family,na.rm=TRUE)+
  ggplot2::scale_fill_gradientn(colours=c(s$low,s$mid,s$high),values=scales::rescale(c(limits[1],0,limits[2]),from=limits),
   limits=limits,breaks=unlist(s$legend_breaks),oob=scales::squish,na.value=if(is.null(s$na_colour))"#E5E5E5"else s$na_colour,name=NULL)+
  ggplot2::scale_y_continuous(breaks=seq_len(n),labels=unname(labels[rev(layout$row_order)]),limits=c(.5,n+.5),expand=c(0,0),position="right")+
  ggplot2::scale_x_continuous(breaks=seq_along(layout$column_order),labels=layout$column_order,limits=c(.5,length(layout$column_order)+.5),expand=c(0,0))+
  ggplot2::guides(fill=ggplot2::guide_colourbar(title=NULL,barheight=grid::unit(24,"mm"),barwidth=grid::unit(3,"mm")))+
  ggplot2::labs(x=NULL,y=NULL,fill=NULL)+theme_round1(style)+
  ggplot2::theme(axis.line=ggplot2::element_blank(),axis.ticks=ggplot2::element_blank(),
   axis.text.y=ggplot2::element_text(face="italic",size=if(nrow(layout$matrix)>30)6.5 else style$font$axis_text_size),
   axis.text.x=ggplot2::element_text(angle=s$column_angle,hjust=1,face="bold"),
   legend.position="right",legend.justification=c(0,1),legend.box.just="top",legend.title=ggplot2::element_blank(),legend.text=ggplot2::element_text(size=7),
   plot.margin=ggplot2::margin(5.5,5.5,5.5,0))
 title<-paste(strwrap(panel$title,width=style$layout$title_wrap_wide),collapse="\n")
 patchwork::wrap_plots(tree,p,ncol=2,widths=c(s$tree_width_fraction,1))+
  patchwork::plot_annotation(title=title,theme=ggplot2::theme(text=ggplot2::element_text(family=style$font$family),
   plot.title=ggplot2::element_text(size=style$font$title_size,face="bold",hjust=.5),plot.background=ggplot2::element_rect(fill="white",colour=NA)))
}

plot_species_heatmap <- function(d,panel,root,paths,style) {
 r2_columns(d,c("feature","comparison","coef","formal_q","tested","label"))
 r1_need(nrow(d)==395L&&length(unique(d$feature))==79L,"Saved species heatmap size differs")
 columns<-unlist(style$round2$comparisons);r1_need(setequal(unique(d$comparison),columns),"Heatmap comparisons differ")
 r1_need(is.logical(d$tested)&&!anyNA(d$tested)&&all(is.na(d$coef[!d$tested]))&&all(is.na(d$formal_q[!d$tested])),"Untested cells must retain NA")
 r1_need(all(is.finite(d$coef[d$tested]))&&all(is.finite(d$formal_q[d$tested]))&&all(d$formal_q[d$tested]>=0&d$formal_q[d$tested]<=1),"Invalid saved tested cell")
 z<-d;z$display_column<-z$comparison;z$display_value<-z$coef;z$display_stars<-r2_stars(z$formal_q)
 s<-style$species_heatmap;r1_need(isTRUE(s$row_clustering)&&identical(s$column_clustering,FALSE),"Species heatmap clustering policy differs")
 layout<-r2_heatmap_layout(z,unique(d$feature),columns,"coef",TRUE,FALSE,style)
 limits<-unlist(s$colour_limits);r1_need(identical(as.numeric(limits),c(-2,2)),"Species heatmap display limits differ")
 rows<-match(layout$row_order,d$feature);labels<-setNames(d$label[rows],d$feature[rows])
 panel$title<-"Differential Abundance Coefficients Across Comparisons"
 p<-r2_heatmap_plot(z,layout,labels,limits,s,style,panel)
 list(plot=p,display_data=d,display_layout=layout,audits=list(
  Fig3E_row_order.csv=data.frame(position=seq_along(layout$row_order),feature=layout$row_order,label=unname(labels)),
  Fig3E_column_order.csv=data.frame(position=seq_along(layout$column_order),comparison=layout$column_order),
  Fig3E_colour_scale.csv=data.frame(lower=limits[1],midpoint=0,upper=limits[2],source_column="coef",underlying_truncated=FALSE),
  Fig3E_dendrogram_segments.csv=r2_dendrogram_segments(layout),
  Fig3E_cell_annotations.csv=data.frame(feature=d$feature,comparison=d$comparison,coef=d$coef,formal_q=d$formal_q,tested=d$tested,stars=z$display_stars)),
  caption="Display colour limits are -2 to 2; colours saturate outside these limits but underlying coefficients are not truncated. The left dendrogram and right-side species names share the same actual clustering leaf order, using the 0705 Euclidean/complete-linkage approach on current RC3 values only. Untested cells remain NA and grey; only the temporary clustering copy uses zero for NA. Columns retain the five-comparison order. Stars use saved full-model formal_q: * <0.05, ** <0.01, *** <0.001.")
}
