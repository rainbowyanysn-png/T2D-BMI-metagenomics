plot_clinical_heatmap <- function(d,panel,root,paths,style) {
 r2_columns(d,c("feature","feature_rank","display_name","clinical_trait","spearman_rho","q_bh_60","stars"))
 r1_need(nrow(d)==60L&&length(unique(d$feature))==20L&&setequal(unique(d$clinical_trait),c("BMI","age","sex")),"Clinical 20-by-3 family differs")
 r1_need(all(is.finite(d$spearman_rho)&abs(d$spearman_rho)<=1)&&all(is.finite(d$q_bh_60)&d$q_bh_60>=0&d$q_bh_60<=1),"Invalid saved correlation/q")
 saved_stars<-as.character(d$stars);saved_stars[is.na(saved_stars)]<-""
 r1_need(identical(saved_stars,r2_stars(d$q_bh_60))&&sum(d$q_bh_60<.05)==24L,"Saved stars or significant count differs")
 z<-d;z$display_column<-unname(c(BMI="BMI",age="Age",sex="Sex")[z$clinical_trait]);z$display_value<-z$spearman_rho;z$display_stars<-saved_stars
 first<-!duplicated(d$feature);features<-d$feature[first][order(d$feature_rank[first])]
 s<-style$clinical_heatmap;layout<-r2_heatmap_layout(z,features,c("BMI","Age","Sex"),"spearman_rho",s$row_clustering,s$column_clustering,style)
 limits<-unlist(s$colour_limits);r1_need(identical(as.numeric(limits),c(-.2,.15)),"Clinical heatmap display limits differ")
 labels<-setNames(d$display_name[match(layout$row_order,d$feature)],layout$row_order)
 panel$title<-"Microbial Associations with Clinical Phenotypes"
 p<-r2_heatmap_plot(z,layout,labels,limits,s,style,panel)
 list(plot=p,display_data=d,display_layout=layout,audits=list(
  Fig4F_row_order.csv=data.frame(position=seq_along(layout$row_order),feature=layout$row_order,display_name=unname(labels)),
  Fig4F_column_order.csv=data.frame(position=seq_along(layout$column_order),clinical_trait=layout$column_order),
  Fig4F_colour_scale.csv=data.frame(lower=limits[1],midpoint=0,upper=limits[2],source_column="spearman_rho",underlying_truncated=FALSE),
  Fig4F_dendrogram_segments.csv=r2_dendrogram_segments(layout),
  Fig4F_cell_annotations.csv=data.frame(feature=d$feature,clinical_trait=d$clinical_trait,q_bh_60=d$q_bh_60,stars=saved_stars)),
  caption="Exploratory Spearman correlations; significance was based on BH-adjusted q values across 60 tests. Saved correlations, q_bh_60 and stars are unchanged. Display colour limits are -0.20 to 0.15, with white anchored at zero; colour saturation does not truncate underlying rho values. The left dendrogram and right-side species names use the same actual Euclidean/complete-linkage leaf order, as in the 0705 display approach applied to current RC3 data. Row/column clustering is for display only.")
}
