format_probability <- function(x,label="P") {
 if(length(x)!=1L||is.na(x))return(paste0(label," = NA"))
 r1_need(is.finite(x)&&x>=0&&x<=1,"Invalid saved P value")
 if(x<.001)paste0(label," < 0.001")else paste0(label," = ",sprintf("%.3f",x))
}
r1_group_vectors <- function(style,total=FALSE) {
 spec<-if(total)style$total_groups else style$groups
 fill<-unlist(spec$fill,use.names=TRUE);order<-unlist(spec$order,use.names=FALSE)
 r1_need(identical(names(fill),order),"Colour names/order differ from factor levels")
 dark<-setNames(colorspace::darken(unname(fill),amount=style$colour_adjustment$darken_amount),names(fill))
 r1_need(identical(names(dark),names(fill)),"Darkened palette lost names")
 list(fill=fill,dark=dark,order=order)
}
r1_factor <- function(x,levels){r1_need(!anyNA(x)&&all(as.character(x)%in%levels),"Unmapped category");factor(as.character(x),levels=levels)}
r1_alpha_annotation <- function(root,paths,metric,family) {
 x<-r1_read(root,paths,"BMI25_alpha_LMM_summary.csv");z<-x[tolower(x$Metric)==tolower(metric),,drop=FALSE]
 r1_need(nrow(z)==1L,"Metric is not uniquely identified in saved LMM summary")
 r1_need(all(is.finite(unlist(z[c("Overall_T2D_vs_Control_p","Four_Group_p","ICC")]))),"Missing LMM summary")
 # Values below are assertions of the stated authority, never plot-label constants.
 expected<-if(metric=="Shannon")c(.5174,.0119,.111)else c(.2167,.0103,.123)
 r1_need(all(abs(unlist(z[c("Overall_T2D_vs_Control_p","Four_Group_p","ICC")])-expected)<1e-12),"Saved LMM checkpoint differs; do not guess")
 icc<-sprintf("ICC = %.3f",z$ICC)
 text<-if(family=="alpha_total")paste0("LMM ",format_probability(z$Overall_T2D_vs_Control_p),"; ",icc)else
  if(family=="alpha_four")paste0("Four-group LMM ",format_probability(z$Four_Group_p),"; ",icc)else
  paste0("Pooled LMM: T2D ",format_probability(z$Overall_T2D_vs_Control_p),"; four-group ",format_probability(z$Four_Group_p),"; ",icc)
 list(text=text,statistics=data.frame(Metric=z$Metric,T2D_P=z$Overall_T2D_vs_Control_p,Four_Group_P=z$Four_Group_p,ICC=z$ICC))
}
r1_beta_annotation <- function(root,paths,total) {
 group<-if(total)"disease"else"analysis_group";suffix<-if(total)"total"else"four_groups"
 a<-r1_read(root,paths,paste0("BMI25_raw_Aitchison_PERMANOVA_",suffix,"_blocked.csv"));a<-a[a$term==group,,drop=FALSE]
 b<-r1_read(root,paths,paste0("BMI25_raw_Aitchison_PERMDISP_",if(total)"disease"else"four_groups",".csv"));b<-b[b$term=="Groups",,drop=FALSE]
 r1_need(nrow(a)==1L&&nrow(b)==1L&&a$target_term==group&&b$group_variable==group,"Wrong PERMANOVA/PERMDISP target")
 r1_need(a$n_samples==1182L&&b$n_samples==1182L&&a$n_permutations==999L&&b$n_permutations==999L,"Saved beta sample/permutation contract differs")
 expected<-if(total).00410923019525066 else .00922200738646509
 r1_need(abs(a$R2-expected)<1e-12&&a[["Pr(>F)"]]==.001&&a$test_method=="Marginal PERMANOVA"&&a$distance=="Uncorrected CZM/CLR Aitchison","Wrong marginal target/distance checkpoint; do not borrow overall or ComBat R2")
 text<-paste0("PERMANOVA\n",if(total)"Disease"else"BMI–T2D group"," marginal R² = ",sprintf("%.3f",a$R2),"\n",format_probability(a[["Pr(>F)"]]))
 list(text=text,statistics=data.frame(target_term=group,R2=a$R2,PERMANOVA_F=a$F,PERMANOVA_P=a[["Pr(>F)"]],PERMDISP_F=b$F,PERMDISP_P=b[["Pr(>F)"]],n_samples=a$n_samples,n_permutations=a$n_permutations))
}

r1_dispersion_caption <- function(root,paths) {
 read<-function(suffix){x<-r1_read(root,paths,paste0("BMI25_raw_Aitchison_PERMDISP_",suffix,".csv"));x[x$term=="Groups",,drop=FALSE]}
 a<-read("disease");b<-read("four_groups")
 r1_need(nrow(a)==1L&&nrow(b)==1L&&a$n_samples==1182L&&b$n_samples==1182L&&a$type=="Spatial median"&&b$type=="Spatial median"&&a$n_permutations==999L&&b$n_permutations==999L,"Frozen spatial-median PERMDISP contract differs")
 r1_need(abs(a$F-22.4781054358386)<1e-12&&abs(b$F-6.65286199089215)<1e-12&&a[["Pr(>F)"]]==.001&&b[["Pr(>F)"]]==.001,"PERMDISP checkpoint differs")
 paste0("Multivariate dispersion differed between disease groups and among the four BMI–T2D groups. The ordination patterns therefore cannot be attributed exclusively to differences in group centroids. Disease PERMDISP: F = ",sprintf("%.3f",a$F),", ",format_probability(a[["Pr(>F)"]]),". Four-group PERMDISP: F = ",sprintf("%.3f",b$F),", ",format_probability(b[["Pr(>F)"]]),". Spatial-median PERMDISP, N = 1182, 999 permutations; results are retained in Supplementary Table S3.")
}

