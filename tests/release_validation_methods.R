# Independent validation-only helpers. No fitting and no statistic replacement.
release_check <- function(label,action,rows) {
 error<-"";ok<-tryCatch({action();TRUE},error=function(e){error<<-conditionMessage(e);FALSE})
 rows$data[[length(rows$data)+1L]]<-data.frame(check=label,pass=ok,detail=error)
 message(if(ok)"PASS "else"FAIL ",label,if(nzchar(error))paste0(": ",error))
 invisible(ok)
}
release_compare_frames <- function(a,b,label,out,absolute,relative=0) {
 absolute<-release_numeric_scalar(absolute,"absolute comparison tolerance");relative<-release_numeric_scalar(relative,"relative comparison tolerance")
 assert_true(identical(names(a),names(b))&&nrow(a)==nrow(b),paste("Schema/row count differs:",label))
 results<-list();differences<-list()
 for(n in names(a)){
  x<-a[[n]];y<-b[[n]]
  if(is.numeric(x)&&is.numeric(y)){
   d<-abs(x-y);bound<-absolute+relative*abs(y)
   same<-(is.na(x)&is.na(y)) | (is.finite(x)&is.finite(y)&d<=bound) | (is.infinite(x)&is.infinite(y)&x==y)
   if(grepl("coef|estimate|log2FC",n))same<-same & ((is.na(x)&is.na(y))|sign(x)==sign(y))
   if(grepl("qval|formal_q|padj|q_bh|q_within|qval_CLR",n,ignore.case=TRUE))same<-same & ((is.na(x)&is.na(y))|((x<.05)==(y<.05)))
  }else{
   d<-bound<-rep(NA_real_,length(x));same<-(is.na(x)&is.na(y))|(!is.na(x)&!is.na(y)&enc2utf8(as.character(x))==enc2utf8(as.character(y)))
  }
  bad<-which(is.na(same)|!same)
  results[[n]]<-data.frame(field=n,pass=!length(bad),n_differences=length(bad))
  if(length(bad))differences[[n]]<-data.frame(row=bad,field=n,current=as.character(x[bad]),reference=as.character(y[bad]),absolute_difference=d[bad],allowed_difference=bound[bad])
 }
 safe_write_csv(do.call(rbind,results),file.path(out,paste0(label,"_fields.csv")))
 if(length(differences))safe_write_csv(do.call(rbind,differences),file.path(out,paste0(label,"_differences.csv")))
 assert_true(!length(differences),paste("Per-field difference:",label));invisible(TRUE)
}
release_internal_checks <- function(root,run,plan,out,rows,export=NULL) {
 plan<-release_normalize_tolerances(plan)
 state<-latest_pipeline_state(run$dir)
 registered<-c("inputs",utils::read.csv(file.path(root,"config/pipeline_modules.csv"))$module)
 jobs<-Filter(function(x)x$module_id%in%registered,plan$modules)
 release_check("required completed modules",function(){
  assert_true(!dir.exists(file.path(run$dir,".running.lock")),"Run is active")
  assert_true(setequal(names(state$completed),vapply(jobs,function(x)x$module_id,character(1))),"Unexpected/missing retained module")
 },rows)
 for(job in jobs)local({j<-job;m<-j$module_id
  release_check(paste("module",m),function(){
   verify_module_checkpoint(run,checkpoint_record(run,state$completed[[m]]))
   base<-run$result_dirs[[m]]
   for(file in j$expected_output_files){
    path<-file.path(base,relative_project_path(file$file));assert_file(path)
    d<-utils::read.csv(path,check.names=FALSE,stringsAsFactors=FALSE)
    assert_true(nrow(d)==file$n_rows,paste("Saved row count differs:",file$file))
    if(!is.null(file$named_numeric_checks))release_validate_pooled_checks(d,file$named_numeric_checks)
    if(!is.null(file$formal_q_column)){
     q<-d[[file$formal_q_column]]
     assert_true((is.numeric(q)||!length(q))&&all(is.na(q)|(q>=0&q<=1)),"Invalid formal q column")
     if(!is.null(file$expected_significant))assert_true(sum(q<.05,na.rm=TRUE)==file$expected_significant,paste("Significant count differs:",file$file))
     else assert_true(identical(file$acceptance_status,"PENDING_CLR_ALIGNMENT_REVIEW"),"Missing significance expectation without explicit review status")
     term_bh<-NULL
     if(m%in%c("interactions","bmi23_interactions")&&file$formal_q_column=="formal_q")term_bh<-p.adjust(d$pval,"BH")
     if(file$formal_q_column=="padj"&&"p_value"%in%names(d))term_bh<-p.adjust(d$p_value,"BH")
     if(file$formal_q_column=="q_bh_60")term_bh<-p.adjust(d$p_value,"BH")
     if(file$formal_q_column=="qval_CLR_BH"&&"pval_CLR"%in%names(d))term_bh<-p.adjust(d$pval_CLR,"BH")
     if(file$formal_q_column=="q_within_exposure_BH"){
      term_bh<-rep(NA_real_,nrow(d));for(term in unique(d$metadata)){keep<-d$metadata==term;term_bh[keep]<-p.adjust(d$pval[keep],"BH")}
     }
     if(!is.null(term_bh)) {
      field<-if("feature"%in%names(d))as.character(d$feature)else if("pathway"%in%names(d))as.character(d$pathway)else as.character(seq_len(nrow(d)))
      same_na<-is.na(q)&is.na(term_bh);delta<-abs(q-term_bh)
      exact_sig<-same_na|(!is.na(q)&!is.na(term_bh)&(q<.05)==(term_bh<.05))
      row_pass<-same_na|(!is.na(delta)&delta<=plan$absolute_tolerance&exact_sig)
      family<-if(file$formal_q_column=="q_within_exposure_BH")as.character(d$metadata)else rep("complete target family",nrow(d))
      size<-as.integer(table(family)[family])
      audit<-data.frame(feature=field,family=family,family_size=size,saved_q=q,recomputed_q_validation_only=term_bh,absolute_difference=delta,relative_difference=delta/pmax(abs(term_bh),.Machine$double.xmin),significance_same=exact_sig,pass=row_pass)
      path<-file.path(out,"q_family_checks",paste0(m,"__",gsub("[^A-Za-z0-9_.-]","_",basename(file$file))))
      safe_write_csv(audit,path)
      assert_true(all(row_pass)&&isTRUE(all.equal(q,term_bh,tolerance=plan$absolute_tolerance,check.attributes=FALSE)),paste("Saved BH family differs:",file$file))
     }
    }
   }
   for(contract in j$model_contract){
    c<-utils::read.csv(file.path(base,relative_project_path(contract$file)),stringsAsFactors=FALSE)
    for(n in names(contract$parameters)){
     x<-c$value[c$parameter==n];assert_true(length(x)==1L && identical(x,as.character(contract$parameters[[n]])),paste("Model contract differs:",m,n))
    }
   }
  },rows)
 })
 release_check("aligned CLR model/sample contracts",function(){pub_clr_alignment_contract(run)},rows)
 release_check("CLR revised-result reference freeze",function(){
   c<-Filter(function(x)x$module_id=="clr",plan$modules)[[1]]
   assert_true(!any(vapply(c$expected_output_files,function(f)identical(f$acceptance_status,"PENDING_CLR_ALIGNMENT_REVIEW"),logical(1))),"Three changed CLR results await author review and a new external reference; do not accept old significance counts")
   accepted<-yaml::read_yaml(file.path(root,"config/clr_review_reference.yml"))
   assert_true(identical(accepted$status,"ACCEPTED_CLR_AUTHOR_REVIEW")&&is.character(accepted$reference_index_sha256)&&nchar(accepted$reference_index_sha256)==64L&&identical(accepted$statistical_contract_sha256,sha256_file(file.path(root,"config/statistical_implementation_sha256.csv"))),"Signed CLR reference freeze is absent or tied to different statistical code")
 },rows)
 source(file.path(root,"tests/validate_species_loco.R"),local=environment())
 release_check("full-panel species LOCO acceptance",function(){validate_species_loco(root,run)},rows)
 release_check("formal raw-Aitchison PERMDISP values",function(){final_beta_sources(run)},rows)
 release_check("input sample/feature contracts",function(){
  base<-run$result_dirs[["inputs"]]
  md<-utils::read.csv(file.path(base,"tables/01_sample_manifest.csv"))
  assert_true(nrow(md)==1186L&&!anyDuplicated(md$sample_id)&&sum(md$disease=="T2D")==760&&sum(md$disease=="control")==426,"Sample composition differs")
  assert_true(identical(as.integer(table(factor(md$BMI25_group,levels=c("Non_Overweight_Control","Overweight_Control","Non_Overweight_T2D","Overweight_T2D")))),c(307L,119L,138L,622L)),"BMI25 groups differ")
  assert_true(nrow(utils::read.csv(file.path(base,"tables/01_species_prevalence.csv")))==200L&&nrow(utils::read.csv(file.path(base,"tables/01_pathway_prevalence.csv")))==5945L,"Feature panels differ")
 },rows)
 release_check("access guard evidence",function(){
  f<-file.path(run$dir,"file_access_guard.csv");assert_file(f)
  d<-utils::read.csv(f);assert_true(nrow(d)>0&&!any(d$denied),"Denied/missing file access evidence")
  source(file.path(root,"tests/validate_no_legacy_dependency.R"),local=environment())
  safe_write_csv(validate_no_legacy_dependency(root),file.path(out,"static_dependency_checks.csv"))
 },rows)
 release_check("statistical warning classification",function(){
  files<-list.files(file.path(run$dir,"logs"),pattern="^warnings[.]csv$",full.names=TRUE,recursive=TRUE)
  warnings<-lapply(files,function(f){x<-utils::read.csv(f,stringsAsFactors=FALSE);if(nrow(x)){x$source<-substring(f,nchar(run$dir)+2L);x}else NULL})
  warnings<-Filter(Negate(is.null),warnings)
  if(length(warnings)){
   x<-do.call(rbind,warnings)
   x$classification<-ifelse(grepl("Column no[.]|Row no[.].*80%|containing >80% zeros",x$warning),"PRESPECIFIED_CZM_FILTERING","UNCLASSIFIED")
   reml<-grepl("Fisher scoring|local maximum|Setting tau\\^2",x$warning,ignore.case=TRUE)
   x$classification[reml]<-"UNRESOLVED_REML_NUMERICAL_WARNING"
   safe_write_csv(x,file.path(out,"warning_classification.csv"))
   source(file.path(root,"tests/reml_warning_review.R"),local=environment())
   release_reml_saved_review(root,run,out,x,export)
   assert_true(!any(x$classification%in%c("UNCLASSIFIED","UNRESOLVED_REML_NUMERICAL_WARNING")),"Unresolved numerical warnings require feature-level diagnosis and author review; release remains blocked")
  }else safe_write_csv(data.frame(warning=character(),classification=character()),file.path(out,"warning_classification.csv"))
 },rows)
}
release_reference_checks <- function(root,run,ref,plan,out,rows) {
 ref<-normalizePath(ref,winslash="/",mustWork=TRUE)
 assert_true(!startsWith(ref,paste0(run$dir,"/"))&&!identical(ref,run$dir),"The current Run cannot serve as its own reference")
 spec<-yaml::read_yaml(file.path(root,"config/reference_results_manifest.yml"))
 assert_true(identical(sha256_file(file.path(ref,spec$reference_index)),spec$reference_index_sha256),"Reference index differs from the approved bundle")
 index<-utils::read.csv(file.path(ref,spec$reference_index),stringsAsFactors=FALSE)
 assert_columns(index,c("module","relative_file","reference_file","sha256"),"External reference index")
 assert_true(!anyDuplicated(paste(index$module,index$relative_file)),"Duplicate reference entry")
 jobs<-Filter(function(x)isTRUE(x$requires_external_reference),plan$modules)
 for(job in jobs)local({j<-job
  release_check(paste("authoritative comparison",j$module_id),function(){
   entries<-index[index$module==j$module_id,,drop=FALSE]
   assert_true(nrow(entries)>0,"Module reference missing")
   expected<-vapply(j$expected_output_files,function(x)x$file,character(1))
   assert_true(all(expected%in%entries$relative_file),"Reference omits an expected output")
   for(i in seq_len(nrow(entries))){
    current<-file.path(run$result_dirs[[j$module_id]],relative_project_path(entries$relative_file[i]))
    frozen<-file.path(ref,relative_project_path(entries$reference_file[i]));assert_file(current);assert_file(frozen)
    assert_true(tolower(sha256_file(frozen))==tolower(entries$sha256[i]),"External reference hash mismatch")
    assert_true(!identical(normalizePath(current,winslash="/"),normalizePath(frozen,winslash="/")),"Self comparison is not independent reproduction")
    sep<-if(grepl("[.]tsv$",current))"\t"else","
    a<-utils::read.table(current,sep=sep,header=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"")
    b<-utils::read.table(frozen,sep=sep,header=TRUE,check.names=FALSE,stringsAsFactors=FALSE,comment.char="",quote="\"")
    # Run paths and serialized-file hashes are auxiliary provenance, not statistics.
    # The model/source/snapshot checks elsewhere protect the referenced numerical data.
    if(basename(current)=="Clinical_correlation_input_manifest.csv"){
     auxiliary<-intersect(c("project_file","sha256"),names(a))
     safe_write_csv(data.frame(field=auxiliary,scope="auxiliary provenance; not a statistical equivalence criterion"),file.path(out,paste0(j$module_id,"_",i,"_provenance_scope.csv")))
     a<-a[,setdiff(names(a),auxiliary),drop=FALSE];b<-b[,setdiff(names(b),auxiliary),drop=FALSE]
    }
    absolute<-if(j$module_id=="pooled_pathways")plan$pooled_absolute_tolerance else plan$absolute_tolerance
    relative<-if(j$module_id=="pooled_pathways")plan$pooled_relative_tolerance else 0
    label<-paste0(j$module_id,"_",i)
    release_compare_frames(a,b,label,out,absolute,relative)
   }
  },rows)
 })
}
