# R-level fail-fast I/O interception, not an operating-system sandbox.
rc3_guard_filesystem_paths <- function(path) {
 # file("") is an anonymous connection, not the current directory.
 # HTTP(S)/FTP connections must not be resolved against getwd().
 if(!is.character(path)||!length(path))return(character())
 path<-path[!is.na(path)&nzchar(path)]
 if(!length(path))return(character())
 # Do not let a local file URL bypass the ordinary filesystem-path checks.
 if(any(grepl("^file:",path,ignore.case=TRUE)))
  stop("Release guard requires explicit filesystem paths instead of local file URLs",call.=FALSE)
 path[!grepl("^(https?|ftps?)://",path,ignore.case=TRUE)]
}

rc3_install_guard <- function(root,resume=NULL) {
 canonical<-function(x)vapply(x,function(p){if(!grepl("^([A-Za-z]:|/|\\\\)",p))p<-file.path(getwd(),p);tolower(normalizePath(p,winslash="/",mustWork=FALSE))},character(1),USE.NAMES=FALSE)
 old<-list.dirs(file.path(root,"output"),recursive=FALSE,full.names=TRUE)
 if(!is.null(resume))old<-old[basename(old)!=resume]
 package<-dirname(root)
 denied<-c(file.path(dirname(root),c("T2D_BMI_RC3_FinalIntegration","T2D_BMI_RC3_IntegrationPreview","T2D_BMI_RC3_CleanTest","T2D_BMI_RC3_FigureRefinement","T2D_BMI_RC3_LOCO_Trial","1015返修","1015返修分析后文件整理","1015返修分析后文件整理V2")),old,file.path(package,c("internal_audit","submission_materials")),strsplit(Sys.getenv("T2DBMI_DENY_LEGACY_ROOTS"),"|",fixed=TRUE)[[1]])
 denied<-canonical(denied[nzchar(denied)])
 events<-new.env(parent=emptyenv());events$paths<-character();events$denied<-character()
 guard<-function(path){
  path<-rc3_guard_filesystem_paths(path)
  if(!length(path))return(invisible(NULL))
  x<-canonical(path)
  legacy<-grepl("/(1015返修分析后文件整理(v2|v3)?|1115返修分析后文件整理v3_medicationrestored|t2d_bmi_rc[23][^/]*|t2d_bmi_publication_(rc[23][^/]*|package_[^/]+)|[^/]*(figurerefinement|integrationpreview|loco_trial))(/|$)",x)
  blocked<-legacy|vapply(x,function(p)any(p==denied|startsWith(p,paste0(denied,"/"))),logical(1))
  events$paths<-union(events$paths,x)
  if(any(blocked)){events$denied<-union(events$denied,x[blocked]);stop("Release guard denied legacy/private path: ",path[which(blocked)[1]],call.=FALSE)}
  invisible(NULL)
 }
 assign("v3_guard",guard,.GlobalEnv)
 assign("v3_file_access_report",function()data.frame(path=events$paths,denied=events$paths%in%events$denied),.GlobalEnv)
 for(n in c("file","gzfile","bzfile","xzfile","unz","readRDS","load","source","sys.source","readLines","scan")){
  arg<-if(n%in%c("file","gzfile","bzfile","xzfile","unz"))"description"else if(n=="readLines")"con"else"file"
  trace(n,tracer=substitute(get("v3_guard",.GlobalEnv)(ARG),list(ARG=as.name(arg))),print=FALSE,where=baseenv())
 }
 Sys.unsetenv(c("T2DBMI_REFERENCE_ROOT","T2DBMI_POOLED_REFERENCE_FILE"))
 invisible(TRUE)
}
