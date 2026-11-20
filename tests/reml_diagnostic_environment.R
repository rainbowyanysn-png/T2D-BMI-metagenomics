# Diagnostic-only environment checks; raw DESCRIPTION versions match frozen run logs.
diagnostic_assert_raw_version <- function(actual, expected, package) {
  valid<-function(x)is.character(x)&&length(x)==1L&&!is.na(x)&&nzchar(x)
  if(!valid(actual)||!valid(expected))stop("Missing or invalid version descriptor for ",package,call.=FALSE)
  if(!identical(actual,expected))stop("Frozen package version mismatch: ",package,"; expected DESCRIPTION Version=",expected,"; loaded DESCRIPTION Version=",actual,". Restore renv.lock before diagnosis; no package is installed automatically.",call.=FALSE)
  invisible(TRUE)
}
diagnostic_metafor_version_check <- function(root) {
  versions<-utils::read.csv(file.path(root,"config/package_versions.csv"),stringsAsFactors=FALSE)
  row<-versions[versions$package=="metafor",,drop=FALSE]
  if(nrow(row)!=1L)stop("Frozen package manifest requires exactly one metafor row",call.=FALSE)
  if(!requireNamespace("metafor",quietly=TRUE))stop("Frozen metafor package missing; restore renv.lock before diagnosis",call.=FALSE)
  loaded_path<-getNamespaceInfo(asNamespace("metafor"),"path")
  actual<-as.character(utils::packageDescription("metafor",lib.loc=dirname(loaded_path),fields="Version"))
  diagnostic_assert_raw_version(actual,as.character(row$version),"metafor")
  data.frame(package="metafor",frozen_DESCRIPTION_version=as.character(row$version),loaded_DESCRIPTION_version=actual,
    packageVersion_display=as.character(utils::packageVersion("metafor",lib.loc=dirname(loaded_path))),
    loaded_package_path=loaded_path,comparison="Exact raw DESCRIPTION Version",pass=TRUE,stringsAsFactors=FALSE)
}
