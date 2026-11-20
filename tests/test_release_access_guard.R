# Regression test only: no ExperimentHub access, download or model execution.
test_release_access_guard <- function(root) {
  scope<-new.env(parent=globalenv())
  source(file.path(root,"tests/release_access_guard.R"),local=scope)
  classify<-scope$rc3_guard_filesystem_paths
  stopifnot(identical(classify(""),character()),identical(classify(NA_character_),character()),
    identical(classify(c("https://experimenthub.bioconductor.org","HTTP://example.invalid","ftp://example.invalid")),character()),
    identical(classify(c("","https://example.invalid","relative.csv")),"relative.csv"))
  blocked<-file.path(dirname(root),"T2D_BMI_RC3_FigureRefinement","output","not_read.csv")
  stopifnot(identical(classify(blocked),blocked),identical(classify("relative.csv"),"relative.csv"))
  error<-tryCatch({classify("file:///not_read.csv");NULL},error=function(e)e)
  stopifnot(inherits(error,"error"))
  # Exercise the installed guard directly, without opening a file or changing cwd.
  scope$rc3_install_guard(root)
  on.exit({
    for(n in c("file","gzfile","bzfile","xzfile","unz","readRDS","load","source","sys.source","readLines","scan"))
      suppressMessages(untrace(n,where=baseenv()))
    for(n in c("v3_guard","v3_file_access_report"))if(exists(n,envir=.GlobalEnv,inherits=FALSE))rm(list=n,envir=.GlobalEnv)
  },add=TRUE)
  guard<-get("v3_guard",.GlobalEnv)
  guard("");guard("https://experimenthub.bioconductor.org")
  err<-tryCatch({guard(blocked);NULL},error=function(e)e)
  stopifnot(inherits(err,"error"),grepl("denied legacy/private path",conditionMessage(err),fixed=TRUE))
  guard(file.path(root,"config","source_resources.csv"))
  invisible(TRUE)
}
