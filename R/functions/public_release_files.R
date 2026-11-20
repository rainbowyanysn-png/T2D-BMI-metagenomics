# Shared policy for public inventory generation and package copying; no analysis calls.
public_release_selection <- function(root) {
  policy<-yaml::read_yaml(file.path(root,"config/public_release_files.yml"))
  files<-policy$include_root_files
  for(tree in policy$include_trees)files<-c(files,file.path(tree,list.files(file.path(root,tree),recursive=TRUE,all.files=TRUE,no..=TRUE,full.names=FALSE)))
  files<-unique(gsub("\\\\","/",unlist(files)))
  # Explicit paths first, then directory/name rules shared by inventory and copy.
  files<-setdiff(files,unlist(policy$exclude_paths))
  for(pattern in policy$exclude_regex)files<-files[!grepl(pattern,files,perl=TRUE,ignore.case=TRUE)]
  stopifnot(!any(grepl("(^/|^[A-Za-z]:|(^|/)[.][.](/|$))",files)),all(file.exists(file.path(root,files))))
  files<-sort(files[!file.info(file.path(root,files))$isdir],method="radix")
  stopifnot(all(unlist(policy$required_entries)%in%files))
  # All R/config/tests support is retained unless explicitly identified as development-only.
  data.frame(file=files,bytes=as.numeric(file.info(file.path(root,files))$size),
    sha256=unname(vapply(file.path(root,files),function(f)digest::digest(file=f,algo="sha256"),character(1))),stringsAsFactors=FALSE)
}
public_release_manifests <- function(root) {
  x<-public_release_selection(root)
  utils::write.csv(x,file.path(root,"public_file_manifest.csv"),row.names=FALSE,fileEncoding="UTF-8")
  utils::write.csv(x[c("file","sha256")],file.path(root,"public_file_sha256.csv"),row.names=FALSE,fileEncoding="UTF-8")
  x
}
public_release_copy <- function(root,destination) {
  admission<-yaml::read_yaml(file.path(root,"config/module_admission.yml"))
  needed<-c("production_ready","clean_run_pass","independent_reference_comparison_pass","manual_review_complete","citation_and_license_complete")
  if(!all(vapply(needed,function(n)isTRUE(admission[[n]]),logical(1))))stop("Publication release gate is closed; complete author acceptance before extracting a final package.",call.=FALSE)
  if(file.exists(destination)||dir.exists(destination))stop("Destination must not already exist",call.=FALSE)
  x<-public_release_manifests(root)
  stopifnot(dir.create(destination,recursive=TRUE))
  for(rel in c(x$file,"public_file_manifest.csv","public_file_sha256.csv")) {
    dst<-file.path(destination,rel);dir.create(dirname(dst),recursive=TRUE,showWarnings=FALSE)
    stopifnot(file.copy(file.path(root,rel),dst,overwrite=FALSE),identical(digest::digest(file=file.path(root,rel),algo="sha256"),digest::digest(file=dst,algo="sha256")))
  }
  invisible(x)
}
