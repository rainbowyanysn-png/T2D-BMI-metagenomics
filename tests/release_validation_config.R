# Validation-only readers; these do not change analysis configuration or saved results.
release_numeric_scalar <- function(value, name) {
  if(length(value)!=1L||is.list(value)||is.logical(value)||!(is.numeric(value)||is.character(value)))stop("Invalid numeric validation setting: ",name,"; require one finite nonnegative number",call.=FALSE)
  converted<-suppressWarnings(as.numeric(value))
  if(length(converted)!=1L||is.na(converted)||!is.finite(converted)||converted<0)stop("Invalid numeric validation setting: ",name,"; missing, non-convertible, negative or non-finite",call.=FALSE)
  converted
}
release_normalize_tolerances <- function(plan) {
  for(n in c("absolute_tolerance","pooled_absolute_tolerance","pooled_relative_tolerance"))plan[[n]]<-release_numeric_scalar(plan[[n]],n)
  plan
}
read_release_validation_plan <- function(root) {
  release_normalize_tolerances(yaml::read_yaml(file.path(root,"config/release_validation_plan.yml")))
}
release_publication_config <- function(root) {
  config<-read_config(root);config$project_root<-root
  file<-file.path(root,"config/publication_panel_config.yml");assert_file(file,"Validation publication config")
  config$publication<-yaml::read_yaml(file)
  for(rank in c("genus","species")){
    name<-paste0("publication.pathway_candidate_limit.",rank)
    val<-release_numeric_scalar(config$publication$pathway_candidate_limit[[rank]],name)
    assert_true(val>0 && val==floor(val),paste("Invalid positive integer config:",name))
  }
  for(n in c("font_family","panel_registry"))assert_true(is.character(config$publication[[n]])&&length(config$publication[[n]])==1L&&nzchar(config$publication[[n]]),paste("Missing publication config:",n))
  config$output$font_family<-config$publication$font_family;config$output$font_size<-config$publication$font_size
  config
}
release_validate_pooled_checks <- function(d, expected) {
  assert_columns(d,c("check","observed","expected","pass"),"Pooled numeric checks")
  names_expected<-names(expected)
  assert_true(nrow(d)==length(expected)&&!anyNA(d$check)&&!anyDuplicated(d$check)&&setequal(d$check,names_expected),"Pooled numeric check names must be complete and unique")
  d<-d[match(names_expected,d$check),,drop=FALSE]
  assert_true(is.numeric(d$observed)&&is.numeric(d$expected)&&all(is.finite(d$observed))&&all(is.finite(d$expected)),"Pooled numeric check values must be finite numeric")
  target<-as.numeric(unlist(expected))
  assert_true(identical(as.numeric(d$observed),target)&&identical(as.numeric(d$expected),target),"Pooled observed/expected named values differ")
  assert_true(is.logical(d$pass)&&!anyNA(d$pass)&&all(d$pass),"All three pooled numeric checks must pass")
  invisible(TRUE)
}
