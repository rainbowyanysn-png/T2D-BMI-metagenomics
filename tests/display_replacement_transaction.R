# File transaction only. Definitions have no source-time side effects.
# Shared by real replacement and temp-fixture tests; never calculates statistics.
replacement_commit <- function(targets, backup, protected_hashes=character(), copy_file=base::file.copy, hash_file=sha256_file) {
 need<-function(ok,msg)if(!isTRUE(ok))stop(msg,call.=FALSE)
 scalar_hash<-function(path){x<-unname(as.character(hash_file(path)));need(length(x)==1L&&!is.na(x)&&nzchar(x),paste("Invalid SHA-256 result for",path));x}
 destinations<-names(targets);staged<-unname(targets)
 need(is.character(targets)&&length(targets)>0L&&length(destinations)==length(staged)&&!anyNA(destinations)&&all(nzchar(destinations)),"Replacement requires named staged paths")
 canonical<-function(x)tolower(normalizePath(x,winslash="/",mustWork=TRUE))
 need(all(file.exists(destinations))&&all(file.exists(staged)),"Replacement source/target file missing")
 need(!anyDuplicated(canonical(destinations)),"Duplicate replacement targets")
 need(!any(canonical(destinations)==canonical(staged)),"Replacement staging must differ from targets")
 need(dir.exists(backup)&&!length(list.files(backup,all.files=TRUE,no..=TRUE)),"Backup directory must be new/empty")
 before<-setNames(unname(vapply(destinations,scalar_hash,character(1))),destinations)
 staged_hashes<-unname(vapply(staged,scalar_hash,character(1)))
 backups<-file.path(backup,paste0(sprintf("%02d",seq_along(destinations)),"_",basename(destinations)))
 protected_ok<-function(){
  if(!length(protected_hashes))return(TRUE)
  need(!is.null(names(protected_hashes))&&all(file.exists(names(protected_hashes))),"Protected file missing")
  identical(unname(vapply(names(protected_hashes),scalar_hash,character(1))),unname(as.character(protected_hashes)))
 }
 need(protected_ok(),"Protected file hashes changed before replacement")
 for(j in seq_along(destinations)){
  copied<-copy_file(destinations[[j]],backups[[j]],overwrite=FALSE)
  need(isTRUE(copied),paste("Backup copy failed. Source:",destinations[[j]],"Destination:",backups[[j]]))
  # [[j]] returns a scalar; [j] would retain the filename attribute and falsely fail.
  need(identical(scalar_hash(backups[[j]]),unname(before[[j]])),paste("Backup byte-hash mismatch:",backups[[j]]))
 }
 utils::write.csv(data.frame(file=destinations,backup=backups,sha256=unname(before)),file.path(backup,"verified_backups.csv"),row.names=FALSE)
 done<-integer()
 tryCatch({
  for(j in seq_along(destinations)){
   need(identical(scalar_hash(staged[[j]]),staged_hashes[[j]]),paste("Staged file changed:",staged[[j]]))
   done<-c(done,j) # Include a partially written destination if copying raises/fails.
   need(isTRUE(copy_file(staged[[j]],destinations[[j]],overwrite=TRUE)),paste("Replacement copy failed; close file viewers. Source:",staged[[j]],"Destination:",destinations[[j]]))
   need(identical(scalar_hash(destinations[[j]]),staged_hashes[[j]]),paste("Replacement byte-hash mismatch:",destinations[[j]]))
  }
  need(protected_ok(),"Unrelated/protected file hashes changed during replacement")
 },error=function(e){
  failed<-character()
  for(j in rev(done)){
   restored<-tryCatch({
     if(file.exists(destinations[[j]])&&identical(scalar_hash(destinations[[j]]),unname(before[[j]])))TRUE
     else isTRUE(copy_file(backups[[j]],destinations[[j]],overwrite=TRUE))&&identical(scalar_hash(destinations[[j]]),unname(before[[j]]))
   },error=function(err)FALSE)
   if(!restored)failed<-c(failed,destinations[[j]])
  }
  status<-if(length(failed))paste("ROLLBACK INCOMPLETE:",paste(failed,collapse="; "))else"All touched targets restored to original byte hashes."
  writeLines(c(conditionMessage(e),status),file.path(backup,"transaction_failure.txt"))
  stop(conditionMessage(e),"\n",status,"\nBackup directory: ",backup,call.=FALSE)
 })
 data.frame(file=destinations,backup=backups,before_sha256=unname(before),after_sha256=unname(vapply(destinations,scalar_hash,character(1))),stringsAsFactors=FALSE)
}
