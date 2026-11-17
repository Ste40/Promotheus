args <- commandArgs(TRUE)
root <- if (length(args)) args[1] else "."
source(file.path(root, "Scripts/Microarray/metadata_retrieval.R"))
local({
  tmp <- tempfile("geo-raw-test-"); dir.create(tmp)
  on.exit(unlink(tmp, recursive=TRUE), add=TRUE)
  fixture <- file.path(tmp, "fixture"); dir.create(fixture)
  writeLines("raw-data", file.path(fixture, "GSM123.txt"))
  old <- setwd(fixture)
  utils::tar(file.path(tmp, "valid.tar"), files="GSM123.txt", tar="internal")
  setwd(old)
  archive <- file.path(tmp, "valid.tar")
  empty_listing <- function(...) NULL
  no_sleep <- function(...) NULL
  calls <- 0L
  direct <- function(url, destfile, ...) {
    stopifnot(url == "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE64nnn/GSE64918/suppl/GSE64918_RAW.tar")
    calls <<- calls + 1L
    file.copy(archive, destfile); 0L
  }
  result <- download_geo_raw("GSE64918", file.path(tmp,"empty"), supp_fn=empty_listing, download_fn=direct)
  stopifnot(geo_raw_archive_valid(result), calls==1L)
  # A valid existing archive avoids both listing and download.
  fail <- function(...) stop("Must not be called")
  stopifnot(download_geo_raw("GSE64918",file.path(tmp,"empty"),supp_fn=fail,download_fn=fail)==result)
  # A blocked listing also falls back successfully.
  result <- download_geo_raw("GSE64918",file.path(tmp,"blocked"),supp_fn=fail,download_fn=direct)
  stopifnot(geo_raw_archive_valid(result))
  # Normal GEOquery result must not invoke the fallback.
  listed <- function(gse, makeDirectory, baseDir) {
    file.copy(archive,file.path(baseDir,gse,paste0(gse,"_RAW.tar")))
  }
  stopifnot(geo_raw_archive_valid(download_geo_raw("GSE64918",file.path(tmp,"listed"),supp_fn=listed,download_fn=fail)))
  # Third-attempt success is retained, partial responses are not cached.
  calls <- 0L
  retry <- function(url,destfile,...) {
    calls <<- calls+1L
    if (calls < 3L) {writeLines("<html>unavailable</html>",destfile); return(0L)}
    file.copy(archive,destfile); 0L
  }
  result <- download_geo_raw("GSE64918",file.path(tmp,"retry"),supp_fn=empty_listing,download_fn=retry,sleep_fn=no_sleep)
  stopifnot(calls==3L,geo_raw_archive_valid(result),!length(list.files(dirname(result),pattern="\\.part$")))
  # Exhausted retries must fail explicitly and clean incomplete downloads.
  calls <- 0L
  broken <- function(url,destfile,...) {calls <<- calls+1L; writeLines("partial",destfile); stop("network error")}
  err <- tryCatch(download_geo_raw("GSE64918",file.path(tmp,"broken"),supp_fn=empty_listing,download_fn=broken,sleep_fn=no_sleep),error=identity)
  stopifnot(inherits(err,"error"),grepl("GSE64918 after 3 attempts",conditionMessage(err)),calls==3L,
            !length(list.files(file.path(tmp,"broken/GSE64918"))))
  # Returned error status cannot promote a file, even if it looks like a tar.
  bad_status <- function(url,destfile,...) {file.copy(archive,destfile); 1L}
  err <- tryCatch(download_geo_raw("GSE64918",file.path(tmp,"status"),max_attempts=1L,supp_fn=empty_listing,download_fn=bad_status),error=identity)
  stopifnot(inherits(err,"error"))
})
cat("PASS: GEO RAW download, fallback, validation, retries and cleanup\n")
