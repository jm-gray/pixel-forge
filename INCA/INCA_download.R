library(parallel)

#-------------------------------------------------------------------------------
DownloadMCD12Q2 <- function(x, out_root_dir, overwrite=FALSE){
  # if(!dir.exists(file.path(out_root_dir, x[1]))) system(paste("mkdir", file.path(out_root_dir, x[1])))
  # x form : c(year, tile)
  out_dir <- file.path(out_root_dir, "MCD12Q2", as.character(x[1]))
  existing_mcd12q2_files <- dir(out_dir, pattern=as.character(x[2]), full=T)
  if(length(existing_mcd12q2_files) == 0 | overwrite){
    wget_cmd <- paste('wget -e robots=off -m -np -nd -R -r -l1 -A "*', as.character(x[2]), '*hdf" "https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/6/MCD12Q2/', x[1], '/001', '/" --header "Authorization: Bearer 6CAF0DF2-48B6-11E8-9950-14FA569DBFBA" -P ', out_dir, sep="")
    system(wget_cmd, ignore.stdout=T, ignore.stderr=T)
    # system2(wget_cmd, stdout=stdout_file, stderr=stderr_file)
    existing_mcd12q2_files <- dir(out_dir, pattern=as.character(x[2]), full=T)
  }
  return(existing_mcd12q2_files)
}

#-------------------------------------------------------------------------------
CheckMCD12Q2 <- function(x, out_root_dir, overwrite=FALSE){
  out_dir <- file.path(out_root_dir, "MCD12Q2", as.character(x[1]))
  files <- dir(out_dir, pattern=as.character(x[2]))
  return(length(files))
}

#-------------------------------------------------------------------------------
# download all CONUS MCD12Q2C6 tiles for 2001-2016
conus_tiles <- c(
	paste(paste("h", sprintf("%02d", 8:13), sep=""), "v04", sep=""),
	paste(paste("h", sprintf("%02d", 8:12), sep=""), "v05", sep=""),
	paste(paste("h", sprintf("%02d", 8:10), sep=""), "v06", sep="")
)
years <- 2001:2016

out_root_dir <- "/rsstu/users/j/jmgray2/SEAL/INCA/MCD12Q2C6"
download_df <- data.frame(expand.grid(years, conus_tiles))
names(download_df) <- c("year", "tile")

cl <- makeCluster(8)
clusterExport(cl, c("DownloadMCD12Q2"))

trash <- parApply(cl, download_df, 1, DownloadMCD12Q2, out_root_dir=out_root_dir)

#-------------------------------------------------------------------------------
# Check for missing files and redownload
file_len <- parApply(cl, download_df, 1, CheckMCD12Q2, out_root_dir=out_root_dir)
redownload_df <- download_df[file_len == 0,]
trash <- parApply(cl, redownload_df, 1, DownloadMCD12Q2, out_root_dir=out_root_dir)