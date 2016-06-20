###############################################
# Script to download the beta C6 NBAR files and construct the GEO-style directory structure
#
# Bash submission script:
# #!/bin/bash
#
# R --vanilla < /projectnb/modislc/users/joshgray/DownloadNewNBAR.R

# Submit:
# qsub -V -pe omp 16 -l h_rt=16:00:00 subDownloadNewNBAR.sh


library(parallel)

DownloadFile <- function(x, mcd43a4_dir="/projectnb/modislc/data/mcd12_in/c6/new_nbar/mcd43a4", mcd43a2_dir="/projectnb/modislc/data/mcd12_in/c6/new_nbar/mcd43a2"){
  # get the variables from data.frame row
  tile <- as.character(unlist(x[1]))
  year <- x[2]
  doy <- formatC(as.integer(x[3]), width=3, flag="0")

  # create the new directory, only happens if it doesn't exist
  out_dir <- file.path(mcd43a4_dir, year, doy)
  created_dir <- dir.create(out_dir, recursive=T)
  wget_cmd <- paste("wget -N --user=landtest --password=STlads -P ", out_dir, " ftp://ladssci.nascom.nasa.gov/allData/6/MCD43A4/", year, "/", doy, "/*", tile, "*", sep="")
  system(wget_cmd)

  # create the new directory, only happens if it doesn't exist
  out_dir <- file.path(mcd43a2_dir, year, doy)
  created_dir <- dir.create(out_dir, recursive=T)
  wget_cmd <- paste("wget -N --user=landtest --password=STlads -P ", out_dir, " ftp://ladssci.nascom.nasa.gov/allData/6/MCD43A2/", year, "/", doy, "/*", tile, "*", sep="")
  system(wget_cmd)
}

CheckFile <- function(x, a4=T, mcd43a4_dir="/projectnb/modislc/data/mcd12_in/c6/new_nbar/mcd43a4", mcd43a2_dir="/projectnb/modislc/data/mcd12_in/c6/new_nbar/mcd43a2"){
  # get the variables from data.frame row
  tile <- as.character(unlist(x[1]))
  year <- x[2]
  doy <- formatC(as.integer(x[3]), width=3, flag="0")

  if(a4){
    files <- dir(file.path(mcd43a4_dir, year, doy), pattern=paste(".*", tile, ".*hdf", sep=""))
  }else{
    files <- dir(file.path(mcd43a2_dir, year, doy), pattern=paste(".*", tile, ".*hdf", sep=""))
  }

  ifelse(length(files) > 0, return(T), return(F))
}


# parameters
mcd43a4_dir <- "/projectnb/modislc/data/mcd12_in/c6/new_nbar/mcd43a4"
mcd43a2_dir <- "/projectnb/modislc/data/mcd12_in/c6/new_nbar/mcd43a2"
tiles <- c("h14v03", "h13v03", "h12v03", "h11v03", "h14v04", "h13v04", "h12v04", "h11v04", "h10v04", "h12v05", "h11v05", "h10v05", "h10v06")
years <- 2002:2004

# create a data.frame of all the
all_tile_years <- data.frame()
for(year in years){
  doys <- 1:365
  if(((year %% 4 == 0) & (year %% 100) > 0) | ((year %% 100 == 0) & (year %% 400 == 0))) doys <- 1:366
  tmp <- expand.grid(tiles, year, doys)
  all_tile_years <- rbind(all_tile_years, tmp)
}

# make the cluster
cl <- makeCluster(16)
clusterExport(cl, c("DownloadFile", "CheckFile"))

# check for file existence
all_tile_years$a4exists <- parApply(cl, all_tile_years, 1, CheckFile)
all_tile_years$a2exists <- parApply(cl, all_tile_years, 1, CheckFile, a4=F)

# download only missing files
download_file_years <- subset(all_tile_years, !a4exists & !a2exists)
# parApply(cl, all_tile_years, 1, DownloadFile)
parApply(cl, download_file_years, 1, DownloadFile)
