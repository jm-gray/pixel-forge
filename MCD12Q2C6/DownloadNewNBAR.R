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


# parameters```Z
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



# on my own computer
mcd43a4_dir <- "/Users/jmgray2/Desktop/MCD43A4"
mcd43a2_dir <- "/Users/jmgray2/Desktop/MCD43A2"
tiles <- c("h12v04")
years <- 2004
doys <- 1:365
tmp_df <- rbind(expand.grid(tiles, 2004, 182:365), expand.grid(tiles, 2005, 1:183))
trash <- parApply(cl, tmp_df, 1, DownloadFile, mcd43a4_dir=mcd43a4_dir, mcd43a2_dir=mcd43a2_dir)


GetSDSName <- function(mcd43a4_file_path, band) return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":", paste("MOD_Grid_BRDF:Nadir_Reflectance_Band", band, sep=""), sep = ""))
GetSDSNameQA <- function(mcd43a2_file_path, band) return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a2_file_path, "\":", paste("MOD_Grid_BRDF:BRDF_Albedo_Band_Quality_Band", band, sep=""), sep = ""))


CalcEVI2 <- function(nir_r, red_r, scale_factor=1){
  # s <- stack(in_file)
  nir_v <- values(nir_r)
  red_v <- values(red_r)
  evi2 <- 2.5 * (((nir_v * scale_factor) - (red_v * scale_factor)) / ((nir_v * scale_factor) + (2.4 * (red_v * scale_factor)) + 1))
	tmp_r <- nir_r
	values(tmp_r) <- evi2
  return(tmp_r)
}

a4_files <- dir(mcd43a4_dir, pattern="h12v04", rec=T, full=T)
a2_files <- dir(mcd43a2_dir, pattern="h12v04", rec=T, full=T)

r <- raster(GetSDSName(a4_files[250], 1)); plot(r)

pixel_num <- (1199 * 2400) + 1
evi2 <- rep(NA, length(a4_files))
qas <- rep(NA, length(a2_files))
dates <- as.Date(gsub("MCD.*A([0-9]{7}).*hdf", "\\1", basename(a4_files)), format="%Y%j")
i <- 1
for(i in 1:length(a4_files)){
  print(paste("Doing date:", dates[i]))
  red_r <- raster(GetSDSName(a4_files[i], 1))
  nir_r <- raster(GetSDSName(a4_files[i], 2))
  qa_r <- raster(GetSDSNameQA(a2_files[i], 1))
  evi2_r <- CalcEVI2(nir_r, red_r)
  evi2[i] <- values(evi2_r)[pixel_num]
  qas[i] <- values(qa_r)[pixel_num]
  i <- i + 1
}

plot(dates, evi2, pch=qas + 1, col=qas + 1, cex=0.5)
legend("bottomleft", legend=0:3, pch=1:4, col=1:4, title="MCD43A32 QA")
