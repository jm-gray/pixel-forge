#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
library(parallel)
library(raster)
library(rgdal)
# bsub -q cnr -n 12 -W 12:00 -o /rsstu/users/j/jmgray2/SEAL/UrmiaCheckTiles.out.%J -e /rsstu/users/j/jmgray2/SEAL/UrmiaCheckTiles.err.%J R CMD BATCH --vanilla /rsstu/users/j/jmgray2/SEAL/UrmiaCheckTiles.R /rsstu/users/j/jmgray2/SEAL/UrmiaCheckTiles.Rlog

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
GetIntersectingTiles <- function(pr, pr_shp, tile_shp){
  tmp_pr_shp <- pr_shp[pr_shp$PR == as.integer(pr), ]
  tmp_tile_shp <- tile_shp[apply(over(tile_shp, tmp_pr_shp), 1, function(x) !all(is.na(x))), ]
  return(paste("h", formatC(as.integer(tmp_tile_shp$grid_h_num), width=2, flag="0"), "v", formatC(as.integer(tmp_tile_shp$grid_v_num), width=2, flag="0"), sep=""))
}

CheckFile <- function(tar_file, tiles_by_pr, tiles_dir="/rsstu/users/j/jmgray2/SEAL/NIP/tiled_data"){
  # checks that all the expected tile-dates for a given tar file have been created and can be opened
  tar_date <- as.Date(gsub("L[T|C|E][0-9]{1}[0-9]{6}([0-9]{7}).*", "\\1", basename(tar_file)), format="%Y%j")
  tar_pathrow <- gsub("L[T|C|E][0-9]{1}([0-9]{6}).*", "\\1", basename(tar_file))
  int_tiles <- unlist(tiles_by_pr[tar_pathrow])
  check_tiles <- sapply(int_tiles, CheckTile, date=tar_date, tiles_dir=tiles_dir)
  tmp_df <- data.frame(tar_file=tar_file, tile=int_tiles, pr=tar_pathrow, date=strftime(tar_date, format="%Y%j"), goodtile=check_tiles)
  return(tmp_df)
}

CheckTile <- function(tile, date, tiles_dir){
  # checks for the existence of a given tile-date, and also that it can be opened
  tile_file <- dir(file.path(tiles_dir, tile), pattern=strftime(date, format="%Y%j"), full=T)
  if(length(tile_file) == 0) return(FALSE)
  tmp_s <- try(stack(tile_file))
  if(inherits(tmp_s, 'try-error')) return(FALSE)
  return(TRUE)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Check that all of the expected files were created
load("/rsstu/users/j/jmgray2/SEAL/NIP/all_tar_files.Rdata")
tar_pathrows <- gsub("L[T|C|E][0-9]{1}([0-9]{6}).*", "\\1", basename(tar_files))
pr_shp <- shapefile("/rsstu/users/j/jmgray2/SEAL/NIP/GIS/UrmiaIntersectingLandsatFootprints.shp")
tile_shp <- shapefile("/rsstu/users/j/jmgray2/SEAL/NIP/GIS/UrmiaLandsatTileGrid.shp")
save_file <- "/rsstu/users/j/jmgray2/SEAL/NIP/UrmiaCheckTile.Rdata"
NUM_CORES <- 12

# get the tiles for each pathrow
tiles_by_pr <- lapply(unique(tar_pathrows), GetIntersectingTiles, pr_shp=pr_shp, tile_shp=tile_shp)
names(tiles_by_pr) <- unique(tar_pathrows)

# setup the cluster
cl <- makeCluster(NUM_CORES)
clusterExport(cl, c("CheckFile", "CheckTile"))
clusterEvalQ(cl, {library(raster); library(rgdal)})

check_tiles <- parLapplyLB(cl, tar_files, CheckFile, tiles_by_pr=tiles_by_pr)
tiles_df <- do.call(rbind, check_tiles)
save(tiles_df, file=save_file)

bad_tiles <- subset(tiles_df, !goodtile)
length(unique(bad_tiles$tar_file))
reproc_tar_files <- unique(bad_tiles$tar_file)


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# NEED TO COPY THE BAD TILES INTO SEAL HPC SHARE FOR REPROCESSING
# R CMD BATCH --vanilla ~/Desktop/test.R &
destination_dir <- "/Volumes/rsstu/users/j/jmgray2/SEAL/UrmiaLandsat"
cl <- makeCluster(detectCores() - 1)
# files_to_copy <- tar_in_files[landsat_dates >= as.Date("2000-1-1")]
files_to_copy <- reproc_tar_files
# CopyFile <- function(file_to_move, destination_dir="/Volumes/rsstu/users/j/jmgray2/SEAL/UrmiaLandsat/") system(paste("cp", file_to_move, destination_dir))
system.time(trash <- parLapply(cl, files_to_copy, function(x) system(paste("cp", x, "/Volumes/rsstu/users/j/jmgray2/SEAL/UrmiaLandsat"))))

pdf(file="/Volumes/rsstu/users/j/jmgray2/SEAL/IMDONE.pdf")
plot(1:10)
title("I'm Done Copying!")
dev.off()
