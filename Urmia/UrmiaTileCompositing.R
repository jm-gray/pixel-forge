library(raster)
library(parallel)

GetTileComposite <- function(tile_dir, output_dir="/Volumes/users/j/jmgray2/SEAL/NIP/maxNDVI_comp"){
  # this_dir <- tile_dirs[17]
  these_tif_files <- dir(tile_dir, pattern=".*.tif$", full=T)
  these_dates <- as.Date(gsub(".*_([0-9]{7})_.*", "\\1", basename(these_tif_files)), format="%Y%j")
  # these_tif_files <- these_tif_files[these_dates >= start_date & these_dates <= end_date]
  # these_dates <- these_dates[these_dates >= start_date & these_dates <= end_date]
  GetMaxNDVIComposite(these_tif_files, these_dates, output_dir)
}

GetMaxNDVIComposite <- function(tif_files, dates, output_dir, max_years=c(2000, 2012)){
  years <- as.integer(strftime(dates, format="%Y"))
  for(max_year in max_years){
    tmp_tif_files <- tif_files[years == max_year]
    ndvi_stack <- do.call(stack, lapply(tmp_tif_files, CalcScreenedNDVI))
    which_max_r <- which.max(ndvi_stack)
    which_max <- values(which_max_r)
    # ndvi_v <- values(ndvi_stack)
    # which_max <- unlist(apply(ndvi_v, 1, function(x) which(x == max(x, na.rm=T))))
    output_v <- array(NA, dim=c(length(which_max), 6))
    for(i in sort(unique(which_max))){
      v <- values(stack(tmp_tif_files[i]))
      output_v[which_max == i, ] <- v[which_max == i, 1:6]
    }
    output_s <- stack(tmp_tif_files[1], bands=1:6)
    values(output_s) <- output_v
    # write the output file
    output_file <- file.path(output_dir, paste("Comp_maxNDVI_", max_year, "_", gsub(".*(h[0-9]{2}v[0-9]{2}).*", "\\1", basename(tif_files[1])), ".tif", sep=""))
    writeRaster(output_s, file=output_file, overwrite=T)
  }
}

CalcScreenedNDVI <- function(tif_file){
  s <- stack(tif_file, bands=c(3, 4, 7))
  v <- values(s)
  v_ndvi <- (v[,2] - v[,1]) / (v[,2] + v[,1])
  v_ndvi[v[,3] != 0] <- NA
  v_ndvi[v_ndvi < 0] <- NA
  s_ndvi <- raster(s, 1)
  values(s_ndvi) <- v_ndvi
  return(s_ndvi)
}

# input_tile_root_dir <- "/Volumes/rsstu/users/j/jmgray2/SEAL/NIP/tiled_data"
# tile_dirs <- dir(input_tile_root_dir, pattern="h[0-9]{2}v[0-9]{2}", full=T)
# max_ndvi_output_dir <- "/Volumes/rsstu/users/j/jmgray2/SEAL/NIP/maxNDVI_comp"
# start_date <- as.Date("2000-1-1")
# end_date <- as.Date("2012-12-31")
# max_ndvi_years <- c(2000, 2012)
trash <- dir("/Volumes/users/j/jmgray2/SEAL/NIP/maxNDVI_comp", pattern=".tif")
trash_df <- data.frame(year=gsub(".*_([0-9]{4})_.*", "\\1", trash), tile=gsub(".*_(h[0-9]{2}v[0-9]{2}).tif*", "\\1", trash))
trash_tbl <- table(trash_df)
tiles_done <- colnames(trash_tbl)[apply(trash_tbl, 2, function(x) all(x==1))]

input_tile_root_dir <- "/Volumes/users/j/jmgray2/SEAL/NIP/tiled_data"
tile_dirs <- dir(input_tile_root_dir, pattern="h[0-9]{2}v[0-9]{2}", full=T)
tile_dirs <- tile_dirs[-match(tiles_done, basename(tile_dirs))]

cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("CalcScreenedNDVI", "GetMaxNDVIComposite"))
clusterEvalQ(cl, {library(raster); library(rgdal)})
system.time(trash <- parLapplyLB(cl, tile_dirs, GetTileComposite))
                                                        
