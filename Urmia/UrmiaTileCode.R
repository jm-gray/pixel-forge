#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Takes input Landsat surface reflectance .tar.gz files (pre-collection, preliminary!)
# and creates tiled output for a given tile shapefile, one tile per day of input data
# only L5-style bands 1, 2, 3, 4, 5, and 7 and the cfmask band are retained
# for the 25x25 km Urmia tiles, each daily Landsat GeoTIFF is about 9.7 MB
# Josh Gray, March 2018 and INTO APRIL!
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# bsub -q cnr -n 12 -W 48:00 -o /rsstu/users/j/jmgray2/SEAL/UrmiaTile.out.%J -e /rsstu/users/j/jmgray2/SEAL/UrmiaTile.err.%J R CMD BATCH --vanilla /rsstu/users/j/jmgray2/SEAL/UrmiaTileCode.R /rsstu/users/j/jmgray2/SEAL/UrmiaTileCode.Rlog
# bsub -q cnr -W 48:00 -o /rsstu/users/j/jmgray2/SEAL/UrmiaTile.out.%J -e /rsstu/users/j/jmgray2/SEAL/UrmiaTile.err.%J R CMD BATCH --vanilla /rsstu/users/j/jmgray2/SEAL/UrmiaTileCode.R /rsstu/users/j/jmgray2/SEAL/UrmiaTileCodeSerial.Rlog

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# prelims
library(raster)
library(rgdal)
library(rgeos)
library(parallel)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
MakeTempDirectory <- function(output_root_dir){
  this_tmp_dir <- file.path(output_root_dir, paste(strftime(Sys.time(), format="%a_%b_%d_%H_%M_%S"), paste(sample(letters, 4), collapse=""), sep="_"))
  while(dir.exists(this_tmp_dir)) this_tmp_dir <- file.path(output_root_dir, paste(strftime(Sys.time(), format="%a_%b_%d_%H_%M_%S"), paste(sample(letters, 4), collapse=""), sep="_"))
  system(paste("mkdir", this_tmp_dir))
  return(this_tmp_dir)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
TileLandsatDate <- function(this_date, tar_files, pr_shp, tile_shp_file, tile_shp, output_dir, tmp_output_dir, cleanup_temp=T, return_files=F, delete_tar=F, gdal_program_path=""){
  #-------------------------------------------
  # Get  the sensor, date, and pathrow of all tar files
  tar_sensors <- gsub("(L[T|C|E][0-9]{1}).*", "\\1", basename(tar_files))
  tar_pathrows <- gsub("L[T|C|E][0-9]{1}([0-9]{6}).*", "\\1", basename(tar_files))
  tar_dates <- as.Date(gsub("L[T|C|E][0-9]{1}[0-9]{6}([0-9]{7}).*", "\\1", basename(tar_files)), format="%Y%j")
  proc_dates <- as.Date(gsub("L[T|C|E][0-9]{1}[0-9]{6}[0-9]{7}-SC([0-9]{14}).*", "\\1", basename(tar_files)), format="%Y%m%d%H%M%S")

  # restrict to given date
  tmp_tar_files <- tar_files[tar_dates == this_date]
  tmp_pathrows <- tar_pathrows[tar_dates == this_date]
  tmp_sensors <- tar_sensors[tar_dates == this_date]
  tmp_procdates <- proc_dates[tar_dates == this_date]

  # check for identical tile-dates and favor the most recently processed one
  tmp_df <- data.frame(tf=tmp_tar_files, pr=tmp_pathrows, pd=tmp_procdates)
  most_recent_indices <- by(tmp_df, tmp_df[,"pr"], function(x) x$tf[x$pd == max(x$pd)])
  delete_indices <- c(1:nrow(tmp_df))[!c(1:nrow(tmp_df) %in% most_recent_indices)]
  delete_tar_files <- tmp_tar_files[delete_indices]
  tmp_tar_files <- tmp_tar_files[most_recent_indices]
  tmp_pathrows <- tmp_pathrows[most_recent_indices]
  tmp_sensor <- tmp_sensors[most_recent_indices]
  tmp_procdates <- tmp_procdates[most_recent_indices]

  # delete old tar files
  if(length(delete_tar_files) > 0){
    for(delete_file in delete_tar_files){
      system(paste("rm -f", delete_file))
    }
  }

  # read in the output tile shapefile
  # tile_shp <- shapefile(tile_shp_file)

  # get path row footprints for this date and intersecting output grid tiles
  tmp_pr_shp <- pr_shp[match(tmp_pathrows, pr_shp$PR), ]
  tmp_tile_shp <- tile_shp[apply(over(tile_shp, tmp_pr_shp), 1, function(x) !all(is.na(x))), ]
  # CheckTilesExist(overlapping_tiles, )
  # gsub("(.*)-.*", "\\1", basename(tmp_tar_files[1]))

  #-------------------------------------------
  # Make a VRT from the VRTs for the untarred, stacked, and maybe reprojected single-scene imagery
  single_scene_vrts <- try(sapply(tmp_tar_files, UntarStackProject, output_root_dir=tmp_output_dir, mtl_dir=file.path(output_dir, "MTL"), gdal_program_path=gdal_program_path))
  # try and catch an error in UntarStackProject
  if(inherits(single_scene_vrts, 'try-error')){
    return(FALSE)
  }
  multi_vrt_output_file <- file.path(tmp_output_dir, paste(gsub("(L[T|C|E][4|5|7|8])[0-9]{6}([0-9]{7}).*([0-9]{2}).VRT$", "\\1_\\2_\\3", basename(single_scene_vrts[1])), ".VRT", sep=""))
  tile_output_name_root <- gsub("(.*).VRT", "\\1", basename(multi_vrt_output_file))
  if(length(single_scene_vrts) > 1){
    gdalvrt_cmd <- paste(paste(gdal_program_path, "gdalbuildvrt", sep=""), multi_vrt_output_file, paste(single_scene_vrts, collapse=" "))
    # gdalvrt_cmd <- paste("/Library/Frameworks/GDAL.framework/Programs/gdalbuildvrt", multi_vrt_output_file, paste(single_scene_vrts, collapse=" "))
    system(gdalvrt_cmd)
  }else{
    # there was only one scene, so no need for a new VRT
    # system(paste("mv", single_scene_vrts[1], multi_vrt_output_file))
    multi_vrt_output_file <- single_scene_vrts[1]
  }

  #-------------------------------------------
  # loop through each overlapping tile and use gdalwarp to create an associated tile
  created_files <- rep(NA, dim(tmp_tile_shp)[1])
  for(i in 1:dim(tmp_tile_shp)[1]){
    vNum <- tmp_tile_shp$grid_v_num[i]
    hNum <- tmp_tile_shp$grid_h_num[i]
    tile_output_directory <- file.path(output_dir, paste("h", formatC(as.integer(hNum), width=2, flag="0"), "v", formatC(as.integer(vNum), width=2, flag="0"), sep=""))
    if(!dir.exists(tile_output_directory)) system(paste("mkdir", tile_output_directory))
    # tile_output_file <- file.path(tile_output_directory, gsub(".VRT", paste("_h", formatC(as.integer(hNum), width=2, flag="0"), "v", formatC(as.integer(vNum), width=2, flag="0"), ".tif", sep=""), basename(multi_vrt_output_file)))
    tile_output_file <- file.path(tile_output_directory, paste(tile_output_name_root, "_h", formatC(as.integer(hNum), width=2, flag="0"), "v", formatC(as.integer(vNum), width=2, flag="0"), ".tif", sep=""))
    cwhere_exp <- paste("'grid_v_num = '", as.integer(vNum), "' AND grid_h_num = '", as.integer(hNum), "''", sep="")
    gdalwarp_cmd <- paste(paste(gdal_program_path, "gdalwarp", sep=""), "-overwrite -cutline", tile_shp_file, "-cwhere", cwhere_exp, "-crop_to_cutline --config GDALWARP_IGNORE_BAD_CUTLINE YES", multi_vrt_output_file, tile_output_file)
    # gdalwarp_cmd <- paste("/Library/Frameworks/GDAL.framework/Programs/gdalwarp -overwrite -cutline", output_tiles_shp_file, "-cwhere", cwhere_exp, "-crop_to_cutline --config GDALWARP_IGNORE_BAD_CUTLINE YES", multi_vrt_output_file, tile_output_file)
    system(gdalwarp_cmd)
    created_files[i] <- tile_output_file
  }

  # check that all the created files actually exists; this is necessary b/c we're shelling out the commands and therefore can't catch errors
  for(created_file in created_files) if(!file.exists(created_file)) return(FALSE)

  #-------------------------------------------
  # clean up temporary files
  if(cleanup_temp){
    temp_dirs_to_delete <- unlist(lapply(strsplit(single_scene_vrts, split="/"), function(x) paste(x[1:(length(x) - 1)], collapse="/")))
    for(temp_dir in temp_dirs_to_delete){
      system(paste("rm -rf", temp_dir))
    }
    system(paste("rm -f", multi_vrt_output_file))
  }
  # delete the original tar files if requested
  if(delete_tar){
    for(temp_tar in tmp_tar_files){
      system(paste("rm -rf", temp_tar))
    }
  }

  #-------------------------------------------
  # return the created files if requested
  # if(return_files) return(created_files)

  #-------------------------------------------
  return(TRUE) # returns true for success
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
UntarStackProject <- function(tar_file, output_root_dir, mtl_dir, proj4="+proj=utm +zone=38 +datum=WGS84 +units=m +no_defs", gdal_program_path=""){
  #-------------------------------------------
  # create a unique temporary output directory
  this_tmp_dir <- MakeTempDirectory(output_root_dir)

  #-------------------------------------------
  # specify bands based on which sensor
  sensor <- gsub("(L[T|C|E][0-9]{1}).*", "\\1", basename(tar_file))
  if(sensor == "LC8"){
    untar_expression <- "'*sr_band[2-7].tif' '*cfmask.tif' '*MTL.txt'"
  }else{
    untar_expression <- "'*sr_band[1-5,7].tif' '*cfmask.tif' '*MTL.txt'"
  }

  #-------------------------------------------
  # untar the directory and get a list of the output tif files
  tar_cmd <- paste("tar -xvf", tar_file, "-C", this_tmp_dir, untar_expression)
  system(tar_cmd)
  tif_files <- dir(this_tmp_dir, pattern="*.tif$", full=T)
  mtl_file <- dir(this_tmp_dir, pattern="*MTL.txt$", full=T)

  # move the MTL file to the output folder
  # mtl_dir <- file.path(output_root_dir, "MTL")
  if(!dir.exists(mtl_dir)) system(paste("mkdir", mtl_dir))
  system(paste("mv", mtl_file, mtl_dir))


  #-------------------------------------------
  # check the projection of the data
  tmp_r <- raster(tif_files[1])
  correct_proj <- identical(CRS(proj4), CRS(projection(tmp_r)))

  #-------------------------------------------
  # if the native projection matches proj4, then we just create a multi-band VRT
  # if not, then we reproject each and then create a multi-band VRT
  output_vrt_file <- file.path(this_tmp_dir, paste(gsub(pattern="(.*)_.*", "\\1", basename(tif_files[1])), ".VRT", sep=""))
  cfmask_index <- c(1:length(tif_files))[sapply(tif_files, function(x) grepl("cfmask", x))]
  # cfmask_int16_file <- gsub("cfmask", "cfmaskINT16", tif_files[cfmask_index])
  if(correct_proj){
    # gdalwarp_cmd <- paste("/Library/Frameworks/GDAL.framework/Programs/gdalwarp -ot Int16", tif_files[cfmask_index], cfmask_int16_file)
    # system(gdalwarp_cmd)
    # tif_files <- c(sort(tif_files[-cfmask_index]), cfmask_int16_file)
    tif_files <- c(sort(tif_files[-cfmask_index]), tif_files[cfmask_index])

    gdalvrt_cmd <- paste(paste(gdal_program_path, "gdalbuildvrt", sep=""), "-separate", output_vrt_file, paste(tif_files, collapse=" "))
    # gdalvrt_cmd <- paste("/Library/Frameworks/GDAL.framework/Programs/gdalbuildvrt -separate", output_vrt_file, paste(tif_files, collapse=" "))
    system(gdalvrt_cmd)
    return(output_vrt_file)
  }else{
    for(tmp_tif_file in tif_files){
      tmp_output_file <- gsub(".tif$","_rp.tif", tmp_tif_file) # append "_rp" to the original filename
      # gdalwarp_cmd <- paste("/Library/Frameworks/GDAL.framework/Programs/gdalwarp -ot Int16 -t_srs '", proj4, "'", tmp_tif_file, tmp_output_file)
      gdalwarp_cmd <- paste(paste(gdal_program_path, "gdalwarp", sep=""), "-t_srs '", proj4, "'", tmp_tif_file, tmp_output_file)
      # gdalwarp_cmd <- paste("/Library/Frameworks/GDAL.framework/Programs/gdalwarp -t_srs '", proj4, "'", tmp_tif_file, tmp_output_file)
      system(gdalwarp_cmd)
    }
    tif_files <- dir(this_tmp_dir, pattern="_rp.tif", full=T)
    cfmask_index <- c(1:length(tif_files))[sapply(tif_files, function(x) grepl("cfmask", x))]
    tif_files <- c(sort(tif_files[-cfmask_index]), tif_files[cfmask_index])
    gdalvrt_cmd <- paste(paste(gdal_program_path, "gdalbuildvrt", sep=""), "-separate", output_vrt_file, paste(tif_files, collapse=" "))
    # gdalvrt_cmd <- paste("/Library/Frameworks/GDAL.framework/Programs/gdalbuildvrt -separate", output_vrt_file, paste(tif_files, collapse=" "))
    system(gdalvrt_cmd)
    return(output_vrt_file)
  }
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Start of execution
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# constants:
# for henry2:
gdal_program_path <- "/usr/local/apps/gdal/gcc483-2.1.1/bin/"
landsat_data_dir <- "/rsstu/users/j/jmgray2/SEAL/NIP/UrmiaLandsat"
output_root_dir <- "/rsstu/users/j/jmgray2/SEAL/NIP/tiled_data"
tmp_output_dir <- "/rsstu/users/j/jmgray2/SEAL/NIP/tmp_output"
int_path_row_shp <- shapefile("/rsstu/users/j/jmgray2/SEAL/NIP/GIS/UrmiaIntersectingLandsatFootprints.shp")
output_tiles_shp_file <- "/rsstu/users/j/jmgray2/SEAL/NIP/GIS/UrmiaLandsatTileGrid.shp"
success_Rdata <- "/rsstu/users/j/jmgray2/SEAL/NIP/TileSuccess.Rdata"
NUM_CORES <- 12

# for iMac
# gdal_program_path <- "/Library/Frameworks/GDAL.framework/Programs/"
# landsat_data_dir <- "/Volumes/rsstu/users/j/jmgray2/SEAL/NIP/UrmiaLandsat"
# output_root_dir <- "/Volumes/rsstu/users/j/jmgray2/SEAL/NIP/tiled_data"
# tmp_output_dir <- "/Volumes/rsstu/users/j/jmgray2/SEAL/NIP/tmp_output"
# int_path_row_shp <- shapefile("/Volumes/rsstu/users/j/jmgray2/SEAL/NIP/GIS/UrmiaIntersectingLandsatFootprints.shp")
# output_tiles_shp_file <- "/Volumes/rsstu/users/j/jmgray2/SEAL/NIP/GIS/UrmiaLandsatTileGrid.shp"
# success_Rdata <- "/Volumes/rsstu/users/j/jmgray2/SEAL/NIP/TileSuccess.Rdata"
# NUM_CORES <- 7

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# start working!
if(length(find("dir.exists")) == 0) dir.exists <- function(d) ifelse(is.na(file.info(d)$isdir), FALSE, TRUE) # this function may be necessary b/c the R version on the cluster is not up to date

output_tiles_shp <- shapefile(output_tiles_shp_file)

# get input tar.gz files and associated: sensor, path-row, and date
tar_in_files <- dir(landsat_data_dir, pattern=".*.tar.gz$", full=T, rec=T)
landsat_dates <- as.Date(gsub("L[T|C|E][0-9]{1}[0-9]{6}([0-9]{7}).*", "\\1", basename(tar_in_files)), format="%Y%j")
unique_landsat_dates <- sort(unique(landsat_dates))

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Example usage:
# TileLandsatDate(this_date=unique_landsat_dates[4], tar_files=tar_in_files, pr_shp=int_path_row_shp, tile_shp_file=output_tiles_shp_file, output_dir=output_root_dir, tmp_output_dir=tmp_output_dir, gdal_program_path=gdal_program_path, cleanup_temp=T)

cl <- makeCluster(NUM_CORES)
clusterExport(cl, c("TileLandsatDate", "UntarStackProject", "MakeTempDirectory"))
clusterEvalQ(cl, {library(raster); library(rgdal)})
system.time(tile_success <- parLapplyLB(cl, unique_landsat_dates, TileLandsatDate, tar_files=tar_in_files, pr_shp=int_path_row_shp, tile_shp=output_tiles_shp, tile_shp_file=output_tiles_shp_file, output_dir=output_root_dir, tmp_output_dir=tmp_output_dir, delete_tar=T, gdal_program_path=gdal_program_path))
save(tile_success, file=success_Rdata)

# tetsting
# cum_sum_in_files <- cumsum(sapply(unique_landsat_dates, function(x) length(tar_in_files[landsat_dates == x])))
# which(cum_sum_in_files == 1892)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# It broke, so we'll use the MTL files to try and determine where to restart
# mtl_dir <- file.path(output_root_dir, "MTL")
# mtl_files <- dir(mtl_dir, pattern=".*MTL.txt$")
# mtl_dates <- as.Date(gsub("L[T|C|E][5|7|8][0-9]{6}([0-9]{7}).*", "\\1", mtl_files), format="%Y%j")
# exclude_dates <- as.Date(c("2015197", "2013102", "2014018", "2009235", "2007053", "2007055"), format="%Y%j") # hard to know if these got done properly, so we'll repeat these dates
# mtl_dates <- mtl_dates[-match(exclude_dates, mtl_dates)]
# unique_landsat_dates <- unique_landsat_dates[-match(mtl_dates, unique_landsat_dates)]

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# copy all files since 2000 to our new share
# destination_dir <- "/Volumes/rsstu/users/j/jmgray2/SEAL/UrmiaLandsat"
# cl <- makeCluster(detectCores() - 1)
# files_to_copy <- tar_in_files[landsat_dates >= as.Date("2000-1-1")]
# # CopyFile <- function(file_to_move, destination_dir="/Volumes/rsstu/users/j/jmgray2/SEAL/UrmiaLandsat/") system(paste("cp", file_to_move, destination_dir))
# system.time(trash <- parLapply(cl, files_to_copy, function(x) system(paste("cp", x, "/Volumes/rsstu/users/j/jmgray2/SEAL/UrmiaLandsat"))))
#
# pdf(file="/Volumes/rsstu/users/j/jmgray2/SEAL/IMDONE.pdf")
# plot(1:10)
# title("I'm Done Copying!")
# dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# scratch space!
# how many of each sensor-year?
# scene_year_tbl <- table(landsat_sensor, landsat_year)
# rowSums(scene_year_tbl)
# colSums(scene_year_tbl)

# # scratch plotting of data availability
# year <- 2014
# path_row <- 169035
# # sensor <- "LC8"
# plot(landsat_date, landsat_date, type="n", ylim=c(0, 1), xlim=c(as.Date(paste(year,"-1-1", sep="")), as.Date(paste(year,"-12-31", sep=""))))
# abline(v=landsat_date[landsat_year == year & landsat_sensor == "LC8" & landsat_pathrow == path_row], lty=2, col=1)
# abline(v=landsat_date[landsat_year == year & landsat_sensor == "LE7" & landsat_pathrow == path_row], lty=2, col=2)
#
# ogr2ogr -t_srs '+proj=utm +zone=38 +datum=WGS84 +units=m +no_defs ' UrmiaBasinHydroshed_utm UrmiaBasinHydroshed.shp
# ogr2ogr -t_srs '+proj=utm +zone=38 +datum=WGS84 +units=m +no_defs ' UrmiaBasinLakesHydroshed_utm UrmiaBasinLakesHydroshed.shp
# ogr2ogr -t_srs '+proj=utm +zone=38 +datum=WGS84 +units=m +no_defs ' UrmiaRiverNetworkHydroshed_utm UrmiaRiverNetworkHydroshed.shp
# ogr2ogr -t_srs '+proj=utm +zone=38 +datum=WGS84 +units=m +no_defs ' UrmiaIntersectingLandsatFootprints_utm UrmiaIntersectingLandsatFootprints.shp
