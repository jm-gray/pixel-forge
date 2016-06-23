#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Preprocesses TM and ETM+ data: cuts to path-row overlap, reprojects, fmasks,
# and calculates EVI2.
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
library(raster)
library(rgdal)
library(tools)
library(parallel)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
LandsatOverlapPreprocess <- function(this_dir, cutline_file, out_proj=NA, suffix="_landsat_overlap", out_dir=NA){
  # gets LPDAAC processed TM and ETM+ files in this_dir, reprojects, crops to cutline, uses Fmask band
  # to mask out not-good data, then calculates EVI2

  # get all the files of interest for this particular directory
  mask_file <- dir(path=this_dir, pattern=".*cfmask.tif", full=T)
  band_files <- dir(path=this_dir, pattern=".*sr_band[1-4].tif", full=T)

  if(length(mask_file) != 1 | length(band_files) != 4){
    stop("mask or band files not proper length")
  }

  if(is.na(out_dir)) out_dir <- this_dir # set output dir
  if(is.na(out_proj)) out_proj <- projection(raster(mask_file)) # define output projection

  # take care of the mask file first
  mask_out_file <- file.path(out_dir, paste(file_path_sans_ext(basename(mask_file)), suffix, ".tif", sep=""))
  gdal_cmd <- paste("gdalwarp -overwrite -t_srs", paste("'", out_proj, "'", sep=""), "-crop_to_cutline -cutline", cutline_file, mask_file, mask_out_file)
  system(gdal_cmd)

  # loop through band files
  for(band_file in band_files){
    out_file <- file.path(out_dir, paste(file_path_sans_ext(basename(band_file)), suffix, ".tif", sep=""))
    tmp_out_file <- file.path(out_dir, "temp_fuse.tif")

    # reproject and cut to overlap line
    gdal_cmd <- paste("gdalwarp -dstnodata -9999 -overwrite -t_srs ", paste("'", out_proj, "'", sep=""), "-crop_to_cutline -cutline", cutline_file, band_file, tmp_out_file)
    system(gdal_cmd)

    # use the cfmask band to screen for good-data
    gdal_cmd <- paste("gdal_calc.py -A ", tmp_out_file, " -B ", mask_out_file, " --calc='(A * (B == 0)) + (-9999 * (B > 0))'  --NoDataValue=-9999 --overwrite --outfile=", out_file, sep="")
    system(gdal_cmd)

    # remove temp file
    system(paste("rm", tmp_out_file))
  }

  # calc evi2
  red_file <- dir(out_dir, paste(".*band3.*", suffix, ".tif", sep=""), full=T)
  nir_file <- dir(out_dir, paste(".*band4.*", suffix, ".tif", sep=""), full=T)
  evi2_out_file <- gsub("band3", "evi2", red_file)
  gdal_cmd <- paste("gdal_calc.py -A ", nir_file, " -B ", red_file, " --calc='2.5 * (((0.0001 * A) - (0.0001 * B)) / ((0.0001 * A) + (2.4 * (B * 0.0001)) + 1.0)) * 10000' --NoDataValue=-9999 --overwrite --outfile=", evi2_out_file, sep="")
  system(gdal_cmd)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# create a cutline file for the overlap region between path-row:29-31 and path-row:28-31
landsat_tile_footprints <- readOGR('/projectnb/modislc/projects/te_phenology/landsat_scenes/wrs2_descending/','wrs2_descending')
landsat_overlap <- intersect(landsat_tile_footprints[landsat_tile_footprints$WRSPR == 29031, ], landsat_tile_footprints[landsat_tile_footprints$WRSPR == 28031, ])
landsat_overlap_utm <- spTransform(landsat_overlap, CRS("+proj=utm +zone=14 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
writeOGR(landsat_overlap_utm, "/projectnb/modislc/users/joshgray/DL_Landsat", "landsat_overlap", driver="ESRI Shapefile")
cutline_file <- "/projectnb/modislc/users/joshgray/DL_Landsat/landsat_overlap.shp"

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# run script over all directories to create Landsat KF input
# get all the ETM+ and TM data directories
data_dir <- "/projectnb/modislc/users/joshgray/DL_Landsat"
dirs <- dir(path=data_dir, pattern="(LE7|LT5).*", full=T)
dirs <- dirs[file.info(dirs)$isdir]

# apply the Stack and Subset function to each directory
cl <- makeCluster(16)
clusterExport(cl, c("LandsatOverlapPreprocess"))
clusterEvalQ(cl, {library(tools); library(raster); library(rgdal)})
trash <- parLapply(cl, dirs, LandsatOverlapPreprocess, cutline_file=cutline_file)
