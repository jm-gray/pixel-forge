library(ncdf4)
library(raster)

ClippedRasterFromNCDF <- function(x, SPgrid_raster, example_raster){
  # tmp.df <- expand.grid(x=met_con$dim$lon$vals, y=met_con$dim$lat$vals)
  # coordinates(tmp.df)<-~x+y
  # proj4string(tmp.df) <- CRS("+init=epsg:4326")
  # gridded(tmp.df) <- TRUE
  # r <- raster(tmp.df)
  values(SPgrid_raster) <- x
  met_tmp <- projectRaster(from=SPgrid_raster, crs=CRS(proj4string(example_raster)))
  met_tmp_sub <- crop(met_tmp, extent(example_raster))
  return(met_tmp_sub)
}

out_dir <- "~/Desktop/NebraskaFusion_met"
ncdf_dir <- "/Volumes/research/fer/jmgray2/UofIMet"

cdl_raster <- raster("/Volumes/research/fer/jmgray2/NebraskaFusion/CDL/CDL_2010_landsat_overlap.tif")
tmp_met_files <- dir(ncdf_dir, full=T)
met_con <- nc_open(tmp_met_files[1])
tmp.df <- expand.grid(x=met_con$dim$lon$vals, y=met_con$dim$lat$vals)
coordinates(tmp.df)<-~x+y
proj4string(tmp.df) <- CRS("+init=epsg:4326")
gridded(tmp.df) <- TRUE
SPgrid_raster <- raster(tmp.df)

# loop through each year, dataset, and DOY and: clip, and reproject to be coincident with other data
years <- 2007:2011
for(year in years){
  print(paste("Working on year", year))
  met_files <- dir(ncdf_dir, pattern=as.character(year), full=T)
  for(met_file in met_files){
    print(paste("Working on file:", met_file))
    met_con <- nc_open(met_file)
    var_base_name <- gsub("(.*)_[0-9]{4}.nc$", "\\1", basename(met_file)) # get the variable name
    dat <- ncvar_get(met_con, names(met_con$var))
    for(i in 1:dim(dat)[3]){
      tmp_raster <- ClippedRasterFromNCDF(dat[,,i], SPgrid_raster, cdl_raster)
      out_file <- file.path(out_dir, paste(var_base_name, "_", year, formatC(i, width=3, flag="0"), ".tif", sep=""))
      print(paste("Writing doy:", i))
      writeRaster(tmp_raster, file=out_file)
    }
  }
}

# now calculate GDD
gdd_out_dir <- "~/Desktop/UofIMet_GDD"
for(year in years){
  max_doy <- as.integer(strftime(as.Date(paste(year, "-12-31", sep="")), format="%j"))
  for(doy in 1:max_doy){
    tmp_tmn <- raster(dir(out_dir, pattern=paste("tmmn_", year, formatC(doy, width=3, flag="0"), sep=""), full=T))
    tmp_tmx <- raster(dir(out_dir, pattern=paste("tmmx_", year, formatC(doy, width=3, flag="0"), sep=""), full=T))
    gdd <- ((tmp_tmn + tmp_tmx) / 2) - 273.3
    out_raster_name <- file.path(gdd_out_dir, paste("gdd_", year, formatC(doy, width=3, flag="0"), ".tif", sep=""))
    writeRaster(gdd, file=out_raster_name)
  }
}

# now make gdalwarp statements to resample the GDD data to identical resolution as other data
resample_out_dir <- "~/Desktop/UofIMet_GDD_resampled"
out_files <- dir(gdd_out_dir, full=T)
ex <- extent(cdl_raster)
for(out_file in out_files){
  tmp_out_file <- file.path(resample_out_dir, basename(out_file))
  gdal_cmd <- paste("/Library/Frameworks/GDAL.framework/Programs/gdalwarp -tr", paste(res(cdl_raster), collapse=" "), "-te", ex[1], ex[3], ex[2], ex[4], out_file, tmp_out_file)
  system(gdal_cmd)
}
