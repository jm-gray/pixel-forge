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

cdl_raster <- raster("/Volumes/research/fer/jmgray2/NebraskaFusion/CDL/CDL_2010_landsat_overlap.tif")
tmp.df <- expand.grid(x=met_con$dim$lon$vals, y=met_con$dim$lat$vals)
coordinates(tmp.df)<-~x+y
proj4string(tmp.df) <- CRS("+init=epsg:4326")
gridded(tmp.df) <- TRUE
SPgrid_raster <- raster(tmp.df)

out_dir <- "~/Desktop/NebraskaFusion_met"
ncdf_dir <- "/Volumes/research/fer/jmgray2/UofIMet"

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
