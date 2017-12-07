#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Get Landsat and MODIS data for KF without resampling MODIS
# Josh Gray
# December, 2017
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Prelims
library(raster)
library(rgdal)
library(tools)
library(parallel)

#--------------------------------------------------------------------------------
GetSDSName <- function(mcd43a4_file_path, band){
  # returns an appropriate SDS name for an MCD43A4 file path of the specified band
  return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":", paste("MOD_Grid_BRDF:Nadir_Reflectance_Band", band, sep=""), sep = ""))
}

#--------------------------------------------------------------------------------
CellNumRaster <- function(modis_file_path){
  # creates and returns a RasterLayer containing the cell numbers associated
  # with the raster given by modis_file_path
  r_modis <- raster(GetSDSName(modis_file_path, 1))
  values(r_modis) <- 1:ncell(r_modis)
  return(r_modis)
}

#--------------------------------------------------------------------------------
TileIndexRaster <- function(modis_file_path){
  # creates and returns a RasterLayer containing the MODIS tile index for
  # the raster given by modis_file_path. Value is the "h" number plus the "v"
  # number plus 10000. E.g. for "h28v09" the value is "12809"
  tile_h <- gsub(".*h([0-9]{2}).*", "\\1", basename(modis_file_path))
  tile_v <- gsub(".*v([0-9]{2}).*", "\\1", basename(modis_file_path))
  tile_num <- as.integer(paste("1", tile_h, tile_v, sep=""))
  r_modis <- raster(GetSDSName(modis_file_path, 1))
  values(r_modis) <- tile_num
  return(r_modis)
}

#--------------------------------------------------------------------------------
MakeCellNumberMap <- function(example_modis_files, example_landsat_file, out_dir){
  # creates a cell number map with the geometry of example_landsat_file over a vector of modis files
  # thus, it handles the case where a landsat geometry overlaps multiple MODIS tiles
  cell_num_rasters <- lapply(as.character(example_modis_files), CellNumRaster)
  if(length(cell_num_rasters) > 1){
    merged_cell_num_rasters <- do.call(merge, cell_num_rasters)
  }else{
    merged_cell_num_rasters <- cell_num_rasters[[1]]
  }
  r_landsat <- raster(example_landsat_file)
  tmp_out_file <- file.path(out_dir, paste("modis_temp_", strftime(Sys.time(), format="%a_%b_%d_%H_%M_%S"), paste(sample(letters, 4), collapse=""), ".tif", sep=""))
  writeRaster(merged_cell_num_rasters, file=tmp_out_file)
  landsat_path_row <- gsub(".*_([0-9]{6})_.*", "\\1", basename(example_landsat_file))
  out_file <- file.path(out_dir, paste("modis_cell_num_", landsat_path_row, ".tif", sep=""))
  gdal_cmd <- paste("/Library/Frameworks/GDAL.framework/Programs/gdalwarp -overwrite -dstnodata -9999 -s_srs '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m' -t_srs", paste("'", projection(r_landsat), "'", sep=""), "-te", paste(extent(r_landsat)[c(1, 3, 2, 4)], collapse=" "), "-ts", ncol(r_landsat), nrow(r_landsat), tmp_out_file, out_file)
  system(gdal_cmd)
  rm_cmd <- paste("rm -f", tmp_out_file)
  system(rm_cmd)
  return(out_file)
}

#--------------------------------------------------------------------------------
MakeTileIndexMap <- function(example_modis_files, example_landsat_file, out_dir){
  # creates a tile index map with the geometry of example_landsat_file over a vector of modis files
  # thus, it handles the case where a landsat geometry overlaps multiple MODIS tiles
  tile_index_rasters <- lapply(as.character(example_modis_files), TileIndexRaster)
  if(length(tile_index_rasters) > 1){
    merged_tile_index_rasters <- do.call(merge, tile_index_rasters)
  }else{
    merged_tile_index_rasters <- tile_index_rasters[[1]]
  }
  r_landsat <- raster(example_landsat_file)
  tmp_out_file <- file.path(out_dir, paste("modis_temp_", strftime(Sys.time(), format="%a_%b_%d_%H_%M_%S"), paste(sample(letters, 4), collapse=""), ".tif", sep=""))
  writeRaster(merged_tile_index_rasters, file=tmp_out_file)
  landsat_path_row <- gsub(".*_([0-9]{6})_.*", "\\1", basename(example_landsat_file))
  out_file <- file.path(out_dir, paste("modis_tile_index_", landsat_path_row, ".tif", sep=""))
  gdal_cmd <- paste("/Library/Frameworks/GDAL.framework/Programs/gdalwarp -overwrite -dstnodata -9999 -s_srs '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m' -t_srs", paste("'", projection(r_landsat), "'", sep=""), "-te", paste(extent(r_landsat)[c(1, 3, 2, 4)], collapse=" "), "-ts", ncol(r_landsat), nrow(r_landsat), tmp_out_file, out_file)
  system(gdal_cmd)
  rm_cmd <- paste("rm -f", tmp_out_file)
  system(rm_cmd)
  return(out_file)
}

#--------------------------------------------------------------------------------
CreateModisMaps <- function(path_row, landsat_path_row_shp, modis_tile_shp, landsat_data_dir="/Volumes/research/fer/jmgray2/EastKalimantan/processedData/envi", modis_data_dir="/Users/jmgray2/Desktop/MCD43A4", modis_map_out_dir="/Users/jmgray2/Desktop/MODIS_cell_maps", landsat_suffix="_EK"){
  # Creates both MODIS tile index and cell number maps for a given Landsat path row
  landsat_in_files <- dir(file.path(landsat_data_dir, path_row), pattern=paste(".*", landsat_suffix, "$", sep=""), full=T)
  this_path_row_shp <- landsat_path_row_shp[which(landsat_path_row_shp$PR == path_row),]
  this_path_row_shp_reproj <- spTransform(this_path_row_shp, CRS(projection(modis_tile_shp)))
  modis_grid_intersect <- over(this_path_row_shp_reproj, modis_tile_shp, returnList=T)
  intersecting_modis_tiles <- apply(modis_grid_intersect[[1]][,c("h", "v")], 1, function(x) paste("h", formatC(as.integer(x[1]), width=2, flag="0"), "v", formatC(as.integer(x[2]), width=2, flag="0"), sep=""))
  modis_in_files <- lapply(intersecting_modis_tiles, function(tile) list.files(modis_data_dir, pattern=tile, rec=T, full=T))
  example_modis_files <- lapply(modis_in_files, function(x) x[1])
  modis_cell_num_map_file <- MakeCellNumberMap(example_modis_files, landsat_in_files[1], modis_map_out_dir)
  modis_tile_index_map_file <- MakeTileIndexMap(example_modis_files, landsat_in_files[1], modis_map_out_dir)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Make MODIS cell number and tile index maps at Landsat resolution for all path-rows
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Paths and constants
modis_cell_maps_output_dir <- "/Users/jmgray2/Desktop/MODIS_cell_maps"
landsat_data_dir <- "/Volumes/research/fer/jmgray2/EastKalimantan/processedData/envi"
landsat_suffix <- "_EK"
mcd43a4_data_dir <- "/Users/jmgray2/Desktop/MCD43A4"
path_row_shp_file <- "/Users/jmgray2/Desktop/wrs2_descending/wrs2_descending.shp"
modis_grid_shp_file <- "/Users/jmgray2/Desktop/modis_grid/modis_sinusoidal_grid_world.shp"

# read in the shapefiles for Landsat and MODIS footprints
path_rows_shp <- shapefile(path_row_shp_file)
modis_grid_shp <- shapefile(modis_grid_shp_file)

# get all Landsat scenes
all_path_rows <- dir(landsat_data_dir)
for(path_row in all_path_rows){
  print(paste("Making maps for:", path_row))
  CreateModisMaps(path_row=path_row, landsat_path_row_shp=path_rows_shp, modis_tile_shp=modis_grid_shp, landsat_data_dir=landsat_data_dir, modis_data_dir=mcd43a4_data_dir, modis_map_out_dir=modis_cell_maps_output_dir, landsat_suffix=landsat_suffix)
}
