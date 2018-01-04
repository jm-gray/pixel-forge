#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Get Landsat and MODIS data for KF without resampling MODIS
# Josh Gray
# December, 2017
# First need to use the function "CreateModisMaps" in KF_functions.R to create
# the cell number and tile index rasters for each Landsat path-row
#
# Landsat scenes are assumed to be preprocessed with EK_preprocess.R and thus
# they are: 1) surface reflectance, 2) bands 1:5,7,6 for TM/ETM+ or 2:7,10 for OLI,
# and 3) have as the 8th layer the surf_ref_pixel_qa transformed to FMASK categories
#
# Band-to-band mapping for Landsat-MODIS
# TM/OLI    MODIS   Name
# 1/2       3       Blue
# 2/3       4       Green
# 3/4       1       Red
# 4/5       2       NIR
# 5/6       6       SWIR
# 7         7       SWIRL
# 6/10      -       Thermal/Brightness-Temp
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# To run:
# bsub -q cnr -W 36:00 -n 24 -R "span[ptile=24]" -o /gpfs_common/share02/jmgray2/EK/EK.out.%J -e -o /gpfs_common/share02/jmgray2/EK/EK.err.%J R CMD BATCH --vanilla /gpfs_common/share02/jmgray2/EK/KF_cluster.R
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Prelims
library(raster)
library(rgdal)
library(tools)
library(parallel)
library(dlm)
source("/Users/jmgray2/Documents/pixel-forge/KalmanFiltering/KF_functions.R")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Paths, parameters, and constants
modis_cell_maps_output_dir <- "/gpfs_common/share02/jmgray2/EK/MODIS_cell_maps"
landsat_data_dir <- "/gpfs_common/share02/jmgray2/EK/"
lwmask_dir <- "/gpfs_common/share02/jmgray2/EK/ClipMasks"
landsat_suffix <- "_EK"
mcd43a4_data_dir <- "/gpfs_common/share02/jmgray2/EK/MCD43A4"
mcd43a2_data_dir <- "/gpfs_common/share02/jmgray2/EK/MCD43A2"
path_row_shp_file <- "/gpfs_common/share02/jmgray2/EK/wrs2_descending.shp"
modis_grid_shp_file <- "/gpfs_common/share02/jmgray2/EK/modis_sinusoidal_grid_world.shp"
output_dir <- "/gpfs_common/share02/jmgray2/EK/KF_output"
start_date <- as.Date("2010-1-1")
end_date <- as.Date("2012-12-31")
path_row_to_process <- "116059"
parallel_cores <- 24
landsat_to_modis_bands <- c(3, 4, 1, 2, 6, NA, 7) # maps post-EK_preprocess landsat band ordering to MCD43A4 band numbers; there is not thermal band in modis, so landsat band 6 returns NA; also QA is contained in MCD43A2 so is not included (handled separately in data acquisition)
MAX_OPEN_FILES <- 2e3
rows_to_do <- 100

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Create a cluster
if(is.null(parallel_cores)){
  cl <- makeCluster(detectCores() - 1)
}else{
  cl <- makeCluster(parallel_cores)
}
clusterEvalQ(cl, {library(raster); library(rgdal); library(dlm); source("/Users/jmgray2/Documents/pixel-forge/KalmanFiltering/KF_functions.R")})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Gather the Landsat files and extract the data
landsat_in_files <- dir(file.path(landsat_data_dir, path_row_to_process), pattern=paste(".*", landsat_suffix, "$", sep=""), full=T)
landsat_dates <- as.Date(gsub(".*([0-9]{8}).*", "\\1", basename(landsat_in_files)), format="%Y%m%d")
landsat_in_files <- landsat_in_files[landsat_dates <= end_date & landsat_dates >= start_date]
landsat_dates <- landsat_dates[landsat_dates <= end_date & landsat_dates >= start_date]
landsat_sensor <- gsub(pattern="(L[E|T|C]0[7|5|8]).*", "\\1", basename(landsat_in_files))

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Get the MODIS cell and tile maps, and the LW mask file
cell_num_map <- dir(modis_cell_maps_output_dir, pattern=paste(".*cell_num.*", path_row_to_process, sep=""), full=T)
tile_index_map <- dir(modis_cell_maps_output_dir, pattern=paste(".*tile_index.*", path_row_to_process, sep=""), full=T)
lwmask_file <- dir(lwmask_dir, pattern=paste(path_row_to_process, ".tif$", sep=""), full=T)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# the main loop
tmp_r <- raster(landsat_in_files[1])
is <- 1:(ceiling(nrow(tmp_r) / rows_to_do))
for(i in is){
  # figure out the starting row and number of rows to read
  print(paste("Block", i, "of", length(is)))
  start_row <- ((i - 1) * rows_to_do) + 1
  nrows <- min(rows_to_do, (nrow(tmp_r) - ((i - 1) * rows_to_do))) # last block may have less rows

  # get the landsat data
  landsat_data <- GetValuesGDAL_multiband(landsat_in_files, start_row, nrows, max_open_datasets=MAX_OPEN_FILES)

  # get the LW mask values
  lwmask_values <- GetValuesGDAL(lwmask_file, start_row, nrows, max_open_datasets=MAX_OPEN_FILES) # only 0 values should be processed

  # get the MODIS cell numbers and tile indices
  modis_cell_nums <- GetValuesGDAL(cell_num_map, start_row, nrows, max_open_datasets=MAX_OPEN_FILES)
  modis_tile_indices <- GetValuesGDAL(tile_index_map, start_row, nrows, max_open_datasets=MAX_OPEN_FILES)

  # get a vector of all the tile indices, then get the MCD43A4/A2 files, and subset to date range
  modis_tiles <- sapply(unique(modis_tile_indices), function(x) paste("h", substr(x, 2, 3), "v", substr(x, 4, 5), sep=""))
  mcd43a4_in_files <- lapply(modis_tiles, function(tile) list.files(mcd43a4_data_dir, pattern=tile, rec=T, full=T))
  mcd43a2_in_files <- lapply(modis_tiles, function(tile) list.files(mcd43a2_data_dir, pattern=tile, rec=T, full=T))
  modis_dates <- lapply(mcd43a4_in_files, function(x) as.Date(gsub(".*A([0-9]{7})", "\\1", basename(x)), format="%Y%j"))
  mcd43a4_in_files <- lapply(1:length(mcd43a4_in_files), function(x) return(mcd43a4_in_files[[x]][modis_dates[[x]] <= end_date & modis_dates[[x]] >= start_date]))
  mcd43a2_in_files <- lapply(1:length(mcd43a2_in_files), function(x) return(mcd43a2_in_files[[x]][modis_dates[[x]] <= end_date & modis_dates[[x]] >= start_date]))
  modis_dates <- lapply(1:length(modis_dates), function(x) return(modis_dates[[x]][modis_dates[[x]] <= end_date & modis_dates[[x]] >= start_date]))

  # get MODIS cell and line number ranges for the intersecting tiles
  modis_cell_num_ranges <- lapply(unique(modis_tile_indices), function(x) range(modis_cell_nums[modis_tile_indices == x]))
  modis_line_ranges <- lapply(1:length(mcd43a4_in_files), function(x){ tmp_r <- raster(GetSDSName(mcd43a4_in_files[[x]][1], 1)); return(range(rowFromCell(tmp_r, modis_cell_num_ranges[[x]])))})

  # get the MODIS data
  modis_data <- lapply(1:length(mcd43a4_in_files), GetMODISLines, mcd43a4_in_files=mcd43a4_in_files, mcd43a2_in_files=mcd43a2_in_files, modis_line_ranges=modis_line_ranges, landsat_to_modis_bands=landsat_to_modis_bands, max_open_datasets=MAX_OPEN_FILES)

  # get the modis cell number offsets
  modis_cell_offsets <- lapply(1:length(modis_cell_num_ranges), GetModisCellOffset)

  # apply the Kalman Filter
  tmp_KF_results <- parLapply(cl, 1:nrow(landsat_data), DoKF, landsat_data=landsat_data, modis_data=modis_data, modis_cell_nums=modis_cell_nums, modis_tile_indices=modis_tile_indices, modis_cell_offsets=modis_cell_offsets, landsat_dates=landsat_dates, modis_dates=modis_dates, landsat_sensor=landsat_sensor, lwmask_data=lwmask_values)

  # write to output files
  print("Writing output files...")
  out_files <- file.path(output_dir, paste("KF_output", 1:7, sep=""))
  example_r <- raster(landsat_in_files[1])
  WriteKFResults(tmp_KF_results, out_files, example_r)
}
