# Band-to-band mapping for Landsat-MODIS
# TM/OLI    MODIS   Name
# 1/2       3       Blue
# 2/3       4       Green
# 3/4       1       Red
# 4/5       2       NIR
# 5/6       6       SWIR
# 7         7       SWIRL
# 6/10      -       Thermal/Brightness-Temp

modis_cell_maps_output_dir <- "/Volumes/research/fer/jmgray2/EastKalimantan/MODIS_cell_maps"
landsat_data_dir <- "/Volumes/research/fer/jmgray2/EastKalimantan/processedData/envi"
lwmask_dir <- "/Volumes/research/fer/jmgray2/EastKalimantan/ClipMasks"
landsat_suffix <- "_EK"
mcd43a4_data_dir <- "/Volumes/research/fer/jmgray2/EastKalimantan/MCD43A4"
mcd43a2_data_dir <- "/Volumes/research/fer/jmgray2/EastKalimantan/MCD43A2"
path_row_shp_file <- "/Volumes/research/fer/jmgray2/EastKalimantan/wrs2_descending.shp"
modis_grid_shp_file <- "/Volumes/research/fer/jmgray2/EastKalimantan/modis_grid/modis_sinusoidal_grid_world.shp"
start_date <- as.Date("2010-1-1")
end_date <- as.Date("2012-12-31")
# path_row_to_process <- "116060"
path_row_to_process <- "116059"
parallel_cores <- NULL
landsat_to_modis_bands <- c(3, 4, 1, 2, 6, NA, 7) # maps post-EK_preprocess landsat band ordering to MCD43A4 band numbers; there is not thermal band in modis, so landsat band 6 returns NA; also QA is contained in MCD43A2 so is not included (handled separately in data acquisition)

landsat_in_files <- dir(file.path(landsat_data_dir, path_row_to_process), pattern=paste(".*", landsat_suffix, "$", sep=""), full=T)
landsat_dates <- as.Date(gsub(".*([0-9]{8}).*", "\\1", basename(landsat_in_files)), format="%Y%m%d")
landsat_in_files <- landsat_in_files[landsat_dates <= end_date & landsat_dates >= start_date]
landsat_dates <- landsat_dates[landsat_dates <= end_date & landsat_dates >= start_date]
landsat_sensor <- gsub(pattern="(L[E|T|C]0[7|5|8]).*", "\\1", basename(landsat_in_files))

mcd43a4_in_files <- list.files(mcd43a4_data_dir, pattern="h29v08", rec=T, full=T)
mcd43a2_in_files <- list.files(mcd43a2_data_dir, pattern="h29v08", rec=T, full=T)
modis_dates <- as.Date(sapply(mcd43a4_in_files, function(x) as.Date(gsub(".*A([0-9]{7})", "\\1", basename(x)), format="%Y%j")), origin="1970-1-1")
mcd43a4_in_files <- mcd43a4_in_files[modis_dates >= start_date & modis_dates <= end_date]
mcd43a2_in_files <- mcd43a2_in_files[modis_dates >= start_date & modis_dates <= end_date]
modis_dates <- modis_dates[modis_dates >= start_date & modis_dates <= end_date]

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Constants, example rasters, and plotting extents
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
tmp_r <- raster(landsat_in_files[1])
kf_dates <- seq.Date(as.Date("2010-1-1"), as.Date("2012-12-31"), by="week")
num_time_steps <- length(kf_dates)
num_cols <- ncol(tmp_r)
# rows_to_read <- 1200
rows_to_read <- 1062
start_row <- 2600
int_size <- 4 # in bytes
start_position <- num_time_steps * num_cols * start_row * int_size
records_to_read <- num_time_steps * num_cols * rows_to_read

red_file <- "~/Desktop/KF_output3"
nir_file <- "~/Desktop/KF_output4"
swir_file <- "~/Desktop/KF_output5"

kf_output_dir <- "/Users/jmgray2/Google Drive/Presentations/CoffeeAndViz/KF_CandV_output"
landsat_output_dir <- "/Users/jmgray2/Google Drive/Presentations/CoffeeAndViz/Landsat_CandV_output"
modis_output_dir <- "/Users/jmgray2/Google Drive/Presentations/CoffeeAndViz/MODIS_CandV_output"
kf_output_root <- "KF_CandV_"
landsat_output_root <- "Landsat_CandV_"
modis_output_root <- "MODIS_CandV_"
tmp_modis_out <- "/Users/jmgray2/Google\\ Drive/Presentations/CoffeeAndViz/tmp"

# Center res: 4608 x 1062
# Left res: 4608 x 1200
# Right res: 3328 x 1200

# on the job that got killed:
# 15175431600 bytes written
# 15175431600 / 4 / 157 / num_cols = 3700 lines!

# Coffe & Viz plotting extent: 2499:3698, 700:5307 (rows, columns)
plot_extent <- extent(tmp_r, start_row, start_row + rows_to_read - 1, 700, 700 + 4608 - 1)
# tmp_r_crop <- crop(tmp_r, extent(tmp_r, start_row, start_row + rows_to_read - 1, 1, num_cols))
tmp_r_crop <- crop(tmp_r, plot_extent)
blank_raster <- raster(matrix(NA, nrow=rows_to_read, ncol=num_cols))
# extent(blank_raster) <- extent(tmp_r_crop)
extent(blank_raster) <- extent(tmp_r, start_row, start_row + rows_to_read - 1, 1, num_cols)

DrawProgressBar <- function(this_date, dates, ypos=178000, xextents=c(528000, 661310), cex_text=2.5, lwd=5, cex_point=2.5){
  par(font=2)
  lines(x=xextents, y=rep(ypos, 2), lwd=lwd, col="white")
  ann_dates <- seq.Date(as.Date(paste(strftime(min(dates), format="%Y"), "-1-1", sep="")), as.Date(paste(as.integer(strftime(max(dates), format="%Y")) + 1, "-1-1", sep="")), by="quarter")
  ann_x_at <- seq(xextents[1], xextents[2], len=length(ann_dates))
  text(x=ann_x_at, y=ypos, labels=ann_dates, col="white", cex=cex_text, adj=c(0.5, -0.3))
  marker_x_pos <- xextents[1] + diff(xextents) * (as.integer(this_date - min(ann_dates)) / as.integer(max(ann_dates) - min(ann_dates)))
  points(x=marker_x_pos, y=ypos, pch=16, col="white", cex=cex_point)
}


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Get the data from the binary files
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# get band 3 (red) data
fcon_red <- file(red_file, open="rb")
seek(fcon_red, where=start_position)
red_data <- readBin(fcon_red, what="integer", n=records_to_read)
close(fcon_red)
# get band 4 (nir) data
fcon_nir <- file(nir_file, open="rb")
seek(fcon_nir, where=start_position)
nir_data <- readBin(fcon_nir, what="integer", n=records_to_read)
close(fcon_nir)
# get band 5 (swir) data
fcon_swir <- file(swir_file, open="rb")
seek(fcon_swir, where=start_position)
swir_data <- readBin(fcon_swir, what="integer", n=records_to_read)
close(fcon_swir)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# make the false color composites (543) from KF output
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
for(time_step in 1:num_time_steps){
  print(paste("Making", time_step, "of", num_time_steps))
  # make a red raster
  tmp_red <- blank_raster
  values(tmp_red) <- red_data[seq(time_step, records_to_read, by=num_time_steps)]
  tmp_red_crop <- crop(tmp_red, plot_extent)
  rm(tmp_red)

  # make a nir raster
  tmp_nir <- blank_raster
  values(tmp_nir) <- nir_data[seq(time_step, records_to_read, by=num_time_steps)]
  tmp_nir_crop <- crop(tmp_nir, plot_extent)
  rm(tmp_nir)

  # make a swir raster
  tmp_swir <- blank_raster
  values(tmp_swir) <- swir_data[seq(time_step, records_to_read, by=num_time_steps)]
  tmp_swir_crop <- crop(tmp_swir, plot_extent)
  rm(tmp_swir)

  kf_s <- stack(tmp_swir_crop, tmp_nir_crop, tmp_red_crop)

  kf_output_file <- file.path(kf_output_dir, paste(kf_output_root, formatC(time_step, width=3, flag="0"), ".jpg", sep=""))
  jpeg(file=kf_output_file, width=4608, height=1062, bg="black")
  plotRGB(kf_s, stretch="lin", maxpixels=5e6, colNA="black")
  DrawProgressBar(kf_dates[time_step], kf_dates)
  dev.off()
}
# ffmpeg -framerate 5 -pattern_type glob -i '*.jpg' -pix_fmt yuv420p KF_CandV.mp4


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# make the false color composites (543) from Landsat input
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
GetClosestDate <- function(this_date, dates, max_lag=NULL){
  date_diff <- abs(this_date - dates)
  min_date_diff <- min(date_diff)
  which_min <- which(date_diff == min_date_diff)
  if(length(which_min) > 1) which_min <- which_min[1]
  if(is.null(max_lag)) return(which_min)
  if(min_date_diff <= max_lag) return(which_min)
  return(NA)
}

# do the plotting
blank_s <- kf_s
values(blank_s) <- NA
for(time_step in 1:num_time_steps){
  print(paste("Making", time_step, "of", num_time_steps))
  matched_landsat_index <- GetClosestDate(kf_dates[time_step], landsat_dates, max_lag=3)
  if(is.na(matched_landsat_index)){
    crop_s <- blank_s
  }else{
    if(landsat_sensor[matched_landsat_index] == "LE07" | landsat_sensor[matched_landsat_index] == "LT05") tmp_landsat_s <- stack(landsat_in_files[matched_landsat_index], bands=c(5, 4, 3))
    if(landsat_sensor[matched_landsat_index] == "LC08") tmp_landsat_s <- stack(landsat_in_files[matched_landsat_index], bands=c(6, 5, 4))
    crop_s <- crop(tmp_landsat_s, plot_extent)
  }

  landsat_output_file <- file.path(landsat_output_dir, paste(landsat_output_root, formatC(time_step, width=3, flag="0"), ".jpg", sep=""))
  jpeg(file=landsat_output_file, width=4608, height=1062, bg="black")
  plotRGB(crop_s, stretch="lin", maxpixels=5e6, colNA="black")
  DrawProgressBar(kf_dates[time_step], kf_dates)
  dev.off()
}
# ffmpeg -framerate 5 -pattern_type glob -i '*.jpg' -pix_fmt yuv420p Landsat_CandV.mp4



#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# make the false color composites (543) from MODIS input
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
for(time_step in 1:num_time_steps){
  modis_index <- which(modis_dates == kf_dates[time_step])
  red_sds <- GetSDSName(mcd43a4_in_files[modis_index], 1)
  nir_sds <- GetSDSName(mcd43a4_in_files[modis_index], 2)
  swir_sds <- GetSDSName(mcd43a4_in_files[modis_index], 6)

  # gdalwarp -s_srs '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs' -t_srs '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs' -ot Int16 CONUSDormancyMedian.tif CONUSDormancyMedian_AEANAD83.tif
  tmp_red_out <- file.path(tmp_modis_out, "tmp_red.tif")
  red_gdal_cmd <- paste("/Library/Frameworks/GDAL.framework/Programs/gdalwarp -s_srs '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs' -t_srs '+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0' -te 525532.5 147380.9 663777 179240 -ts 4608 1062", red_sds, tmp_red_out)
  system(red_gdal_cmd)

  tmp_nir_out <- file.path(tmp_modis_out, "tmp_nir.tif")
  nir_gdal_cmd <- paste("/Library/Frameworks/GDAL.framework/Programs/gdalwarp -s_srs '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs' -t_srs '+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0' -te 525532.5 147380.9 663777 179240 -ts 4608 1062", nir_sds, tmp_nir_out)
  system(nir_gdal_cmd)

  tmp_swir_out <- file.path(tmp_modis_out, "tmp_swir.tif")
  swir_gdal_cmd <- paste("/Library/Frameworks/GDAL.framework/Programs/gdalwarp -s_srs '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs' -t_srs '+proj=utm +zone=50 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0' -te 525532.5 147380.9 663777 179240 -ts 4608 1062", swir_sds, tmp_swir_out)
  system(swir_gdal_cmd)

  modis_s <- stack("/Users/jmgray2/Google Drive/Presentations/CoffeeAndViz/tmp/tmp_swir.tif", "/Users/jmgray2/Google Drive/Presentations/CoffeeAndViz/tmp/tmp_nir.tif", "/Users/jmgray2/Google Drive/Presentations/CoffeeAndViz/tmp/tmp_red.tif")

  modis_output_file <- file.path(modis_output_dir, paste(modis_output_root, formatC(time_step, width=3, flag="0"), ".jpg", sep=""))
  jpeg(file=modis_output_file, width=4608, height=1062, bg="black")
  plotRGB(modis_s, stretch="lin", maxpixels=5e6, colNA="black")
  DrawProgressBar(kf_dates[time_step], kf_dates)
  dev.off()

  file.remove("/Users/jmgray2/Google Drive/Presentations/CoffeeAndViz/tmp/tmp_swir.tif")
  file.remove("/Users/jmgray2/Google Drive/Presentations/CoffeeAndViz/tmp/tmp_nir.tif")
  file.remove("/Users/jmgray2/Google Drive/Presentations/CoffeeAndViz/tmp/tmp_red.tif")
}
# ffmpeg -framerate 5 -pattern_type glob -i '*.jpg' -pix_fmt yuv420p MODIS_CandV.mp4
