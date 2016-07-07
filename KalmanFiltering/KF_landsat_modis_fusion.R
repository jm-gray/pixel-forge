#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Image Kalman Filtering for MODIS-Landsat fusion
# Josh Gray, 2016
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Prelims
library(raster)
library(parallel)
library(rgdal)
library(dlm)
library(tools)
library(caTools)
library(RColorBrewer)
source("https://raw.github.com/jm-gray/pixel-forge/master/KalmanFiltering/KF_functions.R")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# main processing
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#-------------------------------------------------------------------------------
# get the landsat data and dates
landsat_data_dir <- "/projectnb/modislc/users/joshgray/DL_Landsat"
landsat_evi2_files <- dir(landsat_data_dir, pattern="evi2_landsat_overlap.tif", rec=T, full=T)
landsat_evi2_files <- landsat_evi2_files[order(GetLandsatDate(landsat_evi2_files))]
landsat_blue_files <- dir(landsat_data_dir, pattern="band1_landsat_overlap.tif", rec=T, full=T)
landsat_blue_files <- landsat_blue_files[order(GetLandsatDate(landsat_blue_files))]
landsat_green_files <- dir(landsat_data_dir, pattern="band2_landsat_overlap.tif", rec=T, full=T)
landsat_green_files <- landsat_green_files[order(GetLandsatDate(landsat_green_files))]
landsat_red_files <- dir(landsat_data_dir, pattern="band3_landsat_overlap.tif", rec=T, full=T)
landsat_red_files <- landsat_red_files[order(GetLandsatDate(landsat_red_files))]
landsat_nir_files <- dir(landsat_data_dir, pattern="band4_landsat_overlap.tif", rec=T, full=T)
landsat_nir_files <- landsat_nir_files[order(GetLandsatDate(landsat_nir_files))]
landsat_dates <- sort(GetLandsatDate(landsat_nir_files))

#-------------------------------------------------------------------------------
# get the modis data and dates
modis_data_dir <- "/projectnb/modislc/users/joshgray/DL_Landsat/MODISOVERLAP"
modis_evi2_files <- dir(modis_data_dir, pattern="evi2_modis_overlap.tif", full=T)
modis_dates <- as.Date(unlist(lapply(modis_evi2_files, GetModisDate)), origin="1970-1-1")
modis_evi2_files <- modis_evi2_files[order(modis_dates)]
modis_blue_files <- dir(modis_data_dir, pattern="blue_modis_overlap.tif", full=T)
modis_blue_files <- modis_blue_files[order(modis_dates)]
modis_green_files <- dir(modis_data_dir, pattern="green_modis_overlap.tif", full=T)
modis_green_files <- modis_green_files[order(modis_dates)]
modis_red_files <- dir(modis_data_dir, pattern="red_modis_overlap.tif", full=T)
modis_red_files <- modis_red_files[order(modis_dates)]
modis_nir_files <- dir(modis_data_dir, pattern="nir_modis_overlap.tif", full=T)
modis_nir_files <- modis_nir_files[order(modis_dates)]
modis_QA_files <- dir(modis_data_dir, pattern="band_quality_modis_overlap.tif", full=T)
modis_QA_files <- modis_QA_files[order(modis_dates)]
modis_snow_files <- dir(modis_data_dir, pattern="snow_modis_overlap.tif", full=T)
modis_snow_files <- modis_snow_files[order(modis_dates)]
modis_dates <- sort(modis_dates)


#---------------------------------------------------------------------
rows_to_read <- 1e3
landsat_v <- GetValuesGDAL(landsat_evi2_files, 1, rows_to_read)
modis_v <- GetValuesGDAL(modis_evi2_files, 1, rows_to_read)
modis_qa_v <- GetValuesGDAL(modis_QA_files, 1, rows_to_read)
V <- cbind(landsat_v, modis_v, modis_qa_v)


# plot demo...
# plot(c(pred_dates, pred_dates, x_dates, y_dates), c(x_smooth, y_smooth, x_v, y_v), type="n", xlab="", ylab="EVI2")
# points(x_dates, x_v, col=1, pch=1)
# points(pred_dates, x_smooth, type="l", col="grey", lwd=1.5)
# points(y_dates, y_v, col=2, pch=2)
# points(pred_dates, y_smooth, type="l", col="pink", lwd=1.5)

# set up the cluster
cl <- makeCluster(16)
clusterEvalQ(cl, {source("https://raw.github.com/jm-gray/pixel-forge/master/KalmanFiltering/KF_functions.R")})

# data may be too large to keep in memory, read and process blocks of lines instead
rows_to_read <- 1e3
tmp_r <- raster(landsat_evi2_files[1])
is <- 1:(ceiling(nrow(tmp_r) / rows_to_read))
for(i in is){
	print(paste("Working on block", i))
	# determine row read parameters
	start_row <- ((i - 1) * rows_to_read) + 1
	# the last block may not have the same number of rows to read as the first
	if(i == is[length(is)]){
		read_rows <- nrow(tmp_r) - ((i - 1) * rows_to_read)
	}else{
		read_rows <- rows_to_read
	}
	# get the data
	landsat_v <- GetValuesGDAL(landsat_evi2_files, start_row, read_rows)
	modis_v <- GetValuesGDAL(modis_evi2_files, start_row, read_rows)
	modis_snow_v <- GetValuesGDAL(modis_snow_files, start_row, read_rows)
	V <- cbind(landsat_v, modis_v, modis_snow_v)
	# calculate MODIS-Landsat RMSE:
	# tmp <- parApply(cl, V, 1, GetSplineRMSE, x_dates=landsat_dates, y_dates=modis_dates)
	tmp <- parApply(cl, V, 1, GetMatchRMSE, x_dates=landsat_dates, y_dates=modis_dates)
	if(i == 1){
		out_rmse <- tmp
	}else{
		out_rmse <- c(out_rmse, tmp)
	}
}
# save(out_rmse, file="/projectnb/modislc/users/joshgray/DL_Landsat/out_rmse.Rdata")
save(out_rmse, file="/projectnb/modislc/users/joshgray/DL_Landsat/out_rmse_match.Rdata")

values(tmp_r) <- out_rmse
qs <- quantile(out_rmse, c(0, 0.02, 0.98, 1), na.rm=T)
breaks <- c(qs[1], seq(qs[2], qs[3], len=254), qs[4])
pal <- colorRampPalette(brewer.pal(11, "Greys"))
plot(tmp_r, breaks=breaks, col=pal(255))


# system.time(tmp <- apply(V, 1, GetSplineRMSE, x_dates=landsat_dates, y_dates=modis_dates))
# system.time(tmp <- parApply(cl, V, 1, GetSplineRMSE, x_dates=landsat_dates, y_dates=modis_dates))
# 474 seconds per 1e3 rows
