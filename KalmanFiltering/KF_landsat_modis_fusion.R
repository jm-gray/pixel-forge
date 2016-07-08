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

#-------------------------------------------------------------------------------
# get the CDL files
cdl_files <- dir("/projectnb/modislc/users/joshgray/DL_Landsat/CDL_sub", pattern="landsat_overlap.tif", full=T)
cdl_years <- as.numeric(substr(basename(cdl_files), 5, 8))

#---------------------------------------------------------------------
# set up the cluster
cl <- makeCluster(16)
clusterEvalQ(cl, {source("https://raw.github.com/jm-gray/pixel-forge/master/KalmanFiltering/KF_functions.R")})

# data may be too large to keep in memory, read and process blocks of lines instead
rows_to_read <- 1e3
tmp_r <- raster(landsat_evi2_files[1])
is <- 1:(ceiling(nrow(tmp_r) / rows_to_read))
for(i in is){
	print(paste("Working on block", i))
	# get the data
	start_row <- ((i - 1) * rows_to_read) + 1
	nrows <- min(rows_to_read, (nrow(tmp_r) - ((i - 1) * rows_to_read))) # last block may have less rows
	cdl_v <- GetValuesGDAL(cdl_files, start_row, nrows)
	landsat_v <- GetValuesGDAL(landsat_evi2_files, start_row, nrows)
	modis_v <- GetValuesGDAL(modis_evi2_files, start_row, nrows)
	modis_snow_v <- GetValuesGDAL(modis_snow_files, start_row, nrows)
	modis_v[modis_snow_v == 1] <- NA # screen out all snow
	# V <- cbind(landsat_v, modis_v, modis_snow_v)
	V <- cbind(cdl_v, landsat_v, modis_v)
	# calculate MODIS-Landsat RMSE:
	# tmp <- parApply(cl, V, 1, GetSplineRMSE, x_dates=landsat_dates, y_dates=modis_dates)
	tmp <- parApply(cl, V, 1, GetMatchRMSE, x_dates=landsat_dates, y_dates=modis_dates, num_cdl_years=length(cdl_files))
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


# plot demo...
# plot(c(pred_dates, pred_dates, x_dates, y_dates), c(x_smooth, y_smooth, x_v, y_v), type="n", xlab="", ylab="EVI2")
plot(c(x_dates, y_dates), c(x_v, y_v), type="n", xlab="", ylab="EVI2")
points(x_dates, x_v, col=1, pch=1)
# points(pred_dates, x_smooth, type="l", col="grey", lwd=1.5)
points(y_dates, y_v, col=2, pch=2)
points(y_dates[y_snow_v == 1], y_v[y_snow_v == 1], col=4, pch=2)
# points(pred_dates, y_smooth, type="l", col="pink", lwd=1.5)


#-------------------------------------------------------------------------------
# Create a figure demonstrating how a prior is determined from the full landsat time series
# > myc
#      cell     value
# 1 1861199  591.9791
# 2 1775462 1226.9110
# 3 1061973  973.9092
# 4 1362014 1773.0994 # Double crop? Alfalfa?
# 5  811625 1029.8065
# 6 1476564  669.3034
# 7 2507182  879.5399
# 8 2492610 1168.7647

mycols <- unlist(lapply(c("Reds", "Blues", "Greens", "Purples"), function(x) brewer.pal(5, x)[4]))
min_quant <- 0.1
x <- V[myc$cell[4], ]
x_v <- x[1:length(x_dates)] / 1e4
y_v <- x[(length(x_dates) + 1):(length(x_dates) + length(y_dates))] / 1e4
y_snow_v <- x[(length(x_dates) + length(y_dates) + 1):length(x)]
doys <- as.numeric(strftime(landsat_dates, format="%j"))
years <- as.numeric(strftime(landsat_dates, format="%Y"))
x_smooth <- GetDOYSpline(x_v, landsat_dates, min_quant=min_quant)

layout(matrix(c(1,1,2,3), nrow=2, byrow=T))
par(mar=c(2, 4, 1, 1), oma=rep(1, 4))

plot(landsat_dates, x_v, pch=16, col=brewer.pal(5, "Greys")[4], xlab="", ylab="EVI2", cex=1.5)
abline(h=quantile(x_v, min_quant, na.rm=T), lty=2, lwd=2, col=brewer.pal(5, "Greys")[4])
# legend("topright", legend=c("Landsat EVI2", "5% background value"), pch=c(16, NA), lty=c(NA, 2), lwd=c(NA, 2), bty="o", horiz=T)

par(mar=c(4, 4, 1, 1))
plot(doys, x_v, xlab="DOY", ylab="EVI2", type="n")
trash <- lapply(sort(unique(years)), function(x) points(doys[years == x], x_v[years == x], col=mycols[x - min(years) + 1], type="p", pch=16))
points(1:365, x_smooth, type="l", col=brewer.pal(5, "Greys")[4], lwd=3)
legend("topleft", legend=c(sort(unique(years)), "spline"), col=c(mycols[1:length(sort(unique(years)))], brewer.pal(5, "Greys")[4]), pch=c(rep(16, length(sort(unique(years)))), NA), lwd=c(rep(NA, length(sort(unique(years)))), 3), cex=1.5, bty="n")

change_rate <- x_smooth[2:length(x_smooth)] / x_smooth[1:(length(x_smooth) - 1)]
plot(1:365, c(change_rate, change_rate[length(change_rate)]), type="l", col=brewer.pal(5, "Greys")[4], lwd=3, xlab="DOY", ylab="Change Rate")
abline(h=1, lty=2, lwd=1, col=1)
# par(new=T)
# plot(1:365, x_smooth, type="l", lty=3, col=brewer.pal(5, "Greys")[4], lwd=3, xlab="", ylab="", xaxt="n", yaxt="n")
