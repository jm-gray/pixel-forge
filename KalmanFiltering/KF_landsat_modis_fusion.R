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
# landsat_data_dir <- "/projectnb/modislc/users/joshgray/DL_Landsat"
landsat_data_dir <- "/Users/jmgray2/Documents/NebraskaFusion/Landsat"
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
# modis_data_dir <- "/projectnb/modislc/users/joshgray/DL_Landsat/MODISOVERLAP"
modis_data_dir <- "/Users/jmgray2/Documents/NebraskaFusion/MODIS"
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
# cdl_files <- dir("/projectnb/modislc/users/joshgray/DL_Landsat/CDL_sub", pattern="landsat_overlap.tif", full=T)
cdl_files <- dir("/Users/jmgray2/Documents/NebraskaFusion/CDL", pattern="landsat_overlap.tif", full=T)
cdl_years <- as.numeric(substr(basename(cdl_files), 5, 8))

#---------------------------------------------------------------------
# set up the cluster
# cl <- makeCluster(16)
cl <- makeCluster(8)
clusterEvalQ(cl, {source("https://raw.github.com/jm-gray/pixel-forge/master/KalmanFiltering/KF_functions.R")})

#---------------------------------------------------------------------
# Calculate the RMSE between Landsat-MODIS at the pixel scale
#	This can be accomplished by only considering matching dates, or from splines
#
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
# save(out_rmse, file="/projectnb/modislc/users/joshgray/DL_Landsat/out_rmse_match.Rdata")
save(out_rmse, file="/Users/jmgray2/Documents/NebraskaFusion/out_rmse_match.Rdata")
rmse <- tmp_r
values(rmse) <- out_rmse

#-------------------------------------------------------------------------------
# get a block of data including a particular pixel
rows_to_read <- 1e3
plot(rmse); myc <- click(rmse, cell=T)
start_row <- rowFromCell(tmp_r, myc$cell[1]) - round(rows_to_read / 2)
nrows <- min(rows_to_read, (nrow(tmp_r) - ((i - 1) * rows_to_read))) # last block may have less rows
cdl_v <- GetValuesGDAL(cdl_files, start_row, nrows)
landsat_v <- GetValuesGDAL(landsat_evi2_files, start_row, nrows)
rmse_v <- getValues(rmse, start_row, nrows)
modis_v <- GetValuesGDAL(modis_evi2_files, start_row, nrows)
modis_snow_v <- GetValuesGDAL(modis_snow_files, start_row, nrows)
modis_v[modis_snow_v == 1] <- NA # screen out all snow
V <- cbind(cdl_v, landsat_v, modis_v)
tmp_cell <- ((rows_to_read / 2) * ncol(rmse)) + colFromCell(rmse, myc$cell[1]) # determine the row in V of original cell
num_cdl_years <- length(cdl_years)

scale_factor <- 1e4
landsat_measurement_error <- 0.05
min_quant <- 0.1
year_to_do <- 2008
# x <- V[myc$cell[3], ]
x <- V[tmp_cell, ]
x_v <- x[(num_cdl_years + 1):(length(landsat_dates) + num_cdl_years)] / scale_factor
y_v <- x[(length(landsat_dates) + num_cdl_years + 1):(length(landsat_dates) + num_cdl_years + length(modis_dates))] / scale_factor
# x_v <- log(x_v); y_v <- log(y_v)
x_doys <- as.numeric(strftime(landsat_dates, format="%j"))
x_years <- as.numeric(strftime(landsat_dates, format="%Y"))
y_doys <- as.numeric(strftime(modis_dates, format="%j"))
y_years <- as.numeric(strftime(modis_dates, format="%Y"))
# x_smooth <- GetDOYSpline(x_v, landsat_dates, min_quant=min_quant)
x_smooth <- GetDOYSpline(x_v, landsat_dates, min_quant=min_quant, spline_spar=0.6)
change_rate <- x_smooth[2:length(x_smooth)] / x_smooth[1:(length(x_smooth) - 1)]
change_rate <- c(change_rate, change_rate[length(change_rate)])
model_error <- var(x_v - x_smooth[x_doys], na.rm=T) # determine the model error
x_y_error <- GetLandsatMODISError(x / scale_factor, landsat_dates, modis_dates) # get Landsat-MODIS error model
dlm <- MakeMultiDLM(num_states=1, sensors=2, time_varying=F)
X(dlm) <- matrix(change_rate, ncol=1)
JGG(dlm) <- 1 # time varying system evolution in column 1 of X
FF(dlm)[2, 1] <- x_y_error$slope
diag(V(dlm)) <- c(landsat_measurement_error, x_y_error$rmse)
W(dlm) <- model_error
m0(dlm) <- x_smooth[1]
C0(dlm) <- landsat_measurement_error
x_tmp <- rep(NA, 365)
x_tmp[x_doys[x_years == year_to_do]] <- x_v[x_years == year_to_do]
y_tmp <- rep(NA, 365)
y_tmp[y_doys[y_years == year_to_do]] <- y_v[y_years == year_to_do]
# Y <- cbind(x_tmp, y_tmp)
Y <- cbind(x_tmp, y_tmp + x_y_error$intercept)
smooth <- dlmSmooth(Y, dlm)
filt_m <- ExtractSmoothMeans(smooth)
filt_se <- ExtractSmoothSE(smooth)
PlotForecast(filt_m, filt_se, signal=Y)
points(x_tmp)
points(y_tmp, col=2)
points(x_smooth, type="l")


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
num_cdl_years <- length(cdl_files)
x <- V[myc$cell[3], ]
x_v <- x[(num_cdl_years + 1):(length(landsat_dates) + num_cdl_years)]
y_v <- x[(length(landsat_dates) + num_cdl_years + 1):(length(landsat_dates) + num_cdl_years + length(modis_dates))]

doys <- as.numeric(strftime(landsat_dates, format="%j"))
years <- as.numeric(strftime(landsat_dates, format="%Y"))
# x_smooth <- GetDOYSpline(x_v, landsat_dates, min_quant=min_quant)
x_smooth <- GetDOYSpline(x_v, landsat_dates, min_quant=min_quant, spline_spar=0.75)
# x_smooth_fill <- GetDOYSpline(x_v, landsat_dates, min_quant=min_quant, head_fill_len=30, tail_fill_len=30)

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


#-------------------------------------------------------------------------------
x_int <- approx(1:365, x_tmp, xout=1:365, rule=2)
y_int <- approx(1:365, y_tmp, xout=1:365, rule=2)

Y <- cbind(x_tmp, y_tmp)
Y_int <- cbind(x_int$y, y_int$y)

ff <- matrix(c(1, 1), nrow=2)
gg <- matrix(1)
# v <- diag(1e6, dim(ff)[1])
v <- diag(c(0.01, 0.0169))
w <- 0.001
m0 <- 0
C0 <- 1e6
my_dlm <- dlm(list(m0=m0, C0=C0, FF=ff, GG=gg, V=v, W=w))


# tmp <- dlmSmooth(Y, my_dlm)
# tmp <- dlmSmooth(Y_int, my_dlm)
# plot(x_tmp, pch=16, col=2, ylim=c(0, 0.5))
# points(y_tmp, pch=16, col=4)
# points(dropFirst(tmp$s), type="l")
# legend("topleft", legend=c("x", "y"), col=c(2, 4), pch=16)

# smooth <- dlmSmooth(Y, my_dlm)
smooth <- dlmSmooth(Y_int, my_dlm)
filt_m <- ExtractSmoothMeans(smooth)
filt_se <- ExtractSmoothSE(smooth)
PlotForecast(filt_m, filt_se, signal=Y)
points(x_tmp)
points(y_tmp, col=2)
# points(x_smooth, type="l")

KF_rmse_fuse <- function(x, x_dates, y_dates, do_rmse=T, do_interp=T, x_measurement_error=0.05, y_measurement_error=0.05, model_error=0.001, scale_factor=1e4, do_plot=F, ylim=c(0, 0.7)){
	na_return <- rep(NA, 365)
	if(all(is.na(x))) return(na_return) # handle the case where there is no data

	x_v <- x[2:(length(x_dates) + 1)] / scale_factor
	y_v <- x[(length(x_dates) + 2):(length(x_dates) + 1 + length(y_dates))] / scale_factor

	x_doys <- as.numeric(strftime(x_dates, format="%j"))
	x_tmp <- rep(NA, 365)
	x_tmp[x_doys] <- x_v
	y_doys <- as.numeric(strftime(y_dates, format="%j"))
	y_tmp <- rep(NA, 365)
	y_tmp[y_doys] <- y_v

	# retrieve the landsat-modis RMSE value
	if(do_rmse){
		rmse <- x[1] / scale_factor
		if(is.na(rmse)) rmse <- 1e6
	}else{
		rmse <- y_measurement_error
	}

	# interpolate missing values if requested
	if(do_interp){
		x_sig <- try(approx(1:365, x_tmp, xout=1:365, rule=2)$y)
		if(inherits(x_sig, 'try-error')) return(na_return)
		y_sig <- try(approx(1:365, y_tmp, xout=1:365, rule=2)$y)
		if(inherits(y_sig, 'try-error')) return(na_return)
	}else{
		x_sig <- x_tmp
		y_sig <- y_tmp
	}

	my_dlm <- MakeMultiDLM(num_states=1, sensors=2, time_varying=F)
	diag(V(my_dlm)) <- c(x_measurement_error, rmse)
	W(my_dlm) <- model_error
	m0(my_dlm) <- min(x_sig, na.rm=T)
	Y <- cbind(x_sig, y_sig)
	smooth <- try(dlmSmooth(Y, my_dlm))
	if(inherits(smooth, 'try-error')) return(na_return)
	filt_m <- try(ExtractSmoothMeans(smooth))
	if(inherits(filt_m, 'try-error')) return(na_return)
	filt_se <- try(ExtractSmoothSE(smooth))
	if(inherits(filt_se, 'try-error')) return(na_return)
	if(do_plot){
		PlotForecast(filt_m, filt_se, signal=list(x_tmp, y_tmp), ylim=ylim)
	}
	return(filt_m)
}

good_cell <- 1664260
bad_cell <- 1512279
trash <- KF_rmse_fuse(V[good_cell, ], x_dates=landsat_dates, y_dates=modis_dates, cdl_year=cdl_years, year_to_do=2008, do_plot=T)
trash <- KF_rmse_fuse(V[bad_cell, ], x_dates=landsat_dates, y_dates=modis_dates, cdl_year=cdl_years, year_to_do=2008, do_plot=T)
random_cell <- sample(1:dim(V)[1],1)
trash <- KF_rmse_fuse(V[random_cell, ], x_dates=landsat_dates, y_dates=modis_dates, cdl_year=cdl_years, year_to_do=2008, do_plot=T)

#==================================================================================
cl <- makeCluster(detectCores())
clusterExport(cl, c("ModisOveralapPreprocessNBAR", "GetSDSName", "GetSDSNameQA", "GetSDSNameSnow", "ModisOveralapPreprocessQA"))
clusterEvalQ(cl, {source("https://raw.github.com/jm-gray/pixel-forge/master/KalmanFiltering/KF_functions.R");library(dlm)})
# system.time(trash <- parApply(cl, V[1:1e4,], 1, KF_rmse_fuse, x_dates=landsat_dates, y_dates=modis_dates, cdl_year=cdl_years, year_to_do=2008))


landsat_data_dir <- "~/Desktop/KF_fusion_data_new/landsat"
landsat_evi2_files <- dir(landsat_data_dir, pattern="evi2.tif", full=T)
landsat_evi2_files <- landsat_evi2_files[order(GetLandsatDate(landsat_evi2_files))]
landsat_dates <- GetLandsatDate(landsat_evi2_files)
landsat_doys <- as.integer(strftime(landsat_dates, format="%j"))

modis_data_dir <- "~/Desktop/KF_fusion_data_new/modis"
modis_evi2_files <- dir(modis_data_dir, pattern="evi2.tif", full=T)
modis_dates <- as.Date(unlist(lapply(modis_evi2_files, GetModisDate)), origin="1970-1-1")
modis_evi2_files <- modis_evi2_files[order(modis_dates)]
modis_doys <- as.integer(strftime(modis_dates, format="%j"))

rmse_sub <- raster("~/Desktop/KF_fusion_data_new/rmse_sub.tif")

year_to_do <- 2008
rows_to_read <- 1e3
tmp_r <- raster(landsat_evi2_files[1])
is <- 1:(ceiling(nrow(tmp_r) / rows_to_read))
for(i in is){
	print(paste("Working on block", i))
	# get the data
	start_row <- ((i - 1) * rows_to_read) + 1
	nrows <- min(rows_to_read, (nrow(tmp_r) - ((i - 1) * rows_to_read))) # last block may have less rows
	# cdl_v <- GetValuesGDAL(cdl_files, start_row, nrows)
	landsat_v <- GetValuesGDAL(landsat_evi2_files, start_row, nrows)
	modis_v <- GetValuesGDAL(modis_evi2_files, start_row, nrows)
	# modis_snow_v <- GetValuesGDAL(modis_snow_files, start_row, nrows)
	# modis_v[modis_snow_v == 1] <- NA # screen out all snow
	rmse_v <- getValues(rmse_sub, start_row, nrows)
	# V <- cbind(rmse_v, cdl_v, landsat_v, modis_v)
	# V <- cbind(rmse_v, cdl_v, landsat_v, modis_v)
	V <- cbind(rmse_v, landsat_v, modis_v)
	rm(rmse_v, cdl_v, landsat_v, modis_v)
	# system.time(kf_block <- parApply(cl, V, 1, KF_rmse_fuse, x_dates=landsat_dates, y_dates=modis_dates))
	system.time(kf_block <- parApply(cl, V, 1, KF_rmse_fuse, x_dates=landsat_dates, y_dates=modis_dates, do_interp=F, do_rmse=F))
	# system.time(kf_block <- parApply(cl, V[1:1e4,], 1, KF_rmse_fuse, x_dates=landsat_dates, y_dates=modis_dates))
	# out_file <- file.path("~/Desktop", paste("kf_block_", i, ".Rdata", sep=""))
	# save(kf_block, file=out_file)
	if(i == 1){
		kf_results <- kf_block
	}else{
		kf_results <- cbind(kf_results, kf_block)
	}
}

# kf_out_dir <- "~/Desktop/KF_fusion_data_new/rmse_kf_output"
kf_out_dir <- "~/Desktop/KF_fusion_data_new/no_rmse_kf_output"
tmp_r <- rmse_sub
for(i in 1:365){
	values(tmp_r) <- kf_block[i, ]
	# out_file <- file.path(kf_out_dir, paste("kf_rmse_doy", formatC(i, width=3, flag="0"), ".tif", sep=""))
	out_file <- file.path(kf_out_dir, paste("kf_no_rmse_doy", formatC(i, width=3, flag="0"), ".tif", sep=""))
	writeRaster(tmp_r, file=out_file)
}

# do the fusion plotting
# fig_out_dir <- "~/Desktop/KF_fusion_data_new/rmse_kf_figs"
fig_out_dir <- "~/Desktop/KF_fusion_data_new/no_rmse_kf_figs"
# out_prefix <- "rmse_kf_"
out_prefix <- "no_rmse_kf_"
blank_raster <- rmse_sub
values(blank_raster) <- NA
modis_r <- landsat_r <- blank_raster
for(i in 1:365){
	kf_r <- raster(dir(kf_out_dir, pattern=formatC(i, width=3, flag="0"), full=T))
	landsat_ind <- which(landsat_doys == i)
	if(length(landsat_ind)!=0){
		landsat_r <- raster(landsat_evi2_files[landsat_ind]) / scale_factor
	}
	modis_ind <- which(modis_doys == i)
	if(length(modis_ind)!=0){
		modis_r <- raster(modis_evi2_files[modis_ind]) / scale_factor
	}
  print(paste("Plotting day", i))
  out_file <- file.path(fig_out_dir, paste(out_prefix, formatC(i, width=3, flag="0"), ".jpg", sep=""))
  jpeg(file=out_file, width=1214, height=round(772+(772/2)), quality=75)
  PlotFusion(landsat_r, modis_r, kf_r, label=paste("DOY:", i))
  dev.off()
}

# ffmpeg -framerate 10 -i rmse_kf_%03d.jpg -pix_fmt yuv420p -vf scale="720:trunc(ow/a/2)*2" rmse_kf_fast.mp4
# ffmpeg -framerate 10 -i no_rmse_kf_%03d.jpg -pix_fmt yuv420p -vf scale="720:trunc(ow/a/2)*2" no_rmse_kf_fast.mp4

#---------------------------------------------------------
# analyze GDD
gdd_data_dir <- "~/Desktop/KF_fusion_data_new/gdd"
gdd_files <- dir(gdd_data_dir, pattern="tif", full=T)
gdd_v <- GetValuesGDAL(gdd_files, start_row, nrows)

#-------------------------------
# preprocess to new subset area
rmse_sub <- raster("~/Desktop/KF_fusion_data_new/rmse_sub.tif")
year_to_do <- 2008
landsat_data_dir <- "/Users/jmgray2/Documents/NebraskaFusion/Landsat"
landsat_evi2_files <- dir(landsat_data_dir, pattern="evi2_landsat_overlap.tif", rec=T, full=T)
landsat_dates <- GetLandsatDate(landsat_evi2_files)
landsat_years <- as.integer(strftime(landsat_dates, format="%Y"))
landsat_evi2_files <- landsat_evi2_files[landsat_years == year_to_do]
landsat_out_dir <- "~/Desktop/KF_fusion_data_new/landsat"
i <- 1
for(in_file in landsat_evi2_files){
	print(paste("Processing", i, "of", length(landsat_evi2_files)))
	r <- raster(in_file)
	r_sub <- crop(r, extent(rmse_sub))
	out_file <- file.path(landsat_out_dir, paste(unlist(strsplit(basename(in_file), split="_"))[1], "_sub_evi2.tif", sep=""))
	writeRaster(r_sub, file=out_file)
	i <- i + 1
}

modis_data_dir <- "/Users/jmgray2/Documents/NebraskaFusion/MODIS"
modis_evi2_files <- dir(modis_data_dir, pattern="evi2_modis_overlap.tif", full=T)
modis_dates <- as.Date(unlist(lapply(modis_evi2_files, GetModisDate)), origin="1970-1-1")
modis_years <- as.integer(strftime(modis_dates, format="%Y"))
modis_evi2_files <- modis_evi2_files[modis_years == year_to_do]
modis_out_dir <- "~/Desktop/KF_fusion_data_new/modis"
i <- 1
for(in_file in modis_evi2_files){
	print(paste("Processing", i, "of", length(modis_evi2_files)))
	r <- raster(in_file)
	r_sub <- crop(r, extent(rmse_sub))
	out_file <- file.path(modis_out_dir, paste(unlist(strsplit(basename(in_file), split="_"))[1], "_sub_evi2.tif", sep=""))
	writeRaster(r_sub, file=out_file)
	i <- i + 1
}

gdd_data_dir <- "/Users/jmgray2/Desktop/UofIMet_GDD_resampled"
gdd_files <- dir(gdd_data_dir, pattern="tif", full=T)
gdd_dates <- as.Date(basename(gdd_files), format="gdd_%Y%j")
gdd_years <- as.integer(strftime(gdd_dates, format="%Y"))
gdd_files <- gdd_files[gdd_years == year_to_do]
gdd_out_dir <- "~/Desktop/KF_fusion_data_new/gdd"
i <- 1
for(in_file in gdd_files){
	print(paste("Processing", i, "of", length(gdd_files)))
	r <- raster(in_file)
	r_sub <- crop(r, extent(rmse_sub))
	out_file <- file.path(gdd_out_dir, paste(gsub("(.*).tif$", "\\1", basename(in_file)), "_sub.tif", sep=""))
	writeRaster(r_sub, file=out_file)
	i <- i + 1
}




# plot the rmse figure
# par(col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white")
pal <- colorRampPalette(brewer.pal(11, "Spectral"))
v <- values(rmse_sub)
qs <- quantile(v, c(0,0.02, 0.98,1),na.rm=T)
breaks <- c(qs[1], seq(qs[2], qs[3], len=254), qs[4])
LEGENDAXISCEX <- 1
LEGENDMAINCEX <- 1
LEGENDWIDTH <- 2
par(oma=rep(2,4))
# plot(lwmask, breaks=c(-1, 0.5, 1.5, 10), col=c(WATERCOLOR, LANDCOLOR, WATERCOLOR), xaxt="n", yaxt="n", legend=F, bty="n", box=FALSE)
plot(rmse_sub, breaks=breaks, col=pal(255), maxpixels=900e3, legend=F, xaxt="n", yaxt="n", bty="n", box=FALSE)
legend_at <- round(seq(breaks[2], breaks[length(breaks) - 1], len=7))
legend_at_date <- legend_at/scale_factor
legend_labels <- c(paste("<", legend_at_date[1]), as.character(legend_at_date[2:(length(legend_at_date) - 1)]), paste(">", legend_at_date[length(legend_at_date)]))
# par(col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white")
plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=pal(length(breaks)-1), legend.width=LEGENDWIDTH, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=LEGENDAXISCEX), legend.args=list(text="", side=3, font=2, line=0.5, cex=LEGENDMAINCEX, col.axis="red", col.sub="red"))
# title("Medsian")
# title("MODIS-Landsat long-term RMSE")

par(bg="white")
plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=pal(length(breaks)-1), legend.width=1, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=LEGENDAXISCEX), legend.args=list(text="", side=3, font=2, line=0.5, cex=LEGENDMAINCEX, col.axis="red", col.sub="red"))

# Combining two Normal Distributions:
# uinit <- 1.25
# varinit <- 1.5
#
# u0 <- 2
# var0 <- 2
# u1 <- 3
# var1 <- 4
#
# u2 <- u0 + ((var0 * (u1 - u0))/(var1 + var0))
# var2 <- var0 - ((var0^2)/(var0+var1))
#
# pal <- colorRampPalette(c("red", "orange", "yellow"))
# cols <- pal(3)
#
# par(col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white")
# t <- seq(-5, 10, len=1e3)
# plot(t, dnorm(t, uinit, sqrt(varinit)), type="l", col=cols[1], lwd=3, ylim=c(0, 0.4), xlab="value", ylab="density")
# plot(t, dnorm(t, uinit, sqrt(varinit)), type="l", col=rgb(1,0,0,0.25), lwd=3, ylim=c(0, 0.4))
# points(t, dnorm(t, u0, sqrt(var0)), type="l", col=cols[1], lwd=3)
# points(t, dnorm(t, u1, sqrt(var1)), type="l", col=cols[3], lwd=3)
# points(t, dnorm(t, u2, sqrt(var2)), type="l", col=cols[2], lwd=3)
