library(tools)
library(raster)
library(rgdal)
library(parallel)
library(RColorBrewer)
library(dlm)
library(reshape2)

source("~/Documents/pixel-forge/KalmanFiltering/KF_functions.R")
cl <- makeCluster(detectCores())
clusterEvalQ(cl, {source("~/Documents/pixel-forge/KalmanFiltering/KF_functions.R"); library(dlm); library(reshape2)})

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
# Calculate the RMSE between Landsat-MODIS at the pixel scale
#	This can be accomplished by only considering matching dates, or from splines
#
# data may be too large to keep in memory, read and process blocks of lines instead
rows_to_read <- 1e3
tmp_r <- raster(landsat_blue_files[1])
is <- 1:(ceiling(nrow(tmp_r) / rows_to_read))
for(i in is){
	print(paste("Working on block", i))
	# get the data
	start_row <- ((i - 1) * rows_to_read) + 1
	nrows <- min(rows_to_read, (nrow(tmp_r) - ((i - 1) * rows_to_read))) # last block may have less rows
	# cdl_v <- GetValuesGDAL(cdl_files, start_row, nrows)
	landsat_v <- GetValuesGDAL(landsat_blue_files, start_row, nrows)
  landsat_v[landsat_v > 1e4] <- NA
	modis_v <- GetValuesGDAL(modis_blue_files, start_row, nrows)
	modis_snow_v <- GetValuesGDAL(modis_snow_files, start_row, nrows)
	modis_v[modis_snow_v == 1] <- NA # screen out all snow
	# V <- cbind(landsat_v, modis_v, modis_snow_v)
	V <- cbind(landsat_v, modis_v)
	# calculate MODIS-Landsat RMSE:
	# tmp <- parApply(cl, V, 1, GetSplineRMSE, x_dates=landsat_dates, y_dates=modis_dates)
	tmp <- parApply(cl, V, 1, GetMatchRMSE, x_dates=landsat_dates, y_dates=modis_dates, num_cdl_years=0)
	if(i == 1){
		out_rmse_blue <- tmp
	}else{
		out_rmse_blue <- c(out_rmse_blue, tmp)
	}

  landsat_v <- GetValuesGDAL(landsat_green_files, start_row, nrows)
  landsat_v[landsat_v > 1e4] <- NA
	modis_v <- GetValuesGDAL(modis_green_files, start_row, nrows)
	modis_snow_v <- GetValuesGDAL(modis_snow_files, start_row, nrows)
	modis_v[modis_snow_v == 1] <- NA # screen out all snow
	# V <- cbind(landsat_v, modis_v, modis_snow_v)
	V <- cbind(landsat_v, modis_v)
	# calculate MODIS-Landsat RMSE:
	# tmp <- parApply(cl, V, 1, GetSplineRMSE, x_dates=landsat_dates, y_dates=modis_dates)
	tmp <- parApply(cl, V, 1, GetMatchRMSE, x_dates=landsat_dates, y_dates=modis_dates, num_cdl_years=0)
	if(i == 1){
		out_rmse_green <- tmp
	}else{
		out_rmse_green <- c(out_rmse_green, tmp)
	}

  landsat_v <- GetValuesGDAL(landsat_red_files, start_row, nrows)
  landsat_v[landsat_v > 1e4] <- NA
	modis_v <- GetValuesGDAL(modis_red_files, start_row, nrows)
	modis_snow_v <- GetValuesGDAL(modis_snow_files, start_row, nrows)
	modis_v[modis_snow_v == 1] <- NA # screen out all snow
	# V <- cbind(landsat_v, modis_v, modis_snow_v)
	V <- cbind(landsat_v, modis_v)
	# calculate MODIS-Landsat RMSE:
	# tmp <- parApply(cl, V, 1, GetSplineRMSE, x_dates=landsat_dates, y_dates=modis_dates)
	tmp <- parApply(cl, V, 1, GetMatchRMSE, x_dates=landsat_dates, y_dates=modis_dates, num_cdl_years=0)
	if(i == 1){
		out_rmse_red <- tmp
	}else{
		out_rmse_red <- c(out_rmse_red, tmp)
	}

  landsat_v <- GetValuesGDAL(landsat_nir_files, start_row, nrows)
  landsat_v[landsat_v > 1e4] <- NA
	modis_v <- GetValuesGDAL(modis_nir_files, start_row, nrows)
	modis_snow_v <- GetValuesGDAL(modis_snow_files, start_row, nrows)
	modis_v[modis_snow_v == 1] <- NA # screen out all snow
	# V <- cbind(landsat_v, modis_v, modis_snow_v)
	V <- cbind(landsat_v, modis_v)
	# calculate MODIS-Landsat RMSE:
	# tmp <- parApply(cl, V, 1, GetSplineRMSE, x_dates=landsat_dates, y_dates=modis_dates)
	tmp <- parApply(cl, V, 1, GetMatchRMSE, x_dates=landsat_dates, y_dates=modis_dates, num_cdl_years=0)
	if(i == 1){
		out_rmse_nir <- tmp
	}else{
		out_rmse_nir <- c(out_rmse_nir, tmp)
	}
}
rmse_blue <- rmse_green <- rmse_red <- rmse_nir <- tmp_r
values(rmse_blue) <- out_rmse_blue
values(rmse_green) <- out_rmse_green
values(rmse_red) <- out_rmse_red
values(rmse_nir) <- out_rmse_nir
save(rmse_blue, rmse_green, rmse_red, rmse_nir, file="~/Desktop/KF_fusion_data_new/multispectral_rmse.Rdata")

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
# gather the multispectral data
landsat_data_dir <- "~/Desktop/KF_fusion_data_new/landsat_new"
landsat_blue_files <- dir(landsat_data_dir, pattern=".*band1.*.tif$", full=T)
landsat_green_files <- dir(landsat_data_dir, pattern=".*band2.*.tif$", full=T)
landsat_red_files <- dir(landsat_data_dir, pattern=".*band3.*.tif$", full=T)
landsat_nir_files <- dir(landsat_data_dir, pattern=".*band4.*.tif$", full=T)

tmp_r <- raster(landsat_blue_files[1])

landsat_sensor <- gsub(pattern="(L[E|T][7|5]).*", "\\1", basename(landsat_blue_files))
landsat_dates <- as.Date(gsub(pattern="L[E|T][7|5][0-9]{6}([0-9]{7}).*", "\\1", basename(landsat_blue_files)), format="%Y%j")

landsat_blue_files <- landsat_blue_files[order(landsat_dates)]
landsat_green_files <- landsat_green_files[order(landsat_dates)]
landsat_red_files <- landsat_red_files[order(landsat_dates)]
landsat_nir_files <- landsat_nir_files[order(landsat_dates)]

landsat_dates <- sort(landsat_dates)
landsat_doys <- as.integer(strftime(landsat_dates, format="%j"))
landsat_years <- as.integer(strftime(landsat_dates, format="%Y"))

Y_landsat_blue <- GetValuesGDAL(landsat_blue_files, start_row=1, n=nrow(tmp_r))
Y_landsat_blue[Y_landsat_blue > 1e4] <- NA
# Y_landsat_blue <- t(parApply(cl, Y_landsat_blue, 1, QuantileMinReplacement))

Y_landsat_green <- GetValuesGDAL(landsat_green_files, start_row=1, n=nrow(tmp_r))
Y_landsat_green[Y_landsat_green > 1e4] <- NA
# Y_landsat_green <- t(parApply(cl, Y_landsat_green, 1, QuantileMinReplacement))

Y_landsat_red <- GetValuesGDAL(landsat_red_files, start_row=1, n=nrow(tmp_r))
Y_landsat_red[Y_landsat_red > 1e4] <- NA
# Y_landsat_red <- t(parApply(cl, Y_landsat_red, 1, QuantileMinReplacement))

Y_landsat_nir <- GetValuesGDAL(landsat_nir_files, start_row=1, n=nrow(tmp_r))
Y_landsat_nir[Y_landsat_nir > 1e4] <- NA
# Y_landsat_nir <- t(parApply(cl, Y_landsat_nir, 1, QuantileMinReplacement))

# i <- 30
# s <- stack(landsat_blue_files[i], landsat_green_files[i], landsat_red_files[i], landsat_nir_files[i])
# plotRGB(s, 4, 3, 2, stretch="lin")

modis_data_dir <- "~/Desktop/KF_fusion_data_new/modis_new"
modis_blue_files <- dir(modis_data_dir, pattern=".*blue.*.tif$", full=T)
modis_green_files <- dir(modis_data_dir, pattern=".*green.*.tif$", full=T)
modis_red_files <- dir(modis_data_dir, pattern=".*red.*.tif$", full=T)
modis_nir_files <- dir(modis_data_dir, pattern=".*nir.*.tif$", full=T)
modis_snow_files <- dir(modis_data_dir, pattern=".*snow.*.tif$", full=T)

modis_dates <- as.Date(gsub(pattern=".*A([0-9]{7}).*", "\\1", basename(modis_blue_files)), format="%Y%j")

modis_blue_files <- modis_blue_files[order(modis_dates)]
modis_green_files <- modis_green_files[order(modis_dates)]
modis_red_files <- modis_red_files[order(modis_dates)]
modis_nir_files <- modis_nir_files[order(modis_dates)]
modis_snow_files <- modis_snow_files[order(modis_dates)]

modis_dates <- sort(modis_dates)
modis_doys <- as.integer(strftime(modis_dates, format="%j"))
modis_years <- as.integer(strftime(modis_dates, format="%Y"))

Y_modis_snow <- GetValuesGDAL(modis_snow_files, start_row=1, n=nrow(tmp_r))
Y_modis_blue <- GetValuesGDAL(modis_blue_files, start_row=1, n=nrow(tmp_r))
Y_modis_blue[Y_modis_snow != 0] <- NA
Y_modis_green <- GetValuesGDAL(modis_green_files, start_row=1, n=nrow(tmp_r))
Y_modis_green[Y_modis_snow != 0] <- NA
Y_modis_red <- GetValuesGDAL(modis_red_files, start_row=1, n=nrow(tmp_r))
Y_modis_red[Y_modis_snow != 0] <- NA
Y_modis_nir <- GetValuesGDAL(modis_nir_files, start_row=1, n=nrow(tmp_r))
Y_modis_nir[Y_modis_snow != 0] <- NA

# get the RMSE data
load("~/Desktop/KF_fusion_data_new/multispectral_rmse.Rdata")
rmse_blue_v <- values(rmse_blue)
rmse_green_v <- values(rmse_green)
rmse_red_v <- values(rmse_red)
rmse_nir_v <- values(rmse_nir)

# get the CDL data
year_of_interest <- 2008
cdl_data_dir <- "~/Desktop/KF_fusion_data_new/cdl_new"
cdl_files <- dir(cdl_data_dir, pattern=".*tif$", full=T)
cdl_years <- as.integer(gsub(pattern="CDL_([0-9]{4}).*", "\\1", basename(cdl_files)))
cdl_r <- raster(cdl_files[which(cdl_years == year_of_interest)])
cdl_v <- values(cdl_r)
cdl_types <- unique(cdl_v)
cdl_names <- list("1"="corn", "4"="sorghum", "5"="soybeans", "13"="pop/orn corn", "24"="winter wheat", "27"="rye", "28"="oats", "36"="alfalfa", "37"="other hay", "61"="fallow", "69"="grapes", "111"="open water", "121"="dev open", "122"="dev low", "123"="dev med", "124"="dev high", "131"="barren", "141"="deciduous forest", "142"="evergreen forest", "143"="mixed forest", "152"="shrubland", "171"="grasslands/herb", "190"="woody wetlands", "195"="herb wetlands", "229"="pumpkins")

# Y <- cbind(cdl_v, rmse_v, Y_landsat_blue, Y_landsat_green, Y_landsat_red, Y_landsat_nir, Y_modis_blue, Y_modis_green, Y_modis_red, Y_modis_nir)
Y <- cbind(cdl_v, rmse_blue_v, rmse_green_v, rmse_red_v, rmse_nir_v, Y_landsat_blue, Y_landsat_green, Y_landsat_red, Y_landsat_nir, Y_modis_blue, Y_modis_green, Y_modis_red, Y_modis_nir)
rm(Y_landsat_blue, Y_landsat_green, Y_landsat_red, Y_landsat_nir, Y_modis_blue, Y_modis_green, Y_modis_red, Y_modis_nir, rmse_blue_v, rmse_green_v, rmse_red_v, rmse_nir_v)

# save(Y, tmp_r, cdl_types, multispectral_cdl_process_cov, landsat_sensor, landsat_dates, modis_dates, file="~/Desktop/KF_fusion_data_new/nebraska_multispec_workspace.Rdata")

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# Get the time-varying process model std dev by CDL land cover type
# if there is only one example of a pixel, then we adopt a constant, default value
# determine a "default" time-varying del_r,del_g,del_b,del_nir covariance matrix over ALL CDL types
set.seed(42)
max_to_sample <- 1e5
sub_sample <- sample(1:dim(Y_landsat_blue)[1], min(max_to_sample, dim(Y_landsat_blue)[1]))

cdl_process_diff_blue <- parApply(cl, Y_landsat_blue[sub_sample,], 1, GetLandsatMultispecProcessDiff, dates=landsat_dates)
cdl_process_diff_blue <- do.call(rbind, lapply(seq(1, dim(cdl_process_diff_blue)[1], by=364), function(i) t(cdl_process_diff_blue)[,i:(i+364-1)]))

cdl_process_diff_green <- parApply(cl, Y_landsat_green[sub_sample,], 1, GetLandsatMultispecProcessDiff, dates=landsat_dates)
cdl_process_diff_green <- do.call(rbind, lapply(seq(1, dim(cdl_process_diff_green)[1], by=364), function(i) t(cdl_process_diff_green)[,i:(i+364-1)]))

cdl_process_diff_red <- parApply(cl, Y_landsat_red[sub_sample,], 1, GetLandsatMultispecProcessDiff, dates=landsat_dates)
cdl_process_diff_red <- do.call(rbind, lapply(seq(1, dim(cdl_process_diff_red)[1], by=364), function(i) t(cdl_process_diff_red)[,i:(i+364-1)]))

cdl_process_diff_nir <- parApply(cl, Y_landsat_nir[sub_sample,], 1, GetLandsatMultispecProcessDiff, dates=landsat_dates)
cdl_process_diff_nir <- do.call(rbind, lapply(seq(1, dim(cdl_process_diff_nir)[1], by=364), function(i) t(cdl_process_diff_nir)[,i:(i+364-1)]))

default_cov <- array(NA, dim=c(4, 4, dim(cdl_process_diff_nir)[2]))
for(i in 1:dim(cdl_process_diff_nir)[2]){
  # print(paste("Doing", i))
  default_cov[,,i] <- TukeyRestrictedCOV(cdl_process_diff_blue[,i], cdl_process_diff_green[,i], cdl_process_diff_red[,i], cdl_process_diff_nir[,i])
}

set.seed(42)
max_to_sample <- 1e5
multispectral_cdl_process_cov <- list()
n <- 1
tb <- table(cdl_v)
for(cdl_id in cdl_types){
  print(paste("Doing CDL ID:", cdl_id))
  if(tb[names(tb) == cdl_id] > 1){
    sub_sample <- sample(1:dim(Y_landsat_blue[cdl_v == cdl_id,])[1], min(max_to_sample, dim(Y_landsat_blue[cdl_v == cdl_id,])[1]))

    cdl_process_diff_blue <- parApply(cl, Y_landsat_blue[cdl_v == cdl_id, ][sub_sample,], 1, GetLandsatMultispecProcessDiff, dates=landsat_dates)
    cdl_process_diff_blue <- do.call(rbind, lapply(seq(1, dim(cdl_process_diff_blue)[1], by=364), function(i) t(cdl_process_diff_blue)[,i:(i+364-1)]))

    cdl_process_diff_green <- parApply(cl, Y_landsat_green[cdl_v == cdl_id, ][sub_sample,], 1, GetLandsatMultispecProcessDiff, dates=landsat_dates)
    cdl_process_diff_green <- do.call(rbind, lapply(seq(1, dim(cdl_process_diff_green)[1], by=364), function(i) t(cdl_process_diff_green)[,i:(i+364-1)]))

    cdl_process_diff_red <- parApply(cl, Y_landsat_red[cdl_v == cdl_id, ][sub_sample,], 1, GetLandsatMultispecProcessDiff, dates=landsat_dates)
    cdl_process_diff_red <- do.call(rbind, lapply(seq(1, dim(cdl_process_diff_red)[1], by=364), function(i) t(cdl_process_diff_red)[,i:(i+364-1)]))

    cdl_process_diff_nir <- parApply(cl, Y_landsat_nir[cdl_v == cdl_id, ][sub_sample,], 1, GetLandsatMultispecProcessDiff, dates=landsat_dates)
    cdl_process_diff_nir <- do.call(rbind, lapply(seq(1, dim(cdl_process_diff_nir)[1], by=364), function(i) t(cdl_process_diff_nir)[,i:(i+364-1)]))

    this_cov <- array(NA, dim=c(4, 4, dim(cdl_process_diff_nir)[2]))
    for(i in 1:dim(cdl_process_diff_nir)[2]){
      this_cov[,,i] <- TukeyRestrictedCOV(cdl_process_diff_blue[,i], cdl_process_diff_green[,i], cdl_process_diff_red[,i], cdl_process_diff_nir[,i])
    }
  }else{
    this_cov <- default_cov
  }
  multispectral_cdl_process_cov[[n]] <- this_cov
  n <- n + 1
}
save(multispectral_cdl_process_cov, cdl_types, file="~/Desktop/KF_fusion_data_new/multispectral_cdl_process_cov.Rdata")


#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
# playing around with EVI2 error and reflectance
# tmp_error_mat <- outer(seq(0,1,by=0.01), seq(0,1,by=0.01), function(x,y) EVI2_error(x,y,u_red=0.05,u_nir=0.05))
# tmp_evi2_mat <- outer(seq(0,1,by=0.01), seq(0,1,by=0.01), CalcEVI2)
# image(tmp_error_mat)
# contour(seq(0,1,by=0.01), seq(0,1,by=0.01), tmp_evi2_mat, add=T)


#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
# Plot the multispectral covariance for each CDL type
PlotMultiCov <- function(x, lwd=2, dashed_lty=2, scale_factor=1e4, legend=T, ylim=NULL, xlab="", ylab=""){
  mycols <- c(brewer.pal(5, "Blues")[4], brewer.pal(5, "Greens")[4], brewer.pal(5, "Reds")[4], brewer.pal(5, "Purples")[4])
  x <- x / scale_factor
  if(is.null(ylim)) ylim <- c(min(x, na.rm=T), max(x, na.rm=T))
  plot(x[1,1,], type="n", ylim=ylim, xlab=xlab, ylab=ylab)

  points(x[1,1,], type="l", col=mycols[1], lwd=lwd) # blue error variance
  points(x[2,2,], type="l", col=mycols[2], lwd=lwd) # green error variance
  points(x[3,3,], type="l", col=mycols[3], lwd=lwd) # red error variance
  points(x[4,4,], type="l", col=mycols[4], lwd=lwd) # nir error variance

  points(x[2,1,], type="l", col=mycols[1], lwd=lwd) # blue:green
  points(x[2,1,], type="l", col=mycols[2], lwd=lwd, lty=dashed_lty) # blue:green

  points(x[3,1,], type="l", col=mycols[1], lwd=lwd) # blue:red
  points(x[3,1,], type="l", col=mycols[3], lwd=lwd, lty=dashed_lty) # blue:red

  points(x[4,1,], type="l", col=mycols[1], lwd=lwd) # blue:nir
  points(x[4,1,], type="l", col=mycols[4], lwd=lwd, lty=dashed_lty) # blue:nir

  points(x[3,2,], type="l", col=mycols[2], lwd=lwd) # green:red
  points(x[3,2,], type="l", col=mycols[3], lwd=lwd, lty=dashed_lty) # green:red

  points(x[4,2,], type="l", col=mycols[2], lwd=lwd) # green:nir
  points(x[4,2,], type="l", col=mycols[4], lwd=lwd, lty=dashed_lty) # green:nir

  points(x[4,3,], type="l", col=mycols[3], lwd=lwd) # red:nir
  points(x[4,3,], type="l", col=mycols[4], lwd=lwd, lty=dashed_lty) # red:nir

  abline(h=0, lty=3)
  if(legend){
    legend(
      "topleft",
      legend=c(expression("var("~rho[blue]~")"), expression("var("~rho[green]~")"), expression("var("~rho[red]~")"), expression("var("~rho[nir]~")"), expression("cov("~rho[blue]~","~rho[green]~")"), expression("cov("~rho[blue]~","~rho[red]~")"), expression("cov("~rho[blue]~","~rho[nir]~")"), expression("cov("~rho[green]~","~rho[red]~")"), expression("cov("~rho[green]~","~rho[nir]~")"), expression("cov("~rho[red]~","~rho[nir]~")")),
      col=c(mycols[1:4], mycols[1], mycols[1], mycols[1], mycols[2], mycols[2], mycols[3]),
      lty=1,
      lwd=lwd
    )
    legend(
      "topleft",
      legend=rep(NA, 10),
      col=c(NA, NA, NA, NA, mycols[2], mycols[3], mycols[4], mycols[3], mycols[4], mycols[4]),
      lty=dashed_lty,
      lwd=lwd,
      bty="n"
    )

  }
}

layout(matrix(1:25, nrow=5, byrow=T))
par(mar=c(2, 2, 1, 1))
lwd <- 2
ylim <- c(-0.075, 0.3)
for(i in 1:length(multispectral_cdl_process_cov)){
  if(i==1){
    PlotMultiCov(multispectral_cdl_process_cov[[i]], lwd=lwd, dashed_lty=3, legend=T, ylim=ylim)
    title(get(as.character(cdl_types[[i]]), cdl_names))
  }else{
    PlotMultiCov(multispectral_cdl_process_cov[[i]], lwd=lwd, dashed_lty=3, legend=F, ylim=ylim)
    title(get(as.character(cdl_types[[i]]), cdl_names))
  }
}


#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# Maybe too big to store in memory all at once, so we do in chunks
# system.time(trash <- parApply(cl, Y[1:1e4,], 1, FuseLandsatModisMultispec, landsat_dates=landsat_dates, modis_dates=modis_dates, landsat_sensor=landsat_sensor, multispectral_cdl_process_cov=multispectral_cdl_process_cov, cdl_types=cdl_types))
library(tools)
library(raster)
library(rgdal)
library(parallel)
library(RColorBrewer)
library(dlm)
library(reshape2)

fig_out_dir <- "~/Desktop"

source("~/Documents/pixel-forge/KalmanFiltering/KF_functions.R")
cl <- makeCluster(detectCores())
clusterEvalQ(cl, {source("~/Documents/pixel-forge/KalmanFiltering/KF_functions.R"); library(dlm); library(reshape2)})

load("nebraska_multispec_workspace.Rdata")

blank_r <- tmp_r; values(blank_r) <- NA
out_s_blue <- do.call(stack, replicate(length(unique(c(landsat_years, modis_years)))*365, blank_r))
out_s_green <- do.call(stack, replicate(length(unique(c(landsat_years, modis_years)))*365, blank_r))
out_s_red <- do.call(stack, replicate(length(unique(c(landsat_years, modis_years)))*365, blank_r))
out_s_nir <- do.call(stack, replicate(length(unique(c(landsat_years, modis_years)))*365, blank_r))

blue_inds <- 1:1460
green_inds <- 1461:2920
red_inds <- 2921:4380
nir_inds <- 4381:5480

rows_to_do <- 1e3
row_seq <- seq(1, nrow(Y), by=rows_to_do)
# out_dir <- "~/Desktop/KF_fusion_data_new/multispec_output_new"
out_file_blue <- file.path(out_dir, "kf_multispec_nebraska_blue.tif")
out_file_green <- file.path(out_dir, "kf_multispec_nebraska_green.tif")
out_file_red <- file.path(out_dir, "kf_multispec_nebraska_red.tif")
out_file_nir <- file.path(out_dir, "kf_multispec_nebraska_nir.tif")
for(i in row_seq){
	i_end <- min(i + rows_to_do - 1, nrow(Y))
	tmp <- parApply(cl, Y[i:i_end,], 1, FuseLandsatModisMultispec, landsat_dates=landsat_dates, modis_dates=modis_dates, landsat_sensor=landsat_sensor, multispectral_cdl_process_cov=multispectral_cdl_process_cov, cdl_types=cdl_types)
	out_s_blue[i:i_end] <- c(t(tmp[blue_inds,]))
	out_s_green[i:i_end] <- c(t(tmp[green_inds,]))
	out_s_red[i:i_end] <- c(t(tmp[red_inds,]))
	out_s_nir[i:i_end] <- c(t(tmp[nir_inds,]))
}

num_years <- length(unique(c(landsat_years, modis_years)))
plot_dates <- sort(as.Date(paste(rep(sort(unique(c(landsat_years, modis_years))), each=365), rep(1:365, num_years), sep="-"), format="%Y-%j"))
for(i in 1:nlayers(out_s_blue)){
	out_file <- file.path(fig_out_dir, paste(as.character(strftime(plot_dates[i], format="%Y%j")), "_kf_multispec.jpg", sep=""))
	# jpeg(file=out_file, width=1625, height=round(581+(581/2)), quality=75)
	jpeg(file=out_file, width=1214, height=772, quality=75)
	plot_s <- stack(raster(out_s_nir, i), raster(out_s_red, i), raster(out_s_green, i))
	plotRGB(plot_s, stretch="lin")
	dev.off()
}
