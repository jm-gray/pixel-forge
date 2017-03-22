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
Y_landsat_blue <- t(parApply(cl, Y_landsat_blue, 1, QuantileMinReplacement))

Y_landsat_green <- GetValuesGDAL(landsat_green_files, start_row=1, n=nrow(tmp_r))
Y_landsat_green[Y_landsat_green > 1e4] <- NA
Y_landsat_green <- t(parApply(cl, Y_landsat_green, 1, QuantileMinReplacement))

Y_landsat_red <- GetValuesGDAL(landsat_red_files, start_row=1, n=nrow(tmp_r))
Y_landsat_red[Y_landsat_red > 1e4] <- NA
Y_landsat_red <- t(parApply(cl, Y_landsat_red, 1, QuantileMinReplacement))

Y_landsat_nir <- GetValuesGDAL(landsat_nir_files, start_row=1, n=nrow(tmp_r))
Y_landsat_nir[Y_landsat_nir > 1e4] <- NA
Y_landsat_nir <- t(parApply(cl, Y_landsat_nir, 1, QuantileMinReplacement))

# i <- 30
# s <- stack(landsat_blue_files[i], landsat_green_files[i], landsat_red_files[i], landsat_nir_files[i])
# plotRGB(s, 4, 3, 2, stretch="lin")

modis_data_dir <- "~/Desktop/KF_fusion_data_new/modis_new"
modis_blue_files <- dir(modis_data_dir, pattern=".*blue.*.tif$", full=T)
modis_green_files <- dir(modis_data_dir, pattern=".*green.*.tif$", full=T)
modis_red_files <- dir(modis_data_dir, pattern=".*red.*.tif$", full=T)
modis_nir_files <- dir(modis_data_dir, pattern=".*nir.*.tif$", full=T)

modis_dates <- as.Date(gsub(pattern=".*A([0-9]{7}).*", "\\1", basename(modis_blue_files)), format="%Y%j")

modis_blue_files <- modis_blue_files[order(modis_dates)]
modis_green_files <- modis_green_files[order(modis_dates)]
modis_red_files <- modis_red_files[order(modis_dates)]
modis_nir_files <- modis_nir_files[order(modis_dates)]

modis_dates <- sort(modis_dates)
modis_doys <- as.integer(strftime(modis_dates, format="%j"))
modis_years <- as.integer(strftime(modis_dates, format="%Y"))

# get the RMSE data
rmse_sub <- raster("~/Desktop/KF_fusion_data_new/rmse_sub.tif")
rmse_v <- values(rmse_sub)

# get the CDL data
year_of_interest <- 2008
cdl_data_dir <- "~/Desktop/KF_fusion_data_new/cdl_new"
cdl_files <- dir(cdl_data_dir, pattern=".*tif$", full=T)
cdl_years <- as.integer(gsub(pattern="CDL_([0-9]{4}).*", "\\1", basename(cdl_files)))
cdl_r <- raster(cdl_files[which(cdl_years == year_of_interest)])
cdl_v <- values(cdl_r)
cdl_types <- unique(cdl_v)
cdl_names <- list("1"="corn", "4"="sorghum", "5"="soybeans", "13"="pop/orn corn", "24"="winter wheat", "27"="rye", "28"="oats", "36"="alfalfa", "37"="other hay", "61"="fallow", "69"="grapes", "111"="open water", "121"="dev open", "122"="dev low", "123"="dev med", "124"="dev high", "131"="barren", "141"="deciduous forest", "142"="evergreen forest", "143"="mixed forest", "152"="shrubland", "171"="grasslands/herb", "190"="woody wetlands", "195"="herb wetlands", "229"="pumpkins")


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
# default_cov <- array(rep(c(500 * diag(4)), 364), dim=c(4, 4, 364))
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
      # print(paste("Doing", i))
      this_cov[,,i] <- TukeyRestrictedCOV(cdl_process_diff_blue[,i], cdl_process_diff_green[,i], cdl_process_diff_red[,i], cdl_process_diff_nir[,i])
    }

    # mycols <- brewer.pal(10, "Set3")
    # plot(this_cov[1,1,], type="l", col=mycols[1], ylim=c(-150, 2e3))
    # points(this_cov[2,1,], type="l", col=mycols[2])
    # points(this_cov[3,1,], type="l", col=mycols[3])
    # points(this_cov[4,1,], type="l", col=mycols[4])
    # points(this_cov[2,2,], type="l", col=mycols[5])
    # points(this_cov[3,2,], type="l", col=mycols[6])
    # points(this_cov[4,2,], type="l", col=mycols[7])
    # points(this_cov[3,3,], type="l", col=mycols[8])
    # points(this_cov[4,3,], type="l", col=mycols[9])
    # points(this_cov[4,4,], type="l", col=mycols[10])
    # legend("topleft", legend=c("bb", "gb", "rb", "nb", "gg", "rg", "ng", "rr", "nr", "nn"), col=mycols, lty=1)
    # abline(h=0, lty=3)
    # cdl_process_diff_sd <- apply(cdl_process_diff, 2, TukeyRestrictedSD)
    # sd_spline <- predict(smooth.spline(1:364,cdl_process_diff_sd, spar=0.8))$y
    # # the splined values may be negative, so replace with minimum actual SD
    # sd_spline[sd_spline <= 0] <- min(cdl_process_diff_sd, na.rm=T)
  }else{
    this_cov <- default_cov
  }
  multispectral_cdl_process_cov[[n]] <- this_cov
  n <- n + 1
}
save(multispectral_cdl_process_cov, cdl_types, file="~/Desktop/KF_fusion_data_new/multispectral_cdl_process_cov.Rdata")


#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
tmp_error_mat <- outer(seq(0,1,by=0.01), seq(0,1,by=0.01), function(x,y) EVI2_error(x,y,u_red=0.05,u_nir=0.05))
tmp_evi2_mat <- outer(seq(0,1,by=0.01), seq(0,1,by=0.01), CalcEVI2)
image(tmp_error_mat)
contour(seq(0,1,by=0.01), seq(0,1,by=0.01), tmp_evi2_mat, add=T)
