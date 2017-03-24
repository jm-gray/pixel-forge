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

Y <- cbind(cdl_v, rmse_v, Y_landsat_blue, Y_landsat_green, Y_landsat_red, Y_landsat_nir, Y_modis_blue, Y_modis_green, Y_modis_red, Y_modis_nir)
rm(Y_landsat_blue, Y_landsat_green, Y_landsat_red, Y_landsat_nir, Y_modis_blue, Y_modis_green, Y_modis_red, Y_modis_nir)


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
tmp_error_mat <- outer(seq(0,1,by=0.01), seq(0,1,by=0.01), function(x,y) EVI2_error(x,y,u_red=0.05,u_nir=0.05))
tmp_evi2_mat <- outer(seq(0,1,by=0.01), seq(0,1,by=0.01), CalcEVI2)
image(tmp_error_mat)
contour(seq(0,1,by=0.01), seq(0,1,by=0.01), tmp_evi2_mat, add=T)


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



#-------------------------------------------------------------------------------
FuseLandsatModisMultispec <- function(x, landsat_dates, modis_dates, landsat_sensor, cdl_tv_sd, cdl_types, scale_factor=1e4, modis_landsat_slope=1, modis_landsat_bias=0, smooth=T, plot=F, ...){
  # extract components of x
  tmp_cdl <- x[1]
  tmp_rmse <- x[2] / scale_factor
  x <- x[c(-1, -2)]

  # extract individual landsat bands
  landsat_inds <- seq(1, (length(landsat_dates) * 4) + 1, by=length(landsat_dates))
  x_landsat_blue <- x[landsat_inds[1]:(landsat_inds[2] - 1)] / scale_factor
  x_landsat_green <- x[landsat_inds[2]:(landsat_inds[3] - 1)] / scale_factor
  x_landsat_red <- x[landsat_inds[3]:(landsat_inds[4] - 1)] / scale_factor
  x_landsat_nir <- x[landsat_inds[4]:(landsat_inds[5] - 1)] / scale_factor

  modis_inds <- seq(landsat_inds[5], landsat_inds[5] + length(modis_dates) * 4, by=length(modis_dates))
  x_modis_blue <- x[modis_inds[1]:(modis_inds[2] - 1)] / scale_factor
  x_modis_green <- x[modis_inds[2]:(modis_inds[3] - 1)] / scale_factor
  x_modis_red <- x[modis_inds[3]:(modis_inds[4] - 1)] / scale_factor
  x_modis_nir <- x[modis_inds[4]:(modis_inds[5] - 1)] / scale_factor

  x_modis <- x[(length(landsat_dates) + 1):(length(landsat_dates) + length(modis_dates))] / scale_factor
  x_landsat_error <- x[(length(landsat_dates) + length(modis_dates) + 1):length(x)]

  # retrieve the proper time varying process error for this land cover type
  tv_sd <- cdl_tv_sd[[which(cdl_types == tmp_cdl)]]$splined / scale_factor
  tv_sd <- c(tv_sd[1], tv_sd) # append the head value b/c sd is 364 long

  # munge to daily series
  num_years <- length(as.integer(sort(unique(c(strftime(landsat_dates, format="%Y"), strftime(modis_dates, format="%Y"))))))

  # do landsat; eliminate leap year day 366
  tmp_landsat <- rep(NA, num_years * 365)
  tmp_landsat_error <- rep(NA, num_years * 365)
  tmp_landsat_doys <- as.integer(strftime(landsat_dates, format="%j"))
  tmp_landsat_years <- as.integer(strftime(landsat_dates, format="%Y"))
  x_landsat <- x_landsat[tmp_landsat_doys <= 365]
  x_landsat_error <- x_landsat_error[tmp_landsat_doys <= 365]
  tmp_landsat_doys <- tmp_landsat_doys[tmp_landsat_doys <= 365]
  tmp_landsat_years <- tmp_landsat_years[tmp_landsat_doys <= 365]
  daily_landsat_inds <- tmp_landsat_doys + 365 * (tmp_landsat_years - min(tmp_landsat_years))
  tmp_landsat[daily_landsat_inds] <- x_landsat
  tmp_landsat_error[daily_landsat_inds] <- x_landsat_error
  daily_landsat_sensor <- rep(NA, length(tmp_landsat))
  daily_landsat_sensor[daily_landsat_inds] <- landsat_sensor

  # do modis; eliminate leap year day 366
  tmp_modis <- rep(NA, num_years * 365)
  tmp_modis_doys <- as.integer(strftime(modis_dates, format="%j"))
  tmp_modis_years <- as.integer(strftime(modis_dates, format="%Y"))
  x_modis <- x_modis[tmp_modis_doys <= 365]
  tmp_modis_doys <- tmp_modis_doys[tmp_modis_doys <= 365]
  tmp_modis_years <- tmp_modis_years[tmp_modis_doys <= 365]
  daily_modis_inds <- tmp_modis_doys + 365 * (tmp_modis_years - min(tmp_modis_years))
  tmp_modis[daily_modis_inds] <- x_modis
  tmp_modis <- tmp_modis - modis_landsat_bias # remove any bias in the MODIS sensor

  # create obs matrix
  y <- cbind(tmp_landsat, tmp_modis)

  # # replace all missing landsat error values with closest not-NA value
  # trash <- sapply(which(is.na(tmp_landsat_error)), function(a, x) x[which.min(abs(x-a))], x=which(!is.na(tmp_landsat_error)))
  # tmp_landsat_error[which(is.na(tmp_landsat_error))] <- tmp_landsat_error[trash]

  tmp_modis_error <- rep(NA, length(tmp_modis))
  tmp_modis_error[!is.na(tmp_modis)] <- tmp_rmse

  # define dlm components
  GG <- matrix(1) # process transition
  W <- matrix(1) # evolution error covariance
  JW <- matrix(1) # time varying evolution error covariance: in col 1 of X
  # FF <- matrix(c(1, 1), nrow=2) # observation matrix
  FF <- matrix(c(1, modis_landsat_slope), nrow=2) # observation matrix
  V <- matrix(c(1, 0, 0, 1), nrow=2) # obs uncertainty
  JV <- matrix(c(2, 0, 0, 3), nrow=2) # time varying obs uncertainty: in cols 2 (landsat) and 3 (modis)
  X <- matrix(1, nrow=length(tmp_modis), ncol=3)
  X[, 1] <- rep(tv_sd, num_years) # evolution error
  X[, 2] <- tmp_landsat_error # landsat observation error
  # X[, 3] <- rep(tmp_rmse, length(tmp_modis)) # modis observation error
  X[, 3] <- tmp_modis_error # modis observation error
  m0 <- 0 # initial state vector
  C0 <- 1 # initial state uncertainty

  # construct the dlm
  tmp_dlm <- dlm(m0=m0, C0=C0, GG=GG, JW=JW, FF=FF, V=V, JV=JV, W=W, X=X)

  # apply the dlm
  if(smooth){
    kf_result <- dlmSmooth(y, tmp_dlm)
    kf_evi <- dropFirst(kf_result$s)
    kf_evi_error <- sqrt(unlist(dropFirst(dlmSvd2var(kf_result$U.S, kf_result$D.S))))
    ret_value <- list(evi=kf_evi, error=kf_evi_error)
  }else{
    kf_result <- dlmFilter(y, tmp_dlm)
    kf_evi <- dropFirst(kf_result$m)
    kf_evi_error <- sqrt(unlist(dropFirst(dlmSvd2var(kf_result$U.C, kf_result$D.C))))
    ret_value <- list(evi=kf_evi, error=kf_evi_error)
  }

  if(plot){
    # PlotForecast(kf_evi, kf_evi_error, split(y, col(y)), ylab="EVI2", xlab="", sigma=1, ...)
    PlotForecast(kf_evi, kf_evi_error, split(y, col(y)), ylab="EVI2", xlab="", ...)
  }

  return(ret_value)
}
