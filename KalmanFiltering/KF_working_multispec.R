#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
# gather the multispectral data
landsat_data_dir <- "~/Desktop/KF_fusion_data_new/landsat_new"
landsat_blue_files <- dir(landsat_data_dir, pattern=".*band1.*.tif$", full=T)
landsat_green_files <- dir(landsat_data_dir, pattern=".*band2.*.tif$", full=T)
landsat_red_files <- dir(landsat_data_dir, pattern=".*band3.*.tif$", full=T)
landsat_nir_files <- dir(landsat_data_dir, pattern=".*band4.*.tif$", full=T)

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

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# Get the time-varying process model std dev by CDL land cover type
# if there is only one example of a pixel, then we adopt a constant, default value

set.seed(42)
max_to_sample <- 1e5
default_cov <- array(rep(c(500 * diag(4)), 364), dim=c(4, 4, 364))
multispectral_cdl_process_cov <- list()
i <- 1
tb <- table(cdl_v)
for(cdl_id in cdl_types){
  print(paste("Doing CDL ID:", cdl_id))
  if(tb[names(tb) == cdl_id] > 1){
    sub_sample <- sample(1:dim(Y_landsat_blue[cdl_v == cdl_id,])[1], min(max_to_sample, dim(Y_landsat_blue[cdl_v == cdl_id,])[1]))

    cdl_process_diff_blue <- parApply(cl, Y_landsat_blue[cdl_v == cdl_id, ][sub_sample,], 1, GetLandsatProcessSD, dates=landsat_dates)
    cdl_process_diff_blue <- do.call(rbind, lapply(seq(1, dim(cdl_process_diff_blue)[1], by=364), function(i) t(cdl_process_diff_blue)[,i:(i+364-1)]))

    cdl_process_diff_green <- parApply(cl, Y_landsat_green[cdl_v == cdl_id, ][sub_sample,], 1, GetLandsatProcessSD, dates=landsat_dates)
    cdl_process_diff_green <- do.call(rbind, lapply(seq(1, dim(cdl_process_diff_green)[1], by=364), function(i) t(cdl_process_diff_green)[,i:(i+364-1)]))

    cdl_process_diff_red <- parApply(cl, Y_landsat_red[cdl_v == cdl_id, ][sub_sample,], 1, GetLandsatProcessSD, dates=landsat_dates)
    cdl_process_diff_red <- do.call(rbind, lapply(seq(1, dim(cdl_process_diff_red)[1], by=364), function(i) t(cdl_process_diff_red)[,i:(i+364-1)]))

    cdl_process_diff_nir <- parApply(cl, Y_landsat_nir[cdl_v == cdl_id, ][sub_sample,], 1, GetLandsatProcessSD, dates=landsat_dates)
    cdl_process_diff_nir <- do.call(rbind, lapply(seq(1, dim(cdl_process_diff_nir)[1], by=364), function(i) t(cdl_process_diff_nir)[,i:(i+364-1)]))

    this_cov <- array(NA, dim=c(4, 4, dim(cdl_process_diff_nir)[2]))
    for(i in 1:dim(cdl_process_diff_nir)[2]){
      print(paste("Doing", i))
      this_cov[,,i] <- TukeyRestrictedCOV(cdl_process_diff_blue[,i], cdl_process_diff_green[,i], cdl_process_diff_red[,i], cdl_process_diff_nir[,i])
    }

    # mycols <- brewer.pal(10, "Set3")
    # plot(this_cov[1,1,], type="l", col=mycols[1], ylim=c(-0.5e3, 15e3))
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
  multispectral_cdl_process_cov[[i]] <- this_cov
  i <- i + 1
}
save(multispectral_cdl_process_cov, cdl_types, file="~/Desktop/KF_fusion_data_new/multispectral_cdl_process_cov.Rdata")


#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
tmp_error_mat <- outer(seq(0,1,by=0.01), seq(0,1,by=0.01), function(x,y) EVI2_error(x,y,u_red=0.05,u_nir=0.05))
tmp_evi2_mat <- outer(seq(0,1,by=0.01), seq(0,1,by=0.01), CalcEVI2)
image(tmp_error_mat)
contour(seq(0,1,by=0.01), seq(0,1,by=0.01), tmp_evi2_mat, add=T)
