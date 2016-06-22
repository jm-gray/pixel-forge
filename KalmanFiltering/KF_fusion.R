library(raster)
library(parallel)
library(rgdal)
library(dlm)
library(tools)
library(caTools)
library(RColorBrewer)
source("https://raw.github.com/jm-gray/pixel-forge/master/KalmanFiltering/KF_functions.R")

# #-------------------------------------------------------------------
# CalcEVI2 <- function(in_file, scale_factor=1e-4){
#   s <- stack(in_file)
#   nir_v <- values(raster(in_file, 4))
#   red_v <- values(raster(in_file, 3))
#   evi2 <- 2.5 * (((nir_v * scale_factor) - (red_v * scale_factor)) / ((nir_v * scale_factor) + (2.4 * (red_v * scale_factor)) + 1))
#   return(evi2)
# }
#
# #-------------------------------------------------------------------
# PlotEVI <- function(r, breaks=c(-1e9, seq(0, 0.75, len=253), 1e9)){
#   par(mar=rep(0,4), oma=rep(0,4), bg="black")
#   myRamp <- colorRampPalette(brewer.pal(11, "Spectral"))
#   plot(r, breaks=breaks, col=myRamp(255), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
# }
#
# #-------------------------------------------------------------------
# PlotFusion <- function(r1, r2, kf, label, evi_breaks=c(-1e9, seq(0, 0.75, len=254), 1e9)){
#   nf <- layout(matrix(c(1,2,3,3), nrow=2, byrow=T), heights=c(0.5, 1))
#   par(mar=rep(0,4), oma=rep(0,4), bg="black")
#   myRamp <- colorRampPalette(brewer.pal(11, "Spectral"))
#   image(r1, breaks=evi_breaks, col=myRamp(255), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
#   image(r2, breaks=evi_breaks, col=myRamp(255), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
#   image(kf, breaks=evi_breaks, col=myRamp(255), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
#   text(par()$usr[1], par()$usr[3], label=label, adj=c(-0.1, -0.1), col='black', cex=2)
# }
#
# #-------------------------------------------------------------------
# ApplyKFSmooth <- function(x, dlm, prior_rates, prior_vars, prior_means, use_cov=FALSE){
#   cluster_id <- x[1,1]
#   if(use_cov){
#     signal_cov <- x[1,2]
#   }else{
#     signal_cov <- 0
#   }
#   y <- x[-1,]
#   tmp_dlm <- dlm
#   X(tmp_dlm)[,1] <- prior_rates[cluster_id, ]
#   X(tmp_dlm)[,2] <- prior_vars[cluster_id, ]
#   m0(tmp_dlm) <- prior_means[cluster_id, 1]
#   C0(tmp_dlm) <- prior_vars[cluster_id, 1]
#   smooth <- dlmSmooth(y, tmp_dlm)
#   return(smooth)
# }
#
# #-------------------------------------------------------------------
# ExtractSmoothMeans <- function(x) dropFirst(x$s)
# ExtractSmoothSE <- function(x) sqrt(unlist(dropFirst(dlmSvd2var(x$U.S, x$D.S))))

#-------------------------------------------------------------------
# get the input data files
landsat_data_dir <- "/projectnb/modislc/users/joshgray/DL_Landsat"
# landsat_data_dir <- "~/Desktop/KF_data/landsat"
landsat_in_files <- dir(path=landsat_data_dir, pattern="LFstack_fuse.tif", full=T, rec=T)
landsat_dates <- do.call("c", lapply(basename(landsat_in_files), function(x) as.Date(paste(substr(x, 10, 13), substr(x, 14, 16), sep="-"), format="%Y-%j")))
landsat_in_files <- landsat_in_files[order(landsat_dates)]
landsat_dates <- sort(landsat_dates)
# subset files and dates to year 2008
landsat_in_files <- landsat_in_files[landsat_dates >= as.Date("2008-1-1") & landsat_dates <= as.Date("2008-12-31")]
landsat_dates <- landsat_dates[landsat_dates >= as.Date("2008-1-1") & landsat_dates <= as.Date("2008-12-31")]

# get the MODIS data
modis_data_dir <- "/projectnb/modislc/users/joshgray/DL_Landsat/MODIS2008"
# modis_data_dir <- "~/Desktop/KF_data/modis"
modis_in_files <- dir(path=modis_data_dir, pattern="utm_30m.tif", full=T, rec=T)
modis_dates <-  as.Date(substr(basename(modis_in_files), 10, 16), format="%Y%j")

#-------------------------------------------------------------------
# construct a data.frame of weekly mean EVI2 pixel values from Landsat
tmp_r <- raster(landsat_in_files[1])
week_bins <- seq.Date(as.Date("2008-1-1"), as.Date("2008-12-31"), by="week")
pred_dates <- as.Date(round((as.integer(week_bins[-length(week_bins)]) + as.integer(week_bins[-1])) / 2), origin="1970-1-1")
landsat_evi2 <- matrix(NA, ncol=length(pred_dates), nrow=ncell(tmp_r))
for(i in 1:length(pred_dates)){
  print(paste("Working on date:", pred_dates[i]))
  period_files <- landsat_in_files[landsat_dates >= week_bins[i] & landsat_dates < week_bins[i + 1]]
  if(length(period_files) > 1){
    for(period_file in period_files){
      tmp_evi2 <- CalcEVI2(period_file)
      if(period_file == period_files[1]){
        evi2_df <- tmp_evi2
      }else{
        evi2_df <- rbind(evi2_df, tmp_evi2)
      }
    }
    evi2 <- colMeans(evi2_df, na.rm=T)
    landsat_evi2[,i] <- evi2
  }else if(length(period_files) == 1){
    evi2 <- CalcEVI2(period_files)
    landsat_evi2[,i] <- evi2
  }
}
landsat_evi2[is.na(landsat_evi2)] <- NA # because NaN values crop up!

#-------------------------------------------------------------------
# construct a data.frame of weekly mean EVI2 pixel values from MODIS
tmp_r <- raster(modis_in_files[1])
week_bins <- seq.Date(as.Date("2008-1-1"), as.Date("2008-12-31"), by="week")
# pred_dates <- as.Date(round((as.integer(week_bins[-length(week_bins)]) + as.integer(week_bins[-1])) / 2), origin="1970-1-1")
modis_evi2 <- matrix(NA, ncol=length(pred_dates), nrow=ncell(tmp_r))
for(i in 1:length(pred_dates)){
  print(paste("Working on date:", pred_dates[i]))
  period_files <- modis_in_files[modis_dates >= week_bins[i] & modis_dates < week_bins[i + 1]]
  if(length(period_files) > 1){
    for(period_file in period_files){
      tmp_evi2 <- CalcEVI2(period_file)
      if(period_file == period_files[1]){
        evi2_df <- tmp_evi2
      }else{
        evi2_df <- rbind(evi2_df, tmp_evi2)
      }
    }
    evi2 <- colMeans(evi2_df, na.rm=T)
    modis_evi2[,i] <- evi2
  }else if(length(period_files) == 1){
    evi2 <- CalcEVI2(period_files)
    modis_evi2[,i] <- evi2
  }
}
modis_evi2[is.na(modis_evi2)] <- NA # because NaN values crop up!

#-------------------------------------------------------------------
#-------------------------------------------------------------------
# construct a non time-varying dlm, a "dumb" fuser
sensor1_sd <- 0.04 # landsat error
sensor2_sd <- 0.06 # modis error
evi2_fusion_dlm <- MakeMultiDLM(1, 2)
diag(V(evi2_fusion_dlm)) <- c(sensor1_sd^2, sensor2_sd^2)
W(evi2_fusion_dlm) <- 0.05^2
m0(evi2_fusion_dlm) <- 0
C0(evi2_fusion_dlm) <- 1

# construct a big array of observations for KF application
# designed for the KF to apply over rows of Y, or the x-dimension (i.e. x,y,z)
Y <- array(NA, dim=c(dim(landsat_evi2)[1], dim(landsat_evi2)[2], 2))
Y[,,1] <- landsat_evi2
Y[,,2] <- modis_evi2

# apply the filter over all the data using a cluster
# cl <- makeCluster(16)
# cl <- makeCluster(4)
cl <- makeCluster(3)
clusterEvalQ(cl, {library(dlm)})
# system.time(filt_results_dumb <- parApply(cl, Y[1:100000,,], 1, dlmFilter, evi2_fusion_dlm)) # 172 seconds!
system.time(filt_results_dumb <- parApply(cl, Y, 1, dlmFilter, evi2_fusion_dlm)) # 172 seconds!
system.time(filtering_means_dumb <- do.call(rbind, parLapply(cl, filt_results_dumb, function(x) dropFirst(x$m)))) # extract the means of the filtered state distributions
# system.time(filtering_se <- do.call(rbind, parLapply(cl, filt_results, function(x) sqrt(unlist(dropFirst(dlmSvd2var(x$U.C, x$D.C))))))) # extract the means of the filtered


# out_dir <- "/projectnb/modislc/users/joshgray/DL_Landsat/Fusion1"
out_dir <- "~/Desktop/FusionResults"
out_prefix <- "DumbFusion"
for(i in 1:dim(filtering_means)[2]){
  print(paste("Plotting week", i))
  out_file <- file.path(out_dir, paste(out_prefix, "_", formatC(i, width=3, flag="0"), ".jpg", sep=""))
  jpeg(file=out_file, width=1026, height=round(671+(671/2)), quality=75)
  kf <- se <- r1 <- r2 <- tmp_r
  values(kf) <- filtering_means[,i]
  values(se) <- filtering_se[,i]
  values(r1) <- landsat_evi2[,i]
  values(r2) <- modis_evi2[,i]
  PlotFusion(r1, r2, kf, label=paste("Week:", i))
  dev.off()
}
# make into a movie with ffmpeg:
# module load ffmpeg/git-2013-12-18-4a2570f
# ffmpeg -framerate 2 -i Fusion1_%03d.jpg -pix_fmt yuv420p -vf scale="720:trunc(ow/a/2)*2" Fusion1.mp4


#-------------------------------------------------------------------
#-------------------------------------------------------------------
# construct a time-varying dlm based on the mean landsat evolution
numyears <- 1
t <- 1:52
mean_landsat <- colMeans(landsat_evi2, na.rm=T)
var_landsat <- apply(landsat_evi2, 2, var, na.rm=T)
smooth_var_landsat <- predict(smooth.spline(t[!is.na(var_landsat)], var_landsat[!is.na(var_landsat)]),1:52)[[2]]
smooth_mean_landsat <- predict(smooth.spline(t[!is.na(mean_landsat)],mean_landsat[!is.na(mean_landsat)]),1:52)[[2]]
prior_rate <- smooth_mean_landsat[-1] / smooth_mean_landsat[-length(smooth_mean_landsat)]
prior_rate <- c(prior_rate, prior_rate[length(prior_rate)])
prior_var <- smooth_var_landsat

# define the model
tv_evi2_fusion_dlm <- MakeMultiDLM(1, 2)
sensor1_sd <- 0.04 # landsat error
sensor2_sd <- 0.06 # modis error
JGG(tv_evi2_fusion_dlm) <- 1 # store the time-varying values of GG in col 1 of X
JW(tv_evi2_fusion_dlm) <- 2 # store the time-varying values of W in col 2 of X
X(tv_evi2_fusion_dlm) <- matrix(c(rep(prior_rate, numyears), rep(prior_var, numyears)), ncol=2)
diag(V(tv_evi2_fusion_dlm)) <- c(sensor1_sd^2, sensor2_sd^2)
m0(tv_evi2_fusion_dlm) <- smooth_mean_landsat[1]
C0(tv_evi2_fusion_dlm) <- as.numeric(var_landsat[1])

# test on specific pixels
# pixel <- 310313
# # pixel <- 187235
# # pixel <- 279313
# y <- Y[pixel,,]
# filt <- dlmFilter(y, tv_evi2_fusion_dlm)
# filt_m <- dropFirst(filt$m)
# filt_se <- sqrt(unlist(dropFirst(dlmSvd2var(filt$U.C, filt$D.C))))
# PlotForecast(filt_m, filt_se, signal=split(y, col(y)))

# apply over entire image
system.time(filt_results <- parApply(cl, Y, 1, dlmFilter, tv_evi2_fusion_dlm)) # 172 seconds!
system.time(filtering_means <- do.call(rbind, parLapply(cl, filt_results, function(x) dropFirst(x$m)))) # extract the means of the filtered state distributions
# system.time(filtering_se <- do.call(rbind, parLapply(cl, filt_results, function(x) sqrt(unlist(dropFirst(dlmSvd2var(x$U.C, x$D.C))))))) #

out_dir <- "/projectnb/modislc/users/joshgray/DL_Landsat/Fusion2"
out_prefix <- "Fusion2"
for(i in 1:dim(filtering_means)[2]){
  print(paste("Plotting week", i))
  out_file <- file.path(out_dir, paste(out_prefix, "_", formatC(i, width=3, flag="0"), ".jpg", sep=""))
  jpeg(file=out_file, width=1026, height=round(671+(671/2)), quality=75)
  kf <- se <- r1 <- r2 <- tmp_r
  values(kf) <- filtering_means[,i]
  # values(se) <- filtering_se[,i]
  values(r1) <- landsat_evi2[,i]
  values(r2) <- modis_evi2[,i]
  PlotFusion(r1, r2, kf, label=paste("Week:", i))
  dev.off()
}
# module load ffmpeg/git-2013-12-18-4a2570f
# ffmpeg -framerate 2 -i Fusion2_%03d.jpg -pix_fmt yuv420p -vf scale="720:trunc(ow/a/2)*2" Fusion2.mp4

#-------------------------------------------------------------------
#-------------------------------------------------------------------
# construct a time-varying dlm based on the LC-specific landsat evolution
# cdl <- raster(dir("/projectnb/modislc/users/joshgray/DL_Landsat/CDL_sub", pattern="2008", full=T))
# cdl_v <- values(cdl)
# cdl_tab <- sort(table(cdl_v), dec=T)
# # lcs <- as.integer(names(cdl_tab[cumsum(cdl_tab/length(cdl_v)) <= 0.99])) # LC's accounting for 99% of total
# lcs <- as.integer(names(cdl_tab))
# lc_means <- do.call(rbind, lapply(lcs, function(x) colMeans(landsat_evi2[cdl_v == x,], na.rm=T)))
# lc_vars <- do.call(rbind, lapply(lcs, function(x) apply(landsat_evi2[cdl_v == x,], 2, var, na.rm=T)))
# t <- 1:52
# lc_smooth_means <- t(apply(lc_means, 1, function(x) predict(smooth.spline(t[!is.na(x)], x[!is.na(x)]), t)[[2]]))
#
# # 1: corn
# # 5: soybean
# # 36: alfalfa
# # 37: other hay/non-alfalfa
# # 121: dev/open
# # 124: dev/high intensity
# # 141: deciduous forest
# # 171: grassland/herb
# # 176: grass/pasture
# # 190: woody-wetlands
# plot(1:52, type="n", ylim=c(-0.1, 1))
# cols <- apply(lc_means, 1, function(x) {ran_col=colors()[sample(1:250)]; points(1:52, x, col=ran_col, type="l", lwd=2); return(ran_col)})
# legend("topleft", legend=lcs, col=cols, lty=1, pch=NA, lwd=2)

# fit 6 kmean clusters to the min/max/diff(range)/var of each time series
cluster_data <- t(apply(landsat_evi2, 1, function(x) {c(min(x, na.rm=T), max(x, na.rm=T), diff(range(x, na.rm=T)), var(x, na.rm=T))}))
cluster_data <- data.frame(1:dim(cluster_data)[1], cluster_data)
cluster_data_na <- na.omit(cluster_data)
cluster_data_na[,-1] <- scale(cluster_data_na[,-1])
cluster_kmeans <- kmeans(cluster_data_na[,-1], 6)
merge_cluster_df <- data.frame(cluster_data_na[,1], cluster_kmeans$cluster)
names(merge_cluster_df) <- c("pixel", "cluster")
merged_df <- merge(merge_cluster_df, data.frame(pixel=cluster_data[,1]), by="pixel", all=T)
merged_df$cluster[is.na(merged_df$cluster)] <- 0
merged_df$cluster <- merged_df$cluster + 1
cluster_vals <- merged_df$cluster[order(merged_df$pixel)]
cluster_Y <- cbind(cluster_vals, cluster_vals)

# compute cluster means and variances to get model priors for each cluster
t <- 1:52
cluster_means <- do.call(rbind, lapply(sort(unique(cluster_vals)), function(x) colMeans(landsat_evi2[cluster_vals == x,], na.rm=T)))
cluster_smooth_means <- t(apply(cluster_means, 1, function(x) predict(smooth.spline(t[!is.na(x)], x[!is.na(x)]), t)[[2]]))
cluster_vars <- do.call(rbind, lapply(sort(unique(cluster_vals)), function(x) apply(landsat_evi2[cluster_vals == x,], 2, var, na.rm=T)))
cluster_vars <- t(apply(cluster_vars, 1, function(x) {x[is.na(x)] <- mean(x, na.rm=T); return(x)})) # get rid of missing cluster variances

prior_rates <- cluster_smooth_means[,-1] / cluster_smooth_means[,-dim(cluster_smooth_means)[2]]
prior_rates <- cbind(prior_rates, prior_rates[,dim(prior_rates)[2]])
prior_rates[1,] <- 1 # the first were unclustered pixels (missing data), so they get no prior
prior_vars <- cluster_vars
prior_vars[1, ] <- mean(prior_vars[1, ])

# append the clusterval information as the first element of the evi2 time series
Y <- array(NA, dim=c(dim(landsat_evi2)[1], dim(landsat_evi2)[2] + 1, 2))
landsat_with_cluster_val <- cbind(cluster_vals, landsat_evi2)
modis_with_cluster_val <- cbind(cluster_vals, modis_evi2)
Y[,,1] <- landsat_with_cluster_val
Y[,,2] <- modis_with_cluster_val

tv_cluster_evi2_fusion_dlm <- MakeMultiDLM(1, 2)
sensor1_sd <- 0.04 # landsat error
sensor2_sd <- 0.06 # modis error
JGG(tv_cluster_evi2_fusion_dlm) <- 1 # store the time-varying values of GG in col 1 of X
JW(tv_cluster_evi2_fusion_dlm) <- 2 # store the time-varying values of W in col 2 of X
X(tv_cluster_evi2_fusion_dlm) <- matrix(c(prior_rates[1,], prior_vars[1,]), ncol=2)
diag(V(tv_cluster_evi2_fusion_dlm)) <- c(sensor1_sd^2, sensor2_sd^2)
m0(tv_cluster_evi2_fusion_dlm) <- 0
C0(tv_cluster_evi2_fusion_dlm) <- 1

# system.time(filt_results <- parApply(cl, Y, 1, ApplyKF, tv_cluster_evi2_fusion_dlm, prior_rates, prior_vars, cluster_smooth_means))
# system.time(smooth_results <- parApply(cl, Y[1:100,,], 1, ApplyKFSmooth, tv_cluster_evi2_fusion_dlm, prior_rates, prior_vars, cluster_smooth_means))
smooth_results <- parApply(cl, Y, 1, ApplyKFSmooth, tv_cluster_evi2_fusion_dlm, prior_rates, prior_vars, cluster_smooth_means)
smooth_means <- do.call(rbind, parLapply(cl, smooth_results, ExtractSmoothMeans))
save(smooth_results, smooth_means, file="/projectnb/modislc/users/joshgray/KFSmoothResults.Rdata")

#-------------------------------------------------------------------
#-------------------------------------------------------------------
# Do it with a covariance mask
system.time(mod_land_cov <- parApply(cl, Y, 1, function(x) cov(x[,1], x[,2], use="na.or.complete")))
# NOTE: cannot put a value in the off-diagnoal entries of V that is larger than the largest
# value on the diagonal, so...how to encode MODIS-Landsat covariance?



#-------------------------------------------------------------------
# Well, this is curious...
y <- Y[1,,]
sensor1_sd <- 0.04 # landsat error
sensor2_sd <- 0.06 # modis error
evi2_fusion_dlm <- MakeMultiDLM(1, 2)
diag(V(evi2_fusion_dlm)) <- c(sensor1_sd^2, sensor2_sd^2)
W(evi2_fusion_dlm) <- 0.05^2
m0(evi2_fusion_dlm) <- 0
C0(evi2_fusion_dlm) <- 1
V(evi2_fusion_dlm) <- cov(Y[1,,],use="na")
W(evi2_fusion_dlm) <- 0.003
filt <- dlmFilter(y, evi2_fusion_dlm)
filt_m <- dropFirst(filt$m)
filt_se <- sqrt(unlist(dropFirst(dlmSvd2var(filt$U.C, filt$D.C))))
PlotForecast(filt_m, filt_se, split(y, col(y)))

# doesn't appear to happen for random walk data...
y1 <- cumsum(rnorm(100))
y2 <- y1 + rnorm(100, sd=1.5)
y <- cbind(y1, y2)
tmp_dlm <- MakeMultiDLM(1, 2)
V(tmp_dlm) <- cov(cbind(y))
# W(tmp_dlm) <- var(y1)
W(tmp_dlm) <- 1
m0(tmp_dlm) <- 0
C0(tmp_dlm) <- 1e3
filt <- dlmFilter(y, tmp_dlm)
filt_m <- dropFirst(filt$m)
filt_se <- sqrt(unlist(dropFirst(dlmSvd2var(filt$U.C, filt$D.C))))
PlotForecast(filt_m, filt_se, split(y, col(y)))
