library(raster)
library(parallel)
library(rgdal)
library(dlm)
library(tools)
library(caTools)
library(RColorBrewer)

#-------------------------------------------------------------------
ExtractSmoothMeans <- function(x) dropFirst(x$s)
ExtractSmoothSE <- function(x) sqrt(unlist(dropFirst(dlmSvd2var(x$U.S, x$D.S))))

#-------------------------------------------------------------------
GetSentinelDate <- function(x) as.Date(substr(unlist(strsplit(basename(x), "_"))[8], 1, 8), format="%Y%m%d")
GetLandsatDate <- function(x) as.Date(substr(basename(x), 10, 16), format="%Y%j")

#-------------------------------------------------------------------
MakeMultiDLM <- function(num_states=1, sensors=1, time_varying=FALSE){
  # NOTE: time_varying option doesn't work right now, but you would just need to create a dlm
  # of the proper state/sensor configuration, and then copy GG to JGG, FF to JFF, and so on
  # then fill in X column indices in the proper places. That is, the architecture of the KF
  # that results from this function is correct

  m0 <- rep(0, num_states) # initial state vector
  C0 <- diag(1e6, num_states)

  GG <- diag(num_states) # state evolution matrix
  W <- diag(1e6, num_states) # evolution error covariance

  # create observation matrix FF
  if(is.list(sensors)){
    sensors_list <- unlist(sensors)
  }else{
    sensors_list <- rep(1:num_states, sensors)
  }
  Frow <- function(state, num_states) c(rep(0, (state - 1)), 1, rep(0, num_states - state))
  FF <- do.call(rbind, lapply(sensors_list, Frow, num_states=num_states))

  V <- diag(1e6, dim(FF)[1]) # observation error covariance

  if(time_varying){
    JFF <- 0 # column index in X of time varying FF values
    JGG <- 0 # column index in X of time varying GG values
    JW <- 0 # column index in X of time varying W values
    JV <- 0 # column index in X of time varying V values
    X <- 0 # matrix holding time varying parameters in columns
    my_dlm <- dlm(m0=m0, C0=C0, GG=GG, W=W, FF=FF, V=V, JGG=JGG, JFF=JFF, JW=JW, JV=JV, X=X)
    return(my_dlm)
  }else{
    my_dlm <- dlm(m0=m0, C0=C0, GG=GG, W=W, FF=FF, V=V)
    return(my_dlm)
  }
}

#-------------------------------------------------------------------
PlotEVI <- function(r, breaks=c(-1e9, seq(0, 0.75, len=253), 1e9)){
  par(mar=rep(0,4), oma=rep(0,4), bg="black")
  myRamp <- colorRampPalette(brewer.pal(11, "Spectral"))
  plot(r, breaks=breaks, col=myRamp(255), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
}

#-------------------------------------------------------------------
PlotFusion <- function(r1, r2, kf, label, evi_breaks=c(-1e9, seq(0, 0.75, len=254), 1e9)){
  nf <- layout(matrix(c(1,2,3,3), nrow=2, byrow=T), heights=c(0.5, 1))
  par(mar=rep(0,4), oma=rep(0,4), bg="black")
  myRamp <- colorRampPalette(brewer.pal(11, "Spectral"))
  image(r1, breaks=evi_breaks, col=myRamp(255), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
  image(r2, breaks=evi_breaks, col=myRamp(255), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
  image(kf, breaks=evi_breaks, col=myRamp(255), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
  text(par()$usr[1], par()$usr[3], label=label, adj=c(-0.1, -0.1), col='black', cex=2)
}

#-------------------------------------------------------------------
# time-varying version of KF smoother
ApplyKFSmooth_tv <- function(x, dlm, prior_rates, prior_vars, prior_means){
  cluster_id <- x[1,1]
  y <- x[-1,]
  tmp_dlm <- dlm
  X(tmp_dlm)[,1] <- prior_rates[cluster_id, ]
  X(tmp_dlm)[,2] <- prior_vars[cluster_id, ]
  m0(tmp_dlm) <- prior_means[cluster_id, 1]
  C0(tmp_dlm) <- prior_vars[cluster_id, 1]
  smooth <- dlmSmooth(y, tmp_dlm)
  return(smooth)
}

#-------------------------------------------------------------------
# not time-varying version of KF smoother
ApplyKFSmooth <- function(y, dlm){
  smooth <- dlmSmooth(y, dlm)
  return(smooth)
}

#-------------------------------------------------------------------
# construct a data.frame of weekly mean EVI pixel values from Landsat
landsat_evi_dir <- "/projectnb/modislc/users/joshgray/Sentinel/med_bakersfield/LandsatDataClipped/EVI_masked"
landsat_evi_files <- dir(landsat_evi_dir, pattern=".*sr_evi.tif", full=T)
landsat_dates <- GetLandsatDate(landsat_evi_files)
landsat_evi_files <- landsat_evi_files[order(landsat_dates)]
landsat_dates <- landsat_dates[order(landsat_dates)]
# restrict to 2015
landsat_evi_files <- landsat_evi_files[landsat_dates >= "2015-1-1" & landsat_dates <= "2015-12-31"]
landsat_dates <- landsat_dates[landsat_dates >= "2015-1-1" & landsat_dates <= "2015-12-31"]

tmp_r <- raster(landsat_evi_files[1])
week_bins <- seq.Date(as.Date("2015-1-1"), as.Date("2015-12-31"), by="week")
pred_dates <- as.Date(round((as.integer(week_bins[-length(week_bins)]) + as.integer(week_bins[-1])) / 2), origin="1970-1-1")
sample_dates <- pred_dates
landsat_evi <- matrix(NA, ncol=length(pred_dates), nrow=ncell(tmp_r))
for(i in 1:length(pred_dates)){
  print(paste("Working on date:", pred_dates[i]))
  # period_files <- landsat_evi_files[landsat_dates >= week_bins[i] & landsat_dates < week_bins[i + 1]]
  period_files <- landsat_evi_files[landsat_dates > week_bins[i] & landsat_dates <= week_bins[i + 1]]
  if(length(period_files) > 1){
    for(period_file in period_files){
      # tmp_evi2 <- CalcEVI2(period_file)
      tmp_evi <- values(raster(period_file))
      tmp_evi[tmp_evi == 0] <- NA
      if(period_file == period_files[1]){
        evi_df <- tmp_evi
      }else{
        evi_df <- rbind(evi_df, tmp_evi)
      }
    }
    evi <- colMeans(evi_df, na.rm=T)
    landsat_evi[,i] <- evi
  }else if(length(period_files) == 1){
    # evi <- CalcEVI2(period_files)
    evi <- values(raster(period_files))
    evi[evi == 0] <- NA
    landsat_evi[,i] <- evi
  }
}
landsat_evi[is.na(landsat_evi)] <- NA # because NaN values crop up!
landsat_evi <- landsat_evi / 1e4

#-------------------------------------------------------------------
# construct a data.frame of weekly mean EVI pixel values from Landsat
sentinel_evi_dir <- "/projectnb/modislc/users/joshgray/Sentinel/med_bakersfield/SentinelDataClipped/EVI"
sentinel_evi_files <- dir(sentinel_evi_dir, pattern=".*EVI.tif", full=T)
sentinel_evi_files <- sentinel_evi_files[1:2]
sentinel_dates <- as.Date(unlist(lapply(sentinel_evi_files, GetSentinelDate)), origin="1970-1-1")
sentinel_evi_files <- sentinel_evi_files[order(sentinel_dates)]
sentinel_dates <- sentinel_dates[order(sentinel_dates)]

tmp_r <- raster(sentinel_evi_files[1])
week_bins <- seq.Date(as.Date("2015-1-1"), as.Date("2015-12-31"), by="week")
pred_dates <- as.Date(round((as.integer(week_bins[-length(week_bins)]) + as.integer(week_bins[-1])) / 2), origin="1970-1-1")
sample_dates <- pred_dates
sentinel_evi <- matrix(NA, ncol=length(pred_dates), nrow=ncell(tmp_r))
for(i in 1:length(pred_dates)){
  print(paste("Working on date:", pred_dates[i]))
  # period_files <- sentinel_evi_files[sentinel_dates >= week_bins[i] & sentinel_dates < week_bins[i + 1]]
  period_files <- sentinel_evi_files[sentinel_dates > week_bins[i] & sentinel_dates <= week_bins[i + 1]]
  if(length(period_files) > 1){
    for(period_file in period_files){
      # tmp_evi2 <- CalcEVI2(period_file)
      tmp_evi <- values(raster(period_file))
      tmp_evi[tmp_evi == 0] <- NA
      if(period_file == period_files[1]){
        evi_df <- tmp_evi
      }else{
        evi_df <- rbind(evi_df, tmp_evi)
      }
    }
    evi <- colMeans(evi_df, na.rm=T)
    sentinel_evi[,i] <- evi
  }else if(length(period_files) == 1){
    # evi <- CalcEVI2(period_files)
    evi <- values(raster(period_files))
    evi[evi == 0] <- NA
    sentinel_evi[,i] <- evi
  }
}
sentinel_evi[is.na(sentinel_evi)] <- NA # because NaN values crop up!

#-------------------------------------------------------------------
# stack the data together, create a DLM, and go to town
cl <- makeCluster(16)
clusterEvalQ(cl, {library(dlm)})

sensor1_sd <- 0.04 # landsat error
sensor2_sd <- 0.04 # sentinel error

# append the clusterval information as the first element of the evi2 time series
Y <- array(NA, dim=c(dim(landsat_evi)[1], dim(landsat_evi)[2], 2))
Y[,,1] <- landsat_evi
Y[,,2] <- sentinel_evi

my_dlm <- MakeMultiDLM(1, 2)
diag(V(my_dlm)) <- c(sensor1_sd^2, sensor2_sd^2)
diag(W(my_dlm)) <- 1
m0(my_dlm) <- 0
C0(my_dlm) <- 1

system.time(smooth_results <- parApply(cl, Y, 1, ApplyKFSmooth, my_dlm))
smooth_means <- do.call(rbind, parLapply(cl, smooth_results, ExtractSmoothMeans))
# save(smooth_results, smooth_means, file="/projectnb/modislc/users/joshgray/KF_landsatsentinel_smooth_dumb.Rdata")

#-----------------------------------------------
# image fusion animation: smoothing, clustered-mean priors
# load KFSmoothResults.Rdata first, then get a tmp_r, then we're off!
out_dir <- "/projectnb/modislc/users/joshgray/Sentinel/med_bakersfield/KF_smooth_output"
out_prefix <- "KF_smooth_"
for(i in 1:dim(smooth_means)[2]){
  print(paste("Plotting week", i))
  out_file <- file.path(out_dir, paste(out_prefix, "_", formatC(i, width=3, flag="0"), ".jpg", sep=""))
  jpeg(file=out_file, width=1625, height=round(581+(581/2)), quality=75)
  kf <- se <- r1 <- r2 <- tmp_r
  values(kf) <- smooth_means[,i]
  # values(se) <- filtering_se[,i]
  values(r1) <- landsat_evi[,i]
  values(r2) <- sentinel_evi[,i]
  PlotFusion(r1, r2, kf, label=paste("Week:", i))
  dev.off()
}
# module load ffmpeg/git-2013-12-18-4a2570f
# ffmpeg -framerate 2 -i KF_smooth__%03d.jpg -pix_fmt yuv420p -vf scale="720:trunc(ow/a/2)*2" KF_smooth.mp4
