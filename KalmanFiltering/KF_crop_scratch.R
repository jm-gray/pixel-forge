library(rgdal)
library(raster)
library(parallel)
library(tools)
library(pracma)
library(zoo)
library(RColorBrewer)

# read in the fields, reproject to UTM and get buffered extent for subsetting
fields <- shapefile("/projectnb/modislc/users/joshgray/DL_Landsat/field_shapes/allpolys_geo.shp")
# fields <- shapefile("~/Desktop/field_shapes/allpolys_geo.shp")
proj <- "+proj=utm +zone=14 +datum=WGS84 +units=m +no_defs"
fields_utm <- spTransform(fields, CRS(proj))
fields_extent <- extent(fields_utm)
sub_extent <- fields_extent + c(-1000, 1000, -1000, 1000) # 1 km buffer
sub_extent <- sub_extent[c(1, 3, 2, 4)] # gdalwarp wants xmin ymin xmax ymax

# read and join field data

# get all the evi2 files
data_dir <- "/projectnb/modislc/users/joshgray/DL_Landsat"
# data_dir <- "~/Desktop/DL_Landsat"
in_files <- dir(path=data_dir, pattern="evi2.tif", full=T, rec=T)
dates <- do.call("c", lapply(basename(in_files), function(x) as.Date(paste(substr(x, 10, 13), substr(x, 14, 16), sep="-"), format="%Y-%j")))
in_files <- in_files[order(dates)]
dates <- sort(dates)

# create a stack
s <- stack(in_files)

# get values
system.time(v <- values(s))

# get the CDL Stack
cdl_s <- stack(dir("/projectnb/modislc/users/joshgray/DL_Landsat/CDL_sub", full=T))
# cdl_s <- stack(dir("~/Desktop/CDL_sub", full=T))
cdl_v <- values(cdl_s)

# rasterize the fields layer
field_ids <- as.numeric(unlist(lapply(fields_utm$pnames,function(x) unlist(strsplit(x,split="_"))[1])))
fields_utm$pnames <- field_ids
fields_r <- rasterize(fields_utm, raster(s, 1), field="pnames")

#-------------------------------------------
# Fit a DLM for 2008 corn
# get 2008 corn data only
v_2008 <- v[, strftime(dates,format="%Y") == "2008"]
dates_2008 <- dates[strftime(dates,format="%Y") == "2008"]
cdl_v_2008 <- cdl_v[, 2]
# corn_v_2008 <- v_2008[cdl_v_2008 == 61, ]
corn_v_2008 <- v_2008[cdl_v_2008 == 1, ]
mean_corn_2008 <- colMeans(corn_v_2008, na.rm=T)
corn_spl <- smooth.spline(mean_corn_2008[!is.na(mean_corn_2008)] ~ dates_2008[!is.na(mean_corn_2008)])
days_2008 <- seq.Date(as.Date("2008-1-1"), as.Date("2008-12-31"), by="day")
corn_smooth_2008 <- predict(corn_spl, as.integer(days_2008))

plot(dates_2008, mean_corn_2008)
points(days_2008, corn_smooth_2008$y, type="l",col=2)

#-------------------------
# Analyze CDL
year_of_interest <- 2008
v_year <- v[, as.integer(strftime(dates, format="%Y")) == year_of_interest]
dates_year <- dates[as.integer(strftime(dates, format="%Y")) == year_of_interest]
days_year <- seq.Date(as.Date(paste(year_of_interest, "-1-1", sep="")), as.Date(paste(year_of_interest, "-12-31", sep="")), by="day")
cdl_i <- which(as.integer(unlist(lapply(names(cdl_s), function(x) substr(x,2,5)))) == year_of_interest)
cdl_v_year <- cdl_v[, cdl_i]
cdl_h <- hist(cdl_v[, cdl_i], breaks=seq(-0.5, 255.5, by=1), plot=F)
tmp <- data.frame(cdl_code=cdl_h$mids, cdl_count=cdl_h$counts, cdl_frac=cdl_h$counts / sum(cdl_h$counts))
cdl_lcs <- tmp$cdl_code[tmp$cdl_frac > 0.01] # all LC types w/ greater than 1% scene area
plot(days_year, rep(0, length(days_year)), ylim=c(0, 0.8), type="n", xlab="", ylab="EVI")
i <- 1
plot_cols <- c("red", "green", "yellow", "blue", "purple", "orange", "black", "grey", "brown")
for(cdl_lc in cdl_lcs){
  tmp <- v_year[cdl_v_year == cdl_lc,]
  mean_lc <- colMeans(tmp, na.rm=T)
  lc_spl <- smooth.spline(mean_lc[!is.na(mean_lc)] ~ dates_year[!is.na(mean_lc)])
  smooth_year <- predict(lc_spl, as.integer(days_year))
  points(dates_year, mean_lc / 1e4, col=plot_cols[i], pch=16, cex=0.7)
  points(days_year, smooth_year$y / 1e4, type="l", col=plot_cols[i])
  i <- i + 1
}
legend("topleft", legend=cdl_lcs, pch=16, lty=1, lwd=2, col=plot_cols)

core_types <- list(c(1,5), c(141,190), c(121), c(122), c(37, 171), c(111))
all_types <- list(c(1,5,4,24,27,28,29,36,229), c(141,190,61,142,152,195), c(121), c(122,123), c(37, 171), c(111,124,131))

#------------------------
# Need to build the KF matrices for each LC type, we do this on a weekly basis
KF_params <- list()
i <- 1
for(core_type in core_types){
  type_means <- colMeans(v_year[cdl_v_year %in% core_type,], na.rm=T)
  lc_spl <- smooth.spline(type_means[!is.na(type_means)] ~ dates_year[!is.na(type_means)])
  smooth_year <- predict(lc_spl, as.integer(days_year))
  x <- smooth_year$y; x_prev <- x[1:(length(x) - 1)]; x_next <- x[2:length(x)]
  x_week <- by(x / 1e4, as.integer(strftime(days_year, format="%U")), mean, na.rm=T)
  # plot(x_week)
  # brks <- identify(x_week)
  brks <- c(10, 50) # fixed at week 10 and 50
  x_week[1:brks[1]] <- x_week[brks[1]]
  x_week[brks[2]:length(x_week)] <- x_week[brks[2]]
  pred_factor <- x_week[2:length(x_week)] / x_week[1:(length(x_week) - 1)]
  pred_factor_weekly <- c(pred_factor, 1) # repeat the tail value

  A <- array(pred_factor_weekly, dim=c(1, 1, 53))
  B <- array(rep(0, 53), dim=c(1, 1, 53)) # control matrix
  u <- array(rep(0, 53), dim=c(53, 1)) # control signal
  H <- array(rep(1, 53), dim=c(1, 1, 53)) # measurement to state transition matrix
  R <- array(rep(0.05, 53), dim=c(1, 1, 53)) # measurement noise covariance
  Q <- array(rep(0.05, 53), dim=c(1, 1, 53)) # process noise covariance
  x_init <- x_week[brks[1]]
  P_init <- 0.05
  KF_params[[i]] <- list(A=A, B=B, u=u, H=H, R=R, Q=Q, x_init=x_init, P_init=P_init)
  i <- i + 1
}

# create a lookup table for LC code and model
model_df <- data.frame(lc_type=unlist(all_types), model_num=rep(1:length(all_types), unlist(lapply(all_types, length))))

# a function to apply the KF to every row
# assumes that the first element of x is the LC type
ApplyKF <- function(x, dates, KF_parameters, KF_LC_table, do_smooth=TRUE, do_plot=FALSE){
  # retrieve the landcover
  x_lc <- x[1]
  # compute weekly mean values for filtering
  x_tmp <- by(x[2:length(x)] / 1e4, as.integer(strftime(dates, format="%U")), mean, na.rm=T)
  x_week <- rep(NA, 53) # a value for all 53 weeks
  # match x_tmp to proper week in x_week
  x_week[unique(as.integer(strftime(dates, format="%U"))) + 1] <- x_tmp
  # retrieve the proper KF matrices for the given landcover
  model_num <- KF_LC_table$model_num[KF_LC_table$lc_type == x_lc]
  KFmats <- KF_parameters[[model_num]]
  # apply the Kalman filter/smoother
  KFsmooth <- KF(x_init=KFmats$x_init, P_init=KFmats$P_init, measurements=x_week, A=KFmats$A, B=KFmats$B, u=KFmats$u, Q=KFmats$Q, H=KFmats$H, R=KFmats$R, smooth=do_smooth)

  if(do_plot){
    # DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
    KFfilter <- KF(x_init=KFmats$x_init, P_init=KFmats$P_init, measurements=x_week, A=KFmats$A, B=KFmats$B, u=KFmats$u, Q=KFmats$Q, H=KFmats$H, R=KFmats$R, smooth=FALSE)
    # DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG

    weekly_dates <- seq.Date(as.Date(paste(strftime(dates[1], format="%Y"), "-1-1", sep="")), as.Date(paste(strftime(dates[1], format="%Y"), "-12-31", sep="")), by="week")
    # weekly_dates <- seq.Date(as.Date("2008-1-1"), as.Date("2008-12-31"), by="week")

    # create uncertainty polygon
    poly_x <- c(weekly_dates, rev(weekly_dates))
    poly_y <- c((KFsmooth[[1]] + KFsmooth[[2]]), rev(KFsmooth[[1]] - KFsmooth[[2]]))

    # do the plotting
    # plot(dates, x[2:length(x)] / 1e4, xlab="", ylab="EVI", type="n", ylim=c(min(poly_y), max(poly_y)))
    plot(dates, x[2:length(x)] / 1e4, xlab="", ylab="EVI", type="n")
    polygon(poly_x, poly_y, col="lightgrey", border=NA)

    points(dates, x[2:length(x)] / 1e4, xlab="", ylab="EVI", type="p")
    points(weekly_dates, KFsmooth[[1]], type="l", col="firebrick2", lwd=2)
    points(weekly_dates, KFfilter[[1]], type="l", col="seagreen", lwd=2)

    # DEBUG
    # add filtered uncertainty
    points(weekly_dates, KFfilter[[1]] + KFfilter[[2]], type="l", lty=2, col="seagreen", lwd=2)
    points(weekly_dates, KFfilter[[1]] - KFfilter[[2]], type="l", lty=2, col="seagreen", lwd=2)

    legend("topleft", legend=c("Filtered", "Smoothed"), lty=c(1,1), col=c("seagreen", "firebrick2"))
  }
  return(KFsmooth[[1]])
}


#------------------------------
# combine the values and CDL data and apply the filter
kf_df <- cbind(v_year, cdl_v_year)
kf_df <- kf_df[, c(dim(kf_df)[2], 1:(dim(kf_df)[2] - 1))]

# single threaded version
# system.time(tmp <- apply(kf_df[1:1e3,], 1, ApplyKF, dates=dates_year, KF_parameters=KF_params, KF_LC_table=model_df))
# try in parallel
cl <- makeCluster(detectCores())
clusterExport(cl, c("KF", "ApplyKF"))
# system.time(kf_df_filtered <- parApply(cl, kf_df[1:1e3,], 1, ApplyKF, dates=dates_year, KF_parameters=KF_params, KF_LC_table=model_df))
system.time(kf_df_filtered <- parApply(cl, kf_df, 1, ApplyKF, dates=dates_year, KF_parameters=KF_params, KF_LC_table=model_df))
# system.time(tmp_filt <- parApply(cl, kf_df[1:1e3,], 1, ApplyKF, dates=dates_year, KF_parameters=KF_params, KF_LC_table=model_df, do_smooth=FALSE))

#------------------------------
# plot some example results
weekly_dates <- seq.Date(as.Date("2008-1-1"), as.Date("2008-12-31"), by="week")
k <- 25 # this is an example of where it does great for corn
trash <- ApplyKF(kf_df[k, ], dates=dates_year, KF_parameters=KF_params, KF_LC_table=model_df, do_plot=TRUE)

plot(raster(s,27))
cells <- click(raster(s,27), cell=T)
ks <- cells$cell
# ks <- c(25, 500, 3e3)
layout(matrix(1:length(ks), nrow=length(ks)))
par(mar=c(2,4,0,0), oma=rep(1,4))
for(k in ks){
  trash <- ApplyKF(kf_df[k, ], dates=dates_year, KF_parameters=KF_params, KF_LC_table=model_df, do_plot=TRUE)
}

#------------------------------
# plot both the input EVI and the smoothed data
myRamp <- colorRampPalette(brewer.pal(11, "Spectral"))
myBreaks <- c(-1e9, seq(0.1, 0.75, len=(num_cols - 2)), 1e9)
out_dir <- "~/Desktop/KF_output"
out_prefix <- "img"
tmp_r <- raster(s, 1)
s_year <- subset(s, c(1:nlayers(s))[strftime(dates, format="%Y") == 2008])

weekly_dates <- seq.Date(as.Date("2008-1-1"), as.Date("2008-12-31"), by="week")
plot_dates <- sort(unique(c(weekly_dates, dates_year)))
i <- 1
rm("kf_prev", "landsat_prev", "kf_prev_date", "landsat_prev_date")
for(this_date in plot_dates){
  print(paste("doing", i))
  out_file <- file.path(out_dir, paste(out_prefix, formatC(i, width=3, flag="0"), ".jpg", sep=""))
  jpeg(file=out_file, width=1026, height=(671*2), quality=75)
  layout(matrix(1:2, nrow=2))
  par(mar=rep(0,4), oma=rep(0,4), bg="black")

  if(this_date %in% dates_year){
    # plot the
    landsat_index <- which(dates_year == this_date)
    plot(raster(s_year, landsat_index) / 1e4, breaks=myBreaks, col=myRamp(num_cols), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
    # Annotate
    text(par()$usr[1], par()$usr[3], labels=as.Date(this_date, origin="1970-1-1"), cex=3, adj=c(-0.1, -0.2), col="white")

    landsat_prev <- raster(s_year, landsat_index) / 1e4
    landsat_prev_date <- this_date
  }else{
    if(exists("landsat_prev")){
      # plot previous image
      plot(landsat_prev, breaks=myBreaks, col=myRamp(num_cols), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
      # annotate with previous date
      text(par()$usr[1], par()$usr[3], labels=as.Date(landsat_prev_date, origin="1970-1-1"), cex=3, adj=c(-0.1, -0.2), col="white")
    }else{
      plot(1:10, type="n", xlab="", ylab="", xaxt="n", yaxt="n")
    }
  }

  if(this_date %in% weekly_dates){
    # plot the KF smoothed
    kf_index <- which(weekly_dates == this_date)
    values(tmp_r) <- kf_df_filtered[kf_index,]
    plot(tmp_r, breaks=myBreaks, col=myRamp(num_cols), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
    # Annotate
    text(par()$usr[1], par()$usr[3], labels=as.Date(this_date, origin="1970-1-1"), cex=3, adj=c(-0.1, -0.2), col="white")

    kf_prev <- tmp_r
    kf_prev_date <- this_date
  }else{
    if(exists("kf_prev")){
      # plot previous image
      plot(kf_prev, breaks=myBreaks, col=myRamp(num_cols), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
      # annotate with previous date
      text(par()$usr[1], par()$usr[3], labels=as.Date(kf_prev_date, origin="1970-1-1"), cex=3, adj=c(-0.1, -0.2), col="white")
    }else{
      plot(1:10, type="n", xlab="", ylab="", xaxt="n", yaxt="n")
    }
  }

  dev.off()
  i <- i + 1
}

# make into a movie with ffmpeg:
# module load ffmpeg/git-2013-12-18-4a2570f
# ffmpeg -framerate 2 -i img%03d.jpg -pix_fmt yuv420p -vf scale="720:trunc(ow/a/2)*2" DL_KF_output.mp4
s



#------------------------
# What about weekly prediction instead?
core_type <- core_types[[1]]
type_means <- colMeans(v_year[cdl_v_year %in% core_type,], na.rm=T)
lc_spl <- smooth.spline(type_means[!is.na(type_means)] ~ dates_year[!is.na(type_means)])
smooth_year <- predict(lc_spl, as.integer(days_year))
x <- smooth_year$y; x_prev <- x[1:(length(x) - 1)]; x_next <- x[2:length(x)]
x_week <- by(x / 1e4, as.integer(strftime(days_year, format="%U")), mean, na.rm=T)
# repeat the head and tail for flatness of min value
x_week[1:11] <- x_week[11]
x_week[50:length(x_week)] <- x_week[50]
pred_factor <- x_week[2:length(x_week)] / x_week[1:(length(x_week) - 1)]
# pred_factor <- 1 + ((pred_factor - 1) / 7 ) # adjust to daily step
pred_factor_weekly <- c(pred_factor, 1) # repeat the tail value

# build a measurement vector
t_m <- 1:53
# z <- v_year[191808, ]/1e4
z <- v_year[129982, ]/1e4
z_week <- by(z, as.integer(strftime(dates_year, format="%U")), mean, na.rm=T)
z_m <- rep(NA, length(t_m))
z_m[match(as.integer(names(z_week)), (t_m - 1))] <- z_week
# z_m[match(dates_year, days_year)] <- z

# build the DLM
A <- array(pred_factor_weekly, dim=c(1, 1, length(t_m)))
B <- array(rep(0, length(t_m)), dim=c(1, 1, length(t_m))) # control matrix
u <- array(rep(0, length(t_m)), dim=c(length(t_m), 1)) # control signal
H <- array(rep(1, length(t_m)), dim=c(1, 1, length(t_m))) # measurement to state transition matrix
R <- array(rep(0.05, length(t_m)), dim=c(1, 1, length(t_m))) # measurement noise covariance
Q <- array(rep(0.01, length(t_m)), dim=c(1, 1, length(t_m))) # process noise covariance

# use the other KF function from KalmanFilter.R
KFest <- KF(0.1, 0.5, z_m, A, B, u, Q, H, R)
KFsmooth <- KF(0.1, 0.5, z_m, A, B, u, Q, H, R, smooth=T)
plot(t_m, z_m, type="n")
# plot(t_m, KFest[[1]], type="n")
poly_x <- c(t_m, rev(t_m))
# poly_y <- c((KFest[[1]] + KFest[[2]]), rev(KFest[[1]] - KFest[[2]]))
poly_y <- c((KFsmooth[[1]] + KFsmooth[[2]]), rev(KFsmooth[[1]] - KFsmooth[[2]]))
polygon(poly_x, poly_y, col="lightgrey", border=NA)
points(t_m, KFest[[1]], col=2, type="l")
points(t_m, KFsmooth[[1]], col=3, type="l")
points(t_m, z_m, type="p")

#------------------------
# fit the DLMs for each type by plotting the average for the core types, identifying
# breakpoints by hand, and then calculating the DLM slope for each section
core_types <- list(c(1,5), c(141,190), c(121), c(122), c(37, 171), c(111))
all_types <- list(c(1,5,4,24,27,28,29,36,229), c(141,190,61,142,152,195), c(121), c(122,123), c(37, 171), c(111,124,131))
core_breaks <- list(rep(NA, length(core_types)))
core_slopes <- list(rep(NA, length(core_types)))
i <- 1
for(core_type in core_types){
  type_means <- colMeans(v_year[cdl_v_year %in% core_type,], na.rm=T)
  lc_spl <- smooth.spline(type_means[!is.na(type_means)] ~ dates_year[!is.na(type_means)])
  smooth_year <- predict(lc_spl, as.integer(days_year))
  # plot(as.integer(strftime(days_year, format="%j")), smooth_year$y)
  x <- smooth_year$y; x_prev <- x[1:(length(x) - 1)]; x_next <- x[2:length(x)]
  plot(x_next/x_prev)
  brks <- c(1, identify(smooth_year$y), 365)
  core_breaks[[i]] <- brks
  print(core_type)
  tmp_slopes <- rep(NA, length(brks) - 1)
  for(j in 1:(length(brks) - 1)){
    tmp_prev <- smooth_year$y[brks[j]:(brks[j + 1] - 1)]
    tmp_next <- smooth_year$y[(brks[j] + 1):brks[j + 1]]
    tmp_slopes[j] <- mean(tmp_next / tmp_prev)
    print(mean(tmp_next / tmp_prev))
  }
  core_slopes[[i]] <- tmp_slopes
  print("-------------------------------")
  i <- i + 1
}

for(brks in core_breaks){
  for(i in 1:(length(brks) - 1)){
  	tmp_prev <- corn_smooth_2008$y[brks[i]:(brks[i + 1] - 1)]
  	tmp_next <- corn_smooth_2008$y[(brks[i] + 1):brks[i + 1]]
  	print(mean(tmp_next / tmp_prev))
  }
}

# brks <- c(150, 205, 235, 290)
brks <- c(1, 150, 205, 235, 290, 365)
for(i in 1:(length(brks) - 1)){
	tmp_prev <- corn_smooth_2008$y[brks[i]:(brks[i + 1] - 1)]
	tmp_next <- corn_smooth_2008$y[(brks[i] + 1):brks[i + 1]]
	print(mean(tmp_next / tmp_prev))
}

#----------------------------------------------------------
# Kalman Filtering Section

#----------------------------------------------------------
# KF Functions
TimeUpdate <- function(t, x, A, B, P, Q, u){
  x_prior <- (A[,,t] %*% x) + (B[,,t] %*% u[t,])
  P_prior <- (A[,,t] %*% P %*% t(A[,,t])) + Q[,,t]
  return(list(x_prior, P_prior))
}

# KF Measurement Update
# z is the (1 x m) measurement vector, H is the (m x m x t) measurement relation matrix, and
# R is the (m x m x t) measurement noise covariance matrix
MeasurementUpdate <- function(t, z, x_prior, P_prior, H, R){
  if(is.na(z)) return(list(x_prior, P_prior)) # no measurement, no measurement update
  K <- (P_prior %*% t(H[,,t])) %*% solve(H[,,t] %*% P_prior %*% t(H[,,t]) + R[,,t])
  x <- x_prior + (K %*% (z - (H[,,t] %*% x_prior)))
  P <- (diag(dim(K)[1]) - (K %*% H[,,t])) %*% P_prior
  return(list(x, P))
}

#----------------------------------------------------------
# build the DLM
A <- array(c(rep(1, length(1:15)), rep(1.2, length(16:27)), rep(0.85, length(28:36)), rep(1, length(37:46))), dim=c(1, 1, 46))
B <- array(rep(0, length(t_m)), dim=c(1, 1, length(t_m))) # control matrix
u <- array(rep(0, length(t_m)), dim=c(length(t_m), 1)) # control signal
H <- array(rep(1, length(t_m)), dim=c(1, 1, length(t_m))) # measurement to state transition matrix
R <- array(rep(0.03, length(t_m)), dim=c(1, 1, length(t_m))) # measurement noise covariance
Q <- array(rep(0.05, length(t_m)), dim=c(1, 1, length(t_m))) # process noise covariance

# iterate through the data and filter the data
x <- 0 # initial state estimate
P <- 0.02 # initial sd estimate
x_est <- rep(NA, length(t_m)) # output estimate vector
P_est <- rep(NA, length(t_m)) # output uncertainty vector
for(t in t_m){
  t_up <- TimeUpdate(t, x, A, B, P, Q, u)
  # x_prior <- t_up[[1]]
  # P_prior <- t_up[[2]]
  m_up <- MeasurementUpdate(t, z_m[t], t_up[[1]], t_up[[2]], H, R)
  # extract posterior estimates of current state from list
  x <-  m_up[[1]]
  P <-  m_up[[2]]
  # update output vectors
  x_est[t] <- x
  P_est[t] <- P
}

# get the 95% confidence interval
conf_f <- function(x, p) qnorm(p, x[1], x[2])
lower_env <- apply(data.frame(x_est, P_est), 1, conf_f, p=0.05)
upper_env <- apply(data.frame(x_est, P_est), 1, conf_f, p=0.95)
# plot
plot(t_m, z_m, type="n", ylim=c(min(lower_env), max(upper_env)), xlab="t", ylab="NBAR-EVI")
polygon(c(t_m, rev(t_m)), c(lower_env, rev(upper_env)), col="lightgrey", border=NA)
points(t_m, x_est, col=2, type="l", lwd=1.5)
points(t_m, z_m, type="p")
legend("topleft", legend=c("Observations", "KF Estimate", "95% Conf. Interval"), pch=c(1, NA, 15), lty=c(NA, 1, NA), col=c(1, 2, "lightgrey"), bty="n", pt.cex=c(1, 1, 2), lwd=c(NA, 1.5, NA))
