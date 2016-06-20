#############################################
# R Kalman Filter with package "dlm"
# Josh Gray
# Nov, 2015
#############################################

library(dlm)
library(parallel)
library(raster)
library(tools)
library(rgdal)
# NOTE: use KF_crop_scratch.R to get the required variables

# translate the KF_params into the format that dlm expects:
# theta(t) = G*theta(t-1) + w(t)
# y(t) = F*theta(t) + v(t)
# v(t) ~ N(0, V)
# w(t) ~ N(0, W)
# theta(0) ~ N(m0, C0)
# m0: x0
# C0: P0
# FF: observation transition matrix (H)
# V: observation variance (R)
# GG: process transition matrix (A)
# W: process variance (Q)
# JFF: column in X to find time varying FF
# JV: column in X to find time varying V
# JGG: column in X to find time varying GG
# JW: column in X to find time varying W
# X: time varying parameters

params <- KF_params[[1]]
z_daily <- kf_df[25, 2:61] / 1e4
z_weekly <- by(z_daily, strftime(dates_year, "%U"), mean, na.rm=T)
z <- rep(NA, 53)
z[as.integer(names(z_weekly)) + 1] <- z_weekly
z[is.nan(z)] <- NA
# z[unique(as.integer(strftime(dates_year, "%U"))) + 1] <- z_weekly
# plot(dates_year, z)
JGG <- 1 # GG is time varying and the coefficients can be found in column 1 of array X
X <- array(1, dim=c(length(z), 1)) # initialize the time varying parameter array to have only one time varying column: GG
X[,1] <- c(params$A) # fill the first column of X with the time varying GG parameters

# define a DLM with time varying JGG and constant process/measurement variance
my_dlm <- dlm(m0=0.1, C0=0.005, FF=1, V=0.002, GG=1, W=0.005, JFF=0, JGG=1, JV=0, JW=0, X=X)

# Kalman Filtering
my_dlm_filter <- dlmFilter(z, my_dlm) # filter the series with the model
v <- unlist(dlmSvd2var(my_dlm_filter$U.C, my_dlm_filter$D.C))
pl <- dropFirst(my_dlm_filter$m) + qnorm(0.05, sd=sqrt(v[-1]))
pu <- dropFirst(my_dlm_filter$m) - qnorm(0.05, sd=sqrt(v[-1]))
plot(1:53, z, col="seagreen")
lines(dropFirst(my_dlm_filter$m), col="brown")
lines(pl, col="brown", lty=2)
lines(pu, col="brown", lty=2)

# Kalman Smoothing
# my_dlm_smooth <- dlmSmooth(z, my_dlm) # smooth the series with the model
my_dlm_smooth <- dlmSmooth(my_dlm_filter) # smooth the series with the model
v <- unlist(dlmSvd2var(my_dlm_smooth$U.S, my_dlm_smooth$D.S))
pl <- dropFirst(my_dlm_smooth$s) + qnorm(0.05, sd=sqrt(v[-1]))
pu <- dropFirst(my_dlm_smooth$s) - qnorm(0.05, sd=sqrt(v[-1]))
plot(1:53, z, col="black", ylim=c(0, 0.8))
lines(dropFirst(my_dlm_smooth$s), col="darkgreen", lwd=2)
lines(pl, col="darkgreen", lty=2)
lines(pu, col="darkgreen", lty=2)

# compare to an uninformative model...

# use MLE to optimize the process and measurement variance and initial values

######################################################
# multisensor fusion version
# simulated second signal by adding noise to smoothed output and adding random missing values
set.seed(42)
z2 <- dropFirst(my_dlm_smooth$s) + rnorm(length(dropFirst(my_dlm_smooth$s)), mean=0, sd=0.025)
z2[sample(1:length(z2), 20)] <- NA
z1_var <- 0.002
z2_var <- 0.005
FF <- array(1, dim=c(2, 1)) # construct observation matrix, now 2x2
JFF <- array(0, dim=c(2, 1)) # we don't allow the observation matrices to be time varying
V <- array(0, dim=c(2, 2))
diag(V) <- c(z1_var, z2_var) # setup observation covariance matrix
JV <- array(0, dim=c(2, 2)) # we don't allow the observation covariance to be time varying
my_multi_dlm <- dlm(m0=0.1, C0=0.005, FF=FF, V=V, GG=1, W=0.005, JFF=JFF, JGG=1, JV=JV, JW=0, X=X)
z_multi <- cbind(z, z2) # combine observations together in array
my_multi_dlm_filter <- dlmFilter(z_multi, my_multi_dlm) # smooth the series with the model
my_multi_dlm_smooth <- dlmSmooth(my_multi_dlm_filter)
v <- unlist(dlmSvd2var(my_multi_dlm_smooth$U.S, my_multi_dlm_smooth$D.S))
pl <- dropFirst(my_multi_dlm_smooth$s) + qnorm(0.05, sd=sqrt(v[-1]))
pu <- dropFirst(my_multi_dlm_smooth$s) - qnorm(0.05, sd=sqrt(v[-1]))
plot(1:53, z, col="black", ylim=c(0, 0.8))
points(1:53, z2, col="black", pch=4)
lines(dropFirst(my_multi_dlm_smooth$s), col="darkgreen", lwd=2)
lines(pl, col="darkgreen", lty=2)
lines(pu, col="darkgreen", lty=2)
legend("topleft", legend=c("z1", "z2", "KF Smoothed", "95% CI"), pch=c(1, 4, NA, NA), lty=c(NA, NA, 1, 2), lwd=c(NA, NA, 2, 1), col=c(1, 1, "darkgreen", "darkgreen"), bty="n")



######################################################
# multispectral, multisensor
# here we observe n=3 states, each with two different sensors (p=3*2=6), and T=100 (N)
######################################################

MakeMultisensorDLM <- function(n, p, N=NA){
  # prototypes a DLM appropriate for multisensor fusion. Specifically: a model for n states observed with p sensors, where each state is observed by every sensor and a is expected to be: [a1, a2, ..., an, b1, b2, ..., bn, c1, c2, ..., cn] for sensors a, b, c

  # initialize
  m0 <- array(0, dim=c(n, 1))
  C0 <- array(0, dim=c(n, n))
  diag(C0) <- 1

  # process model parameters
  GG <- array(0, dim=c(n, n)) # coefficients
  diag(GG) <- 1 # initialize to 1
  W <- array(0, dim=c(n, n)) # state covariance
  diag(W) <- 1 # initialize diagonals to 1; no covariance between states

  # observation model parameters
  FF <- array(0, dim=c(p, n)) # observation coefficient matrix
  FF[col(FF) == (row(FF) - (n * ((row(FF) - 1) %/% n)))] <- 1 # stupid way to set the "diagonals" to 1
  V <- array(0, dim=c(p, p)) # observation variance matrix
  diag(V) <- 1 # intialize observation variance to 1 w/ no sensor covariance

  ## time varying parameters
  if(!is.na(N)){
    JGG <- array(0, dim=c(n, n)) # time varying coefficients
    diag(JGG) <- 1:n # diagonals have column # in X where time varying GG is found
    iX <- (n + 1) # next starting column # in X
    JW <- array(0, dim=c(n, n)) # time varying variance
    diag(JW) <- iX:(iX + n - 1)
    iX <- (iX + n) # next starting column # in X

    JFF <- array(0, dim=c(p, n)) # to allow or time varying coefficients
    # all of this chicanery is to get the proper column number in X in the proper place here
    # JFF[col(JFF) == (row(JFF) - (n * ((row(JFF) - 1) %/% n)))] <- iX:(iX + p - 1)
    tmp <- t(JFF)
    tmp[t(col(JFF)) == (t(row(JFF)) - (n * ((t(row(JFF)) - 1) %/% n)))] <- iX:(iX + p - 1)
    JFF <- t(tmp)
    iX <- iX + p # increment iX
    JV <- array(0, dim=c(p, p)) # to allow time varying sensor variances
    diag(JV) <- iX:(iX + p - 1) # the column in X where we find the time varying variances for each sensor
    iX <- iX + p

    # NOTE: need to finish this part out
    # set up the time varying coefficients
    X <- array(1, dim=c(N, (iX - 1)))
    X[, 1] <- rep(1, N) # time varying GG coefficients for state 1
    X[, 2] <- rep(1, N) # time varying GG coefficients for state 2
    # ...  and so on...

    # create the dlm
    msmstv_dlm <- dlm(m0=m0, C0=C0, FF=FF, V=V, GG=GG, W=W, JFF=JFF, JGG=JGG, JV=JV, X=X) # multistate, multisensor, time varying parameter DLM

    return(msmstv_dlm)
  }else{
    # create the dlm
    msms_dlm <- dlm(m0=m0, C0=C0, FF=FF, V=V, GG=GG, W=W) # multistate, multisensor DLM
    return(msms_dlm)
  }
}

#-----------------------------------------------
N <- 100 # number of time steps
n <- 3 # number of states
p <- 6 # number of sensors

# create some random walk data: 3 states observed by two sensors (n=3, p=6)
set.seed(42)
x1 <- rep(0, N)
x2 <- rep(0, N)
x3 <- rep(0, N)
x1var <- x2var <- x3var <- 1
s1var <- s2var <- 0.5
for(i in 2:N){
  x1[i] <- x1[i - 1] + rnorm(1, 0, sqrt(x1var))
  x2[i] <- x2[i - 1] + rnorm(1, 0, sqrt(x2var))
  x3[i] <- x3[i - 1] + rnorm(1, 0, sqrt(x3var))
}
# add noise
s11 <- x1 + rnorm(N, 0, sqrt(s1var))
s12 <- x2 + rnorm(N, 0, sqrt(s1var))
s13 <- x3 + rnorm(N, 0, sqrt(s1var))
s21 <- x1 + rnorm(N, 0, sqrt(s2var))
s22 <- x2 + rnorm(N, 0, sqrt(s2var))
s23 <- x3 + rnorm(N, 0, sqrt(s2var))

# add missing data
s1_num_missing <- 35
s2_num_missing <- 50
s1_missing_sample <- sample(1:N, s1_num_missing)
s2_missing_sample <- sample(1:N, s2_num_missing)
s11[s1_missing_sample] <- s12[s1_missing_sample] <- s13[s1_missing_sample] <- NA
s21[s2_missing_sample] <- s22[s2_missing_sample] <- s23[s2_missing_sample] <- NA

# bundle into observation array
z <- rbind(s11, s12, s13, s21, s22, s23)

#-----------------------------------------------
# create the dlm
msms_dlm <- MakeMultisensorDLM(n, p)

# use dlmFilter and dlmSmooth
# msms_dlm_filter <- dlmFilter(z, msms_dlm) # this doesn't work...
msms_dlm_filter <- dlmFilter(t(z), msms_dlm) # but this does
msms_dlm_smooth <- dlmSmooth(t(z), msms_dlm) # but this does

#-----------------------------------------------
# extract the standard errors, sd is Nxn and s[i,j] is the sd of the jth state at time i
sd_filter <- sqrt(do.call(rbind, lapply(dropFirst(dlmSvd2var(msms_dlm_filter$U.C, msms_dlm_filter$D.C)), diag))) # for filtering version
sd_smooth <- sqrt(do.call(rbind, lapply(dropFirst(dlmSvd2var(msms_dlm_smooth$U.S, msms_dlm_smooth$D.S)), diag))) # for smoothing version
CI_level <- 0.05
CI_filter <- qnorm(CI_level, sd=sd_filter) # add/sub to m/s for upper/lower CI
CI_smooth <- qnorm(CI_level, sd=sd_smooth) # add/sub to m/s for upper/lower CI

#-----------------------------------------------
# plot the results
smooth <- TRUE # if FALSE, then the filtering distributions are plotted, otherwise the smoothing ones
layout(matrix(1:3, nrow=3))
par(mar=c(1, 4, 1, 1))

# state 1
plot(x1, type="l", lwd=2, col="darkgrey")
points(s11, col=2, pch=2); points(s21, col=4, pch=4)
if(!smooth){
  points(dropFirst(msms_dlm_filter$m)[,1], type="l", col="black", lty=2, lwd=2)
  points(dropFirst(msms_dlm_filter$m)[,1] + CI_filter[,1], type="l", col="black", lty=3, lwd=1)
  points(dropFirst(msms_dlm_filter$m)[,1] - CI_filter[,1], type="l", col="black", lty=3, lwd=1)
}else{
  points(dropFirst(msms_dlm_smooth$s)[,1], type="l", col="black", lty=2, lwd=2)
  points(dropFirst(msms_dlm_smooth$s)[,1] + CI_smooth[,1], type="l", col="black", lty=3, lwd=1)
  points(dropFirst(msms_dlm_smooth$s)[,1] - CI_smooth[,1], type="l", col="black", lty=3, lwd=1)
}

# state 2
plot(x2, type="l", lwd=2, col="darkgrey")
points(s12, col=2, pch=2); points(s22, col=4, pch=4)
if(!smooth){
  points(dropFirst(msms_dlm_filter$m)[,2], type="l", col="black", lty=2, lwd=2)
  points(dropFirst(msms_dlm_filter$m)[,2] + CI_filter[,2], type="l", col="black", lty=3, lwd=1)
  points(dropFirst(msms_dlm_filter$m)[,2] - CI_filter[,2], type="l", col="black", lty=3, lwd=1)
}else{
  points(dropFirst(msms_dlm_smooth$s)[,2], type="l", col="black", lty=2, lwd=2)
  points(dropFirst(msms_dlm_smooth$s)[,2] + CI_smooth[,2], type="l", col="black", lty=3, lwd=1)
  points(dropFirst(msms_dlm_smooth$s)[,2] - CI_smooth[,2], type="l", col="black", lty=3, lwd=1)
}

# state 3
plot(x3, type="l", lwd=2, col="darkgrey")
points(s13, col=2, pch=2); points(s23, col=4, pch=4)
if(!smooth){
  points(dropFirst(msms_dlm_filter$m)[,3], type="l", col="black", lty=2, lwd=2)
  points(dropFirst(msms_dlm_filter$m)[,3] + CI_filter[,3], type="l", col="black", lty=3, lwd=1)
  points(dropFirst(msms_dlm_filter$m)[,3] - CI_filter[,3], type="l", col="black", lty=3, lwd=1)
}else{
  points(dropFirst(msms_dlm_smooth$s)[,3], type="l", col="black", lty=2, lwd=2)
  points(dropFirst(msms_dlm_smooth$s)[,3] + CI_smooth[,3], type="l", col="black", lty=3, lwd=1)
  points(dropFirst(msms_dlm_smooth$s)[,3] - CI_smooth[,3], type="l", col="black", lty=3, lwd=1)
}

# make a legend
legend("bottomright", legend=c("signal", "sensor 1", "sensor 2", "KF", "KF 95% CI"), lwd=c(2, NA, NA, 2, 1), col=c("darkgrey", 2, 4, "black", "black"), pch=c(NA, 2, 4, NA, NA), lty=c(1, NA, NA, 2, 3))


######################################################
######################################################
######################################################
# Apply to real data: Landsat/MODIS fusion problem
# we fuse R, G, B, and NIR bands utilizing a Kalman smoother
# based on an non-innovative, temporally static process model

# get the landsat data
landsat_data_dir <- "/projectnb/modislc/users/joshgray/DL_Landsat"
landsat_in_files <- dir(path=landsat_data_dir, pattern="LFstack_fuse.tif", full=T, rec=T)
landsat_dates <- do.call("c", lapply(basename(landsat_in_files), function(x) as.Date(paste(substr(x, 10, 13), substr(x, 14, 16), sep="-"), format="%Y-%j")))
landsat_in_files <- landsat_in_files[order(landsat_dates)]
landsat_dates <- sort(landsat_dates)
# subset files and dates to year 2008
landsat_in_files <- landsat_in_files[landsat_dates >= as.Date("2008-1-1") & landsat_dates <= as.Date("2008-12-31")]
landsat_dates <- landsat_dates[landsat_dates >= as.Date("2008-1-1") & landsat_dates <= as.Date("2008-12-31")]

# get the MODIS data
modis_data_dir <- "/projectnb/modislc/users/joshgray/DL_Landsat/MODIS2008"
modis_in_files <- dir(path=modis_data_dir, pattern="utm_30m.tif", full=T, rec=T)
modis_dates <-  as.Date(substr(basename(modis_in_files), 10, 16), format="%Y%j")

# loop through all dates to build the data cube
# NOTE: this needs to be heavily modified in order to scale
pred_dates <- seq.Date(as.Date("2008-1-1"), as.Date("2008-12-31"), by="day")
# the observation array should be NxpxM where N is the # of time steps, p is the number of sensors, and M is the number of pixels. Thus, we apply over the 3rd dimension, feeding dlmSmooth an Nxp array each time
M <- ncell(raster(landsat_in_files[1]))
N <- length(pred_dates)
p <- 8 # the number of sensors
Z <- array(NA, dim=c(N, p, M))

landsat_sensors <- 1:4 # which of the second dimension indices belong to landsat sensors?
modis_sensors <- 5:8 # which of the second dimension indices belong to modis sensors?

i <- 1
for(this_date in pred_dates){
  print(paste("Working on", as.Date(this_date, origin="1970-1-1")))
  # check for landsat data
  if(this_date %in% landsat_dates){
    print("Getting Landsat data")
    tmp_s <- stack(landsat_in_files[which(landsat_dates == this_date)])
    tmp_Z <- values(tmp_s)
    Z[i, landsat_sensors, ] <- t(tmp_Z)
  }else{
    # there is no landsat data to grab
    print("No Landsat data")
  }

  # check for modis data
  if(this_date %in% modis_dates){
    print("Getting MODIS data")
    tmp_s <- stack(modis_in_files[which(modis_dates == this_date)])
    tmp_Z <- values(tmp_s)
    Z[i, modis_sensors, ] <- t(tmp_Z)
  }else{
    # there is no modis data to grab
    print("No MODIS data")
  }
  i <- i + 1
}

# create a dlm for 4 states and 8 sensors
landsat_modis_dlm <- MakeMultisensorDLM(n=4, p=8)
m0(landsat_modis_dlm) <- rep(5e3, 4) # initialize to 0.5 reflectance
diag(C0(landsat_modis_dlm)) <- rep(1e6, 4) # intialize variance
# assuming that the max change from day to day is 300 (0.03), then the process variance is less than 300^2
diag(W(landsat_modis_dlm)) <- rep(2.5e4, 4) # process noise
# diag(V(landsat_modis_dlm)) <- rep(1e3, 8) # measurement noise all constant at 1e3
landsat_var <- 62.5e3 # assume that s.d. is 0.025 units (0.025*1e4)^2
diag(V(landsat_modis_dlm)) <- c(rep(landsat_var, 4), rep((3 * landsat_var), 4)) # inflate MODIS variance by 3

# function to plot the kalman filter for the landsat_modis fusion model
PlotKF <- function(obs, dlm_model, pred_dates, CI_level=0.05, modis_color="dodgerblue", landsat_color="firebrick2", kf_color="black", kf_uncertainty_col="black", pt_cex=1.5, legend_pos="topright"){
  # smooth the data and plot
  tmp <- dlmSmooth(obs, dlm_model) # do the KF smoothing
  sd_smooth <- sqrt(do.call(rbind, lapply(dropFirst(dlmSvd2var(tmp$U.S, tmp$D.S)), diag)))
  CI_smooth <- qnorm(CI_level, sd=sd_smooth) # add/sub to m/s for upper/lower CI

  layout(matrix(1:4, nrow=2))
  par(mar=c(3, 4, 1, 1), oma=rep(0.5, 4))
  ylim <- c(min(min(dropFirst(tmp$s[, 1]) - CI_smooth[, 1], na.rm=T), min(obs[,1], na.rm=T), min(obs[,5], na.rm=T)), max(max(dropFirst(tmp$s[, 1]) + CI_smooth[, 1], na.rm=T), max(obs[,1], na.rm=T), max(obs[,5], na.rm=T)))
  plot(pred_dates, obs[,1], col=landsat_color, pch=1, ylim=ylim, cex=pt_cex, xlab="", ylab="Surface Reflectance x 1e4")
  points(pred_dates, obs[,5], col=modis_color, pch=2, cex=pt_cex)
  points(pred_dates, dropFirst(tmp$s[, 1]), type="l", lwd=2, col=kf_color, lty=2)
  points(pred_dates, dropFirst(tmp$s[, 1]) + CI_smooth[, 1], type="l", lwd=1, col=kf_uncertainty_col, lty=3)
  points(pred_dates, dropFirst(tmp$s[, 1]) - CI_smooth[, 1], type="l", lwd=1, col=kf_uncertainty_col, lty=3)
  title("Blue reflectance")

  ylim <- c(min(min(dropFirst(tmp$s[, 2]) - CI_smooth[, 2], na.rm=T), min(obs[,2], na.rm=T), min(obs[,6], na.rm=T)), max(max(dropFirst(tmp$s[, 2]) + CI_smooth[, 2], na.rm=T), max(obs[,2], na.rm=T), max(obs[,6], na.rm=T)))
  plot(pred_dates, obs[,2], col=landsat_color, pch=1, ylim=ylim, cex=pt_cex, xlab="", ylab="Surface Reflectance x 1e4")
  points(pred_dates, obs[,6], col=modis_color, pch=2, cex=pt_cex)
  points(pred_dates, dropFirst(tmp$s[, 2]), type="l", lwd=2, col=kf_color, lty=2)
  points(pred_dates, dropFirst(tmp$s[, 2]) + CI_smooth[, 2], type="l", lwd=1, col=kf_uncertainty_col, lty=3)
  points(pred_dates, dropFirst(tmp$s[, 2]) - CI_smooth[, 2], type="l", lwd=1, col=kf_uncertainty_col, lty=3)
  title("Green reflectance")

  ylim <- c(min(min(dropFirst(tmp$s[, 3]) - CI_smooth[, 3], na.rm=T), min(obs[,3], na.rm=T), min(obs[,7], na.rm=T)), max(max(dropFirst(tmp$s[, 3]) + CI_smooth[, 3], na.rm=T), max(obs[,3], na.rm=T), max(obs[,7], na.rm=T)))
  plot(pred_dates, obs[,3], col=landsat_color, pch=1, ylim=ylim, cex=pt_cex, xlab="", ylab="Surface Reflectance x 1e4")
  points(pred_dates, obs[,7], col=modis_color, pch=2, cex=pt_cex)
  points(pred_dates, dropFirst(tmp$s[, 3]), type="l", lwd=2, col=kf_color, lty=2)
  points(pred_dates, dropFirst(tmp$s[, 3]) + CI_smooth[, 3], type="l", lwd=1, col=kf_uncertainty_col, lty=3)
  points(pred_dates, dropFirst(tmp$s[, 3]) - CI_smooth[, 3], type="l", lwd=1, col=kf_uncertainty_col, lty=3)
  title("Red reflectance")

  ylim <- c(min(min(dropFirst(tmp$s[, 4]) - CI_smooth[, 4], na.rm=T), min(obs[,4], na.rm=T), min(obs[,8], na.rm=T)), max(max(dropFirst(tmp$s[, 4]) + CI_smooth[, 4], na.rm=T), max(obs[,4], na.rm=T), max(obs[,8], na.rm=T)))
  plot(pred_dates, obs[,4], col=landsat_color, pch=1, ylim=ylim, cex=pt_cex, xlab="", ylab="Surface Reflectance x 1e4")
  points(pred_dates, obs[,8], col=modis_color, pch=2, cex=pt_cex)
  points(pred_dates, dropFirst(tmp$s[, 4]), type="l", lwd=2, col=kf_color, lty=2)
  points(pred_dates, dropFirst(tmp$s[, 4]) + CI_smooth[, 4], type="l", lwd=1, col=kf_uncertainty_col, lty=3)
  points(pred_dates, dropFirst(tmp$s[, 4]) - CI_smooth[, 4], type="l", lwd=1, col=kf_uncertainty_col, lty=3)
  title("NIR reflectance")

  # make a legend
  legend(legend_pos, legend=c("MODIS", "Landsat", "KF", "KF CI"), pch=c(1, 2, NA, NA), lty=c(NA, NA, 2, 3), col=c(modis_color, landsat_color, kf_color, kf_uncertainty_col), lwd=c(NA, NA, 2, 1), cex=pt_cex)
}

# plot a few series
i <- 1e3
PlotKF(Z[,,i], landsat_modis_dlm, pred_dates, legend_pos="bottomleft")

# apply to the entire Z matrix
ApplyKF <- function(x, dlm_model){
  # function to compute the KF smoother on the matrix x and return the state and uncertainty estimates
  tmp <- dlmSmooth(x, dlm_model) # do the KF smoothing
  tmp_sd <- sqrt(do.call(rbind, lapply(dropFirst(dlmSvd2var(tmp$U.S, tmp$D.S)), diag)))
  # return(cbind(dropFirst(tmp$s), tmp_sd))
  return(list(dropFirst(tmp$s), tmp_sd))
}

cl <- makeCluster(16)
clusterEvalQ(cl, {library(dlm)})
clusterExport(cl, c("landsat_modis_dlm", "ApplyKF"))
# system.time(test <- parApply(cl, Z[,,1:1e3], 3, ApplyKF, dlm_model=landsat_modis_dlm))
system.time(landsat_modis_KF_results <- parApply(cl, Z, 3, ApplyKF, dlm_model=landsat_modis_dlm))
# test[[100]][[1]][,4] # the fourth smoothed system state
save(landsat_modis_KF_results, file="/projectnb/modislc/users/joshgray/DL_Landsat/MODIS_Landsat_Fusion/KF_out.Rdata")

# extract the results
extract_f <- function(x, i) unlist(x[[1]][, i])
system.time(nir_v <- do.call(rbind, parLapply(cl, landsat_modis_KF_results, extract_f, i=4)))
system.time(red_v <- do.call(rbind, parLapply(cl, landsat_modis_KF_results, extract_f, i=3)))
system.time(green_v <- do.call(rbind, parLapply(cl, landsat_modis_KF_results, extract_f, i=2)))


# break all the results out, make rasters, and plot it!
out_img_dir <- "/projectnb/modislc/users/joshgray/DL_Landsat/MODIS_Landsat_Fusion"
img_prefix <- "FusionImage"
tmp_r <- raster(landsat_in_files[1])
NA_stack <- stack(landsat_in_files[1])
NA_stack <- setValues(NA_stack, NA)

# loop through all doys and make the picture
i <- 1
for(this_date in pred_dates){
  print(paste("Doing", this_date))
  # create the outfile name
  out_file <- file.path(out_img_dir, paste(img_prefix, "_", formatC(i, width=3, flag="0"), ".jpg", sep=""))
  jpeg(file=out_file, width=(1026*2), height=(671*2), quality=75)
  layout(matrix(c(1, 2, 3, 3), nrow=2, byrow=T))
  par(mar=rep(1, 4), oma=rep(0.5, 4))

  # check for and plot landsat data
  if(this_date %in% landsat_dates){
    landsat_s <- stack(landsat_in_files[which(landsat_dates == this_date)])
    plotRGB(landsat_s, 4, 3, 2, scale=1e4, colNA="black", stretch="lin")
  }else{
    plotRGB(NA_stack, 4, 3, 2, scale=1e4, colNA="black", stretch="lin")
  }

  # check for and plot modis data
  if(this_date %in% modis_dates){
    modis_s <- stack(modis_in_files[which(modis_dates == this_date)])
    plotRGB(modis_s, 4, 3, 2, scale=1e4, colNA="black", stretch="lin")
  }else{
    plotRGB(NA_stack, 4, 3, 2, scale=1e4, colNA="black", stretch="lin")
  }

  # create a raster of the KF smoothed output
  nir_tmp <- tmp_r
  values(nir_tmp) <- nir_v[,i]
  red_tmp <- tmp_r
  values(red_tmp) <- red_v[,i]
  green_tmp <- tmp_r
  values(green_tmp) <- green_v[,i]
  tmp_s <- stack(nir_tmp, red_tmp, green_tmp)
  plotRGB(tmp_s, 1, 2, 3, scale=1e4, colNA="black", stretch="lin")
  text(mean(par()$usr[1:2]), par()$usr[3], labels=this_date, col="yellow", cex=3, adj=c(0.5, -0.5))

  dev.off()
  i <- i + 1
}




MakeRGBPicture <- function(this_date, nir_v, red_v, green_v, landsat_dates, landsat_in_files, modis_dates, modis_in_files, pred_dates, NA_stack, tmp_r, out_img_dir="/projectnb/modislc/users/joshgray/DL_Landsat/MODIS_Landsat_Fusion2", img_prefix="Fusion"){
  print(paste("Doing", this_date))
  # create the outfile name
  i <- which(pred_dates == this_date)
  out_file <- file.path(out_img_dir, paste(img_prefix, "_", formatC(i, width=3, flag="0"), ".jpg", sep=""))
  jpeg(file=out_file, width=(1026*2), height=(671*2), quality=75)
  print(paste("Writing to:", out_file))

  layout(matrix(c(1, 2, 3, 3), nrow=2, byrow=T))
  par(mar=rep(1, 4), oma=rep(0.5, 4))

  # check for and plot landsat data
  if(this_date %in% landsat_dates){
    landsat_s <- stack(landsat_in_files[which(landsat_dates == this_date)])
    plotRGB(landsat_s, 4, 3, 2, scale=1e4, colNA="black", stretch="lin")
  }else{
    plotRGB(NA_stack, 4, 3, 2, scale=1e4, colNA="black", stretch="lin")
  }

  # check for and plot modis data
  if(this_date %in% modis_dates){
    modis_s <- stack(modis_in_files[which(modis_dates == this_date)])
    plotRGB(modis_s, 4, 3, 2, scale=1e4, colNA="black", stretch="lin")
  }else{
    plotRGB(NA_stack, 4, 3, 2, scale=1e4, colNA="black", stretch="lin")
  }

  # create a raster of the KF smoothed output
  nir_tmp <- tmp_r
  values(nir_tmp) <- nir_v[,i]
  red_tmp <- tmp_r
  values(red_tmp) <- red_v[,i]
  green_tmp <- tmp_r
  values(green_tmp) <- green_v[,i]
  tmp_s <- stack(nir_tmp, red_tmp, green_tmp)
  plotRGB(tmp_s, 1, 2, 3, scale=1e4, colNA="black", stretch="lin")
  text(mean(par()$usr[1:2]), par()$usr[3], labels=this_date, col="yellow", cex=3, adj=c(0.5, -0.5))

  dev.off()
}

# this doesn't work for some reason...
cl <- makeCluster(16)
clusterExport(cl, c("MakeRGBPicture"))
tmp_r <- raster(landsat_in_files[1])
NA_stack <- stack(landsat_in_files[1])
NA_stack <- setValues(NA_stack, NA)

trash <- parLapply(cl, pred_dates, MakeRGBPicture, nir_v=nir_v, red_v=red_v, green_v=green_v, landsat_date=landsat_dates, landsat_in_files=landsat_in_files, modis_dates=modis_dates, modis_in_files=modis_in_files, pred_dates=pred_dates, NAstack=NA_stack, tmp_r=tmp_r)

for(this_date in pred_dates){
  print(paste("Working on:", this_date))
  MakeRGBPicture(this_date, nir_v, red_v, green_v, landsat_dates, landsat_in_files, modis_dates, modis_in_files, pred_dates, out_img_dir="/projectnb/modislc/users/joshgray/DL_Landsat/MODIS_Landsat_Fusion")
}

# MakeRGBPicture(pred_dates[125], nir_v, red_v, green_v, landsat_dates, landsat_in_files, modis_dates, modis_in_files, pred_dates, out_img_dir="/projectnb/modislc/users/joshgray/DL_Landsat/MODIS_Landsat_Fusion")


# make into a movie with ffmpeg:
# module load ffmpeg/git-2013-12-18-4a2570f
# ffmpeg -framerate 4 -i FusionImage_%03d.jpg -pix_fmt yuv420p -vf scale="720:trunc(ow/a/2)*2" multisensor_output.mp4
