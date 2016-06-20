#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Functions to aid in implementing image Kalman Filters in R
# Copyright Josh Gray, 2016
# https://raw.github.com/jm-gray/pixel-forge/master/KalmanFiltering/KF_functions.R
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Utility Functions
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#-------------------------------------------------------------------------------
GetSentinelDate <- function(x) as.Date(substr(unlist(strsplit(basename(x), "_"))[8], 1, 8), format="%Y%m%d")

#-------------------------------------------------------------------------------
GetLandsatDate <- function(x) as.Date(substr(basename(x), 10, 16), format="%Y%j")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Functions for creating, implementing, and extracting results from image KF's
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#-------------------------------------------------------------------------------
MakeMultiDLM <- function(num_states=1, sensors=1, time_varying=FALSE){
  # creates a prototype, non-time varying KF (dlm package) for multiple states and
  # multiple sensors. If all system states are observed by all sensors, then the
  # "sensors" parameter should be the integer number of sensors. If not, "sensors"
  # can be a list with length(num_sensors) where the nth element is a vector of
  # state indices observed by the nth sensor. Example:
  # kf1 <- MakeMultiDLM() # single state, single sensor
  # kf2 <- MakeMultiDLM(5, 3) # 5 states all observed by 3 sensors
  # # 3 states, 2 sensors, sensor A sees states 1 & 2, and B sees states 2 & 3:
  # kf3 <- MakeMultiDLM(3, list(1:2, 2:3))
  # 1 state, 1 sensor, but the evolution coefficient (GG) varies randomly in time:
  # kf4 <- MakeMultiDLM(time_varying=TRUE)
  # X(kf4) <- t(rnorm(10))
  # GG(kf4) <- 1

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
    # time-varying parameters not populated here, but stubbed for future expansion
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

#-------------------------------------------------------------------------------
# this seems stupid
ExtractSmoothMeans <- function(x) dropFirst(x$s)

#-------------------------------------------------------------------------------
ExtractSmoothMeans_multi <- function(x, band=1){
  tmp <- dropFirst(x$s)
  return(tmp[, band])
}

#-------------------------------------------------------------------------------
ExtractSmoothSE <- function(x) sqrt(unlist(dropFirst(dlmSvd2var(x$U.S, x$D.S))))

#-------------------------------------------------------------------------------
ApplyKFSmooth_tv <- function(x, dlm, prior_rates, prior_vars, prior_means){
  # time-varying version of KF smoother
  # the first element of x specifies to which "cluster" (e.g. landcover type) the
  # t.s. belongs, and then gets the t.v. parameters from prior_rates and prior_vars
  # which are expected to be lists of sufficient length and size

  cluster_id <- x[1,1] # retrieve the "cluster" ID to find
  y <- x[-1,]
  tmp_dlm <- dlm
  X(tmp_dlm)[,1] <- prior_rates[cluster_id, ]
  X(tmp_dlm)[,2] <- prior_vars[cluster_id, ]
  m0(tmp_dlm) <- prior_means[cluster_id, 1]
  C0(tmp_dlm) <- prior_vars[cluster_id, 1]
  smooth <- dlmSmooth(y, tmp_dlm)
  return(smooth)
}

#-------------------------------------------------------------------------------
ApplyKFSmooth <- function(y, dlm){
  # not time-varying version of KF smoother
  # NOTE: this is pretty stupid...

  smooth <- dlmSmooth(y, dlm)
  return(smooth)
}

#-------------------------------------------------------------------------------
ApplyAnnualKF <- function(y, dlm, period=52){
  tmp_dlm <- dlm
  filt_m <- c()
  filt_se <- c()
  for(i in seq(1, length(y), by=period)){
    tmp_y <- y[i:(i + period - 1)]
    filt <- dlmFilter(tmp_y, tmp_dlm)
    filt_m <- c(filt_m, dropFirst(filt$m))
    filt_se <- c(filt_se, sqrt(unlist(dropFirst(dlmSvd2var(filt$U.C, filt$D.C)))))
    m0(tmp_dlm) <- filt_m[length(filt_m)]
    C0(tmp_dlm) <- filt_se[length(filt_m)]^2
  }
  return(list(filt_m, filt_se))
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Error quantification functions
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#-------------------------------------------------------------------------------
GetKFSmoothError <- function(x, kf, error_ind){
  og_obs <- x[error_ind]
  x[error_ind] <- NA
  smooth <- dlmSmooth(x, kf)
  tmp_error <- og_obs - dropFirst(smooth$s)[error_ind]
  return(tmp_error)
}

#-------------------------------------------------------------------------------
GetAllDayKFError <- function(x, kf){
  tmp_fun <- function(ind, x, kf) GetKFSmoothError(x, kf, ind)
  i <- 1:length(x)
  results <- unlist(lapply(i, tmp_fun, x, kf))
  return(results)
}

#-------------------------------------------------------------------------------
ApplyGetError <- function(ind, Y, kf){
  if(length(dim(Y)) == 3){
    x <- Y[ind,,]
  }else{
    x <- Y[ind,]
  }
  return(GetAllDayKFError(x, kf))
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Functions to simulate SVI time series w/ double-logistic functions
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#------------------------------------------------
doublesigmoid <- function(t, params) return(0.5 * (tanh(((t - params$M) - params$c1) / params$w1) - tanh(((t - params$M) - params$c2) / params$w2)))

#------------------------------------------------
SimSVI <- function(numyears, params, seed=42, quiet=T){
  set.seed(seed)
  svi <- c()
  t <- 1:365
  for(i in 1:numyears){
    M <- rnorm(1, mean=params$midpoint[1], sd=params$midpoint[2])
    c1 <- -0.5 * rnorm(1, mean=params$gup_duration[1], sd=params$gup_duration[2])
    c2 <- 0.5 * rnorm(1, mean=params$gdown_duration[1], sd=params$gdown_duration[2])
    w1 <- rnorm(1, mean=params$gup_rate[1], sd=params$gup_rate[2])
    w2 <- rnorm(1, mean=params$gdown_rate[1], sd=params$gdown_rate[2])
    svi <- c(svi, doublesigmoid(t, list(M=M, c1=c1, c2=c2, w1=abs(w1), w2=abs(w2))))
    if(!quiet) print(c(M=M, c1=c1, c2=c2, w1=w1, w2=w2))
  }
  return(svi)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Functions for visualization of results
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#-------------------------------------------------------------------------------
PlotEVI <- function(r, breaks=c(-1e9, seq(0, 0.75, len=253), 1e9)){
  par(mar=rep(0,4), oma=rep(0,4), bg="black")
  myRamp <- colorRampPalette(brewer.pal(11, "Spectral"))
  plot(r, breaks=breaks, col=myRamp(255), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
}

#-------------------------------------------------------------------------------
PlotFusion <- function(r1, r2, kf, label, evi_breaks=c(-1e9, seq(0, 0.75, len=254), 1e9)){
  nf <- layout(matrix(c(1,2,3,3), nrow=2, byrow=T), heights=c(0.5, 1))
  par(mar=rep(0,4), oma=rep(0,4), bg="black")
  myRamp <- colorRampPalette(brewer.pal(11, "Spectral"))
  image(r1, breaks=evi_breaks, col=myRamp(255), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
  image(r2, breaks=evi_breaks, col=myRamp(255), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
  image(kf, breaks=evi_breaks, col=myRamp(255), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
  text(par()$usr[1], par()$usr[3], label=label, adj=c(-0.1, -0.1), col='black', cex=2)
}

#-------------------------------------------------------------------------------
PlotFusionRGB <- function(r1, r2, kf, label, evi_breaks=c(-1e9, seq(0, 0.75, len=254), 1e9)){
  nf <- layout(matrix(c(1,2,3,3), nrow=2, byrow=T), heights=c(0.5, 1))
  par(mar=rep(0,4), oma=rep(0,4), bg="black")
  # myRamp <- colorRampPalette(brewer.pal(11, "Spectral"))
  plotRGB(r1, 3, 2, 1, colNA="black", stretch="lin")
  plotRGB(r2, 3, 2, 1, colNA="black", stretch="lin")
  plotRGB(kf, 3, 2, 1, colNA="black", stretch="lin")
  # image(r1, breaks=evi_breaks, col=myRamp(255), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
  # image(r2, breaks=evi_breaks, col=myRamp(255), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
  # image(kf, breaks=evi_breaks, col=myRamp(255), legend=FALSE, xaxt="n", yaxt="n", colNA="black")
  text(par()$usr[1], par()$usr[3], label=label, adj=c(-0.1, -0.1), col='black', cex=2)
}

#-------------------------------------------------------------------------------
PlotForecast <- function(filt_m, filt_se, signal, t=NULL, conf_level=0.05, ylim=NULL, ...){
  if(is.null(t)) t <- 1:length(filt_m)
  colmain <- "#3182BD"; colerr <- "#BDD7E7"; colsignal <- "#636363"
  lower <- filt_m + qnorm(conf_level, sd=filt_se)
  upper <- filt_m + qnorm(1 - conf_level, sd=filt_se)

  if(is.null(ylim)){
    ylim <- range(c(upper, lower), na.rm=T) * c(0.9, 1.1)
    if(diff(range(c(upper, lower), na.rm=T)) > (2 * diff(range(filt_m, na.rm=T)))){
      # ylim <- range(filt_m, na.rm=T) * c(0.55, 1.5)
      # ylim <- c(quantile(lower,c(0.1)), quantile(upper,c(0.9)))
      ylim <- range(filt_m, na.rm=T) + sd(filt_m)*c(-1, 1)
    }
  }

  # ylim <- c(min(lower, na.rm=T), max(upper, na.rm=T)) * c(0.9, 1.1)
  plot(t, filt_m, type="n", lwd=2, col=2, lty=2, ylim=ylim, ...)
  polygon(x=c(t, rev(t)), y=c(upper, rev(lower)), border=NA, col=colerr)
  points(t, filt_m, type="l", lwd=2, col=colmain, lty=1)
  if(is.list(signal)){
    for(i in 1:length(signal)){
      tmp_y <- signal[[i]]
      points(t, tmp_y, type="p", pch=i, cex=0.75, col=colsignal)
    }
  }else{
    points(t, signal, type="p", pch=1, cex=0.75, col=colsignal)
  }
}
