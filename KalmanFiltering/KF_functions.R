#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Functions to aid in implementing image Kalman Filters in R
# Josh Gray, 2016
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

#-------------------------------------------------------------------------------
GetModisDate <- function(x) as.Date(unlist(strsplit(basename(x), split="\\."))[2], format="A%Y%j")

#-------------------------------------------------------------------------------
CalcEVI2 <- function(in_file, scale_factor=1e-4){
  s <- stack(in_file)
  nir_v <- values(raster(in_file, 4))
  red_v <- values(raster(in_file, 3))
  evi2 <- 2.5 * (((nir_v * scale_factor) - (red_v * scale_factor)) / ((nir_v * scale_factor) + (2.4 * (red_v * scale_factor)) + 1))
  return(evi2)
}

#-------------------------------------------------------------------------------
GetValuesGDAL <- function(dsets, start_row, n, max_open_datasets=2.75e3) {
	# extracts n rows from all dsets starting at start_row

	# determine how many blocks of datasets to open simultaneously and dataset size
	num_blocks <- ceiling(length(dsets) / max_open_datasets)
	nrows <- nrow(raster(dsets[1]))
	ncols <- ncol(raster(dsets[1]))

	# initialize output matrix
	out_vals <- matrix(NA, nrow=ncols * n, ncol=length(dsets))

	# loop through all datasets in a block and extract data
	for(dset_block in 1:num_blocks){
		# determine which dataset to start and end on
		dset_start <- ((dset_block - 1) * max_open_datasets) + 1
		dset_end <- min((dset_block * max_open_datasets), length(dsets))
		gds <- list()

		# open all datasets
		for (i in dset_start:dset_end) { gds <- c(gds, GDAL.open(dsets[i])) }

		# loop through each dataset and extract rows
		for (j in 1:length(gds)) {
      val <- getRasterData(gds[[j]], offset = c((start_row - 1), 0), region.dim = c(n, ncols), as.is = FALSE)
			out_vals[, (dset_start + j - 1)] <- c(val)
		}
		# close all files
		for (i in 1:length(gds)) { GDAL.close(gds[[i]]) }
	}
	return(out_vals)
}

#-------------------------------------------------------------------------------
GetSplineRMSE <- function(x, x_dates, y_dates, pred_dates=NULL){
	# fits a daily smoothing spline to time series x and y (landsat and modis)
	# and then calculates the RMSE
	x_v <- x[1:length(x_dates)]
	y_v <- x[(length(x_dates) + 1):(length(x_dates) + length(y_dates))]
	y_snow_v <- x[(length(x_dates) + length(y_dates) + 1):length(x)]
  y_v[y_snow_v == 1] <- NA # screen snow obs out of y (MODIS)

	if(is.null(pred_dates)) pred_dates <- seq(min(x_dates, y_dates, na.rm=T), max(x_dates, y_dates, na.rm=T), by=1)

  x_smooth <- try(predict(smooth.spline(x_dates[!is.na(x_dates) & !is.na(x_v)], x_v[!is.na(x_dates) & !is.na(x_v)]), as.numeric(pred_dates))$y, silent=T)
  if(inherits(x_smooth, 'try-error')) return(NA)

	y_smooth <- try(predict(smooth.spline(y_dates[!is.na(y_dates) & !is.na(y_v)], y_v[!is.na(y_dates) & !is.na(y_v)]), as.numeric(pred_dates))$y, silent=T)
  if(inherits(y_smooth, 'try-error')) return(NA)

	return(sqrt(mean((x_smooth - y_smooth)^2, na.rm=T)))
}

#-------------------------------------------------------------------------------
GetMatchRMSE <- function(x, x_dates, y_dates){
	# calculates RMSE of matching dates in x/y
  x_match_inds <- x_dates %in% y_dates
  y_match_inds <- y_dates %in% x_dates
	x_v <- x[1:length(x_dates)]
	y_v <- x[(length(x_dates) + 1):(length(x_dates) + length(y_dates))]
	y_snow_v <- x[(length(x_dates) + length(y_dates) + 1):length(x)]
  y_v[y_snow_v == 1] <- NA # screen snow obs out of y (MODIS)
  x_match <- x_v[x_match_inds]
  y_match <- y_v[y_match_inds]
	return(sqrt(mean((x_match - y_match)^2, na.rm=T)))
}

#-------------------------------------------------------------------------------
GetBandBRDFQualities <- function(x){
	# represent the binary conversion of x as a 32 x N matrix of values
	M <- matrix(as.integer(intToBits(t(x))), ncol=32, byrow=T)
	M <- M[, 32:1] # reverse column order

	# get the individual bit-words and convert to decimal
	QA_Fill <- M[1]
	NOT_USED <- M[2:4]
	B7QA <- sum(M[5:8] * t(c(8,4,2,1)))
	B6QA <- sum(M[9:12] * t(c(8,4,2,1)))
	B5QA <- sum(M[13:16] * t(c(8,4,2,1)))
	B4QA <- sum(M[17:20] * t(c(8,4,2,1)))
	B3QA <- sum(M[21:24] * t(c(8,4,2,1)))
	B2QA <- sum(M[25:28] * t(c(8,4,2,1)))
	B1QA <- sum(M[29:32] * t(c(8,4,2,1)))

	return(c(QA_Fill, B1QA, B2QA, B3QA, B4QA, B5QA, B6QA, B7QA))
}

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
