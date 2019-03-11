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
GetValuesGDAL_multiband <- function(dsets, start_row, n, max_open_datasets=2.75e3) {
	# extracts n rows from all dsets starting at start_row

	# determine how many blocks of datasets to open simultaneously and dataset size
	num_blocks <- ceiling(length(dsets) / max_open_datasets)
  tmp_r <- stack(dsets[1])
	nrows <- nrow(tmp_r)
	ncols <- ncol(tmp_r)
  nbands <- nlayers(tmp_r)

	# initialize output matrix
  out_vals <- array(NA, dim=c(ncols * n, length(dsets), nbands))

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
      val <- getRasterData(gds[[j]], offset = c((start_row - 1), 0), region.dim = c(n, ncols), as.is = TRUE)
      val <- matrix(val, nrow=ncols * n, ncol=nbands)
      out_vals[, (dset_start + j - 1), ] <- c(val)
		}
		# close all files
		for (i in 1:length(gds)) { GDAL.close(gds[[i]]) }
	}
	return(out_vals)
}


#-------------------------------------------------------------------------------
GetSplineRMSE <- function(x, x_dates, y_dates, pred_dates=NULL, num_cdl_years=4){
	# fits a daily smoothing spline to time series x and y (landsat and modis)
	# and then calculates the RMSE
  x_v <- x[(num_cdl_years + 1):(length(x_dates) + num_cdl_years)]
  y_v <- x[(length(x_dates) + num_cdl_years + 1):(length(x_dates) + num_cdl_years + length(y_dates))]
	if(is.null(pred_dates)) pred_dates <- seq(min(x_dates, y_dates, na.rm=T), max(x_dates, y_dates, na.rm=T), by=1)
  x_smooth <- try(predict(smooth.spline(x_dates[!is.na(x_dates) & !is.na(x_v)], x_v[!is.na(x_dates) & !is.na(x_v)]), as.numeric(pred_dates))$y, silent=T)
  if(inherits(x_smooth, 'try-error')) return(NA)
	y_smooth <- try(predict(smooth.spline(y_dates[!is.na(y_dates) & !is.na(y_v)], y_v[!is.na(y_dates) & !is.na(y_v)]), as.numeric(pred_dates))$y, silent=T)
  if(inherits(y_smooth, 'try-error')) return(NA)
	return(sqrt(mean((x_smooth - y_smooth)^2, na.rm=T)))
}

#-------------------------------------------------------------------------------
GetMatchRMSE <- function(x, x_dates, y_dates, num_cdl_years=4){
	# calculates RMSE of matching dates in x/y
  if(all(is.na(x))) return(NA) # check for special case of no data
  x_match_inds <- x_dates %in% y_dates
  y_match_inds <- y_dates %in% x_dates
  x_v <- x[(num_cdl_years + 1):(length(x_dates) + num_cdl_years)]
  y_v <- x[(length(x_dates) + num_cdl_years + 1):(length(x_dates) + num_cdl_years + length(y_dates))]
  x_match <- x_v[x_match_inds]
  y_match <- y_v[y_match_inds]
  rmse <- sqrt(mean((x_match - y_match)^2, na.rm=T))
  if(is.nan(rmse)){
    return(NA)
  }else{
    return(rmse)
  }
}

#-------------------------------------------------------------------------------
EVI2_error <- function(ref_red, ref_nir, u_red, u_nir){
  # uses the law of prop of uncertainty to calculate the EVI2 error
  # assumes independent errors in red and nir reflectance
  # Markham and Helder 2012 give TM and ETM absolute calibration errors as 7 and 5%, rep (one sigma)
  G <- 2.5; L <- 1; C <- 2.4
  del_EVI2_red <- ((-1 * G) * (L + (ref_nir * (1 + C)))) / ((L + ref_nir + (C * ref_red))^2)
  del_EVI2_nir <- (G * (L + (ref_red * (1 + C)))) / ((L + ref_nir + (C * ref_red))^2)
  return(sqrt(((del_EVI2_nir^2) * ((u_nir * ref_nir)^2)) + ((del_EVI2_red^2) * ((u_red * ref_red)^2))))
}

#-------------------------------------------------------------------------------
GetLandsatMODISError <- function(x, x_dates, y_dates, num_cdl_years=4){
	# calculates RMSE of matching dates in x/y
  return_NA <- list(intercept=NA,slope=NA,rmse=NA)
  if(all(is.na(x))) return(return_NA) # check for special case of no data
  x_match_inds <- x_dates %in% y_dates
  y_match_inds <- y_dates %in% x_dates
  x_v <- x[(num_cdl_years + 1):(length(x_dates) + num_cdl_years)]
  y_v <- x[(length(x_dates) + num_cdl_years + 1):(length(x_dates) + num_cdl_years + length(y_dates))]
  x_match <- x_v[x_match_inds]
  y_match <- y_v[y_match_inds]
  rmse <- sqrt(mean((x_match - y_match)^2, na.rm=T))
  lm1 <- try(lm(x_match ~ y_match))
  if(inherits(lm1, 'try-error')) return(return_NA)
  if(is.nan(rmse)){
    return(return_NA)
  }else{
    return(list(intercept=lm1$coefficients[1], slope=lm1$coefficients[2], rmse=rmse))
  }
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

#-------------------------------------------------------------------------------
GetDOYSpline <- function(x, dates, min_quant=0.05, spline_spar=NULL, head_fill_len=NA, tail_fill_len=NA){
  if(all(is.na(x))) return(rep(NA, 365))
	# calculates an annual spline fit to the data in x
	# values below the quantile specified by min_quant are filled
	doys <- as.numeric(strftime(dates, format="%j"))
	years <- as.numeric(strftime(dates, format="%Y"))
	# calculate background minimum value as quantile
	qmin <- quantile(x, min_quant, na.rm=T)
	# x[x < qmin] <- NA
  x[x < qmin] <- qmin
	# add background value to head/tail from first/last non-NA value to constrain spline fit
	start_doy <- min(doys[!is.na(x)], na.rm=T)
	end_doy <- max(doys[!is.na(x)], na.rm=T)
	# doys <- c(doys, 1:(start_doy - 1), (end_doy + 1):365)
  doys <- c(doys, seq(from=1, to=start_doy - 1, len=start_doy - 1), seq(from=end_doy + 1, to=365, len=365 - end_doy))
	x <- c(x, rep(qmin, start_doy - 1), rep(qmin, 365 - end_doy))
  # if head/tail fills lengths are supplied, fill w/ background value:
  if(!is.na(head_fill_len)){
    doys <- c(doys, 1:head_fill_len)
    x <- c(x, rep(qmin, head_fill_len))
  }
  if(!is.na(tail_fill_len)){
    doys <- c(doys, (366 - tail_fill_len):365)
    x <- c(x, rep(qmin, tail_fill_len))
  }
	# fit smoothing spline
	x_smooth <- predict(smooth.spline(doys[!is.na(x)], x[!is.na(x)], spar=spline_spar), 1:365)$y
	return(x_smooth)
}

#-------------------------------------------------------------------------------
GetLandsatProcessSD <- function(x, dates){
  # fits a spline to each year's worth of daily Landsat data, returns a matrix of the diffs
  # with as many rows as there are years
  tmp_years <- as.integer(strftime(dates, format="%Y"))
  unique_years <- unique(tmp_years)
  for(i in 1:length(unique_years)){
    tmp_spline <- GetDOYSpline(x[tmp_years == unique_years[i]], dates[tmp_years == unique_years[i]])
    if(i == 1){
      out_val <- diff(tmp_spline)
    }else{
      out_val <- rbind(out_val, diff(tmp_spline))
    }
  }
  return(c(t(out_val)))
}

#-------------------------------------------------------------------------------
GetLandsatMultispecProcessDiff <- function(x, dates, k=1.5){
  # fits a spline to each year's worth of daily Landsat data, returns a matrix of the diffs
  # with as many rows as there are years
  tmp_years <- as.integer(strftime(dates, format="%Y"))
  unique_years <- unique(tmp_years)
  x <- TukeyOutlierRemoval(x, k=k)
  for(i in 1:length(unique_years)){
    tmp_x <- x[tmp_years == unique_years[i]]
    tmp_dates <- dates[tmp_years == unique_years[i]]
    tmp_doys <- as.integer(strftime(tmp_dates, format="%j"))
    # # eliminate outliers
    # qs <- quantile(tmp_x, na.rm=T)
    # tukey_range <- c(qs[2] - k * (qs[4] - qs[2]), qs[4] + k * (qs[4] - qs[2]))
    # tmp_x[tmp_x < tukey_range[1] | tmp_x > tukey_range[2]] <- NA

    if(length(tmp_x[!is.na(tmp_x)]) < 4){
      # if there's not enough data to fit the spline we just use NA
      tmp_spline <- rep(NA, 365)
    }else{
      tmp_sp <- smooth.spline(tmp_doys[!is.na(tmp_x)], tmp_x[!is.na(tmp_x)], spar=NULL)
      tmp_spline <- predict(tmp_sp, 1:365)$y
    }
    # tmp_spline <- GetDOYSpline(x[tmp_years == unique_years[i]], dates[tmp_years == unique_years[i]])
    if(i == 1){
      out_val <- diff(tmp_spline)
    }else{
      out_val <- rbind(out_val, diff(tmp_spline))
    }
  }
  return(c(t(out_val)))
}

#-------------------------------------------------------------------------------
TukeyRestrictedSD <- function(x, k=1.5){
  # returns the std dev of x, with Tukey outliers screened
  qs <- quantile(x, na.rm=T)
  tukey_range <- c(qs[2] - k * (qs[4] - qs[2]), qs[4] + k * (qs[4] - qs[2]))
  return(sd(x[x >= tukey_range[1] & x <= tukey_range[2]], na.rm=T))
}

#-------------------------------------------------------------------------------
TukeyRestrictedCOV <- function(a, b, c, d, k=1.5, retCOV=T){
  # returns cov(a, b, c, d), with Tukey outliers screened
  # cov(cbind(cdl_process_diff_blue[1,], cdl_process_diff_green[1,], cdl_process_diff_red[1,], cdl_process_diff_nir[1,]))/1e4
  qs_a <- quantile(a, na.rm=T)
  tukey_range_a <- c(qs_a[2] - k * (qs_a[4] - qs_a[2]), qs_a[4] + k * (qs_a[4] - qs_a[2]))
  a_good_inds <- a >= tukey_range_a[1] & a <= tukey_range_a[2]

  qs_b <- quantile(b, na.rm=T)
  tukey_range_b <- c(qs_b[2] - k * (qs_b[4] - qs_b[2]), qs_b[4] + k * (qs_b[4] - qs_b[2]))
  b_good_inds <- b >= tukey_range_b[1] & b <= tukey_range_b[2]

  qs_c <- quantile(c, na.rm=T)
  tukey_range_c <- c(qs_c[2] - k * (qs_c[4] - qs_c[2]), qs_c[4] + k * (qs_c[4] - qs_c[2]))
  c_good_inds <- c >= tukey_range_c[1] & c <= tukey_range_c[2]

  qs_d <- quantile(d, na.rm=T)
  tukey_range_d <- c(qs_d[2] - k * (qs_d[4] - qs_d[2]), qs_d[4] + k * (qs_d[4] - qs_d[2]))
  d_good_inds <- d >= tukey_range_d[1] & d <= tukey_range_d[2]

  all_good_inds <- a_good_inds & b_good_inds & c_good_inds & d_good_inds
  if(retCOV){
    return(cov(cbind(a[all_good_inds], b[all_good_inds], c[all_good_inds], d[all_good_inds]), use="complete"))
  }else{
    return(cor(cbind(a[all_good_inds], b[all_good_inds], c[all_good_inds], d[all_good_inds]), use="complete"))
  }
  # i <- 34
  # a <- cdl_process_diff_blue[,i]
  # b <- cdl_process_diff_green[,i]
  # c <- cdl_process_diff_red[,i]
  # d <- cdl_process_diff_nir[,i]
  # cov(cbind(a, b, c, d))
  # TukeyRestrictedCOV(a,b,c,d)
  # cor(cbind(a, b, c, d))
  # TukeyRestrictedCOV(a,b,c,d, retCOV=F)
}


#-------------------------------------------------------------------------------
QuantileMinReplacement <- function(x, qval=0.05){
  replacement_value <- quantile(x[x >= 0], qval, na.rm=T)
  x[x < 0] <- replacement_value
  return(x)
}

#-------------------------------------------------------------------------------
TukeyOutlierRemoval <- function(x, k=1.5){
  qs <- quantile(x, na.rm=T)
  tukey_range <- c(qs[2] - k * (qs[4] - qs[2]), qs[4] + k * (qs[4] - qs[2]))
  x[x < tukey_range[1] | x > tukey_range[2]] <- NA
  return(x)
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

# #-------------------------------------------------------------------------------
# FuseLandsatModisEVI <- function(x, landsat_dates, modis_dates, landsat_sensor, cdl_tv_sd, cdl_types, scale_factor=1e4, modis_landsat_slope=1, modis_landsat_bias=0, smooth=T, plot=F, ...){
#   # extract components of x
#   tmp_cdl <- x[1]
#   tmp_rmse <- x[2] / scale_factor
#   x <- x[c(-1, -2)] / scale_factor
#   x_landsat <- x[1:length(landsat_dates)]
#   x_modis <- x[(length(landsat_dates) + 1):length(x)]
#
#   # retrieve the proper time varying process error for this land cover type
#   tv_sd <- cdl_tv_sd[[which(cdl_types == tmp_cdl)]]$splined / scale_factor
#   tv_sd <- c(tv_sd[1], tv_sd) # append the head value b/c sd is 364 long
#
#   # munge to daily series
#   num_years <- length(as.integer(sort(unique(c(strftime(landsat_dates, format="%Y"), strftime(modis_dates, format="%Y"))))))
#
#   # do landsat; eliminate leap year day 366
#   tmp_landsat <- rep(NA, num_years * 365)
#   tmp_landsat_doys <- as.integer(strftime(landsat_dates, format="%j"))
#   tmp_landsat_years <- as.integer(strftime(landsat_dates, format="%Y"))
#   x_landsat <- x_landsat[tmp_landsat_doys <= 365]
#   tmp_landsat_doys <- tmp_landsat_doys[tmp_landsat_doys <= 365]
#   tmp_landsat_years <- tmp_landsat_years[tmp_landsat_doys <= 365]
#   daily_landsat_inds <- tmp_landsat_doys + 365 * (tmp_landsat_years - min(tmp_landsat_years))
#   tmp_landsat[daily_landsat_inds] <- x_landsat
#   daily_landsat_sensor <- rep(NA, length(tmp_landsat))
#   daily_landsat_sensor[daily_landsat_inds] <- landsat_sensor
#
#   # do modis; eliminate leap year day 366
#   tmp_modis <- rep(NA, num_years * 365)
#   tmp_modis_doys <- as.integer(strftime(modis_dates, format="%j"))
#   tmp_modis_years <- as.integer(strftime(modis_dates, format="%Y"))
#   x_modis <- x_modis[tmp_modis_doys <= 365]
#   tmp_modis_doys <- tmp_modis_doys[tmp_modis_doys <= 365]
#   tmp_modis_years <- tmp_modis_years[tmp_modis_doys <= 365]
#   daily_modis_inds <- tmp_modis_doys + 365 * (tmp_modis_years - min(tmp_modis_years))
#   tmp_modis[daily_modis_inds] <- x_modis
#   tmp_modis <- tmp_modis - modis_landsat_bias # remove any bias in the MODIS sensor
#
#   # create obs matrix
#   y <- cbind(tmp_landsat, tmp_modis)
#
#   # specify landsat TM and ETM+ errors (7% and 5%, resp)
#   tmp_landsat_error <- rep(NA, length(tmp_landsat))
#   tmp_landsat_error[daily_landsat_sensor == "LE7" & !is.na(daily_landsat_sensor)] <- tmp_landsat[daily_landsat_sensor == "LE7" & !is.na(daily_landsat_sensor)] * 0.05
#   tmp_landsat_error[daily_landsat_sensor == "LT5" & !is.na(daily_landsat_sensor)] <- tmp_landsat[daily_landsat_sensor == "LT5" & !is.na(daily_landsat_sensor)] * 0.07
#
#   # # replace all missing landsat error values with closest not-NA value
#   # trash <- sapply(which(is.na(tmp_landsat_error)), function(a, x) x[which.min(abs(x-a))], x=which(!is.na(tmp_landsat_error)))
#   # tmp_landsat_error[which(is.na(tmp_landsat_error))] <- tmp_landsat_error[trash]
#
#   tmp_modis_error <- rep(NA, length(tmp_modis))
#   tmp_modis_error[!is.na(tmp_modis)] <- tmp_rmse
#
#   # define dlm components
#   GG <- matrix(1) # process transition
#   W <- matrix(1) # evolution error covariance
#   JW <- matrix(1) # time varying evolution error covariance: in col 1 of X
#   # FF <- matrix(c(1, 1), nrow=2) # observation matrix
#   FF <- matrix(c(1, modis_landsat_slope), nrow=2) # observation matrix
#   V <- matrix(c(1, 0, 0, 1), nrow=2) # obs uncertainty
#   JV <- matrix(c(2, 0, 0, 3), nrow=2) # time varying obs uncertainty: in cols 2 (landsat) and 3 (modis)
#   X <- matrix(1, nrow=length(tmp_modis), ncol=3)
#   X[, 1] <- rep(tv_sd, num_years) # evolution error
#   X[, 2] <- tmp_landsat_error # landsat observation error
#   # X[, 3] <- rep(tmp_rmse, length(tmp_modis)) # modis observation error
#   X[, 3] <- tmp_modis_error # modis observation error
#   m0 <- 0 # initial state vector
#   C0 <- 1 # initial state uncertainty
#
#   # construct the dlm
#   tmp_dlm <- dlm(m0=m0, C0=C0, GG=GG, JW=JW, FF=FF, V=V, JV=JV, W=W, X=X)
#
#   # apply the dlm
#   if(smooth){
#     kf_result <- dlmSmooth(y, tmp_dlm)
#     kf_evi <- dropFirst(kf_result$s)
#     kf_evi_error <- sqrt(unlist(dropFirst(dlmSvd2var(kf_result$U.S, kf_result$D.S))))
#     ret_value <- list(evi=kf_evi, error=kf_evi_error)
#   }else{
#     kf_result <- dlmFilter(y, tmp_dlm)
#     kf_evi <- dropFirst(kf_result$m)
#     kf_evi_error <- sqrt(unlist(dropFirst(dlmSvd2var(kf_result$U.C, kf_result$D.C))))
#     ret_value <- list(evi=kf_evi, error=kf_evi_error)
#   }
#
#   if(plot){
#     # PlotForecast(kf_evi, kf_evi_error, split(y, col(y)), ylab="EVI2", xlab="", sigma=1, ...)
#     PlotForecast(kf_evi, kf_evi_error, split(y, col(y)), ylab="EVI2", xlab="", ...)
#   }
#
#   return(ret_value)
# }

#-------------------------------------------------------------------------------
FuseLandsatModisEVI <- function(x, landsat_dates, modis_dates, landsat_sensor, cdl_tv_sd, cdl_types, scale_factor=1e4, modis_landsat_slope=1, modis_landsat_bias=0, smooth=T, plot=F, ...){
  # extract components of x
  tmp_cdl <- x[1]
  tmp_rmse <- x[2] / scale_factor
  x <- x[c(-1, -2)]
  x_landsat <- x[1:length(landsat_dates)] / scale_factor
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

#-------------------------------------------------------------------------------
FuseLandsatModisMultispec <- function(x, landsat_dates, modis_dates, landsat_sensor, multispectral_cdl_process_cov, cdl_types, ret_error=F, scale_factor=1e4, modis_landsat_slope=1, modis_landsat_bias=0, smooth=T, plot=F, ...){
  # extract components of x
  tmp_cdl <- x[1]
  tmp_rmse_blue <- x[2] / scale_factor
  tmp_rmse_green <- x[3] / scale_factor
  tmp_rmse_red <- x[4] / scale_factor
  tmp_rmse_nir <- x[5] / scale_factor
  x <- x[c(-1, -2, -3, -4, -5)]

  # extract individual landsat bands
  landsat_inds <- seq(1, (length(landsat_dates) * 4) + 1, by=length(landsat_dates))
  x_landsat_blue <- x[landsat_inds[1]:(landsat_inds[2] - 1)] / scale_factor
  x_landsat_green <- x[landsat_inds[2]:(landsat_inds[3] - 1)] / scale_factor
  x_landsat_red <- x[landsat_inds[3]:(landsat_inds[4] - 1)] / scale_factor
  x_landsat_nir <- x[landsat_inds[4]:(landsat_inds[5] - 1)] / scale_factor

  x_error_landsat_blue <- rep(NA, length(x_landsat_blue))
  x_error_landsat_blue[landsat_sensor == "LT5"] <- x_landsat_blue[landsat_sensor == "LT5"] * 0.07
  x_error_landsat_blue[landsat_sensor == "LE7"] <- x_landsat_blue[landsat_sensor == "LE7"] * 0.05
  x_error_landsat_green <- rep(NA, length(x_landsat_green))
  x_error_landsat_green[landsat_sensor == "LT5"] <- x_landsat_green[landsat_sensor == "LT5"] * 0.07
  x_error_landsat_green[landsat_sensor == "LE7"] <- x_landsat_green[landsat_sensor == "LE7"] * 0.05
  x_error_landsat_red <- rep(NA, length(x_landsat_red))
  x_error_landsat_red[landsat_sensor == "LT5"] <- x_landsat_red[landsat_sensor == "LT5"] * 0.07
  x_error_landsat_red[landsat_sensor == "LE7"] <- x_landsat_red[landsat_sensor == "LE7"] * 0.05
  x_error_landsat_nir <- rep(NA, length(x_landsat_nir))
  x_error_landsat_nir[landsat_sensor == "LT5"] <- x_landsat_nir[landsat_sensor == "LT5"] * 0.07
  x_error_landsat_nir[landsat_sensor == "LE7"] <- x_landsat_nir[landsat_sensor == "LE7"] * 0.05

  modis_inds <- seq(landsat_inds[5], landsat_inds[5] + length(modis_dates) * 4, by=length(modis_dates))
  x_modis_blue <- x[modis_inds[1]:(modis_inds[2] - 1)] / scale_factor
  x_modis_green <- x[modis_inds[2]:(modis_inds[3] - 1)] / scale_factor
  x_modis_red <- x[modis_inds[3]:(modis_inds[4] - 1)] / scale_factor
  x_modis_nir <- x[modis_inds[4]:(modis_inds[5] - 1)] / scale_factor

  # retrieve the proper time varying process error for this land cover type
  # tmp_cov <- multispectral_cdl_process_cov[[which(cdl_types == tmp_cdl)]] / scale_factor
  this_cov <- array(NA, dim=c(4, 4, 365))
  # this_cov[,,2:365] <- tmp_cov
	this_cov[,,2:365] <- multispectral_cdl_process_cov[[which(cdl_types == tmp_cdl)]] / scale_factor
  this_cov[,,1] <- this_cov[,,2]


  # munge to daily series
  num_years <- length(as.integer(sort(unique(c(strftime(landsat_dates, format="%Y"), strftime(modis_dates, format="%Y"))))))
  tmp_landsat_doys <- as.integer(strftime(landsat_dates, format="%j"))
  tmp_landsat_years <- as.integer(strftime(landsat_dates, format="%Y"))

  # do landsat; eliminate leap year day 366
  tmp_landsat_blue <- tmp_landsat_green <- tmp_landsat_red <- tmp_landsat_nir <- tmp_error_landsat_blue <- tmp_error_landsat_green <- tmp_error_landsat_red <- tmp_error_landsat_nir <- rep(NA, num_years * 365)
  x_landsat_blue <- x_landsat_blue[tmp_landsat_doys <= 365]
  x_error_landsat_blue <- x_error_landsat_blue[tmp_landsat_doys <= 365]
  x_landsat_green <- x_landsat_green[tmp_landsat_doys <= 365]
  x_error_landsat_green <- x_error_landsat_green[tmp_landsat_doys <= 365]
  x_landsat_red <- x_landsat_red[tmp_landsat_doys <= 365]
  x_error_landsat_red <- x_error_landsat_red[tmp_landsat_doys <= 365]
  x_landsat_nir <- x_landsat_nir[tmp_landsat_doys <= 365]
  x_error_landsat_nir <- x_error_landsat_nir[tmp_landsat_doys <= 365]
  tmp_landsat_doys <- tmp_landsat_doys[tmp_landsat_doys <= 365]
  tmp_landsat_years <- tmp_landsat_years[tmp_landsat_doys <= 365]

  # make into daily, gappy time series
  daily_landsat_inds <- tmp_landsat_doys + 365 * (tmp_landsat_years - min(tmp_landsat_years))
  tmp_landsat_blue[daily_landsat_inds] <- x_landsat_blue
  tmp_error_landsat_blue[daily_landsat_inds] <- x_error_landsat_blue
  tmp_landsat_green[daily_landsat_inds] <- x_landsat_green
  tmp_error_landsat_green[daily_landsat_inds] <- x_error_landsat_green
  tmp_landsat_red[daily_landsat_inds] <- x_landsat_red
  tmp_error_landsat_red[daily_landsat_inds] <- x_error_landsat_red
  tmp_landsat_nir[daily_landsat_inds] <- x_landsat_nir
  tmp_error_landsat_nir[daily_landsat_inds] <- x_error_landsat_nir

  daily_landsat_sensor <- rep(NA, length(tmp_landsat_blue))
  daily_landsat_sensor[daily_landsat_inds] <- landsat_sensor

  # do modis; eliminate leap year day 366
  tmp_modis_blue <- tmp_modis_green <- tmp_modis_red <- tmp_modis_nir <- rep(NA, num_years * 365)
  tmp_modis_doys <- as.integer(strftime(modis_dates, format="%j"))
  tmp_modis_years <- as.integer(strftime(modis_dates, format="%Y"))

  x_modis_blue <- x_modis_blue[tmp_modis_doys <= 365]
  x_modis_green <- x_modis_green[tmp_modis_doys <= 365]
  x_modis_red <- x_modis_red[tmp_modis_doys <= 365]
  x_modis_nir <- x_modis_nir[tmp_modis_doys <= 365]
  tmp_modis_doys <- tmp_modis_doys[tmp_modis_doys <= 365]
  tmp_modis_years <- tmp_modis_years[tmp_modis_doys <= 365]

	# munge daily MODIS reflectance and error vectors
  daily_modis_inds <- tmp_modis_doys + 365 * (tmp_modis_years - min(tmp_modis_years))
  tmp_modis_blue[daily_modis_inds] <- x_modis_blue
  tmp_modis_green[daily_modis_inds] <- x_modis_green
  tmp_modis_red[daily_modis_inds] <- x_modis_red
  tmp_modis_nir[daily_modis_inds] <- x_modis_nir
	tmp_error_modis_blue <- rep(NA, length(tmp_modis_blue))
	tmp_error_modis_green <- rep(NA, length(tmp_modis_green))
	tmp_error_modis_red <- rep(NA, length(tmp_modis_red))
	tmp_error_modis_nir <- rep(NA, length(tmp_modis_nir))
	tmp_error_modis_blue[daily_modis_inds] <- tmp_rmse_blue
	tmp_error_modis_green[daily_modis_inds] <- tmp_rmse_green
	tmp_error_modis_red[daily_modis_inds] <- tmp_rmse_red
	tmp_error_modis_nir[daily_modis_inds] <- tmp_rmse_nir

  # create obs matrix
  # y <- cbind(tmp_landsat, tmp_modis)
  y <- cbind(tmp_landsat_blue, tmp_landsat_green, tmp_landsat_green, tmp_landsat_nir, tmp_modis_blue, tmp_modis_green, tmp_modis_red, tmp_modis_nir)

  # # replace all missing landsat error values with closest not-NA value
  # trash <- sapply(which(is.na(tmp_landsat_error)), function(a, x) x[which.min(abs(x-a))], x=which(!is.na(tmp_landsat_error)))
  # tmp_landsat_error[which(is.na(tmp_landsat_error))] <- tmp_landsat_error[trash]

  # define dlm components
  GG <- diag(4) # process transition
  W <- diag(4) # evolution error covariance
	JW <- diag(4) # time varying evolution error covariance: in cols 1:10 of X
	JW[lower.tri(JW, diag=T)] <- 1:10
	JW[upper.tri(JW)] <- t(JW)[upper.tri(JW)]

	# NOTE: fix X dims here!
	X <- matrix(1, nrow=length(tmp_modis_blue), ncol=18)
	X[, 1] <- rep(this_cov[1,1,], num_years)
	X[, 2] <- rep(this_cov[2,1,], num_years)
	X[, 3] <- rep(this_cov[3,1,], num_years)
	X[, 4] <- rep(this_cov[4,1,], num_years)
	X[, 5] <- rep(this_cov[2,2,], num_years)
	X[, 6] <- rep(this_cov[2,3,], num_years)
	X[, 7] <- rep(this_cov[2,4,], num_years)
	X[, 8] <- rep(this_cov[3,3,], num_years)
	X[, 9] <- rep(this_cov[3,4,], num_years)
	X[, 10] <- rep(this_cov[4,4,], num_years)
	# 1    2    3    4
	# 2    5    6    7
	# 3    6    8    9
	# 4    7    9   10
	FF <- rbind(diag(4), diag(4))
  V <- matrix(c(1, 0, 0, 1), nrow=2) # obs uncertainty
	V <- diag(8)
	JV <- diag(8) * 11:18 # time varying obs uncertainty: in cols 11 through 18
	X[, 11] <- tmp_error_landsat_blue # landsat blue error
	X[, 12] <- tmp_error_landsat_green # landsat green error
	X[, 13] <- tmp_error_landsat_red # landsat red error
	X[, 14] <- tmp_error_landsat_nir # landsat nir error
	X[, 15] <- tmp_error_modis_blue # modis blue error
	X[, 16] <- tmp_error_modis_green # modis green error
	X[, 17] <- tmp_error_modis_red # modis red error
	X[, 18] <- tmp_error_modis_nir # modis nir error

  m0 <- matrix(rep(0, 4), nrow=4) # initial state vector
  C0 <- diag(4) * 0.1 # initial state uncertainty

  # construct the dlm
  tmp_dlm <- dlm(m0=m0, C0=C0, GG=GG, JW=JW, FF=FF, V=V, JV=JV, W=W, X=X)

  # apply the dlm
  if(smooth){
    kf_result <- dlmSmooth(y, tmp_dlm)
    kf_surf_ref <- dropFirst(kf_result$s)
		if(ret_error){
			kf_surf_ref_error <- dropFirst(dlmSvd2var(kf_result$U.S, kf_result$D.S))
			kf_surf_ref_error <- sapply(kf_surf_ref_error, diag)
			ret_value <- list(evi=c(kf_surf_ref), error=c(t(kf_surf_ref_error)))
		}else{
			ret_value <- c(kf_surf_ref)
		}
  }else{
    kf_result <- dlmFilter(y, tmp_dlm)
    kf_surf_ref <- dropFirst(kf_result$m)
    kf_surf_ref_error <- sqrt(unlist(dropFirst(dlmSvd2var(kf_result$U.C, kf_result$D.C))))
    ret_value <- list(evi=kf_surf_ref, error=kf_surf_ref_error)
  }

  if(plot){
    # PlotForecast(kf_evi, kf_evi_error, split(y, col(y)), ylab="EVI2", xlab="", sigma=1, ...)
    # PlotForecast(kf_evi, kf_evi_error, split(y, col(y)), ylab="EVI2", xlab="", ...)
  }

  return(ret_value)
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

#-------------------------------------------------------------------------------
# Quantify prediction/interpolation error for any number of consecutive missing Landsat values
# Note: they are not consecutive, but consecutive not missing. That is, for a series: 1, NA, NA, NA, 2, NA, 5
# if we wanted two "consecutive" not missing values, the series would be: NA, NA, NA, NA, NA, NA, 5 with error
# quantified for indices 1 and 5
GetErrorLandsatMODISFusion <- function(x, landsat_dates, modis_dates, landsat_sensor, cdl_tv_sd, cdl_types, scale_factor=1e4, smooth=T, consecutive_missing=1, missing_fraction=NULL, missing_iterations=10){
    not_na_landsat_indices <- which(!is.na(x[3:(2 + length(landsat_dates))])) + 2
    max_landsat_index <- 2 + length(landsat_dates)
    output_errors <- NULL
    kf_dates <- sort(as.Date(paste(rep(unique(strftime(c(landsat_dates, modis_dates), format="%Y")), 365), 1:365, sep="-"), format="%Y-%j"))
    if(!is.null(missing_fraction)){
      for(i in 1:missing_iterations){
        tmp_errors <- rep(NA, length(landsat_dates))
        missing_indices <- sort(sample(not_na_landsat_indices, round(length(not_na_landsat_indices) * missing_fraction)))
        missing_dates <- landsat_dates[missing_indices - 2]
        tmp_x <- x
        # tmp_x[2:length(tmp_x)] <- tmp_x[2:length(tmp_x)] / scale_factor
        og_evi <- tmp_x[missing_indices] / scale_factor
        tmp_x[missing_indices] <- NA
        kf_result <- FuseLandsatModisEVI(tmp_x, landsat_dates, modis_dates, landsat_sensor, cdl_tv_sd, cdl_types, scale_factor=1e4, smooth=T, plot=F)
        kf_evi <- kf_result$evi[sapply(kf_dates, function(x) any(x==missing_dates))]
        tmp_error <- og_evi - kf_evi
        tmp_errors[missing_indices - 2] <- tmp_error
        # append the cdl flag and RMSE
        tmp_errors <- c(x[1:2], tmp_errors)
        if(is.null(output_errors)){
          output_errors <- tmp_errors
        }else{
          output_errors <- rbind(output_errors, tmp_errors)
        }
      }
      return(output_errors)
    }

    # doing consecutive missing instead
    for(start_miss_index in not_na_landsat_indices){
      # print(start_miss_index) # debug
      tmp_errors <- rep(NA, length(landsat_dates))
      missing_indices <- sort(not_na_landsat_indices[which(not_na_landsat_indices == start_miss_index):(which(not_na_landsat_indices == start_miss_index) + consecutive_missing - 1)])
      missing_dates <- landsat_dates[missing_indices - 2]
      # check that we're not over the maximum landsat index, if so break the loop b/c we're done
      if(any(is.na(missing_indices))) break
      # if(any(is.na(missing_indices))) next
      # if(max(missing_indices) > max_landsat_index) next # this probably never happens...
      tmp_x <- x
      # tmp_x[2:length(tmp_x)] <- tmp_x[2:length(tmp_x)] / scale_factor
      og_evi <- tmp_x[missing_indices] / scale_factor
      tmp_x[missing_indices] <- NA
      kf_result <- FuseLandsatModisEVI(tmp_x, landsat_dates, modis_dates, landsat_sensor, cdl_tv_sd, cdl_types, scale_factor=1e4, smooth=T, plot=F)
      # kf_dates <- sort(as.Date(paste(rep(unique(strftime(c(landsat_dates, modis_dates), format="%Y")), 365), 1:365, sep="-"), format="%Y-%j"))
      kf_evi <- kf_result$evi[sapply(kf_dates, function(x) any(x==missing_dates))]
      tmp_error <- og_evi - kf_evi
      tmp_errors[missing_indices - 2] <- tmp_error
      # append the cdl flag
      tmp_errors <- c(x[1], tmp_errors)
      if(is.null(output_errors)){
        output_errors <- tmp_errors
      }else{
        output_errors <- rbind(output_errors, tmp_errors)
      }
    }
    return(output_errors)
}

#-------------------------------------------------------------------------------
# Loops over 10% to 90% random missing data in increments of 10% and quantifies
# error at omitted Landsat observations. Each missing fraction is resampled miss_iter times
ProgressiveMissingFraction <- function(x, landsat_dates, modis_dates, landsat_sensor, cdl_process_sds, cdl_types, miss_iter=10){
  total_errors <- NULL
  x_lo <- x
  # x_lo[(length(landsat_dates) + 3):length(x_lo)] <- NA # eliminate all MODIS measurements
  x_lo[(length(landsat_dates) + 3):(length(landsat_dates) + 2 + length(modis_dates))] <- NA # eliminate all MODIS measurements

  for(miss_frac in seq(0.1, 0.9, by=0.1)){
    fusion_errors <- GetErrorLandsatMODISFusion(x, landsat_dates=landsat_dates, modis_dates=modis_dates, landsat_sensor=landsat_sensor, cdl_tv_sd=cdl_process_sds, cdl_types=cdl_types, consecutive_missing=1, missing_fraction=miss_frac, missing_iterations=miss_iter)
    fusion_errors <- data.frame(fusion_errors)
    names(fusion_errors) <- c("cdl_code", "rmse", strftime(landsat_dates, format="%j"))
    long_fusion_errors <- melt(fusion_errors, measure.vars=3:223, variable.name="doy")
    long_fusion_errors$type <- "fusion"

    lo_errors <- GetErrorLandsatMODISFusion(x_lo, landsat_dates=landsat_dates, modis_dates=modis_dates, landsat_sensor=landsat_sensor, cdl_tv_sd=cdl_process_sds, cdl_types=cdl_types, consecutive_missing=1, missing_fraction=miss_frac, missing_iterations=miss_iter)
    lo_errors <- data.frame(lo_errors)
    names(lo_errors) <- c("cdl_code", "rmse", strftime(landsat_dates, format="%j"))
    long_lo_errors <- melt(lo_errors, measure.vars=3:223, variable.name="doy")
    long_lo_errors$type <- "landsat"

    tmp_total_errors <- rbind(long_lo_errors, long_fusion_errors)
    tmp_total_errors$miss_frac <- miss_frac

    if(is.null(total_errors)){
      total_errors <- tmp_total_errors
    }else{
      total_errors <- rbind(total_errors, tmp_total_errors)
    }
  }
  total_errors <- total_errors[!is.na(total_errors$value), ]
  total_errors$doy <- as.numeric(as.character(total_errors$doy)) # convert DOY to numeric
  return(total_errors)
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
PlotForecast <- function(filt_m, filt_se, signal, t=NULL, conf_level=0.95, sigma=NULL, ylim=NULL, pt_cex=1, colmain=rgb(0.42, 0.68, 0.84), colerr=rgb(0.42, 0.68, 0.84, 0.5), colsignal="#636363", ...){
  if(is.null(t)) t <- 1:length(filt_m)
  # colmain <- "#3182BD"; colerr <- "#BDD7E7"; colsignal <- "#636363"
  #6BAED6 = rgb(0.42, 0.68, 0.84)
  # rgb(0.42, 0.68, 0.84, 0.5)


  if(is.null(sigma)) sigma <- qnorm(1 - (1 - conf_level) / 2)
  lower <- filt_m + (sigma * filt_se)
  upper <- filt_m - (sigma * filt_se)

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
      points(t, tmp_y, type="p", pch=i, cex=pt_cex, col=colsignal)
    }
  }else{
    points(t, signal, type="p", pch=1, cex=pt_cex, col=colsignal)
  }
}


#--------------------------------------------------------------------------------
GetSDSName <- function(mcd43a4_file_path, band){
  # returns an appropriate SDS name for an MCD43A4 file path of the specified band
  return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":", paste("MOD_Grid_BRDF:Nadir_Reflectance_Band", band, sep=""), sep = ""))
}

#--------------------------------------------------------------------------------
GetSDSName_SnowBRDF <- function(mcd43a2_file_path){
  # returns the SDS name for the Snow BRDF Albedo layer of an MCD43A2 file
  return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a2_file_path, "\":MOD_Grid_BRDF:Snow_BRDF_Albedo", sep = ""))
}

#--------------------------------------------------------------------------------
GetSDSName_AlbedoBandQA <- function(mcd43a2_file_path, band){
  # returns the SDS name for the BRDF Albedo Band Quality for a particular band of an MCD43A2 file
  return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a2_file_path, "\":", paste("MOD_Grid_BRDF:BRDF_Albedo_Band_Quality_Band", band, sep=""), sep = ""))
}

#--------------------------------------------------------------------------------
CellNumRaster <- function(modis_file_path){
  # creates and returns a RasterLayer containing the cell numbers associated
  # with the raster given by modis_file_path
  r_modis <- raster(GetSDSName(modis_file_path, 1))
  values(r_modis) <- 1:ncell(r_modis)
  return(r_modis)
}

#--------------------------------------------------------------------------------
TileIndexRaster <- function(modis_file_path){
  # creates and returns a RasterLayer containing the MODIS tile index for
  # the raster given by modis_file_path. Value is the "h" number plus the "v"
  # number plus 10000. E.g. for "h28v09" the value is "12809"
  tile_h <- gsub(".*h([0-9]{2}).*", "\\1", basename(modis_file_path))
  tile_v <- gsub(".*v([0-9]{2}).*", "\\1", basename(modis_file_path))
  tile_num <- as.integer(paste("1", tile_h, tile_v, sep=""))
  r_modis <- raster(GetSDSName(modis_file_path, 1))
  values(r_modis) <- tile_num
  return(r_modis)
}

#--------------------------------------------------------------------------------
MakeCellNumberMap <- function(example_modis_files, example_landsat_file, out_dir){
  # creates a cell number map with the geometry of example_landsat_file over a vector of modis files
  # thus, it handles the case where a landsat geometry overlaps multiple MODIS tiles
  cell_num_rasters <- lapply(as.character(example_modis_files), CellNumRaster)
  if(length(cell_num_rasters) > 1){
    merged_cell_num_rasters <- do.call(merge, cell_num_rasters)
  }else{
    merged_cell_num_rasters <- cell_num_rasters[[1]]
  }
  r_landsat <- raster(example_landsat_file)
  tmp_out_file <- file.path(out_dir, paste("modis_temp_", strftime(Sys.time(), format="%a_%b_%d_%H_%M_%S"), paste(sample(letters, 4), collapse=""), ".tif", sep=""))
  writeRaster(merged_cell_num_rasters, file=tmp_out_file)
  landsat_path_row <- gsub(".*_([0-9]{6})_.*", "\\1", basename(example_landsat_file))
  out_file <- file.path(out_dir, paste("modis_cell_num_", landsat_path_row, ".tif", sep=""))
  gdal_cmd <- paste("/Library/Frameworks/GDAL.framework/Programs/gdalwarp -overwrite -dstnodata -9999 -s_srs '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m' -t_srs", paste("'", projection(r_landsat), "'", sep=""), "-te", paste(extent(r_landsat)[c(1, 3, 2, 4)], collapse=" "), "-ts", ncol(r_landsat), nrow(r_landsat), tmp_out_file, out_file)
  system(gdal_cmd)
  rm_cmd <- paste("rm -f", tmp_out_file)
  system(rm_cmd)
  return(out_file)
}

#--------------------------------------------------------------------------------
MakeTileIndexMap <- function(example_modis_files, example_landsat_file, out_dir){
  # creates a tile index map with the geometry of example_landsat_file over a vector of modis files
  # thus, it handles the case where a landsat geometry overlaps multiple MODIS tiles
  tile_index_rasters <- lapply(as.character(example_modis_files), TileIndexRaster)
  if(length(tile_index_rasters) > 1){
    merged_tile_index_rasters <- do.call(merge, tile_index_rasters)
  }else{
    merged_tile_index_rasters <- tile_index_rasters[[1]]
  }
  r_landsat <- raster(example_landsat_file)
  tmp_out_file <- file.path(out_dir, paste("modis_temp_", strftime(Sys.time(), format="%a_%b_%d_%H_%M_%S"), paste(sample(letters, 4), collapse=""), ".tif", sep=""))
  writeRaster(merged_tile_index_rasters, file=tmp_out_file)
  landsat_path_row <- gsub(".*_([0-9]{6})_.*", "\\1", basename(example_landsat_file))
  out_file <- file.path(out_dir, paste("modis_tile_index_", landsat_path_row, ".tif", sep=""))
  gdal_cmd <- paste("/Library/Frameworks/GDAL.framework/Programs/gdalwarp -overwrite -dstnodata -9999 -s_srs '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m' -t_srs", paste("'", projection(r_landsat), "'", sep=""), "-te", paste(extent(r_landsat)[c(1, 3, 2, 4)], collapse=" "), "-ts", ncol(r_landsat), nrow(r_landsat), tmp_out_file, out_file)
  system(gdal_cmd)
  rm_cmd <- paste("rm -f", tmp_out_file)
  system(rm_cmd)
  return(out_file)
}

#--------------------------------------------------------------------------------
CreateModisMaps <- function(path_row, landsat_path_row_shp, modis_tile_shp, landsat_data_dir="/Volumes/research/fer/jmgray2/EastKalimantan/processedData/envi", modis_data_dir="/Users/jmgray2/Desktop/MCD43A4", modis_map_out_dir="/Users/jmgray2/Desktop/MODIS_cell_maps", landsat_suffix="_EK"){
  # Creates both MODIS tile index and cell number maps for a given Landsat path row
  landsat_in_files <- dir(file.path(landsat_data_dir, path_row), pattern=paste(".*", landsat_suffix, "$", sep=""), full=T)
  this_path_row_shp <- landsat_path_row_shp[which(landsat_path_row_shp$PR == path_row),]
  this_path_row_shp_reproj <- spTransform(this_path_row_shp, CRS(projection(modis_tile_shp)))
  modis_grid_intersect <- over(this_path_row_shp_reproj, modis_tile_shp, returnList=T)
  intersecting_modis_tiles <- apply(modis_grid_intersect[[1]][,c("h", "v")], 1, function(x) paste("h", formatC(as.integer(x[1]), width=2, flag="0"), "v", formatC(as.integer(x[2]), width=2, flag="0"), sep=""))
  modis_in_files <- lapply(intersecting_modis_tiles, function(tile) list.files(modis_data_dir, pattern=tile, rec=T, full=T))
  example_modis_files <- lapply(modis_in_files, function(x) x[1])
  modis_cell_num_map_file <- MakeCellNumberMap(example_modis_files, landsat_in_files[1], modis_map_out_dir)
  modis_tile_index_map_file <- MakeTileIndexMap(example_modis_files, landsat_in_files[1], modis_map_out_dir)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Make MODIS cell number and tile index maps at Landsat resolution for all path-rows:
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Paths and constants
# modis_cell_maps_output_dir <- "/Users/jmgray2/Desktop/MODIS_cell_maps"
# landsat_data_dir <- "/Volumes/research/fer/jmgray2/EastKalimantan/processedData/envi"
# landsat_suffix <- "_EK"
# mcd43a4_data_dir <- "/Users/jmgray2/Desktop/MCD43A4"
# path_row_shp_file <- "/Users/jmgray2/Desktop/wrs2_descending/wrs2_descending.shp"
# modis_grid_shp_file <- "/Users/jmgray2/Desktop/modis_grid/modis_sinusoidal_grid_world.shp"
#
# # read in the shapefiles for Landsat and MODIS footprints
# path_rows_shp <- shapefile(path_row_shp_file)
# modis_grid_shp <- shapefile(modis_grid_shp_file)
#
# # get all Landsat scenes
# all_path_rows <- dir(landsat_data_dir)
# for(path_row in all_path_rows){
#   print(paste("Making maps for:", path_row))
#   CreateModisMaps(path_row=path_row, landsat_path_row_shp=path_rows_shp, modis_tile_shp=modis_grid_shp, landsat_data_dir=landsat_data_dir, modis_data_dir=mcd43a4_data_dir, modis_map_out_dir=modis_cell_maps_output_dir, landsat_suffix=landsat_suffix)
# }

#--------------------------------------------------------------------------------
GetMODISLines <- function(x, mcd43a4_in_files, mcd43a2_in_files, modis_line_ranges, landsat_to_modis_bands, cast_as_int=T, scale_factor=1e4, max_open_datasets=2e3){
  # for use in an lapply expression only, argument "x" is the index into
  # mcd43a4_in_files, mcd43a2_in_files, and modis_line_ranges which MUST exist
  # within the current scope. It is a list index used to get the proper MCD43A4/A2
  # file names and MODIS line ranges

  # get the MODIS data for the specified lines
  tmp_modis_r <- raster(GetSDSName(mcd43a4_in_files[[x]][1], 1)) # read in a temporary MODIS file for geometry
  modis_tmp_data <- array(NA, dim=c(ncol(tmp_modis_r) * (diff(modis_line_ranges[[x]]) + 1), length(mcd43a4_in_files[[x]]), length(landsat_to_modis_bands) + 1))
  modis_tmp_data_dim <- dim(modis_tmp_data)
  i <- 1
  for(modis_band in landsat_to_modis_bands){
    if(!is.na(modis_band)){
      mcd43a4_tmp_names <- GetSDSName(mcd43a4_in_files[[x]], modis_band)
      mcd43a4_tmp_data <- GetValuesGDAL(mcd43a4_tmp_names, start_row=modis_line_ranges[[x]][1], n=diff(modis_line_ranges[[x]]) + 1, max_open_datasets=max_open_datasets)
      modis_tmp_data[, , i] <- mcd43a4_tmp_data
      rm(mcd43a4_tmp_data)
    }else{
      # no MODIS data for this band (thermal) so we input a fill matrix instead
      modis_tmp_data[, , i] <- matrix(NA, nrow=dim(modis_tmp_data)[1], ncol=dim(modis_tmp_data)[2])
    }
    i <- i + 1
  }
  # get A2 data; we assume that all band Albedo Band QA are the same and use only band 1; we fill with NA where Snow BRDF Albedo is 1
  mcd43a2_snow_tmp_names <- GetSDSName_SnowBRDF(mcd43a2_in_files[[x]])
  mcd43a4_snow_tmp_data <- GetValuesGDAL(mcd43a2_snow_tmp_names, start_row=modis_line_ranges[[x]][1], n=diff(modis_line_ranges[[x]]) + 1, max_open_datasets=max_open_datasets)
  mcd43a2_albedo_qa_tmp_names <- GetSDSName_AlbedoBandQA(mcd43a2_in_files[[x]], 1)
  mcd43a4_albedo_qa_tmp_data <- GetValuesGDAL(mcd43a2_albedo_qa_tmp_names, start_row=modis_line_ranges[[x]][1], n=diff(modis_line_ranges[[x]]) + 1, max_open_datasets=max_open_datasets)

  # make Albedo Band QA "4" where it was a snow retrieval
  # NOTE: should this be made "NA" instead?
  mcd43a4_albedo_qa_tmp_data[mcd43a4_snow_tmp_data == 1] <- 4
  modis_tmp_data[, , i] <- mcd43a4_albedo_qa_tmp_data

  # scale and cast as integer if requested
  if(cast_as_int){
    # only multiply bands 1:7, leave the QA band as-is
    modis_tmp_data[,, 1:7] <- modis_tmp_data[,,1:7] * scale_factor
    modis_tmp_data <- as.integer(modis_tmp_data) # cast as integer
    modis_tmp_data <- array(modis_tmp_data, dim=modis_tmp_data_dim)
  }

  return(modis_tmp_data)
}

#--------------------------------------------------------------------------------
GetModisCellOffset <- function(x){
  # gets the appropriate offset number to adjust original MODIS cell numbers to the extracted lines
  # original numbers are adjusted by subtracting the minimum cell number of the minimum extracted line number
  # verified to work:
  # get the MODIS cell coordinate for a particular landsat_cell_number
  # tmp_modis_cell_num <- cellFromXY(tmp_modis_r, spTransform(SpatialPoints(xyFromCell(tmp_r, landsat_cell_num), CRS(projection(tmp_r))), CRS(projection(tmp_modis_r))))
  tmp_modis_r <- raster(GetSDSName(mcd43a4_in_files[[x]][1], 1)) # read in a temporary MODIS file for geometry
  tmp_offset <- (modis_line_ranges[[x]][1] - 1) * ncol(tmp_modis_r)
  return(tmp_offset)
}

#--------------------------------------------------------------------------------
AssembleSinglePixelTimeSeries <- function(landsat_cell_num, landsat_data, modis_data, modis_cell_nums, modis_tile_indices, modis_cell_offsets, landsat_dates, modis_dates, landsat_sensor, lwmask_data, qa_band_index=8, default_rmse=1e6){
  # returns a list object containing the multispectral landsat and modis time series,
  # their associated dates, and the landsat sensor identifier for the landsat pixel number provided
  # as the first argument. The appropriate modis_cell_nums, modis_tile_indices, and modis_cell_offsets
  # must be created and provided, in addition to the landsat_data and modis_data for the give lines chunk

  # find the appropriate MODIS data set
  the_modis_dataset <- which(unique(modis_tile_indices) == modis_tile_indices[landsat_cell_num])

  # check if this is a land pixel, if not don't bother doing anything else
  if(lwmask_data[landsat_cell_num,] != 0){
    return(list(land_pixel=F, nbands=dim(landsat_data)[3] - 1, landsat_data=NA, modis_data=NA, landsat_dates=landsat_dates, modis_dates=modis_dates[[the_modis_dataset]], landsat_sensor=NA, landsat_modis_rmse=NA))
  }

  tmp_landsat_data <- landsat_data[landsat_cell_num,,]
  modis_cell_num <- modis_cell_nums[landsat_cell_num]
  modis_tile_index <- modis_tile_indices[landsat_cell_num]

  # get the appropriate MODIS cell number
  modis_cell_num <- modis_cell_num - modis_cell_offsets[[the_modis_dataset]]
  tmp_modis_data <- modis_data[[the_modis_dataset]][modis_cell_num,,]

  # calculate the band-by-band RMSE for matched dates
  good_landsat_obs <- tmp_landsat_data[, qa_band_index] == 0 & !is.na(tmp_landsat_data[, qa_band_index])
  # if there are 1 or less good_landsat_obs, we can't find RMSE, so use default
  if(length(sum(good_landsat_obs)) <= 1){
    matched_rmse <- rep(default_rmse, dim(tmp_landsat_data)[2] - 1)
  }else{
    matched_modis_dates <- match(landsat_dates[good_landsat_obs], modis_dates[[the_modis_dataset]])
    matched_rmse <- apply((tmp_landsat_data[good_landsat_obs,] - tmp_modis_data[matched_modis_dates,])^2, 2, function(x) sqrt(mean(x, na.rm=T)))[-qa_band_index]
    matched_rmse[is.nan(matched_rmse) | is.na(matched_rmse)] <- default_rmse
  }

  # construct and return the single-pixel data object
  ret_object <- list(land_pixel=T, nbands=dim(landsat_data)[3] - 1, landsat_data=tmp_landsat_data, modis_data=tmp_modis_data, landsat_dates=landsat_dates, modis_dates=modis_dates[[the_modis_dataset]], landsat_sensor=landsat_sensor, landsat_modis_rmse=matched_rmse)
  return(ret_object)
}

#--------------------------------------------------------------------------------
MakeLandsatMODISKFData <- function(data_object, temp_res="daily", agg_func="median", qa_band_index=8, oli_mult_error=0.05, etm_mult_error=0.05, tm_mult_error=0.07){
  # this expands/aggregates landsat and modis data to a regularly sampled matrix of time series values (daily or weekly supported)
  # it does QA/QC screening (just removes snow values from MODIS, and anyting where FMASK != 0 from landsat)
  # returns the data matrix suitable for KF filtering, a vector of the matched date MODIS-Landsat RMSE,
  # and the Landsat time-varying uncertainty matrix

  min_date <- as.Date(paste(strftime(min(c(data_object$landsat_dates, data_object$modis_dates)), format="%Y"), "-1-1", sep=""))
  max_date <- as.Date(paste(strftime(max(c(data_object$landsat_dates, data_object$modis_dates)), format="%Y"), "-12-31", sep=""))
  if(temp_res == "daily"){
    # do daily
    pred_dates <- seq.Date(min_date, max_date, by="day")
  }else{
    # do weekly
    pred_dates <- seq.Date(min_date, max_date, by="week")
  }

  # check if this is a land pixel, if not don't bother doing anything else
  if(!data_object$land_pixel){
    return(list(land_pixel=F, nbands=data_object$nbands, Y=NA, P_landsat_tv=NA, dates=pred_dates, landsat_modis_rmse=NA))
  }

  # determine if we are aggregating using "mean" or "median"
  if(agg_func == "mean"){
    by_func <- function(x, by_inds) by(x, by_inds, FUN=mean, na.rm=T)
  }else{
    by_func <- function(x, by_inds) by(x, by_inds, FUN=median, na.rm=T)
  }

  # get the mean of all modis data by the pred_dates
  # first, filter out missing observations using the MODIS QA (NOTE: only filters out snow!)
  tmp_modis_qa <- rbind(replicate(dim(data_object$modis_data[, -qa_band_index])[2], data_object$modis_data[, qa_band_index]))
  tmp_modis_data <- data_object$modis_data[, -qa_band_index]
  tmp_modis_data[is.na(tmp_modis_qa)] <- NA
  modis_by_indices <- findInterval(data_object$modis_dates, pred_dates, all.inside=F)
  unique_modis_by_indices <- sort(unique(modis_by_indices))
  tmp_modis_data <- apply(tmp_modis_data, 2, by_func, by_inds=modis_by_indices)

  # get the mean of all Landsat data by the pred_dates
  tmp_landsat_qa <- rbind(replicate(dim(data_object$landsat_data[, -qa_band_index])[2], data_object$landsat_data[, qa_band_index]))
  tmp_landsat_data <- data_object$landsat_data[, -qa_band_index]
  tmp_landsat_data[tmp_landsat_qa != 0 & !is.na(tmp_landsat_qa)] <- NA
  landsat_by_indices <- findInterval(data_object$landsat_dates, pred_dates, all.inside=F)
  unique_landsat_by_indices <- sort(unique(landsat_by_indices))
  tmp_landsat_data <- apply(tmp_landsat_data, 2, by_func, by_inds=landsat_by_indices)

  # get the time-varying Landsat errors
  landsat_error_multiplier <- rep(NA, nrow(tmp_landsat_data))
  landsat_error_multiplier[data_object$landsat_sensor[order(unique(landsat_by_indices))] == "LC08"] <- oli_mult_error
  landsat_error_multiplier[data_object$landsat_sensor[order(unique(landsat_by_indices))] == "LE07"] <- etm_mult_error
  landsat_error_multiplier[data_object$landsat_sensor[order(unique(landsat_by_indices))] == "LT05"] <- tm_mult_error
  landsat_error <- t(t(tmp_landsat_data) * landsat_error_multiplier)

  # prototype the output matrix for KF
  Y <- matrix(NA, nrow=dim(tmp_landsat_data)[2] * 2, ncol=length(pred_dates))
  Y_landsat_p <- matrix(NA, nrow=dim(tmp_landsat_data)[2], ncol=length(pred_dates))
  nbands <- dim(tmp_landsat_data)[2]
  # NOTE: this is probably bad juju b/c functions shouldn't have side effects. However, this avoids a for-loop and so is a
  # bit more readable. We also can ensure that this never runs in parallel, so it _should_ be ok
  trash <- lapply(unique_landsat_by_indices, function(x) Y[1:nbands, x] <<- tmp_landsat_data[which(unique_landsat_by_indices == x), ])
  trash <- lapply(unique_modis_by_indices, function(x) Y[(nbands + 1):(nbands * 2), x] <<- tmp_modis_data[which(unique_modis_by_indices == x), ])
  trash <- lapply(unique_landsat_by_indices, function(x) Y_landsat_p[1:nbands, x] <<- landsat_error[which(unique_landsat_by_indices == x), ])

  return(list(land_pixel=T, nbands=data_object$nbands, Y=Y, P_landsat_tv=Y_landsat_p, dates=pred_dates, landsat_modis_rmse=data_object$landsat_modis_rmse))
}

#--------------------------------------------------------------------------------
DoKF <- function(landsat_cell_num, landsat_data, modis_data, modis_cell_nums, modis_tile_indices, modis_cell_offsets, landsat_dates, modis_dates, landsat_sensor, lwmask_data, qa_band_index=8){
  data_object <- AssembleSinglePixelTimeSeries(landsat_cell_num, landsat_data=landsat_data, modis_data=modis_data, modis_cell_nums=modis_cell_nums, modis_tile_indices=modis_tile_indices, modis_cell_offsets=modis_cell_offsets, landsat_dates=landsat_dates, modis_dates=modis_dates, landsat_sensor=landsat_sensor, lwmask_data=lwmask_data)
  KF_data_object <- MakeLandsatMODISKFData(data_object, temp_res="weekly")

  if(!KF_data_object$land_pixel){
    # if this is not a land pixel, then return a matrix of NAs
    return(matrix(NA, nrow=length(KF_data_object$dates), ncol=KF_data_object$nbands))
  }

  # assemble the dlm
  mydlm <- MakeMultiDLM(KF_data_object$nbands, 2)
  # mydlm$m0 <- KF_data_object$Y[1:nbands, min(which(!is.na(KF_data_object$Y[1, ])))] # set initial condition to first non-missing Landsat observation (fragile: perhaps not all Landsat band obs are non-missing!)
  mydlm$m0 <- apply(KF_data_object$Y[1:KF_data_object$nbands,], 1, function(x) x[min(which(!is.na(x)))])

  # if there are NA m0 values we fill 0.2
  mydlm$m0[is.na(mydlm$m0)] <- 2000

  diag(mydlm$V) <- c(rep(1e6, KF_data_object$nbands), KF_data_object$landsat_modis_rmse) # modis-landsat long-term RMSE is in columns 2:5 of Y
  mydlm$JV <- diag(c(1:KF_data_object$nbands, rep(0, KF_data_object$nbands))) # allow the Landsat observation error to be time-varying, use long-term MODIS-landsat RMSE for MODIS obs error
  X <- t(KF_data_object$P_landsat_tv)
  mydlm$X <- X

  # and the process error matrix
  process_error <- 10  # assume arbitrary constant process error
  diag(mydlm$W) <- rep(process_error, KF_data_object$nbands)

  # run the filter
  m_smooth <- dlmSmooth(t(KF_data_object$Y), mydlm) # smoothing
  return(dropFirst(m_smooth$s))
  # plotting
  # layout(matrix(1:8, nrow=2, byrow=T))
  # par(mar=rep(1, 4))
  # for(i in 1:KF_data_object$nbands) plotdlmresults(Y=KF_data_object$Y, dlm_result=m_smooth, dates=KF_data_object$dates, plot_band=i)

}

#--------------------------------------------------------------------------------
WriteKFResults <- function(data_to_write, out_files, example_r, cast_as_int=T){
  # writes multiband KF output lists to the specified output files as BIP ENVI raster files

  # convert data_to_write into a matrix
  data_mat <- array(unlist(data_to_write), dim=c(dim(data_to_write[[1]])[1], dim(data_to_write[[1]])[2], length(data_to_write)))

  # loop through each band and write the data
  for(i in 1:dim(data_mat)[2]){
    out_file <- out_files[i]
    # create or open the file for appending
    if(!file.exists(out_file)){
      ff <- file(out_file, 'wb') # create the file and open for writing
    }else{
      ff <- file(out_file, 'ab') # file exists, append to the end
    }
    # write the data
    if(cast_as_int){
      envi_data_type <- 3
      writeBin(as.integer(round(data_mat[,i,])), ff)
    }else{
      envi_data_type <- 5
      writeBin(c(data_mat[,i,]), ff)
    }
    # writeBin(c(data_mat[,i,]), ff)
    close(ff)
    # create a header if it doesn't yet exist
    out_hdr <- paste(out_file, ".hdr", sep="")
    if(!file.exists(out_hdr)){
      temp_txt = paste(
        "ENVI description = { KF Output }",
        "\nsamples = ", ncol(example_r),
        "\nlines =", nrow(example_r),
        "\nbands = ", dim(data_mat)[1],
        # "\npixel size = {", paste(res(example_r), collapse=", "), "}",
        "\nheader offset = 0",
        "\nfile type = ENVI Standard",
        "\ndata type =", envi_data_type,
        "\ninterleave = bip",
        "\nbyte order = 0",
        paste("\nmap info = {UTM, 1, 1", xmin(example_r), ymax(example_r), res(example_r)[1], res(example_r)[2], gsub(".*\\+zone=([0-9]*).*", "\\1", projection(example_r)), "North, WGS-84}", sep=", "),
        "\ncoordinate system string = {}",
        sep="")
      sink(out_hdr)
      cat(temp_txt)
      sink()
    }
  }
}

#--------------------------------------------------------------------------------
WriteENVIMultiband <- function(data_to_write, out_file, example_r, cast_as_int=T){
    if(!is.matrix(data_to_write)){
        print("Provided data is not a matrix, aborting")
        return(NA)
    }

    # create or open the file for appending
    if(!file.exists(out_file)){
        ff <- file(out_file, 'wb') # create the file and open for writing
    }else{
        ff <- file(out_file, 'ab') # file exists, append to the end
    }

    # write the data
    if(cast_as_int){
        envi_data_type <- 3
        # writeBin(as.integer(round(data_mat[,i,])), ff)
        writeBin(as.integer(round(c(t(data_to_write)))), ff)
    }else{
        envi_data_type <- 5
        # writeBin(c(data_mat[,i,]), ff)
        writeBin(c(t(data_to_write)), ff)
    }
    close(ff)
    # create a header if it doesn't yet exist
    out_hdr <- paste(out_file, ".hdr", sep="")
    if(!file.exists(out_hdr)){
        temp_txt = paste(
        "ENVI description = { KF Output }",
        "\nsamples = ", ncol(example_r),
        "\nlines =", nrow(example_r),
        "\nbands = ", dim(data_to_write)[2],
        # "\npixel size = {", paste(res(example_r), collapse=", "), "}",
        "\nheader offset = 0",
        "\nfile type = ENVI Standard",
        "\ndata type =", envi_data_type,
        "\ninterleave = bip",
        "\nbyte order = 0",
        paste("\nmap info = {UTM, 1, 1", xmin(example_r), ymax(example_r), res(example_r)[1], res(example_r)[2], gsub(".*\\+zone=([0-9]*).*", "\\1", projection(example_r)), "North, WGS-84}", sep=", "),
        "\ncoordinate system string = {}",
        sep="")
        sink(out_hdr)
        cat(temp_txt)
        sink()
    }
}

#--------------------------------------------------------------------------------
plotdlmresults <- function(Y, dlm_result, dates, plot_band){
  ylim <- range(c(Y[c(plot_band, plot_band + nbands), ], dropFirst(dlm_result$s)[, plot_band]), na.rm=T)
  plot(dates, Y[plot_band, ], col=1, ylim=ylim, xlab="", ylab="Surface Reflectance")
  points(dates, Y[plot_band + nbands, ], col=2)
  points(KF_data_object$dates, dropFirst(dlm_result$s)[, plot_band], type="l", col="grey")
  title(paste("Band:", plot_band))
}
