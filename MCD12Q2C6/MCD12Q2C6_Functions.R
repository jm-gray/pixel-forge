#==============================================================================
#==============================================================================
# Functions for making a stack of NBAR data, calculating VIs, spline smoothing,
# and screening/filling snow contaminated and outlier values
#==============================================================================
#==============================================================================

#---------------------------------------------------------------------
# these functions return an SDS filename string so the raster package can handle the HDF inputs
# Get_SDS_name_red <- function(mcd43a4_file_path){
# 	paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":MOD_Grid_BRDF:Nadir_Reflectance_Band1", sep = "")
# }
#
# Get_SDS_name_nir <- function(mcd43a4_file_path){
# 	paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":MOD_Grid_BRDF:Nadir_Reflectance_Band2", sep = "")
# }
#
# Get_SDS_name_blue <- function(mcd43a4_file_path){
# 	paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":MOD_Grid_BRDF:Nadir_Reflectance_Band3", sep = "")
# }
#
# Get_SDS_name_band4 <- function(mcd43a4_file_path){
# 	paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":MOD_Grid_BRDF:Nadir_Reflectance_Band4", sep = "")
# }
#
# Get_SDS_name_band6 <- function(mcd43a4_file_path){
# 	paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":MOD_Grid_BRDF:Nadir_Reflectance_Band6", sep = "")
# }

Get_SDS_name_band1 <- function(mcd43a4_file_path, fill_tile){
	ifelse(
		is.na(mcd43a4_file_path),
		return(fill_tile),
		return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":MOD_Grid_BRDF:Nadir_Reflectance_Band1", sep = ""))
	)
}

Get_SDS_name_band2 <- function(mcd43a4_file_path, fill_tile){
	ifelse(
		is.na(mcd43a4_file_path),
		return(fill_tile),
		return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":MOD_Grid_BRDF:Nadir_Reflectance_Band2", sep = ""))
	)
}

Get_SDS_name_band3 <- function(mcd43a4_file_path, fill_tile){
	ifelse(
		is.na(mcd43a4_file_path),
		return(fill_tile),
		return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":MOD_Grid_BRDF:Nadir_Reflectance_Band3", sep = ""))
	)
}

Get_SDS_name_band4 <- function(mcd43a4_file_path, fill_tile){
	ifelse(
		is.na(mcd43a4_file_path),
		return(fill_tile),
		return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":MOD_Grid_BRDF:Nadir_Reflectance_Band4", sep = ""))
	)
}

Get_SDS_name_band6 <- function(mcd43a4_file_path, fill_tile){
	ifelse(
		is.na(mcd43a4_file_path),
		return(fill_tile),
		return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":MOD_Grid_BRDF:Nadir_Reflectance_Band6", sep = ""))
	)
}

Get_SDS_name_NBAR_QA <- function(mcd43a2_file_path, fill_tile){
	ifelse(
		is.na(mcd43a2_file_path),
		return(fill_tile),
		return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a2_file_path, "\":MOD_Grid_BRDF:BRDF_Albedo_Quality", sep = ""))
	)
}

Get_SDS_name_NBAR_snow <- function(mcd43a2_file_path, fill_tile){
	ifelse(
		is.na(mcd43a2_file_path),
		return(fill_tile),
		return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a2_file_path, "\":MOD_Grid_BRDF:Snow_BRDF_Albedo", sep = ""))
	)
}

Get_SDS_name_NBAR_bandquality <- function(mcd43a2_file_path, fill_tile){
	ifelse(
		is.na(mcd43a2_file_path),
		return(fill_tile),
		return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a2_file_path, "\":MOD_Grid_BRDF:BRDF_Albedo_Band_Quality", sep = ""))
	)
}

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


#---------------------------------------------------------------------
# returns NBAR bands 1, 2, 3, 4, and 6, as well as QA, snow, and bandquality datasets for a particular tile, within a date range
GetTileDataSets <- function(tile, start_date=as.Date("2013-1-1"), end_date=as.Date("2014-12-31"), nbar_data_dir="/projectnb/modislc/data/mcd12_in/c5/mcd43a4", nbarqa_data_dir="/projectnb/modislc/data/mcd12_in/c5/mcd43a2", tmp_out_dir="/projectnb/modislc/users/joshgray/MCD12Q2C6/filltiles"){
  # expected_dates[!expected_dates %in% actual_dates]
	# make filepaths for expected data
	years <- as.numeric(strftime(start_date,format="%Y")):as.numeric(strftime(end_date,format="%Y"))

	nbar_file_path_df <- data.frame(nbar_data_dir, expand.grid(years, formatC(seq(1, 361, by=8), width=3, flag="0"), paste("*", as.character(tile), "*.hdf", sep="")))
	dates <- as.Date(apply(nbar_file_path_df, 1, function(x) paste(x[2:3], collapse="-")), format="%Y-%j")
	nbar_file_path_df[] <- lapply(nbar_file_path_df, as.character) # Thanks, Hadley
	nbar_files <- apply(nbar_file_path_df, 1, function(x) {Sys.glob(paste(x, collapse="/"))})
	nbar_files[unlist(lapply(nbar_files, length)) == 0] <- NA
	nbar_files <- unlist(nbar_files)[order(dates)]

	nbarqa_file_path_df <- data.frame(nbarqa_data_dir, expand.grid(years, formatC(seq(1, 361, by=8), width=3, flag="0"), paste("*", as.character(tile), "*.hdf", sep="")))
	nbarqa_file_path_df[] <- lapply(nbarqa_file_path_df, as.character) # Thanks, Hadley
	nbarqa_files <- apply(nbarqa_file_path_df, 1, function(x) {Sys.glob(paste(x, collapse="/"))})
	nbarqa_files[unlist(lapply(nbarqa_files, length)) == 0] <- NA
	nbarqa_files <- unlist(nbarqa_files)[order(dates)]

	dates <- sort(dates)

	# Create a fill file if any expected files are missing
	if(any(is.na(c(nbar_files, nbarqa_files)))){
		fill_tile <- CreateFillTile(tile, tmp_out_dir)
	}else{
		fill_tile <- NA
	}

	# create data sets
	data_sets <- c(
		unlist(lapply(nbar_files, Get_SDS_name_band1, fill_tile)),
		unlist(lapply(nbar_files, Get_SDS_name_band2, fill_tile)),
		unlist(lapply(nbar_files, Get_SDS_name_band3, fill_tile)),
		unlist(lapply(nbar_files, Get_SDS_name_band4, fill_tile)),
		unlist(lapply(nbar_files, Get_SDS_name_band6, fill_tile)),
		unlist(lapply(nbarqa_files, Get_SDS_name_NBAR_QA, fill_tile)),
		unlist(lapply(nbarqa_files, Get_SDS_name_NBAR_snow, fill_tile)),
		unlist(lapply(nbarqa_files, Get_SDS_name_NBAR_bandquality, fill_tile))
	)

	return(list(data_sets, dates))
}

#---------------------------------------------------------------------
# Creates a tile datset filled with NA's, returns the path
CreateFillTile <- function(tile, tmp_dir, nbar_data_dir="/projectnb/modislc/data/mcd12_in/c5/mcd43a4"){
	# check to see if there's already a fill tile
	fill_tile_path <- file.path(tmp_dir, paste(tile, "fill.tif", sep=""))
	if(!file.exists(fill_tile_path)){
		# get an example file to copy by finding a random tile date
		set.seed(42)
		keep_going <- T
		while(keep_going){
			tmp_file <- Sys.glob(paste(nbar_data_dir, sample(2002:2014, 1), sample(formatC(seq(1, 361, by=8), width=3, flag="0"), 1), paste("*", tile, "*.hdf", sep=""), sep="/"))
			if(length(tmp_file) == 0){
				keep_going <- T
			}else{
				keep_going <- F
			}
		}

		# copy the example raster and set all values to NA
		fill_tile <- raster(Get_SDS_name_band1(tmp_file))
		fill_tile <- setValues(fill_tile, NA)

		writeRaster(fill_tile, fill_tile_path)
	}

	return(fill_tile_path)
}

#---------------------------------------------------------------------
# Keith Ma's benchmarking found this to be much faster than raster's getValues()
# this version doesn't know anything about a maximum number of simultaneously open files
getValuesGDAL_old <- function(dsets, start_row, n) {
	# open all files
	gds <- list()
	for (i in 1:length(dsets)) { gds <- c(gds, GDAL.open(dsets[i])) }

	nrows <- dim(gds[[1]])[1]
	ncols <- dim(gds[[1]])[2]
	# sumval <- 0.0
	out_vals <- matrix(NA, nrow=nrows * n, ncol=length(dsets))

	for (j in 1:length(gds)) {
		val <- getRasterData(gds[[j]], offset = c((start_row - 1)*n, 0), region.dim = c(n, ncols), as.is = FALSE)
		out_vals[, j] <- c(val)
		# if (dbg == TRUE ) sumval <- sumval+sum(val, na.rm = TRUE)
	}

	# close all files
	for (i in 1:length(gds)) { GDAL.close(gds[[i]]) }

	# return
	return(out_vals)
}

#---------------------------------------------------------------------
# extract n rows from all dsets starting at start_row
getValuesGDAL <- function(dsets, start_row, n, max_open_datasets=2.75e3) {
	# determine how many blocks of datasets to open simultaneously and dataset size
	num_blocks <- ceiling(length(dsets) / max_open_datasets)
	nrows <- nrow(raster(dsets[1]))
	ncols <- ncol(raster(dsets[1]))

	# initialize output matrix
	out_vals <- matrix(NA, nrow=nrows * n, ncol=length(dsets))

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
			val <- getRasterData(gds[[j]], offset = c((start_row - 1)*n, 0), region.dim = c(n, ncols), as.is = FALSE)
			out_vals[, (dset_start + j - 1)] <- c(val)
		}
		# close all files
		for (i in 1:length(gds)) { GDAL.close(gds[[i]]) }
	}
	return(out_vals)
}

#---------------------------------------------------------------------
# For reading from the binary, splined output files which contain 365 bands (days) of spline-smoothed EVI2 data/resids/QA
read_int_file <- function(in_file, lines_to_read=2400, start_line=1, samples=2400, dsize=2, nbands=365, s_flag=T, bsq_flag=F){
  # start_line begins at 1!!!
  where_to_start <- (start_line - 1) * samples * dsize * nbands
  size_to_read <- lines_to_read * nbands * samples
  f <- file(in_file, "rb") # open the file
  seek(f, where_to_start) # position the file connection
  temp <- readBin(f, integer(), n=size_to_read, endian= "little", size=dsize, signed=s_flag) # read all of the data into an array
  close(f) # close the file
  # re-order the temp array into a matrix; if bsq_flag then read by row, otherwise by column
  ifelse(bsq_flag, byr_flag<-FALSE, byr_flag<-TRUE)
  temp <- matrix(temp, ncol=nbands, byrow=byr_flag)
  return(temp) # return the data
}

#---------------------------------------------------------------------
# takes a vector of MCD43A4 filenames and returns a vector of their dates
DateFunc <- function(x){
	as.Date(strsplit(basename(x),split="\\.")[[1]][2], format="A%Y%j")
}

#---------------------------------------------------------------------
# calculates time series of evi2 and ndsi, then replaces evi2 values where
# ndsi > ndsi_thresh with the ev2_snow_quant quantile value of snow-free evi2
# expects a vector of time series, one after the other and of equal length:
# red, nir, blue, band4, and band6
SnowFilterEVI2 <- function(x, evi2_snow_quant=0.02, ndsi_thresh=0, throttle_min=T, hampel_thresh=NULL, hampel_window=3){
	# get indices in x of each component
	band_length <- length(x) / 5
	red_inds <- 1:band_length
	nir_inds <- (band_length + 1):(2 * band_length)
	band4_inds <- ((3 * band_length) + 1):(4 * band_length)
	band6_inds <- ((4 * band_length) + 1):(5 * band_length)

	# calculate EVI2 and NDSI
	evi2 <- 2.5 * ((x[nir_inds] - x[red_inds]) / (x[nir_inds] + (2.4 * x[red_inds]) + 1))
	ndsi <- (x[band4_inds] - x[band6_inds]) / (x[band4_inds] + x[band6_inds])

	# make a vector to hold the snow fill information
	snow_fills <- rep(0, length(evi2))

	# fill snow-corrupted EVI2 values
	min_snow_free_evi2 <- quantile(evi2[ndsi <= ndsi_thresh], evi2_snow_quant, na.rm=T)
	# snow_fill_count <- sum(ndsi > ndsi_thresh, na.rm=T) # number of filled evi2 values due to snow
	evi2[ndsi > ndsi_thresh] <- min_snow_free_evi2
	snow_fills[ndsi > ndsi_thresh] <- 1 # flag snow filled values

	# replace any values below min_snow_free_evi2 with min_snow_free_evi2
	if(throttle_min){
		evi2[evi2 < min_snow_free_evi2] <- min_snow_free_evi2
		# snow_fill_count <- snow_fill_count + sum(evi2 < min_snow_free_evi2, na.rm=T)
		snow_fills[evi2 < min_snow_free_evi2] <- 2 # flag minimum filled values
	}

	# do Hampel filtering for spike removal
	if(!is.null(hampel_thresh)){
		tmp_hampel <- hampel(evi2[!is.na(evi2)], k=hampel_window, t0=hampel_thresh)
		evi2[!is.na(evi2)] <- tmp_hampel$y
	}

	# return(evi2)
	return(list(evi2, snow_fills))
}

#---------------------------------------------------------------------
SnowFilterEVI2_mv <- function(x, dates, evi2_snow_quant=0.02, ndsi_thresh=0, throttle_min=T, hampel_thresh=NULL, hampel_window=1, quantile_method=3, max_snow_fill_ratio=1.25){
		# get indices in x of each component
		band_length <- length(x) / 5
		red_inds <- 1:band_length
		nir_inds <- (band_length + 1):(2 * band_length)
		band4_inds <- ((3 * band_length) + 1):(4 * band_length)
		band6_inds <- ((4 * band_length) + 1):(5 * band_length)

		# calculate EVI2 and NDSI
		evi2 <- 2.5 * ((x[nir_inds] - x[red_inds]) / (x[nir_inds] + (2.4 * x[red_inds]) + 1))
		ndsi <- (x[band4_inds] - x[band6_inds]) / (x[band4_inds] + x[band6_inds])

		# make a vector to hold the snow fill information
		snow_fills <- rep(0, length(evi2))

		# fill snow values with NA (initially, later we fill with a min snow free value)
		evi2[ndsi >= ndsi_thresh] <- NA
		snow_fills[ndsi >= ndsi_thresh] <- 1

		# do Hampel filtering for spike removal
		if(!is.null(hampel_thresh)){
			evi2_hampel <- hampel(evi2[!is.na(evi2)], k=hampel_window, t0=hampel_thresh)
			evi2[!is.na(evi2)] <- evi2_hampel$y

			# mark outlier eliminated data
			snow_fills[!is.na(evi2)][evi2_hampel$ind] <- 3
		}

		# apply a moving window to fill snow values
		# gets the previous and next year, computes min_evi2_snow, and fills the year of interest
		years <- as.numeric(strftime(dates, format="%Y"))
		unique_years <- sort(unique(years))
		for(year_of_interest in unique_years){
			# Debug
			# print(year_of_interest)

			tmp_inds <- (year_of_interest - 1) <= years & years <= (year_of_interest + 1)

			tmp_evi2 <- evi2[tmp_inds]
			tmp_ndsi <- ndsi[tmp_inds]
			tmp_snow_fills <- snow_fills[tmp_inds]

			# min_snow_free_evi2 <- quantile(tmp_evi2[tmp_ndsi <= ndsi_thresh], evi2_snow_quant, na.rm=T)
			# min_snow_free_evi2 <- quantile(tmp_evi2[tmp_ndsi <= ndsi_thresh], evi2_snow_quant, na.rm=T, type=quantile_method)
			min_snow_free_evi2_3year <- quantile(tmp_evi2[tmp_ndsi <= ndsi_thresh], evi2_snow_quant, na.rm=T, type=quantile_method)
			min_snow_free_evi2_1year <- quantile(tmp_evi2[years[tmp_inds] == year_of_interest & tmp_ndsi <= ndsi_thresh], evi2_snow_quant, na.rm=T, type=quantile_method)

			# check if the 3-year snow fill is substantially different than the 1 year value
			# if it is, we fill with the 1-year value, if not, the 3-year
			if((max(c(min_snow_free_evi2_3year, min_snow_free_evi2_1year)) / min(c(min_snow_free_evi2_3year, min_snow_free_evi2_1year))) <= max_snow_fill_ratio){
				min_snow_free_evi2 <- min_snow_free_evi2_3year
			}else{
				min_snow_free_evi2 <- min_snow_free_evi2_1year
			}

			tmp_evi2[tmp_ndsi > ndsi_thresh] <- min_snow_free_evi2
			tmp_snow_fills[tmp_ndsi > ndsi_thresh] <- 1 # flag snow filled values

			# replace any values below min_snow_free_evi2 with the min_snow_free_evi2
			if(throttle_min){
				tmp_evi2[tmp_evi2 < min_snow_free_evi2] <- min_snow_free_evi2
				tmp_snow_fills[tmp_evi2 < min_snow_free_evi2] <- 2
			}

			# fill year_of_interest with values
			evi2[years == year_of_interest] <- tmp_evi2[years[tmp_inds] == year_of_interest]
			ndsi[years == year_of_interest] <- tmp_ndsi[years[tmp_inds] == year_of_interest]
			snow_fills[years == year_of_interest] <- tmp_snow_fills[years[tmp_inds] == year_of_interest]
		}

		return(list(evi2, snow_fills))
	}



#---------------------------------------------------------------------
SnowFilterEVI2_mv_test <- function(evi2, ndsi, dates, evi2_snow_quant=0.02, ndsi_thresh=0, throttle_min=T, hampel_thresh=NULL, hampel_window=1, quantile_method=3, max_snow_fill_ratio=1.25){
	# get indices in x of each component
	# band_length <- length(x) / 5
	# red_inds <- 1:band_length
	# nir_inds <- (band_length + 1):(2 * band_length)
	# band4_inds <- ((3 * band_length) + 1):(4 * band_length)
	# band6_inds <- ((4 * band_length) + 1):(5 * band_length)
	#
	# # calculate EVI2 and NDSI
	# evi2 <- 2.5 * ((x[nir_inds] - x[red_inds]) / (x[nir_inds] + (2.4 * x[red_inds]) + 1))
	# ndsi <- (x[band4_inds] - x[band6_inds]) / (x[band4_inds] + x[band6_inds])

	# make a vector to hold the snow fill information
	snow_fills <- rep(0, length(evi2))

	# fill snow values with NA (initially, later we fill with a min snow free value)
	#DEBUG
	print(paste("NDSI_THRESH:", ndsi_thresh))

	evi2[ndsi >= ndsi_thresh] <- NA
	snow_fills[ndsi >= ndsi_thresh] <- 1

	# do Hampel filtering for spike removal
	if(!is.null(hampel_thresh)){
		evi2_hampel <- hampel(evi2[!is.na(evi2)], k=hampel_window, t0=hampel_thresh)
		evi2[!is.na(evi2)] <- evi2_hampel$y

		# mark outlier eliminated data
		snow_fills[!is.na(evi2)][evi2_hampel$ind] <- 3
	}

	# apply a moving window to fill snow values
	# gets the previous and next year, computes min_evi2_snow, and fills the year of interest
	years <- as.numeric(strftime(dates, format="%Y"))
	unique_years <- sort(unique(years))
	for(year_of_interest in unique_years){
		tmp_inds <- (year_of_interest - 1) <= years & years <= (year_of_interest + 1)

		tmp_evi2 <- evi2[tmp_inds]
		tmp_ndsi <- ndsi[tmp_inds]
		tmp_snow_fills <- snow_fills[tmp_inds]

		# min_snow_free_evi2 <- quantile(tmp_evi2[tmp_ndsi <= ndsi_thresh], evi2_snow_quant, na.rm=T)
		# min_snow_free_evi2 <- quantile(tmp_evi2[tmp_ndsi <= ndsi_thresh], evi2_snow_quant, na.rm=T, type=quantile_method)
		min_snow_free_evi2_3year <- quantile(tmp_evi2[tmp_ndsi <= ndsi_thresh], evi2_snow_quant, na.rm=T, type=quantile_method)
		min_snow_free_evi2_1year <- quantile(tmp_evi2[years[tmp_inds] == year_of_interest & tmp_ndsi <= ndsi_thresh], evi2_snow_quant, na.rm=T, type=quantile_method)

		# check if the 3-year snow fill is substantially different than the 1 year value
		# if it is, we fill with the 1-year value, if not, the 3-year
		if((max(c(min_snow_free_evi2_3year, min_snow_free_evi2_1year)) / min(c(min_snow_free_evi2_3year, min_snow_free_evi2_1year))) <= max_snow_fill_ratio){
			min_snow_free_evi2 <- min_snow_free_evi2_3year
		}else{
			min_snow_free_evi2 <- min_snow_free_evi2_1year
		}

		tmp_evi2[tmp_ndsi > ndsi_thresh] <- min_snow_free_evi2
		tmp_snow_fills[tmp_ndsi > ndsi_thresh] <- 1 # flag snow filled values

		# replace any values below min_snow_free_evi2 with the min_snow_free_evi2
		if(throttle_min){
			tmp_evi2[tmp_evi2 < min_snow_free_evi2] <- min_snow_free_evi2
			tmp_snow_fills[tmp_evi2 < min_snow_free_evi2] <- 2
		}

		# fill year_of_interest with values
		evi2[years == year_of_interest] <- tmp_evi2[years[tmp_inds] == year_of_interest]
		ndsi[years == year_of_interest] <- tmp_ndsi[years[tmp_inds] == year_of_interest]
		snow_fills[years == year_of_interest] <- tmp_snow_fills[years[tmp_inds] == year_of_interest]
	}

	return(list(evi2, snow_fills))
}

#----------------------------------------------------------
#----------------------------------------------------------
# Functions for smoothing Landsat time series
#----------------------------------------------------------
IsSpike <- function(x, thresh=2, minResid=0){
  # this function detects negative spikes with an ad-hoc, 3-point method.
	# the deviation of the middle point from the line between the first and last
	# point is calculated. If it is negative, the absolute value of the
	# ratio of this deviation and the change between the first and last points is
	# calculated. If the ratio is greater than "thresh", and the deviation is greater
	# than minResid, then the point is considered a negative spike and T is returned
  # y <- x
  # x <- 1:3
  # dev <- y[2] - ((((y[3] - y[1]) * (x[2] - x[1])) / (x[3] - x[1])) + y[1])
	dev <- x[2] - ((x[3] - x[1]) / 2) - x[1]
  # devRatio <- dev / (y[3] - y[1])
	devRatio <- dev / (x[3] - x[1])
  ifelse((abs(dev) > minResid) & (dev < 0) & (abs(devRatio) > thresh), T, F)
}

CheckSpike <- function(x, thresh=2, minResid=0){
	# applies the IsSpike function to a time series
	x_og <- x # preserve original vector
	x_outs <- rep(F, length(x_og)) # create the outlier output vector
	x <- x_og[!is.na(x_og)] # subset to non missing values
	# apply the IsSpike function on a rolling 3-point window
	outs <- rollapply(x, 3, IsSpike, thresh=thresh, minResid=minResid, partial=T)
	outs[is.na(outs)] <- F # replace NAs with FALSE
	x_outs[!is.na(x_og)] <- outs # expand to size of original x, accounting for missing values
	return(x_outs)
}

LandsatSmoother <- function(x, dates, fill_quant=0.05, x_min=0, spike_thresh=2, min_resid=0.1, spline_spar=NULL, fill_doy=NULL, doAnnual=T, padHeadTail=T){
	# special function to screen/smooth/fill the Landsat time series

	# if there are less than 4 valid values, we can't fit a spline
	if(sum(!is.na(x)) < 4) return(NA)

	# create daily dates and VI vector
	if(doAnnual){
		# if doAnnual is TRUE, then the pred dates are for full calendar years, regardless of the min/max of dates
		pred_dates <- seq(as.Date(paste(strftime(min(dates), format="%Y"), "-1-1", sep="")), as.Date(paste(strftime(max(dates), format="%Y"), "-12-31", sep="")), by="day")
	}else{
		pred_dates <- seq(min(dates), max(dates), by="day")
	}
	x_doy <- rep(NA, length(pred_dates))
	x_doy[pred_dates %in% dates] <- x # assign original data to proper dates

  # despike
	x_tmp <- x_doy[!is.na(x_doy)] # remove missing values for despiking
	outliers <- CheckSpike(x_tmp, thresh=spike_thresh, minResid=min_resid)
  x_tmp[outliers] <- NA
  x_doy[!is.na(x_doy)] <- x_tmp # replace the original values with outlier screened values

  # fill min values
  minVI <- quantile(x_doy[x_doy > x_min], fill_quant, na.rm=T)
  x_doy[x_doy < minVI] <- minVI

	# pad the head and tail with minVI, if requested
	if(padHeadTail){
		x_doy[1:min(which(!is.na(x_doy)))] <- minVI
		x_doy[max(which(!is.na(x_doy))):length(x_doy)] <- minVI
	}

	# fill dormant period if requested, each period is filled with 3 values only
	num_fill_values <- 3
	if(!is.null(fill_doy)){
		doy <- as.numeric(strftime(pred_dates, format="%j"))
		for(fill_seg in fill_doy){
			fill_doys <- c(fill_seg[1], sample((fill_seg[1] + 1):(fill_seg[2] - 1), (num_fill_values - 2)), fill_seg[2])
			x_doy[doy %in% fill_doys] <- minVI
		}
	}

	# experiment with weighting high values the most
	# xe <- ecdf(x_doy)
	# w <- xe(x_doy)
	# w[w < 0.7 & !is.na(w)] <- 1
	# w[w >= 0.7 & !is.na(w)] <- w[w >= 0.7 & !is.na(w)] * 10

	# check once more for too few non-missing values
	if(sum(!is.na(x_doy)) < 4) return(NA)

  # smooth with a spline to get continuous daily series
  # pred_dates <- seq(min(dates), max(dates), by="day")
  spl <- smooth.spline(pred_dates[!is.na(x_doy)], x_doy[!is.na(x_doy)], spar=spline_spar)
	# weighted version
	# spl <- smooth.spline(pred_dates[!is.na(x_doy)], x_doy[!is.na(x_doy)], w=w[!is.na(x_doy)], spar=spline_spar)
  x_smooth <- predict(spl, as.numeric(pred_dates))$y

  # screen again for low values
  x_smooth[x_smooth < minVI] <- minVI

  return(list(pred_dates, x_smooth))
}

#---------------------------------------------------------------------
# NOTE: there's no outlier identification or removal in SmoothSnowFill!
SmoothSnowFill <- function(evi2, snow_fills, dates, pred_dates, pheno_pars, plot=F){
	# check for less than two values necessary for interpolation
	if(sum(!is.na(evi2)) < 2) return(NA)

  # Plotting
	if(plot){
		plot(dates, evi2, type="n")
		points(dates[!snow_fills], evi2[!snow_fills], pch=1)
		points(dates[snow_fills], evi2[snow_fills], pch=4, col=4)
	}

  # initialize the snow_fill value
  snow_fill <- NA
  # get the "minimum" value among snow-free observations
  snow_free_min <- quantile(evi2[!snow_fills], pheno_pars$snow_free_min_quant, na.rm=T)
  # get the "maximum" value among snow-present observations
  snow_max <- quantile(evi2[snow_fills], pheno_pars$snow_max_quant, na.rm=T)
  # the fill value is the maximum among snow_free_min, and snow_max
  snow_fill <- max(snow_free_min, snow_max, na.rm=T)
  # everyting less than the snow_fill value is filled with the snow_fill value
  if(!is.na(snow_fill)) evi2[evi2 < snow_fill] <- snow_fill

  # Plotting
	if(plot){
		abline(h=snow_free_min, lty=3, col=4)
		abline(h=snow_max, lty=4, col=4)
	}

  # fit a smoothing spline and predict over the pred_dates
  smooth_evi2 <- predict(smooth.spline(dates[!is.na(evi2)], evi2[!is.na(evi2)]), pred_dates)[[2]]
	# smooth_evi2 <- predict(smooth.spline(dates[!is.na(evi2)], evi2[!is.na(evi2)]), as.numeric(pred_dates))[[2]]

	# fill any post-spline values less than snow_fill to snow_fill
	smooth_evi2[smooth_evi2 < snow_fill] <- snow_fill

  # Plotting
	if(plot) points(pred_dates, smooth_evi2, col=2, type="l", lty=2)

  return(smooth_evi2)
}


#---------------------------------------------------------------------
# This function takes a time series w/ dates (x, dates) and returns a spline smoothed time series with outliers removed.
# Outliers are identified as points with absolute value more than out_sigma * sd, where sd is the residual
# standard deviation between the input data and the initial spline fit, and out_sigma is a variable
# coefficient. The spline smoothing parameter spline_spar controls the smoothness of the fit (see spline.smooth help)
# and out_iterations controls the number of times that outliers are checked and removed w/ subsequent spline refit
# pred_dates is a vector of dates where spline smoothed predictions of x are desired. If NA, then a daily series spanning
# min(dates)-max(dates) is returned
# SplineAndOutlierRemoval <- function(x, dates, out_sigma=3, spline_spar=0.3, out_iterations=1, pred_dates=NA){
SplineAndOutlierRemoval <- function(x, dates, out_quant=0.99, spline_spar=0.3, out_iterations=1, pred_dates=NA, interp=T, throttle_min=T){
	# check for all NA
	# if(all(is.na(x))) return(NA)
	# check for less than two values necessary for interpolation
	if(sum(!is.na(x)) < 2) return(NA)

	# get the minimum value for later throttling
	x_min <- min(x, na.rm=T)

	# first, we linearly interplate missing values if requested
	if(interp){
		# we interpolate missing values w/ linear zoo function
		x <- coredata(na.approx(zoo(as.numeric(x), order.by=dates), na.rm=F))
	}


	# convert dates
	dates <- as.integer(dates) # spline doesn't work with dates

	# if prediction dates aren't provided, we assume we want daily ones
	if(all(is.na(pred_dates)))
		pred_dates <- min(dates, na.rm=T):max(dates, na.rm=T)

	if(out_iterations == 0){
		spl <- try(smooth.spline(dates[!is.na(x)], x[!is.na(x)], spar=spline_spar), silent=T)
		if(inherits(spl, 'try-error')){
			# print("Failed to fit smoothing spline")
			return(NA)
		}

		smooth_x <- try(predict(spl, pred_dates)$y, silent=T) # calculate spline smoothed values
		if(inherits(smooth_x, 'try-error')){
			# print("Failed to predict with spline")
			return(NA)
		}

		if(throttle_min){
			smooth_x[smooth_x < x_min] <- x_min
		}

		return(smooth_x)
	}

	# eliminate outliers and respline
	for(i in 1:out_iterations){
		# fit a smoothing spline to non-missing data
		spl <- try(smooth.spline(dates[!is.na(x)], x[!is.na(x)], spar=spline_spar), silent=T)
		if(inherits(spl, 'try-error')){
			# print("Failed to fit smoothing spline")
			return(NA)
		}

		smooth_x <- try(predict(spl, dates)$y, silent=T) # calculate spline smoothed values
		if(inherits(smooth_x, 'try-error')){
			# print("Failed to predict with spline")
			return(NA)
		}

		# smooth_x_resid <- x - smooth_x # calculate residuals from spline
		# smooth_x_resid_sd <- try(sd(smooth_x_resid, na.rm=T), silent=T) # standard dev of absolute value of residuals
		# if(inherits(smooth_x_resid_sd, 'try-error')){
		# 	print("Failed to get sd of residuals")
		# 	return(NA)
		# }
		#
		# outliers <- abs(smooth_x_resid) > out_sigma * smooth_x_resid_sd
		# outliers[is.na(outliers)] <- F

		smooth_x_resid <- abs(x - smooth_x) # calculate absolute value of residuals from spline
		outliers <- smooth_x_resid >= quantile(smooth_x_resid, out_quant, na.rm=T)
		outliers[is.na(outliers)] <- F

		if(sum(outliers) > 0){
			# if we found outliers, eliminate them in x and refit up to iterations
			x[outliers] <- NA
		}else{
			# if we didn't find any outliers, we abandon the iteration and return the smoothed values
			smooth_x_return <- try(predict(spl, pred_dates)$y, silent=T)
			if(inherits(smooth_x_return, 'try-error')){
				# print("No outliers, but failed to predict with final spline")
				return(NA)
			}else{
				# if we need to throttle, do that now
				if(throttle_min){
					smooth_x_return[smooth_x_return < x_min] <- x_min
				}
				return(smooth_x_return)
			}
		}
	}

	# fit the spline to the outlier screened data, then return the predicted series
	spl <- try(smooth.spline(dates[!is.na(x)], x[!is.na(x)], spar=spline_spar), silent=T)
	if(inherits(spl, 'try-error')){
		# print("Failed to predict with final spline")
		return(NA)
	}else{
		smooth_x_return <- try(predict(spl, pred_dates)$y, silent=T)
		if(inherits(smooth_x_return, 'try-error')){
			return(NA)
		}else{
			# if we need to throttle, do that now
			if(throttle_min){
				smooth_x_return[smooth_x_return < x_min] <- x_min
			}
			return(smooth_x_return)
		}
	}
}

#==============================================================================
#==============================================================================
# segmentation functions
#==============================================================================
#==============================================================================

#---------------------------------------------------------------------
# Example of using the pracma package's "findpeaks"
# plotSeg <- function(findpeak_result, t, gup_col=rgb(0,1,0,0.5), gdown_col=rgb(1,0,0,0.5)){
# 	gup_x <- c(rep(t[findpeak_result[3]], 2), rep(t[findpeak_result[2]], 2))
# 	gup_y <- c(par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3])
# 	polygon(gup_x, gup_y, col=gup_col, border=NA)
#
# 	gdown_x <- c(rep(t[findpeak_result[2]], 2), rep(t[findpeak_result[4]], 2))
# 	gdown_y <- c(par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3])
# 	polygon(gdown_x, gdown_y, col=gdown_col, border=NA)
#
# }
#
# tmp_peaks <- findpeaks(evi2_smoothed/1e3, nups=14, ndowns=14, minpeakheight=0.3, threshold=0.05, minpeakdistance=30)
#
# plot(pred_dates,evi2_smoothed/1e3,type="l")
# apply(tmp_peaks, 1, plotSeg, t=pred_dates)

#---------------------------------------------------------------------
# returns ascending segments in x
GetSegs_snow_troughs <- function(x, peaks, max_seg_length, min_seg_amplitude, agg_amp_frac=0.15){
	# check if peaks is NA and return NA if so
	if(all(is.na(peaks))) return(NA)

	# find ascending nodes by locating the minimum value before the
	# peak and within search window, checking the amplitude
	ascending_segs <- NULL
	for(i in 1:length(peaks)){
		# define search window
		search_window_start <- max(
			1,
			peaks[i - 1],
			peaks[i] - max_seg_length
		)

		# find troughs in search window, return one closest to peak that meets min_seg_amplitude requirement,
		# if there are none, find minimum value
		troughs <- GetLocalMaxMin(x[search_window_start:peaks[i]], minima=T)
		if(length(troughs) == 0){
			# no troughs, use min value in search window
			potential_trough <- search_window_start + which(x[search_window_start:peaks[i]] == min(x[search_window_start:peaks[i]], na.rm=T))[1] - 1

			# check amplitude, add segment to list if it covers minimum amplitude. Otherwise, do nothing and move on
			if(x[peaks[i]] - x[potential_trough] > min_seg_amplitude){
				if(is.null(ascending_segs)){
					ascending_segs <- list(c(potential_trough, peaks[i]))
				}else{
					ascending_segs <- c(ascending_segs, list(c(potential_trough, peaks[i])))
				}
			}
		}else{
			# so there are troughs, we want to find the one furthest from the peak that yields a segment that is both
			# greater than min_seg_amplitude, and increases the amplitude more than the agg_amp_frac.
			seg_amp <- 0
			true_trough <- NULL
			for(trough in rev(troughs)){
				potential_trough <- trough + search_window_start
				tmp_amp <- x[peaks[i]] - x[potential_trough]
				if((tmp_amp >= min_seg_amplitude) & (tmp_amp >= ((1 + agg_amp_frac) * seg_amp))){
					true_trough <- potential_trough
					seg_amp <- tmp_amp
				}
			}
			if(!is.null(true_trough)){
				if(is.null(ascending_segs)){
					ascending_segs <- list(c(true_trough, peaks[i]))
				}else{
					ascending_segs <- c(ascending_segs, list(c(true_trough, peaks[i])))
				}
			}
		}
	}
	if(!is.null(ascending_segs)){
		return(ascending_segs)
	}else{
		return(NA)
	}
}

##################################################
# returns ascending segments in x
# must return 1) integrated segment SVI, 2) min segment SVI, 3) max segment SVI, 4) spline fit RMSE, 5) missing data in segment
GetFullSegs_new <- function(x, peaks, max_seg_length, min_seg_amplitude, agg_amp_frac=0.15){
	# check if peaks is NA and return NA if so
	if(all(is.na(peaks))) return(NA)

	# find ascending nodes by locating the minimum/trough value before the
	# peak and within search window, checking the amplitude
	full_segs <- NULL
	tmp_seg <- NA

	for(i in 1:length(peaks)){
		#################################
		# First, find ascending segments
		# define search window
		search_window_start <- max(
			1,
			peaks[i - 1],
			peaks[i] - max_seg_length
		)

		# find troughs in search window, return one closest to peak that meets min_seg_amplitude requirement,
		# if there are none, find minimum value
		troughs <- GetLocalMaxMin(x[search_window_start:peaks[i]], minima=T)
		if(length(troughs) == 0){
			# no troughs, use min value in search window
			potential_trough <- search_window_start + which(x[search_window_start:peaks[i]] == min(x[search_window_start:peaks[i]], na.rm=T))[1] - 1

			# check amplitude, add segment to list if it covers minimum amplitude. Otherwise, do nothing and move on
			tmp_amp <- x[peaks[i]] - x[potential_trough]
			# if(x[peaks[i]] - x[potential_trough] > min_seg_amplitude){
			if(!is.na(tmp_amp)){
				if(x[peaks[i]] - x[potential_trough] > min_seg_amplitude){
					tmp_seg <- list(c(potential_trough, peaks[i]))
				}
			}
		}else{
			# so there are troughs, we want to find the one furthest from the peak that yields a segment that is both
			# greater than min_seg_amplitude, and increases the amplitude more than the agg_amp_frac.
			seg_amp <- 0
			true_trough <- NULL
			for(trough in rev(troughs)){
				potential_trough <- trough + search_window_start
				tmp_amp <- x[peaks[i]] - x[potential_trough]

				if(!is.na(tmp_amp)){
					if((tmp_amp >= min_seg_amplitude) & (tmp_amp >= ((1 + agg_amp_frac) * seg_amp))){
						true_trough <- potential_trough
						seg_amp <- tmp_amp
					}
				}
			}
			if(!is.null(true_trough)){
				tmp_seg <- c(true_trough, peaks[i])
			}
		}# end check for troughs in search window

		#################################
		# Find end of segment
		# if we found an ascending segment, then we need to find an associated trough/min point after the peak
		if(!all(is.na(tmp_seg))){
			search_window_end <- min(
				length(x),
				peaks[i + 1], # this can result in NA when it reads off the end
				peaks[i] + max_seg_length,
				na.rm=T
			)

			# find troughs in search window, return one closest to peak that meets min_seg_amplitude requirement,
			# if there are none, find minimum value
			troughs <- GetLocalMaxMin(x[peaks[i]:search_window_end], minima=T)

			if(length(troughs) == 0){
				# no troughs, use min value in search window
				potential_trough <- peaks[i] + which(x[peaks[i]:search_window_end] == min(x[peaks[i]:search_window_end], na.rm=T))[1] - 1

				# check amplitude, add segment to list if it covers minimum amplitude. Otherwise, do nothing and move on
				if(x[peaks[i]] - x[potential_trough] > min_seg_amplitude){
					if(is.null(full_segs)){
						full_segs <- list(c(tmp_seg, potential_trough))
					}else{
						full_segs <- c(full_segs, list(c(tmp_seg, potential_trough)))
					}
				}
			}else{
				# so there are troughs, we want to find the one furthest from the peak that yields a segment that is both
				# greater than min_seg_amplitude, and increases the amplitude more than the agg_amp_frac.
				seg_amp <- 0
				true_trough <- NULL
				for(trough in troughs){
					potential_trough <- peaks[i] + trough - 1
					tmp_amp <- x[peaks[i]] - x[potential_trough]
					if(!is.na(tmp_amp)){
						if((tmp_amp >= min_seg_amplitude) & (tmp_amp >= ((1 + agg_amp_frac) * seg_amp))){
							true_trough <- potential_trough
							seg_amp <- tmp_amp
						}
					}
				}
				if(!is.null(true_trough)){
					if(is.null(full_segs)){
						full_segs <- list(c(tmp_seg, true_trough))
					}else{
						full_segs <- c(full_segs, list(c(tmp_seg, true_trough)))
					}
				}
			}# end check for troughs in search window
		}# end check for valid ascending segment

		tmp_seg <- NA # reset tmp_seg to NA

	} # end for
	if(!is.null(full_segs)){
		return(full_segs)
	}else{
		return(NA)
	}
}



# returns ascending segments in x
GetFullSegs <- function(x, peaks, max_seg_length, min_seg_amplitude, agg_amp_frac=0.15){
	# check if peaks is NA and return NA if so
	if(all(is.na(peaks))) return(NA)

	# find ascending nodes by locating the minimum/trough value before the
	# peak and within search window, checking the amplitude
	full_segs <- NULL
	tmp_seg <- NA

	for(i in 1:length(peaks)){
		#################################
		# First, find ascending segments
		# define search window
		search_window_start <- max(
			1,
			peaks[i - 1],
			peaks[i] - max_seg_length
		)

		# find troughs in search window, return one closest to peak that meets min_seg_amplitude requirement,
		# if there are none, find minimum value
		troughs <- GetLocalMaxMin(x[search_window_start:peaks[i]], minima=T)
		if(length(troughs) == 0){
			# no troughs, use min value in search window
			potential_trough <- search_window_start + which(x[search_window_start:peaks[i]] == min(x[search_window_start:peaks[i]], na.rm=T))[1] - 1

			# check amplitude, add segment to list if it covers minimum amplitude. Otherwise, do nothing and move on
			tmp_amp <- x[peaks[i]] - x[potential_trough]
			# if(x[peaks[i]] - x[potential_trough] > min_seg_amplitude){
			if(!is.na(tmp_amp)){
				if(x[peaks[i]] - x[potential_trough] > min_seg_amplitude){
					tmp_seg <- list(c(potential_trough, peaks[i]))
				}
			}
		}else{
			# so there are troughs, we want to find the one furthest from the peak that yields a segment that is both
			# greater than min_seg_amplitude, and increases the amplitude more than the agg_amp_frac.
			seg_amp <- 0
			true_trough <- NULL
			for(trough in rev(troughs)){
				potential_trough <- trough + search_window_start
				tmp_amp <- x[peaks[i]] - x[potential_trough]

				if(!is.na(tmp_amp)){
					if((tmp_amp >= min_seg_amplitude) & (tmp_amp >= ((1 + agg_amp_frac) * seg_amp))){
						true_trough <- potential_trough
						seg_amp <- tmp_amp
					}
				}
			}
			if(!is.null(true_trough)){
				tmp_seg <- c(true_trough, peaks[i])
			}
		}# end check for troughs in search window

		#################################
		# Find end of segment
		# if we found an ascending segment, then we need to find an associated trough/min point after the peak
		if(!all(is.na(tmp_seg))){
			search_window_end <- min(
				length(x),
				peaks[i + 1], # this can result in NA when it reads off the end
				peaks[i] + max_seg_length,
				na.rm=T
			)

			# find troughs in search window, return one closest to peak that meets min_seg_amplitude requirement,
			# if there are none, find minimum value
			troughs <- GetLocalMaxMin(x[peaks[i]:search_window_end], minima=T)

			if(length(troughs) == 0){
				# no troughs, use min value in search window
				potential_trough <- peaks[i] + which(x[peaks[i]:search_window_end] == min(x[peaks[i]:search_window_end], na.rm=T))[1] - 1

				# check amplitude, add segment to list if it covers minimum amplitude. Otherwise, do nothing and move on
				if(x[peaks[i]] - x[potential_trough] > min_seg_amplitude){
					if(is.null(full_segs)){
						full_segs <- list(c(tmp_seg, potential_trough))
					}else{
						full_segs <- c(full_segs, list(c(tmp_seg, potential_trough)))
					}
				}
			}else{
				# so there are troughs, we want to find the one furthest from the peak that yields a segment that is both
				# greater than min_seg_amplitude, and increases the amplitude more than the agg_amp_frac.
				seg_amp <- 0
				true_trough <- NULL
				for(trough in troughs){
					potential_trough <- peaks[i] + trough - 1
					tmp_amp <- x[peaks[i]] - x[potential_trough]
					if(!is.na(tmp_amp)){
						if((tmp_amp >= min_seg_amplitude) & (tmp_amp >= ((1 + agg_amp_frac) * seg_amp))){
							true_trough <- potential_trough
							seg_amp <- tmp_amp
						}
					}
				}
				if(!is.null(true_trough)){
					if(is.null(full_segs)){
						full_segs <- list(c(tmp_seg, true_trough))
					}else{
						full_segs <- c(full_segs, list(c(tmp_seg, true_trough)))
					}
				}
			}# end check for troughs in search window
		}# end check for valid ascending segment

		tmp_seg <- NA # reset tmp_seg to NA

	} # end for
	if(!is.null(full_segs)){
		return(full_segs)
	}else{
		return(NA)
	}
}

# returns ascending segments in x
GetFullSegs_old <- function(x, peaks, max_seg_length, min_seg_amplitude, agg_amp_frac=0.15){
	# check if peaks is NA and return NA if so
	if(all(is.na(peaks))) return(NA)

	# find ascending nodes by locating the minimum/trough value before the
	# peak and within search window, checking the amplitude
	full_segs <- NULL
	tmp_seg <- NA

	for(i in 1:length(peaks)){
		#################################
		# First, find ascending segments
		# define search window
		search_window_start <- max(
			1,
			peaks[i - 1],
			peaks[i] - max_seg_length
		)

		# find troughs in search window, return one closest to peak that meets min_seg_amplitude requirement,
		# if there are none, find minimum value
		troughs <- GetLocalMaxMin(x[search_window_start:peaks[i]], minima=T)
		if(length(troughs) == 0){
			# no troughs, use min value in search window
			potential_trough <- search_window_start + which(x[search_window_start:peaks[i]] == min(x[search_window_start:peaks[i]], na.rm=T))[1] - 1

			# check amplitude, add segment to list if it covers minimum amplitude. Otherwise, do nothing and move on
			if(x[peaks[i]] - x[potential_trough] > min_seg_amplitude){
				# if(is.null(full_segs)){
				# 	full_segs <- list(c(potential_trough, peaks[i]))
				# }else{
				# 	full_segs <- c(full_segs, list(c(potential_trough, peaks[i])))
				# }
				tmp_seg <- list(c(potential_trough, peaks[i]))
			}
		}else{
			# so there are troughs, we want to find the one furthest from the peak that yields a segment that is both
			# greater than min_seg_amplitude, and increases the amplitude more than the agg_amp_frac.
			seg_amp <- 0
			true_trough <- NULL
			for(trough in rev(troughs)){
				potential_trough <- trough + search_window_start
				tmp_amp <- x[peaks[i]] - x[potential_trough]

				if((tmp_amp >= min_seg_amplitude) & (tmp_amp >= ((1 + agg_amp_frac) * seg_amp))){
					true_trough <- potential_trough
					seg_amp <- tmp_amp
				}
			}
			if(!is.null(true_trough)){
				tmp_seg <- c(true_trough, peaks[i])
			}
		}# end check for troughs in search window

		#################################
		# Find end of segment
		# if we found an ascending segment, then we need to find an associated trough/min point after the peak
		if(!is.na(tmp_seg)){
			search_window_end <- min(
				length(x),
				peaks[i + 1], # this can result in NA when it reads off the end
				peaks[i] + max_seg_length,
				na.rm=T
			)

			# find troughs in search window, return one closest to peak that meets min_seg_amplitude requirement,
			# if there are none, find minimum value
			troughs <- GetLocalMaxMin(x[peaks[i]:search_window_end], minima=T)

			if(length(troughs) == 0){
				# no troughs, use min value in search window
				potential_trough <- peaks[i] + which(x[peaks[i]:search_window_end] == min(x[peaks[i]:search_window_end], na.rm=T))[1] - 1

				# check amplitude, add segment to list if it covers minimum amplitude. Otherwise, do nothing and move on
				if(x[peaks[i]] - x[potential_trough] > min_seg_amplitude){
					if(is.null(full_segs)){
						full_segs <- list(c(tmp_seg, potential_trough))
					}else{
						full_segs <- c(full_segs, list(c(tmp_seg, potential_trough)))
					}
				}
			}else{
				# so there are troughs, we want to find the one furthest from the peak that yields a segment that is both
				# greater than min_seg_amplitude, and increases the amplitude more than the agg_amp_frac.
				seg_amp <- 0
				true_trough <- NULL
				for(trough in troughs){
					potential_trough <- peaks[i] + trough - 1
					tmp_amp <- x[peaks[i]] - x[potential_trough]

					if((tmp_amp >= min_seg_amplitude) & (tmp_amp >= ((1 + agg_amp_frac) * seg_amp))){
						true_trough <- potential_trough
						seg_amp <- tmp_amp
					}
				}
				if(!is.null(true_trough)){
					if(is.null(full_segs)){
						full_segs <- list(c(tmp_seg, true_trough))
					}else{
						full_segs <- c(full_segs, list(c(tmp_seg, true_trough)))
					}
				}
			}# end check for troughs in search window
		}# end check for valid ascending segment

		tmp_seg <- NA # reset tmp_seg to NA

	} # end for
	if(!is.null(full_segs)){
		return(full_segs)
	}else{
		return(NA)
	}
}

#------------------------------------------------------------------
# returns the index of the first/last value  of x that is
# greater/less than the value of thresh. If gup is False (greendown)
# then it returns the first/last value of x that is less/greater than
# the value of thresh. first/last and greater/less determined by first_greater
# NOTE: if thresh is 1 or 0, rounding error can be a problem. Now we round
# the threshold and each of the evi values to 6 decimal places to compensate
GetThresh <- function(thresh_value, x, first_greater=T, gup=T){
	if(gup){
		if(first_greater){
			return(min(which(round(x, 6) >= round(thresh_value, 6))))
		}else{
			return(max(which(round(x, 6) <= round(thresh_value, 6))))
		}
	}else{
		if(first_greater){
			return(min(which(round(x, 6) <= round(thresh_value, 6))))
		}else{
			return(max(which(round(x, 6) >= round(thresh_value, 6))))
		}
	}
}

#------------------------------------------------------------------
GetSegThresh <- function(seg, x, thresh, gup=T){
	if(gup){
		# check for valid greenup segment
		if(!is.na(seg[1]) & !is.na(seg[2])){
			gup_thresh <- x[seg[1]] + ((x[seg[2]] - x[seg[1]]) * thresh)
			gup_thresh_index <- GetThresh(gup_thresh, x[seg[1]:seg[2]], first_greater=T, gup=T)
			return(gup_thresh_index + seg[1] - 1)
		}else{
			return(NA)
		}
	}else{
		# check for valid greendown segment
		if(!is.na(seg[2]) & !is.na(seg[3])){
			gdown_thresh <- x[seg[3]] + ((x[seg[2]] - x[seg[3]]) * thresh)
			gdown_thresh_index <- GetThresh(gdown_thresh, x[seg[2]:seg[3]], first_greater=T, gup=F)
			return(gdown_thresh_index + seg[2] - 1)
		}else{
			return(NA)
		}
	}
}
# example application
# gup_dates <- pred_dates[unlist(lapply(full_segs, GetSegThresh, x, 0.5, T))]

#------------------------------------------------------------------
GetPhenoDates <- function(segs, x, dates, gup_threshes, gdown_threshes){
	pheno_dates <- list()
	for(gup_thresh in gup_threshes){
		tmp_dates <- list(dates[unlist(lapply(segs, GetSegThresh, x, gup_thresh, gup=T), use.names=F)])
		if(all(is.na(unlist(tmp_dates, use.names=F)))){
			tmp_dates <- NA
		}
		pheno_dates <- c(pheno_dates, tmp_dates)
	}

	for(gdown_thresh in gdown_threshes){
		tmp_dates <- list(dates[unlist(lapply(segs, GetSegThresh, x, gdown_thresh, gup=F), use.names=F)])
		if(all(is.na(unlist(tmp_dates, use.names=F)))){
			tmp_dates <- NA
		}
		pheno_dates <- c(pheno_dates, tmp_dates)
	}
	return(list(pheno_dates))
}

#------------------------------------------------------------------
R2 <- function(x, resids){
	if(all(is.na(x)) | all(is.na(x + resids))) return(NA)
	return(cor(x, x + resids, use="complete")^2)
}

#------------------------------------------------------------------
SegMet <- function(seg, x, dates, pheno_pars){
	# prototype the return list (should probably be done as an object instead)
	seg_metrics=list(
		# phenometrics as indices of x
		ogi=NA,
		midgup=NA,
		mat=NA,
		peak=NA,
		sen=NA,
		midgdown=NA,
		dor=NA,

		# segment properties
		evi2_area=NA, #
		evi2_min=NA, # GUP min (not max amp when min(gdown)>min(gup), but consistent in the case of LCLUC)
		evi2_amp=NA, # GUP amplitude (not max amp when amp(gdown)>amp(gup), but consistent in the case of LCLUC)
		frac_filled_gup=NA, # fraction missing, snow, interpolated, or below-min filled
		frac_filled_gdown=NA,
		length_gup=NA,
		length_gdown=NA,
		spline_r2_gup=NA,
		spline_r2_gdown=NA,

		# phenometric QA scores
		ogi_qual=NA,
		midgup_qual=NA,
		mat_qual=NA,
		peak_qual=NA,
		sen_qual=NA,
		midgdown_qual=NA,
		dor_qual=NA
	)

	# unpack the x vector into component parts
	x[x==pheno_pars$nbar_NA_value] <- NA
	data_length <- length(x) / 3
	evi <- x[1:data_length] / pheno_pars$nbar_scale_factor
	evi_gup <- evi[seg[1]:seg[2]]
	evi_gdown <- evi[seg[2]:seg[3]]
	# residuals
	resids <- x[(data_length + 1):(2 * data_length)] / pheno_pars$nbar_scale_factor
	resids_gup <- resids[seg[1]:seg[2]]
	resids_gdown <- resids[seg[2]:seg[3]]
	# snow/fill flags
	# 0 is non-snow observed
	# 1 is snow observed
	# 2 is non-snow interpolated.  Is used to initialize array so values between stride are going to be 2 even if they are snow.
	# 3 is snow interpolated
	# 4 is observed but evi is below minimum evi so that it is filled.
	snowflags <- x[(2 * data_length + 1):(3 * data_length)]
	snowflags[snowflags == 2] <- 0 # NOTE: I think flag=2 only happens for stride-fills, which are irrelevant for QA purposes
	snowflags_gup <- snowflags[seg[1]:seg[2]]
	snowflags_gdown <- snowflags[seg[2]:seg[3]]

	# calculate GUP minimum, amplitude, and thresholds
	evi2_min_gup <- evi_gup[1]
	seg_metrics$evi2_min <- evi2_min_gup # set output values
	evi2_amp_gup <- evi_gup[length(evi_gup)] - evi2_min_gup
	seg_metrics$evi2_amp <- evi2_amp_gup
	ogi_thresh <- evi2_min_gup + evi2_amp_gup * pheno_pars$ogi_thresh
	midgup_thresh <- evi2_min_gup + evi2_amp_gup * pheno_pars$midgup_thresh
	mat_thresh <- evi2_min_gup + evi2_amp_gup * pheno_pars$mat_thresh

	# calculate GDOWN minimum, amplitude, and thresholds
	evi2_min_gdown <- evi_gdown[length(evi_gdown)]
	evi2_amp_gdown <- evi_gdown[1] - evi2_min_gdown
	sen_thresh <- evi2_min_gdown + evi2_amp_gdown * pheno_pars$sen_thresh
	midgdown_thresh <- evi2_min_gdown + evi2_amp_gdown * pheno_pars$midgdown_thresh
	dor_thresh <- evi2_min_gdown + evi2_amp_gdown * pheno_pars$dor_thresh

	# calculate evi area: segment area above segment minimum
	tmp_area <- evi[seg[1]:seg[2]] - evi2_min_gup
	tmp_area[tmp_area < 0] <- 0
	evi_area <- sum(tmp_area, na.rm=T)
	seg_metrics$evi2_area <- evi_area # set output value

	# calculate gup dates
	# NOTE: GetThresh could fail here?
	ogi_seg_ind <- GetThresh(ogi_thresh, evi_gup, first_greater=T, gup=T) + seg[1] - 1
	seg_metrics$ogi <- as.numeric(dates[ogi_seg_ind] - as.Date("1970-1-1"))
	midgup_seg_ind <- GetThresh(midgup_thresh, evi_gup, first_greater=T, gup=T) + seg[1] - 1
	seg_metrics$midgup <- as.numeric(dates[midgup_seg_ind] - as.Date("1970-1-1"))
	mat_seg_ind <- GetThresh(mat_thresh, evi_gup, first_greater=T, gup=T) + seg[1] - 1
	seg_metrics$mat <- as.numeric(dates[mat_seg_ind] - as.Date("1970-1-1"))

	# get peak date
	peak_seg_ind <- seg[2]
	seg_metrics$peak <- as.numeric(dates[peak_seg_ind] - as.Date("1970-1-1"))

	# get gdown dates
	sen_seg_ind <- GetThresh(sen_thresh, evi_gdown, first_greater=T, gup=F) + seg[2] - 1
	seg_metrics$sen <- as.numeric(dates[sen_seg_ind] - as.Date("1970-1-1"))
	midgdown_seg_ind <- GetThresh(midgdown_thresh, evi_gdown, first_greater=T, gup=F) + seg[2] - 1
	seg_metrics$midgdown <- as.numeric(dates[midgdown_seg_ind] - as.Date("1970-1-1"))
	dor_seg_ind <- GetThresh(dor_thresh, evi_gdown, first_greater=T, gup=F) + seg[2] - 1
	seg_metrics$dor <- as.numeric(dates[dor_seg_ind] - as.Date("1970-1-1"))

	# calculate GUP and GDOWN spline R^2
	r2_gup <- R2(evi_gup, resids_gup)
	seg_metrics$spline_r2_gup <- r2_gup
	r2_gdown <- R2(evi_gdown, resids_gdown)
	seg_metrics$spline_r2_gdown <- r2_gdown

	# calculate GUP/GDOWN length, and filled fractions
	length_gup <- as.numeric(dates[seg[2]] - dates[seg[1]])
	seg_metrics$length_gup <- length_gup
	frac_filled_gup <- sum(snowflags_gup != 0) / length_gup
	seg_metrics$frac_filled_gup <- frac_filled_gup
	length_gdown <- as.numeric(dates[seg[3]] - dates[seg[2]])
	seg_metrics$length_gdown <- length_gdown
	frac_filled_gdown <- sum(snowflags_gdown != 0 & snowflags_gdown != 2) / length_gdown
	seg_metrics$frac_filled_gdown <- frac_filled_gdown

	# calculate phenometric-vicinity filled/missing/snow fraction
	ogi_frac_filled <- sum(snowflags[max(0, (ogi_seg_ind - pheno_pars$qual_buffer_days)):min(length(snowflags), (ogi_seg_ind + pheno_pars$qual_buffer_days))] != 0) / (min(length(snowflags), (ogi_seg_ind + pheno_pars$qual_buffer_days)) - max(0, (ogi_seg_ind - pheno_pars$qual_buffer_days)) + 1)
	midgup_frac_filled <- sum(snowflags[max(0, (midgup_seg_ind - pheno_pars$qual_buffer_days)):min(length(snowflags), (midgup_seg_ind + pheno_pars$qual_buffer_days))] != 0) / (min(length(snowflags), (midgup_seg_ind + pheno_pars$qual_buffer_days)) - max(0, (midgup_seg_ind - pheno_pars$qual_buffer_days)) + 1)
	mat_frac_filled <- sum(snowflags[max(0, (mat_seg_ind - pheno_pars$qual_buffer_days)):min(length(snowflags), (mat_seg_ind + pheno_pars$qual_buffer_days))] != 0) / (min(length(snowflags), (mat_seg_ind + pheno_pars$qual_buffer_days)) - max(0, (mat_seg_ind - pheno_pars$qual_buffer_days)) + 1)
	peak_frac_filled <- sum(snowflags[max(0, (peak_seg_ind - pheno_pars$qual_buffer_days)):min(length(snowflags), (peak_seg_ind + pheno_pars$qual_buffer_days))] != 0) / (min(length(snowflags), (peak_seg_ind + pheno_pars$qual_buffer_days)) - max(0, (peak_seg_ind - pheno_pars$qual_buffer_days)) + 1)
	sen_frac_filled <- sum(snowflags[max(0, (sen_seg_ind - pheno_pars$qual_buffer_days)):min(length(snowflags), (sen_seg_ind + pheno_pars$qual_buffer_days))] != 0) / (min(length(snowflags), (sen_seg_ind + pheno_pars$qual_buffer_days)) - max(0, (sen_seg_ind - pheno_pars$qual_buffer_days)) + 1)
	midgdown_frac_filled <- sum(snowflags[max(0, (midgdown_seg_ind - pheno_pars$qual_buffer_days)):min(length(snowflags), (midgdown_seg_ind + pheno_pars$qual_buffer_days))] != 0) / (min(length(snowflags), (midgdown_seg_ind + pheno_pars$qual_buffer_days)) - max(0, (midgdown_seg_ind - pheno_pars$qual_buffer_days)) + 1)
	dor_frac_filled <- sum(snowflags[max(0, (dor_seg_ind - pheno_pars$qual_buffer_days)):min(length(snowflags), (dor_seg_ind + pheno_pars$qual_buffer_days))] != 0) / (min(length(snowflags), (dor_seg_ind + pheno_pars$qual_buffer_days)) - max(0, (dor_seg_ind - pheno_pars$qual_buffer_days)) + 1)

	# calculate phenometric-vicinity R2
	ogi_r2 <- R2(evi[max(0, (ogi_seg_ind - pheno_pars$qual_buffer_days)):min(length(evi), (ogi_seg_ind + pheno_pars$qual_buffer_days))], resids[max(0, (ogi_seg_ind - pheno_pars$qual_buffer_days)):min(length(evi), (ogi_seg_ind + pheno_pars$qual_buffer_days))])
	midgup_r2 <- R2(evi[max(0, (midgup_seg_ind - pheno_pars$qual_buffer_days)):min(length(evi), (midgup_seg_ind + pheno_pars$qual_buffer_days))], resids[max(0, (midgup_seg_ind - pheno_pars$qual_buffer_days)):min(length(evi), (midgup_seg_ind + pheno_pars$qual_buffer_days))])
	mat_r2 <- R2(evi[max(0, (mat_seg_ind - pheno_pars$qual_buffer_days)):min(length(evi), (mat_seg_ind + pheno_pars$qual_buffer_days))], resids[max(0, (mat_seg_ind - pheno_pars$qual_buffer_days)):min(length(evi), (mat_seg_ind + pheno_pars$qual_buffer_days))])
	peak_r2 <- R2(evi[max(0, (peak_seg_ind - pheno_pars$qual_buffer_days)):min(length(evi), (peak_seg_ind + pheno_pars$qual_buffer_days))], resids[max(0, (peak_seg_ind - pheno_pars$qual_buffer_days)):min(length(evi), (peak_seg_ind + pheno_pars$qual_buffer_days))])
	sen_r2 <- R2(evi[max(0, (sen_seg_ind - pheno_pars$qual_buffer_days)):min(length(evi), (sen_seg_ind + pheno_pars$qual_buffer_days))], resids[max(0, (sen_seg_ind - pheno_pars$qual_buffer_days)):min(length(evi), (sen_seg_ind + pheno_pars$qual_buffer_days))])
	midgdown_r2 <- R2(evi[max(0, (midgdown_seg_ind - pheno_pars$qual_buffer_days)):min(length(evi), (midgdown_seg_ind + pheno_pars$qual_buffer_days))], resids[max(0, (midgdown_seg_ind - pheno_pars$qual_buffer_days)):min(length(evi), (midgdown_seg_ind + pheno_pars$qual_buffer_days))])
	dor_r2 <- R2(evi[max(0, (dor_seg_ind - pheno_pars$qual_buffer_days)):min(length(evi), (dor_seg_ind + pheno_pars$qual_buffer_days))], resids[max(0, (dor_seg_ind - pheno_pars$qual_buffer_days)):min(length(evi), (dor_seg_ind + pheno_pars$qual_buffer_days))])

	# calculate phenometric quality scores
	# NOTE: what to do if the R2 is NA?, currently the qual is NA, perhaps it should be 0?
	ogi_qual <- ((pheno_pars$qual_fill_weight * (1 - ogi_frac_filled)) + (pheno_pars$qual_r2_weight * ogi_r2)) / (pheno_pars$qual_fill_weight + pheno_pars$qual_r2_weight)
	seg_metrics$ogi_qual <- ogi_qual
	midgup_qual <- ((pheno_pars$qual_fill_weight * (1 - midgup_frac_filled)) + (pheno_pars$qual_r2_weight * midgup_r2)) / (pheno_pars$qual_fill_weight + pheno_pars$qual_r2_weight)
	seg_metrics$midgup_qual <- midgup_qual
	mat_qual <- ((pheno_pars$qual_fill_weight * (1 - mat_frac_filled)) + (pheno_pars$qual_r2_weight * mat_r2)) / (pheno_pars$qual_fill_weight + pheno_pars$qual_r2_weight)
	seg_metrics$mat_qual <- mat_qual
	peak_qual <- ((pheno_pars$qual_fill_weight * (1 - peak_frac_filled)) + (pheno_pars$qual_r2_weight * peak_r2)) / (pheno_pars$qual_fill_weight + pheno_pars$qual_r2_weight)
	seg_metrics$peak_qual <- peak_qual
	sen_qual <- ((pheno_pars$qual_fill_weight * (1 - sen_frac_filled)) + (pheno_pars$qual_r2_weight * sen_r2)) / (pheno_pars$qual_fill_weight + pheno_pars$qual_r2_weight)
	seg_metrics$sen_qual <- sen_qual
	midgdown_qual <- ((pheno_pars$qual_fill_weight * (1 - midgdown_frac_filled)) + (pheno_pars$qual_r2_weight * midgdown_r2)) / (pheno_pars$qual_fill_weight + pheno_pars$qual_r2_weight)
	seg_metrics$midgdown_qual <- midgdown_qual
	dor_qual <- ((pheno_pars$qual_fill_weight * (1 - dor_frac_filled)) + (pheno_pars$qual_r2_weight * dor_r2)) / (pheno_pars$qual_fill_weight + pheno_pars$qual_r2_weight)
	seg_metrics$dor_qual <- dor_qual

	return(seg_metrics)
}

#------------------------------------------------------------------
# Jamie Olson's local min/max method from here: http://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
GetLocalMaxMin <- function(x, partial=TRUE, minima=FALSE){
	if(minima){
		if(partial){
			which(diff(c(FALSE,diff(x)>0,TRUE))>0)
		}else{
			which(diff(diff(x)>0)>0)+1
		}
	}else{
		if(partial){
			which(diff(c(TRUE,diff(x)>=0,FALSE))<0)
		}else{
			which(diff(diff(x)>=0)<0)+1
		}
	}
}

#------------------------------------------------------------------
# function to get gcc amplitude
# GetGCCAmplitude <- function(x, dates, out_sigma=3, spline_spar=0.3, out_iterations=1, min_peak_to_peak_distance=180, min_seg_amplitude=0.15, max_seg_length=180){
GetGCCAmplitude <- function(x, dates, pheno_pars){

	pred_dates <- min(dates, na.rm=T):max(dates, na.rm=T)

	# smooth and eliminate outliers
	smooth_evi2 <- try(SplineAndOutlierRemoval(x, dates, out_quant=pheno_pars$out_quant, spline_spar=pheno_pars$spline_spar, out_iterations=pheno_pars$out_iterations, pred_dates=pred_dates), silent=T)
	if(inherits(smooth_evi2, 'try-error'))
		return(NA)

	# find valid peaks in the time series
	valid_peaks <- try(FilterPeaks(smooth_evi2, FindPotentialPeaks_minmax(smooth_evi2), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance), silent=T)
	if(inherits(valid_peaks, 'try-error'))
		return(NA)

	# find full segments
	full_segs <- try(GetFullSegs(smooth_evi2, valid_peaks, max_seg_length=pheno_pars$max_seg_length, min_seg_amplitude=pheno_pars$min_seg_amplitude), silent=T)
	if(inherits(full_segs, 'try-error')){
		return(NA)
	}

	gcc_amp <- list()
	for(seg in full_segs){
		if(!is.na(seg[1]) & !is.na(seg[2])){
			gcc_amp <- c(gcc_amp, smooth_evi2[seg[2]] - smooth_evi2[seg[1]])
		}else{
			gcc_amp <- c(gcc_amp, NA)
		}
	}
	return(gcc_amp)
}

#---------------------------------------------------------------------
# Returns indices in x of potential peak values using the GetLocalMaxMin function
FindPotentialPeaks_minmax <- function(x){
	# check for all NA
	if(all(is.na(x))) return(NA)

	# check for constant x
	if((min(x, na.rm=T) - max(x, na.rm=T)) == 0) return(NA)

	potential_peaks <- try(GetLocalMaxMin(x), silent=T)
	if(inherits(potential_peaks, 'try-error')){
		return(NA)
	}else{
		# return the list of potential peaks
		return(potential_peaks)
	}
}

#---------------------------------------------------------------------
# screens peaks to make sure that they are separated by at least min_peak_to_peak_distance
# if they're not, the lesser of the two peaks is eliminated.
# valid peaks must also be greater than the min_peak_quantile of x
FilterPeaks <- function(x, potential_peaks, min_peak_to_peak_distance, min_peak_quantile=0.2){
	# check that potential_peaks is not NA
	if(all(is.na(potential_peaks))) return(NA)

	# first eliminate all peaks that have an absolute magnitude below min_peak_quantile
	# potential_peaks <- potential_peaks[!(x[potential_peaks] < quantile(x[potential_peaks], min_peak_quantile))] # among the peaks
	potential_peaks <- potential_peaks[!(x[potential_peaks] < quantile(x, min_peak_quantile))] # among the EVI2 time series

	# we loop through all the potential peaks, only moving on when we have satisfied
	# the minimum distance between peaks requirement
	i <- 1

	while(i < length(potential_peaks)){
		# find the distance to next peak
		peak_to_peak_dist <- potential_peaks[i + 1] - potential_peaks[i]

		if(peak_to_peak_dist < min_peak_to_peak_distance){
			# eliminate the smaller of the two potential peaks
			if(x[potential_peaks[i + 1]] >= x[potential_peaks[i]]){
				potential_peaks <- potential_peaks[-1 * i]
			}else{
				potential_peaks <- potential_peaks[-1 * (i + 1)]
			}
		}else{
			# distance is fine, move on
			i <- i + 1
		}
		# if we've eliminated the last peak, return NA
		if(length(potential_peaks) == 0) return(NA)
	}
	return(potential_peaks)
}

#==============================================================================
#==============================================================================
# high-level functions
#==============================================================================
#==============================================================================

#---------------------------------------------------------------------

 DoPhenology <- function(x, dates, pheno_pars){
	# set up pred_dates as a daily vector from min/max of dates
	pred_dates <- min(dates, na.rm=T):max(dates, na.rm=T)

	# get a time series of snow-filtered EVI2
	# COMMENT OUT FOR PHENOCAM
	# filt_evi2 <- x
	# filt_return <- try(SnowFilterEVI2(x, evi2_snow_quant=pheno_pars$evi2_snow_quant, ndsi_thresh=pheno_pars$ndsi_thresh), silent=T)
	# if(inherits(filt_return, 'try-error'))
	# 	return(NA)

	filt_return <- try(SnowFilterEVI2_mv(x, dates, evi2_snow_quant=pheno_pars$evi2_snow_quant, ndsi_thresh=pheno_pars$ndsi_thresh, max_snow_fill_ratio=pheno_pars$max_snow_fill_ratio), silent=T)
	if(inherits(filt_return, 'try-error'))
		return(NA)

	filt_evi2 <- filt_return[[1]] # The EVI2 time series
	snow_fills <- filt_return[[2]] # Boolean vector indicating snow fills

	# smooth and eliminate outliers

	# USE THIS ONE FOR PHENOCAM
	# smooth_evi2 <- try(SplineAndOutlierRemoval(x, dates, out_quant=pheno_pars$out_quant, spline_spar=pheno_pars$spline_spar, out_iterations=pheno_pars$out_iterations, pred_dates=pred_dates), silent=T)

	smooth_evi2 <- try(SplineAndOutlierRemoval(filt_evi2, dates, out_quant=pheno_pars$out_quant, spline_spar=pheno_pars$spline_spar, out_iterations=pheno_pars$out_iterations, pred_dates=pred_dates), silent=T)
	if(inherits(smooth_evi2, 'try-error'))
		return(NA)

	# find valid peaks in the time series
	valid_peaks <- try(FilterPeaks(smooth_evi2, FindPotentialPeaks_minmax(smooth_evi2), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance), silent=T)
	if(inherits(valid_peaks, 'try-error'))
		return(NA)

	# find full segments
	full_segs <- try(GetFullSegs(smooth_evi2, valid_peaks, max_seg_length=pheno_pars$max_seg_length, min_seg_amplitude=pheno_pars$min_seg_amplitude, agg_amp_frac=pheno_pars$agg_amp_frac), silent=T)
	if(inherits(full_segs, 'try-error')){
		return(NA)
	}

	# get PhenoDates
	pheno_dates <- try(GetPhenoDates(full_segs, smooth_evi2, pred_dates, pheno_pars$gup_threshes, pheno_pars$gdown_threshes), silent=T)
	if(inherits(pheno_dates, "try-error")){
		return(NA)
	}

	# get the segment metrics
	if(!all(is.na(unlist(pheno_dates, use.names=F)))){
		# seg_metrics <- lapply(full_segs, GetSegMetrics, smooth_evi2, x, pred_dates, dates)
		# seg_metrics <- lapply(full_segs, GetSegMetrics, smooth_evi2, filt_evi2, pred_dates, dates)
		seg_metrics <- lapply(full_segs, GetSegMetrics, smooth_evi2, filt_evi2, pred_dates, dates, snow_fills)
		seg_min <- list(unlist(seg_metrics, use.names=F)[seq(1, length(unlist(seg_metrics, use.names=F)), by=9)])
		seg_max <- list(unlist(seg_metrics, use.names=F)[seq(2, length(unlist(seg_metrics, use.names=F)), by=9)])
		seg_int <- list(unlist(seg_metrics, use.names=F)[seq(3, length(unlist(seg_metrics, use.names=F)), by=9)])
		gup_rsq <- list(unlist(seg_metrics, use.names=F)[seq(4, length(unlist(seg_metrics, use.names=F)), by=9)])
		gup_missing <- list(unlist(seg_metrics, use.names=F)[seq(5, length(unlist(seg_metrics, use.names=F)), by=9)])
		gup_snowfills <- list(unlist(seg_metrics, use.names=F)[seq(6, length(unlist(seg_metrics, use.names=F)), by=9)])
		gdown_rsq <- list(unlist(seg_metrics, use.names=F)[seq(7, length(unlist(seg_metrics, use.names=F)), by=9)])
		gdown_missing <- list(unlist(seg_metrics, use.names=F)[seq(8, length(unlist(seg_metrics, use.names=F)), by=9)])
		gdown_snowfills <- list(unlist(seg_metrics, use.names=F)[seq(9, length(unlist(seg_metrics, use.names=F)), by=9)])

		# append this new information to the output
		# pheno_dates <- c(unlist(pheno_dates, rec=F), seg_min, seg_max, seg_int, gup_rsq, gup_missing, gdown_rsq, gdown_missing)
		pheno_dates <- c(unlist(pheno_dates, rec=F, use.names=F), seg_min, seg_max, seg_int, gup_rsq, gup_missing, gup_snowfills, gdown_rsq, gdown_missing, gdown_snowfills)
	}

	return(pheno_dates)
}

#---------------------------------------------------------------------
DoPhenologyVIIRS <- function(x, dates, pheno_pars){
  # set up pred_dates as a daily vector from min/max of dates
  pred_dates <- min(dates, na.rm=T):max(dates, na.rm=T)

  # get a time series of snow-filtered EVI2
  # COMMENT OUT FOR PHENOCAM
  # filt_evi2 <- x
  # filt_return <- try(SnowFilterEVI2(x, evi2_snow_quant=pheno_pars$evi2_snow_quant, ndsi_thresh=pheno_pars$ndsi_thresh), silent=T)
  # if(inherits(filt_return, 'try-error'))
  # 	return(NA)

  # filt_return <- try(SnowFilterEVI2_mv(x, dates, evi2_snow_quant=pheno_pars$evi2_snow_quant, ndsi_thresh=pheno_pars$ndsi_thresh, max_snow_fill_ratio=pheno_pars$max_snow_fill_ratio), silent=T)
  # if(inherits(filt_return, 'try-error'))
  # return(NA)

  filt_evi2 <- x[1:(length(x) / 2)] # The EVI2 time series
  qa <- x[((length(x) / 2) + 1):length(x)] # Vector of QA values
  snow_fills <- rep(FALSE, length(filt_evi2)) # Boolean vector indicating snow fills
  snow_fills[qa == 1] <- TRUE

  # smooth and eliminate outliers

  # USE THIS ONE FOR PHENOCAM
  # smooth_evi2 <- try(SplineAndOutlierRemoval(x, dates, out_quant=pheno_pars$out_quant, spline_spar=pheno_pars$spline_spar, out_iterations=pheno_pars$out_iterations, pred_dates=pred_dates), silent=T)

  # NOTE: there's no outlier identification or removal in SmoothSnowFill!
  smooth_evi2 <- try(SmoothSnowFill(filt_evi2, snow_fills, dates, pred_dates, pheno_pars), silent=T)
  if(inherits(smooth_evi2, 'try-error'))
  return(NA)

  # smooth_evi2 <- try(SplineAndOutlierRemoval(filt_evi2, dates, out_quant=pheno_pars$out_quant, spline_spar=pheno_pars$spline_spar, out_iterations=pheno_pars$out_iterations, pred_dates=pred_dates), silent=T)
  # if(inherits(smooth_evi2, 'try-error'))
  # return(NA)

  # find valid peaks in the time series
  valid_peaks <- try(FilterPeaks(smooth_evi2, FindPotentialPeaks_minmax(smooth_evi2), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance), silent=T)
  if(inherits(valid_peaks, 'try-error'))
  return(NA)

  # find full segments
  full_segs <- try(GetFullSegs(smooth_evi2, valid_peaks, max_seg_length=pheno_pars$max_seg_length, min_seg_amplitude=pheno_pars$min_seg_amplitude, agg_amp_frac=pheno_pars$agg_amp_frac), silent=T)
  if(inherits(full_segs, 'try-error')){
  return(NA)
  }

  # get PhenoDates
  pheno_dates <- try(GetPhenoDates(full_segs, smooth_evi2, pred_dates, pheno_pars$gup_threshes, pheno_pars$gdown_threshes), silent=T)
  if(inherits(pheno_dates, "try-error")){
  return(NA)
  }

  # get the segment metrics
  if(!all(is.na(unlist(pheno_dates, use.names=F)))){
  # seg_metrics <- lapply(full_segs, GetSegMetrics, smooth_evi2, x, pred_dates, dates)
  # seg_metrics <- lapply(full_segs, GetSegMetrics, smooth_evi2, filt_evi2, pred_dates, dates)
  seg_metrics <- lapply(full_segs, GetSegMetrics, smooth_evi2, filt_evi2, pred_dates, dates, snow_fills)
  seg_min <- list(unlist(seg_metrics, use.names=F)[seq(1, length(unlist(seg_metrics, use.names=F)), by=9)])
  seg_max <- list(unlist(seg_metrics, use.names=F)[seq(2, length(unlist(seg_metrics, use.names=F)), by=9)])
  seg_int <- list(unlist(seg_metrics, use.names=F)[seq(3, length(unlist(seg_metrics, use.names=F)), by=9)])
  gup_rsq <- list(unlist(seg_metrics, use.names=F)[seq(4, length(unlist(seg_metrics, use.names=F)), by=9)])
  gup_missing <- list(unlist(seg_metrics, use.names=F)[seq(5, length(unlist(seg_metrics, use.names=F)), by=9)])
  gup_snowfills <- list(unlist(seg_metrics, use.names=F)[seq(6, length(unlist(seg_metrics, use.names=F)), by=9)])
  gdown_rsq <- list(unlist(seg_metrics, use.names=F)[seq(7, length(unlist(seg_metrics, use.names=F)), by=9)])
  gdown_missing <- list(unlist(seg_metrics, use.names=F)[seq(8, length(unlist(seg_metrics, use.names=F)), by=9)])
  gdown_snowfills <- list(unlist(seg_metrics, use.names=F)[seq(9, length(unlist(seg_metrics, use.names=F)), by=9)])

  # append this new information to the output
  # pheno_dates <- c(unlist(pheno_dates, rec=F), seg_min, seg_max, seg_int, gup_rsq, gup_missing, gdown_rsq, gdown_missing)
  pheno_dates <- c(unlist(pheno_dates, rec=F, use.names=F), seg_min, seg_max, seg_int, gup_rsq, gup_missing, gup_snowfills, gdown_rsq, gdown_missing, gdown_snowfills)
  }

  return(pheno_dates)
}

#---------------------------------------------------------------------
DoPhenologyC6 <- function(x, dates, pheno_pars){
 # set up pred_dates as a daily vector from min/max of dates
 pred_dates <- min(dates, na.rm=T):max(dates, na.rm=T)

 # get a time series of snow-filtered EVI2
 # COMMENT OUT FOR PHENOCAM
 # filt_evi2 <- x
 # filt_return <- try(SnowFilterEVI2(x, evi2_snow_quant=pheno_pars$evi2_snow_quant, ndsi_thresh=pheno_pars$ndsi_thresh), silent=T)
 # if(inherits(filt_return, 'try-error'))
 # 	return(NA)

 # filt_return <- try(SnowFilterEVI2_mv(x, dates, evi2_snow_quant=pheno_pars$evi2_snow_quant, ndsi_thresh=pheno_pars$ndsi_thresh, max_snow_fill_ratio=pheno_pars$max_snow_fill_ratio), silent=T)
 # if(inherits(filt_return, 'try-error'))
 #  return(NA)
 #
 # filt_evi2 <- filt_return[[1]] # The EVI2 time series
 # snow_fills <- filt_return[[2]] # Boolean vector indicating snow fills
 #
 # # smooth and eliminate outliers
 #
 # # USE THIS ONE FOR PHENOCAM
 # # smooth_evi2 <- try(SplineAndOutlierRemoval(x, dates, out_quant=pheno_pars$out_quant, spline_spar=pheno_pars$spline_spar, out_iterations=pheno_pars$out_iterations, pred_dates=pred_dates), silent=T)
 #
 # smooth_evi2 <- try(SplineAndOutlierRemoval(filt_evi2, dates, out_quant=pheno_pars$out_quant, spline_spar=pheno_pars$spline_spar, out_iterations=pheno_pars$out_iterations, pred_dates=pred_dates), silent=T)
 # if(inherits(smooth_evi2, 'try-error'))
 #  return(NA)

 # find valid peaks in the time series
 valid_peaks <- try(FilterPeaks(x, FindPotentialPeaks_minmax(x), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance), silent=T)
 if(inherits(valid_peaks, 'try-error'))
	 return(NA)

 # find full segments
 full_segs <- try(GetFullSegs(x, valid_peaks, max_seg_length=pheno_pars$max_seg_length, min_seg_amplitude=pheno_pars$min_seg_amplitude, agg_amp_frac=pheno_pars$agg_amp_frac), silent=T)
 if(inherits(full_segs, 'try-error')){
	 return(NA)
 }

 # get PhenoDates
 pheno_dates <- try(GetPhenoDates(full_segs, x, pred_dates, pheno_pars$gup_threshes, pheno_pars$gdown_threshes), silent=T)
 if(inherits(pheno_dates, "try-error")){
	 return(NA)
 }

 # get the segment metrics
 if(!all(is.na(unlist(pheno_dates, use.names=F)))){
	 # seg_metrics <- lapply(full_segs, GetSegMetrics, smooth_evi2, x, pred_dates, dates)
	 # seg_metrics <- lapply(full_segs, GetSegMetrics, smooth_evi2, filt_evi2, pred_dates, dates)

	#  seg_metrics <- lapply(full_segs, GetSegMetrics, smooth_evi2, filt_evi2, pred_dates, dates, snow_fills)
	#  seg_min <- list(unlist(seg_metrics, use.names=F)[seq(1, length(unlist(seg_metrics, use.names=F)), by=9)])
	#  seg_max <- list(unlist(seg_metrics, use.names=F)[seq(2, length(unlist(seg_metrics, use.names=F)), by=9)])
	#  seg_int <- list(unlist(seg_metrics, use.names=F)[seq(3, length(unlist(seg_metrics, use.names=F)), by=9)])
	#  gup_rsq <- list(unlist(seg_metrics, use.names=F)[seq(4, length(unlist(seg_metrics, use.names=F)), by=9)])
	#  gup_missing <- list(unlist(seg_metrics, use.names=F)[seq(5, length(unlist(seg_metrics, use.names=F)), by=9)])
	#  gup_snowfills <- list(unlist(seg_metrics, use.names=F)[seq(6, length(unlist(seg_metrics, use.names=F)), by=9)])
	#  gdown_rsq <- list(unlist(seg_metrics, use.names=F)[seq(7, length(unlist(seg_metrics, use.names=F)), by=9)])
	#  gdown_missing <- list(unlist(seg_metrics, use.names=F)[seq(8, length(unlist(seg_metrics, use.names=F)), by=9)])
	#  gdown_snowfills <- list(unlist(seg_metrics, use.names=F)[seq(9, length(unlist(seg_metrics, use.names=F)), by=9)])

	 # append this new information to the output
	 # pheno_dates <- c(unlist(pheno_dates, rec=F), seg_min, seg_max, seg_int, gup_rsq, gup_missing, gdown_rsq, gdown_missing)
	#  pheno_dates <- c(unlist(pheno_dates, rec=F, use.names=F), seg_min, seg_max, seg_int, gup_rsq, gup_missing, gup_snowfills, gdown_rsq, gdown_missing, gdown_snowfills)
	pheno_dates <- c(unlist(pheno_dates, rec=F, use.names=F))
 }

 return(pheno_dates)
}

#---------------------------------------------------------------------
WritePhenologyData <- function(out_file, pheno_values_bip, tile, nbands=44){
	# write the data to a file as binary integers (4 bytes!)
	if(!file.exists(out_file)){
		ff <- file(out_file, 'wb') # create the file and open for writing
	}else{
		ff <- file(out_file, 'ab') # file exists, append to the end
	}

	# where_to_start <- (start_line - 1) * samples * natural_size * nbands
	# seek(ff, where_to_start) # position the file connection
	writeBin(pheno_values_bip, ff)
	close(ff)

	# this code courtesy of Damien Sulla Menashe
	out_hdr <- paste(out_file, ".hdr", sep="")
	if(!file.exists(out_hdr)){
		# variables specifying the upper left corner of the modis grid and the tile/pixel lengths
		uly_map = 10007554.677
		ulx_map = -20015109.354
		lry_map = -10007554.677
		lrx_map = 20015109.354
		pix = 463.312716525
		dims = 2400

		# calculate the upper left corner of the current tile
		tile_h <- as.integer(substr(tile, 2, 3))
		tile_v <- as.integer(substr(tile, 5, 6))
		ulx = ulx_map + (tile_h * pix * dims)
		uly = uly_map - (tile_v * pix * dims)

		# create the header text
		temp_txt = paste("ENVI description = { MCD12Q2 Collection 6 }\nlines = ", dims, "\nsamples = ", dims, "\nbands = ", nbands, "\nheader offset = 0\nfile type = ENVI Standard\ndata type = 3\ninterleave = bip\nbyte order = 0\nmap info = {Sinusoidal, 1, 1,", ulx, ", ", uly, ", ", pix, ", ", pix, "}", "\ncoordinate system string = {PROJCS[\"Sinusoidal\",GEOGCS[\"GCS_unnamed ellipse\",DATUM[\"D_unknown\",SPHEROID[\"Unknown\",6371007.181,0]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Sinusoidal\"],PARAMETER[\"central_meridian\",0],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"Meter\",1]]}", sep="")

		# write the header to a file
		sink(out_hdr)
		cat(temp_txt)
		sink()
	}
}

#---------------------------------------------------------------------
AnnualPhenologyC6 <- function(x, dates, pheno_pars, plot=F){
	# get a template return value
	annual_pheno_metrics <- PhenoReturnValue()

	# split the data up: first third are smoothed nbar-evi2, second third are residuals, and last third are snowflags
	# also scale, and fill for NA
	x[x==pheno_pars$nbar_NA_value] <- NA
	data_length <- length(x) / 3
	evi <- x[1:data_length]
	evi <- evi / pheno_pars$nbar_scale_factor
	resids <- x[(data_length + 1):(2 * data_length)]
	resids <- resids / pheno_pars$nbar_scale_factor
	snowflags <- x[(2 * data_length + 1):(3 * data_length)]

	# find valid peaks in the time series
	valid_peaks <- try(FilterPeaks(evi, FindPotentialPeaks_minmax(evi), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance), silent=T)
	if(inherits(valid_peaks, 'try-error') | is.na(valid_peaks)){
		annual_pheno_metrics$fill_code <- 1 # no valid peaks
		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
		return(annual_pheno_metrics)
	}

	# find full segments
	full_segs <- try(GetFullSegs(evi, valid_peaks, max_seg_length=pheno_pars$max_seg_length, min_seg_amplitude=pheno_pars$min_seg_amplitude, agg_amp_frac=pheno_pars$agg_amp_frac), silent=T)
	if(inherits(full_segs, 'try-error') | is.na(full_segs)){
		annual_pheno_metrics$fill_code <- 2 # no valid segments
		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
		return(annual_pheno_metrics)
	}

	# eliminate segs that end before, or start after the period of interest
	seg_overlaps <- unlist(lapply(full_segs, SegOverlapsPeriod, dates, pheno_pars$pheno_period_start, pheno_pars$pheno_period_end))
	if(all(!seg_overlaps)){
		# no segs within period of interest
		annual_pheno_metrics$fill_code <- 2 # no valid segments
		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
		return(annual_pheno_metrics)
	}
	full_segs <- full_segs[seg_overlaps]

	# get the segment metrics
	seg_metrics <- lapply(full_segs, SegMet, x=x, dates=dates, pheno_pars=pheno_pars)

	# limit to segments where the peak is in the period of interest
	is_in_period <- unlist(lapply(seg_metrics, PeakInPeriod, start_date=pheno_pars$pheno_period_start, end_date=pheno_pars$pheno_period_end))
	segs_in_period <- seg_metrics[is_in_period]

	# assign the seg metric values to the output return value
	num_cycles <- length(segs_in_period) # total number of valid peaks within the period of interest
	if(num_cycles == 0){
		# no segments in period of interest, return default annual pheno metrics
		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
		return(annual_pheno_metrics)
	}else{
		# set first cycle return metrics
		cycle1_seg_met <- segs_in_period[[1]]
		annual_pheno_metrics <- SetReturnValues(annual_pheno_metrics, cycle1_seg_met, cycle=1, num_cycles=num_cycles, fill_code=0)
		if(num_cycles > 1){
			# set the second cycle
			cycle2_seg_met <- segs_in_period[[2]]
			annual_pheno_metrics <- SetReturnValues(annual_pheno_metrics, cycle2_seg_met, cycle=2)
		}
		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
		return(annual_pheno_metrics)
	}
}

#---------------------------------------------------------------------
ScaleToIntegerAndSetNA <- function(annual_pheno_metrics, scale=1e4, outNA=32767){
	annual_pheno_metrics$evi_area_cycle1=round(annual_pheno_metrics$evi_area_cycle1 / 1e3 * scale) # so we don't exceed 32767!
	annual_pheno_metrics$evi_amp_cycle1=round(annual_pheno_metrics$evi_amp_cycle1 * scale)
	annual_pheno_metrics$evi_min_cycle1=round(annual_pheno_metrics$evi_min_cycle1 * scale)
	annual_pheno_metrics$frac_filled_gup_cycle1=round(annual_pheno_metrics$frac_filled_gup_cycle1 * scale)
	annual_pheno_metrics$frac_filled_gdown_cycle1=round(annual_pheno_metrics$frac_filled_gdown_cycle1 * scale)
	annual_pheno_metrics$ogi_qual_cycle1=round(annual_pheno_metrics$ogi_qual_cycle1 * scale)
	annual_pheno_metrics$midgup_qual_cycle1=round(annual_pheno_metrics$midgup_qual_cycle1 * scale)
	annual_pheno_metrics$mat_qual_cycle1=round(annual_pheno_metrics$mat_qual_cycle1 * scale)
	annual_pheno_metrics$peak_qual_cycle1=round(annual_pheno_metrics$peak_qual_cycle1 * scale)
	annual_pheno_metrics$sen_qual_cycle1=round(annual_pheno_metrics$sen_qual_cycle1 * scale)
	annual_pheno_metrics$midgdown_qual_cycle1=round(annual_pheno_metrics$midgdown_qual_cycle1 * scale)
	annual_pheno_metrics$dor_qual_cycle1=round(annual_pheno_metrics$dor_qual_cycle1 * scale)

	annual_pheno_metrics$evi_area_cycle2=round(annual_pheno_metrics$evi_area_cycle2 / 1e3 * scale) # so we don't exceed 32767!
	annual_pheno_metrics$evi_amp_cycle2=round(annual_pheno_metrics$evi_amp_cycle2 * scale)
	annual_pheno_metrics$evi_min_cycle2=round(annual_pheno_metrics$evi_min_cycle2 * scale)
	annual_pheno_metrics$frac_filled_gup_cycle2=round(annual_pheno_metrics$frac_filled_gup_cycle2 * scale)
	annual_pheno_metrics$frac_filled_gdown_cycle2=round(annual_pheno_metrics$frac_filled_gdown_cycle2 * scale)
	annual_pheno_metrics$ogi_qual_cycle2=round(annual_pheno_metrics$ogi_qual_cycle2 * scale)
	annual_pheno_metrics$midgup_qual_cycle2=round(annual_pheno_metrics$midgup_qual_cycle2 * scale)
	annual_pheno_metrics$mat_qual_cycle2=round(annual_pheno_metrics$mat_qual_cycle2 * scale)
	annual_pheno_metrics$peak_qual_cycle2=round(annual_pheno_metrics$peak_qual_cycle2 * scale)
	annual_pheno_metrics$sen_qual_cycle2=round(annual_pheno_metrics$sen_qual_cycle2 * scale)
	annual_pheno_metrics$midgdown_qual_cycle2=round(annual_pheno_metrics$midgdown_qual_cycle2 * scale)
	annual_pheno_metrics$dor_qual_cycle2=round(annual_pheno_metrics$dor_qual_cycle2 * scale)
	annual_pheno_metrics[is.na(annual_pheno_metrics)] <- outNA
	return(as.integer(annual_pheno_metrics))
}

#---------------------------------------------------------------------
SetReturnValues <- function(annual_pheno_metrics, seg_met, cycle=1, num_cycles=NA, fill_code=NA){
	if(!is.na(num_cycles)){
		annual_pheno_metrics$num_cycles <- num_cycles
	}

	if(!is.na(fill_code)){
		annual_pheno_metrics$fill_code <- fill_code
	}

	if(cycle==1){
		# cycle 1 metrics
		annual_pheno_metrics$evi_area_cycle1=seg_met$evi2_area
		annual_pheno_metrics$evi_amp_cycle1=seg_met$evi2_amp
		annual_pheno_metrics$evi_min_cycle1=seg_met$evi2_min
		annual_pheno_metrics$frac_filled_gup_cycle1=seg_met$frac_filled_gup
		annual_pheno_metrics$frac_filled_gdown_cycle1=seg_met$frac_filled_gdown
		annual_pheno_metrics$length_gup_cycle1=seg_met$length_gup
		annual_pheno_metrics$length_gdown_cycle1=seg_met$length_gdown
		annual_pheno_metrics$ogi_cycle1=seg_met$ogi
		annual_pheno_metrics$midgup_cycle1=seg_met$midgup
		annual_pheno_metrics$mat_cycle1=seg_met$mat
		annual_pheno_metrics$peak_cycle1=seg_met$peak
		annual_pheno_metrics$sen_cycle1=seg_met$sen
		annual_pheno_metrics$midgdown_cycle1=seg_met$midgdown
		annual_pheno_metrics$dor_cycle1=seg_met$dor
		annual_pheno_metrics$ogi_qual_cycle1=seg_met$ogi_qual
		annual_pheno_metrics$midgup_qual_cycle1=seg_met$midgup_qual
		annual_pheno_metrics$mat_qual_cycle1=seg_met$mat_qual
		annual_pheno_metrics$peak_qual_cycle1=seg_met$peak_qual
		annual_pheno_metrics$sen_qual_cycle1=seg_met$sen_qual
		annual_pheno_metrics$midgdown_qual_cycle1=seg_met$midgdown_qual
		annual_pheno_metrics$dor_qual_cycle1=seg_met$dor_qual
	}else{
		# cycle 2 metrics
		annual_pheno_metrics$evi_area_cycle2=seg_met$evi2_area
		annual_pheno_metrics$evi_amp_cycle2=seg_met$evi2_amp
		annual_pheno_metrics$evi_min_cycle2=seg_met$evi2_min
		annual_pheno_metrics$frac_filled_gup_cycle2=seg_met$frac_filled_gup
		annual_pheno_metrics$frac_filled_gdown_cycle2=seg_met$frac_filled_gdown
		annual_pheno_metrics$length_gup_cycle2=seg_met$length_gup
		annual_pheno_metrics$length_gdown_cycle2=seg_met$length_gdown
		annual_pheno_metrics$ogi_cycle2=seg_met$ogi
		annual_pheno_metrics$midgup_cycle2=seg_met$midgup
		annual_pheno_metrics$mat_cycle2=seg_met$mat
		annual_pheno_metrics$peak_cycle2=seg_met$peak
		annual_pheno_metrics$sen_cycle2=seg_met$sen
		annual_pheno_metrics$midgdown_cycle2=seg_met$midgdown
		annual_pheno_metrics$dor_cycle2=seg_met$dor
		annual_pheno_metrics$ogi_qual_cycle2=seg_met$ogi_qual
		annual_pheno_metrics$midgup_qual_cycle2=seg_met$midgup_qual
		annual_pheno_metrics$mat_qual_cycle2=seg_met$mat_qual
		annual_pheno_metrics$peak_qual_cycle2=seg_met$peak_qual
		annual_pheno_metrics$sen_qual_cycle2=seg_met$sen_qual
		annual_pheno_metrics$midgdown_qual_cycle2=seg_met$midgdown_qual
		annual_pheno_metrics$dor_qual_cycle2=seg_met$dor_qual
	}

	return(annual_pheno_metrics)
}

#---------------------------------------------------------------------
DefaultPhenoParameters <- function(year_of_interest){
	pheno_pars <- list(
		min_peak_to_peak_distance=120,
		min_seg_amplitude=0.15,
		agg_amp_frac=0.15,
		max_seg_length=200,
	  ogi_thresh=0.15,
	  midgup_thresh=0.5,
	  mat_thresh=0.95, #???
	  sen_thresh=0.8, #???
	  midgdown_thresh=0.5,
	  dor_thresh=0.15,
	  nbar_scale_factor=1e4,
	  nbar_NA_value=32767,
	  qual_buffer_days=14,
	  qual_r2_weight=1,
	  qual_fill_weight=4,
		out_float_scale=1e4,
		out_NA_value=32767,
	  pheno_period_start=as.numeric(as.Date(paste(year_of_interest, "-1-1", sep=""))),
	  pheno_period_end=as.numeric(as.Date(paste(year_of_interest, "-12-31", sep="")))
	)
	return(pheno_pars)
}

#---------------------------------------------------------------------
ReadPhenoParameters <- function(parameter_file, year_of_interest){
	return(DefaultPhenoParameters(year_of_interest))
}

#---------------------------------------------------------------------
# Checks if the peak of the segment with metrics "segmet" are between "start_date" and "end_date"
PeakInPeriod <- function(segmet, start_date, end_date){
	if(segmet$peak >= start_date & segmet$peak <= end_date){
		return(TRUE)
	}else{
		return(FALSE)
	}
}

PhenoReturnValue <- function(default_value=NA){
	ret_value <- list(
		# segment metrics
		num_cycles=default_value, # 1
		fill_code=default_value, # 2

		# cycle 1 metrics
		evi_area_cycle1=default_value, # 3
		evi_amp_cycle1=default_value, # 4
		evi_min_cycle1=default_value, # 5
		frac_filled_gup_cycle1=default_value, # 6
		frac_filled_gdown_cycle1=default_value, # 7
		length_gup_cycle1=default_value, # 8
		length_gdown_cycle1=default_value, # 9
		ogi_cycle1=default_value, # 10
		midgup_cycle1=default_value, # 11
		mat_cycle1=default_value, # 12
		peak_cycle1=default_value, # 13
		sen_cycle1=default_value, # 14
		midgdown_cycle1=default_value, # 15
		dor_cycle1=default_value, # 16
		ogi_qual_cycle1=default_value, # 17
		midgup_qual_cycle1=default_value, # 18
		mat_qual_cycle1=default_value, # 19
		peak_qual_cycle1=default_value, # 20
		sen_qual_cycle1=default_value, # 21
		midgdown_qual_cycle1=default_value, # 22
		dor_qual_cycle1=default_value, # 23

		# cycle 2 metrics
		evi_area_cycle2=default_value, # 24
		evi_amp_cycle2=default_value, # 25
		evi_min_cycle2=default_value, # 26
		frac_filled_gup_cycle2=default_value, # 27
		frac_filled_gdown_cycle2=default_value, # 28
		length_gup_cycle2=default_value, # 29
		length_gdown_cycle2=default_value, # 30
		ogi_cycle2=default_value, # 31
		midgup_cycle2=default_value, # 32
		mat_cycle2=default_value, # 33
		peak_cycle2=default_value, # 34
		sen_cycle2=default_value, # 35
		midgdown_cycle2=default_value, # 36
		dor_cycle2=default_value, # 37
		ogi_qual_cycle2=default_value, # 38
		midgup_qual_cycle2=default_value, # 39
		mat_qual_cycle2=default_value, # 40
		peak_qual_cycle2=default_value, # 41
		sen_qual_cycle2=default_value, # 42
		midgdown_qual_cycle2=default_value, # 43
		dor_qual_cycle2=default_value # 44
	)
	return(ret_value)
}

#---------------------------------------------------------------------
# retrieves 3 full years of splined C6 NBAR-EVI2 data, residuals, and snow flags as a data.frame
# where each row is a time series of evi, resid, and snow: evi001, ..., evi365, resid001, ..., resid365, flag001, ..., flag365
Get3YearDataChunk <- function(evi2_files, resid_files, snow_files, year_of_interest, start_line, lines_to_read){
	# small function to retrieve the year from the nbar file names
	# NOTE: this is FRAGILE! A more robust solution would be better
	year_function <- function(x) unlist(strsplit(basename(x), split="\\."))[3]
	evi2_years <- as.integer(unlist(lapply(evi2_files, year_function)))
	resid_years <- as.integer(unlist(lapply(resid_files, year_function)))
	snow_years <- as.integer(unlist(lapply(snow_files, year_function)))
	years_of_interest <- c(year_of_interest - 1, year_of_interest, year_of_interest + 1)

	# get indices for year before, year of interest, and year after in the files
	evi2_file_inds <- match(years_of_interest, evi2_years)
	resid_file_inds <- match(years_of_interest, resid_years)
	snow_file_inds <- match(years_of_interest, snow_years)

	# check that all needed data is available
	if(any(is.na(c(evi2_file_inds, resid_file_inds, snow_file_inds)))){
		print("Some data are missing, aborting")
		return(NA)
	}

	tmp <- cbind(
		read_int_file(evi2_files[evi2_file_inds[1]], lines_to_read=lines_to_read, start_line=start_line),
		read_int_file(evi2_files[evi2_file_inds[2]], lines_to_read=lines_to_read, start_line=start_line),
		read_int_file(evi2_files[evi2_file_inds[3]], lines_to_read=lines_to_read, start_line=start_line),
		read_int_file(resid_files[resid_file_inds[1]], lines_to_read=lines_to_read, start_line=start_line),
		read_int_file(resid_files[resid_file_inds[2]], lines_to_read=lines_to_read, start_line=start_line),
		read_int_file(resid_files[resid_file_inds[3]], lines_to_read=lines_to_read, start_line=start_line),
		read_int_file(snow_files[snow_file_inds[1]], lines_to_read=lines_to_read, start_line=start_line),
		read_int_file(snow_files[snow_file_inds[2]], lines_to_read=lines_to_read, start_line=start_line),
		read_int_file(snow_files[snow_file_inds[3]], lines_to_read=lines_to_read, start_line=start_line)
	) # end data read/cbind data

	return(tmp)
}

#---------------------------------------------------------------------
# Checks if a single seg overlaps the period interest at all, if not then we don't want to bother processing
SegOverlapsPeriod <- function(seg, dates, period_start, period_end){
	if(dates[seg[3]] < period_start | dates[seg[1]] > period_end){
		return(FALSE)
	}else{
		return(TRUE)
	}
}

#---------------------------------------------------------------------
DoPhenologyLandsat <- function(x, dates, pheno_pars){
	# run the landsat filter/smoother
	filt_evi2 <- x
	snow_fills <- rep(0, length(filt_evi2))

	# LandsatSmoother <- function(x, dates, fill_quant=0.05, x_min=0, spike_thresh=2, min_resid=0.1, spline_spar=NULL, fill_doy=NULL){
	tmp <- LandsatSmoother(x, dates, fill_quant=pheno_pars$LandsatFillQuant, spike_thresh=pheno_pars$LandsatSpikeThresh, min_resid=pheno_pars$LandsatMinResid, fill_doy=pheno_pars$LandsatFillDOY, doAnnual=pheno_pars$LandsatDoAnnual, padHeadTail=pheno_pars$LandsatPadHeadTail)
	if(is.na(tmp)){
		# all values were missing
		# NOTE: probably better way to handle this condition than in LandsatSmoother...
		return(NA)
	}
	pred_dates <- tmp[[1]]
	smooth_evi2 <- tmp[[2]]

	# find valid peaks in the time series
	valid_peaks <- try(FilterPeaks(smooth_evi2, FindPotentialPeaks_minmax(smooth_evi2), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance), silent=T)
	if(inherits(valid_peaks, 'try-error'))
		return(NA)

	# find full segments
	full_segs <- try(GetFullSegs(smooth_evi2, valid_peaks, max_seg_length=pheno_pars$max_seg_length, min_seg_amplitude=pheno_pars$min_seg_amplitude, agg_amp_frac=pheno_pars$agg_amp_frac), silent=T)
	if(inherits(full_segs, 'try-error')){
		return(NA)
	}

	# get PhenoDates
	pheno_dates <- try(GetPhenoDates(full_segs, smooth_evi2, pred_dates, pheno_pars$gup_threshes, pheno_pars$gdown_threshes), silent=T)
	if(inherits(pheno_dates, "try-error")){
		return(NA)
	}

	# get the segment metrics
	if(!all(is.na(unlist(pheno_dates, use.names=F)))){
		# seg_metrics <- lapply(full_segs, GetSegMetrics, smooth_evi2, x, pred_dates, dates)
		# seg_metrics <- lapply(full_segs, GetSegMetrics, smooth_evi2, filt_evi2, pred_dates, dates)
		seg_metrics <- lapply(full_segs, GetSegMetrics, smooth_evi2, filt_evi2, pred_dates, dates, snow_fills)
		seg_min <- list(unlist(seg_metrics, use.names=F)[seq(1, length(unlist(seg_metrics, use.names=F)), by=9)])
		seg_max <- list(unlist(seg_metrics, use.names=F)[seq(2, length(unlist(seg_metrics, use.names=F)), by=9)])
		seg_int <- list(unlist(seg_metrics, use.names=F)[seq(3, length(unlist(seg_metrics, use.names=F)), by=9)])
		gup_rsq <- list(unlist(seg_metrics, use.names=F)[seq(4, length(unlist(seg_metrics, use.names=F)), by=9)])
		gup_missing <- list(unlist(seg_metrics, use.names=F)[seq(5, length(unlist(seg_metrics, use.names=F)), by=9)])
		gup_snowfills <- list(unlist(seg_metrics, use.names=F)[seq(6, length(unlist(seg_metrics, use.names=F)), by=9)])
		gdown_rsq <- list(unlist(seg_metrics, use.names=F)[seq(7, length(unlist(seg_metrics, use.names=F)), by=9)])
		gdown_missing <- list(unlist(seg_metrics, use.names=F)[seq(8, length(unlist(seg_metrics, use.names=F)), by=9)])
		gdown_snowfills <- list(unlist(seg_metrics, use.names=F)[seq(9, length(unlist(seg_metrics, use.names=F)), by=9)])

		# append this new information to the output
		# pheno_dates <- c(unlist(pheno_dates, rec=F), seg_min, seg_max, seg_int, gup_rsq, gup_missing, gdown_rsq, gdown_missing)
		pheno_dates <- c(unlist(pheno_dates, rec=F, use.names=F), seg_min, seg_max, seg_int, gup_rsq, gup_missing, gup_snowfills, gdown_rsq, gdown_missing, gdown_snowfills)
	}

	return(pheno_dates)
}


GetSegMetrics <- function(seg, x_smooth, x_raw, smooth_dates, raw_dates, snow_fills){
	if(is.na(seg)){
		return(NA)
	}
	# get the subset of the smoothed and original time series
	tmp_seg_smooth <- x_smooth[seg[1]:seg[3]]
	tmp_gup_smooth <- x_smooth[seg[1]:seg[2]]
	tmp_gdown_smooth <- x_smooth[seg[2]:seg[3]]

	# # get the greenup segment minimum/maximum SVI
	# seg_min <- min(tmp_gup_smooth, na.rm=T)
	# seg_max <- max(tmp_gup_smooth, na.rm=T)

	# get the full segment minimum/maximum SVI
	seg_min <- min(tmp_seg_smooth, na.rm=T)
	seg_max <- max(tmp_seg_smooth, na.rm=T)


	# get the segment integrated SVI: the sum of values above the segment minimum
	# seg_int <- sum(tmp_seg_smooth)
	seg_int <- sum(tmp_seg_smooth - seg_min)

	# get greenup segment spline R^2
	gup_raw_date_inds <- which(raw_dates >= smooth_dates[seg[1]] & raw_dates <= smooth_dates[seg[2]]) # indices in raw data of gup segment
	gup_smooth_date_inds <- match(raw_dates[gup_raw_date_inds], smooth_dates) # indices of raw dates in smooth dates

	gup_raw_data <- x_raw[gup_raw_date_inds] # get the raw data associated with the gup segment
	gup_smooth_data <- x_smooth[gup_smooth_date_inds] # get the smoothed values associated with each raw data value
	gup_snow_data <- snow_fills[gup_raw_date_inds] # get the snow data associated with the gup segment

	# calculate the coeff of determination for the spline fit to the greenup segment
	gup_seg_rsquared <- 1 - (sum((gup_raw_data - gup_smooth_data)^2, na.rm=T) / sum((gup_raw_data - mean(gup_raw_data, na.rm=T))^2, na.rm=T))

	# get greenup segment spline R^2
	gdown_raw_date_inds <- which(raw_dates >= smooth_dates[seg[2]] & raw_dates <= smooth_dates[seg[3]]) # indices in raw data of gdown segment
	gdown_smooth_date_inds <- match(raw_dates[gdown_raw_date_inds], smooth_dates) # indices of raw dates in smooth dates

	gdown_raw_data <- x_raw[gdown_raw_date_inds] # get the raw data associated with the gdown segment
	gdown_smooth_data <- x_smooth[gdown_smooth_date_inds] # get the smoothed values associated with each raw data value
	gdown_snow_data <- snow_fills[gdown_raw_date_inds] # get the raw data associated with the gdown segment

	# calculate the coeff of determination for the spline fit to the greenup segment
	gdown_seg_rsquared <- 1 - (sum((gdown_raw_data - gdown_smooth_data)^2, na.rm=T) / sum((gdown_raw_data - mean(gdown_raw_data, na.rm=T))^2, na.rm=T))

	# get the segment missing raw data
	gup_missing <- sum(is.na(gup_raw_data)) / length(gup_raw_data)
	gdown_missing <- sum(is.na(gdown_raw_data)) / length(gdown_raw_data)

	# get the segment snow fill information
	gup_snowfilled <- sum(gup_snow_data != 0, na.rm=T)
	gdown_snowfilled <- sum(gdown_snow_data != 0, na.rm=T)

	# return the metrics as: seg_min, seg_max, seg_int, gup_r2, gup_missing_percentage, gdown_missing_percentage
	# return(c(seg_min, seg_max, seg_int, gup_seg_rsquared, gup_missing, gdown_seg_rsquared, gdown_missing))
	return(c(seg_min, seg_max, seg_int, gup_seg_rsquared, gup_missing, gup_snowfilled, gdown_seg_rsquared, gdown_missing, gdown_snowfilled))
}

#---------------------------------------------------------------------
DoPhenology_old <- function(x, dates, pheno_pars){
	# set up pred_dates as a daily vector from min/max of dates
	pred_dates <- min(dates, na.rm=T):max(dates, na.rm=T)

	# get a time series of snow-filtered EVI2
	# COMMENT OUT FOR PHENOCAM
	filt_evi2 <- try(SnowFilterEVI2(x, evi2_snow_quant=pheno_pars$evi2_snow_quant, ndsi_thresh=pheno_pars$ndsi_thresh), silent=T)
	if(inherits(filt_evi2, 'try-error'))
		return(NA)

	# smooth and eliminate outliers

	# USE THIS ONE FOR PHENOCAM
	# smooth_evi2 <- try(SplineAndOutlierRemoval(x, dates, out_quant=pheno_pars$out_quant, spline_spar=pheno_pars$spline_spar, out_iterations=pheno_pars$out_iterations, pred_dates=pred_dates), silent=T)

	smooth_evi2 <- try(SplineAndOutlierRemoval(filt_evi2, dates, out_quant=pheno_pars$out_quant, spline_spar=pheno_pars$spline_spar, out_iterations=pheno_pars$out_iterations, pred_dates=pred_dates), silent=T)
	if(inherits(smooth_evi2, 'try-error'))
		return(NA)

	# find valid peaks in the time series
	valid_peaks <- try(FilterPeaks(smooth_evi2, FindPotentialPeaks_minmax(smooth_evi2), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance), silent=T)
	if(inherits(valid_peaks, 'try-error'))
		return(NA)

	# find full segments
	full_segs <- try(GetFullSegs(smooth_evi2, valid_peaks, max_seg_length=pheno_pars$max_seg_length, min_seg_amplitude=pheno_pars$min_seg_amplitude, agg_amp_frac=pheno_pars$agg_amp_frac), silent=T)
	if(inherits(full_segs, 'try-error')){
		return(NA)
	}

	# get PhenoDates
	pheno_dates <- try(GetPhenoDates(full_segs, smooth_evi2, pred_dates, pheno_pars$gup_threshes, pheno_pars$gdown_threshes), silent=T)
	if(inherits(pheno_dates, "try-error")){
		return(NA)
	}

	return(pheno_dates)
}


#==============================================================================
#==============================================================================
# extraction functions
#==============================================================================
#==============================================================================

#-----------------------------------------------
# works on a "result" object from DoPhenology, returns
# DOY for the ith phenometric for the specified year and cycle number
ExtractDOYYear <- function(x, i, year, cycle=1){
	if(is.na(x)){
		return(NA)
	}

	tmp <- try(x[[1]][[i]], silent=T)
	if(inherits(tmp, "try-error")){
		return(NA)
	}

	my_dates <- try(
		as.Date(
			tmp,
			origin=as.Date("1970-1-1")
		),
		silent=T
	)
	if(inherits(my_dates, "try-error")){
		return(NA)
	}

	years <- try(
		as.integer(
			strftime(
				my_dates,
				format="%Y"
			)
		),
		silent=T
	)
	if(inherits(years, "try-error")){
		return(NA)
	}

	matched_year_date <- my_dates[which(years == year)][cycle]

	if(length(matched_year_date > 0)){
		return(as.integer(strftime(matched_year_date, format="%j")))
	}else{
		return(NA)
	}
}



#-----------------------------------------------
# works on a "result" object from DoPhenology, and retrieves
# the ith threshold phenometric of the nth vegetation cycle
# with i_key phenometric between start_date and end_date
ExtractSingleDate <- function(x, i, i_key, start_date, end_date, n=1){
	if(is.na(x)){
		return(NA)
	}

	# something funky with certain row blocks and the number of nested lists
	# this corrects for those problems during the extraction process
	tmp <- try(x[[1]], silent=T)
	if(inherits(tmp, "try-error")){
		return(NA)
	}else if(typeof(tmp) == "list"){
		# tmp <- tmp[[i]]
		# x <- unlist(unlist(x, rec=F), rec=F)
		x <- unlist(x, rec=F, use.names=F)
	}else{
		# tmp <- x[[i]]
		# x <- unlist(x, rec=F)
		x <- x
	}

	# convert to dates
	key_dates <- try(as.Date(x[[i_key]], origin="1970-1-1"), silent=T)
	if(inherits(key_dates, "try-error")){
		return(NA)
	}

	# find indices for i_key phenometrics in the valid window
	valid_indices <- try(which(key_dates >= start_date & key_dates <= end_date), silent=T)
	if(inherits(valid_indices, "try-error")){
		return(NA)
	}
	# subset_dates <- my_dates[my_dates >= start_date & my_dates <= end_date]

	# get the ith phenometric(s)
	valid_dates <- try(x[[i]][valid_indices], silent=T)
	if(inherits(valid_dates, "try-error")){
		return(NA)
	}

	# get the specified cycle, checking if it exists first
	if(length(valid_dates) < n){
		return(NA)
	}else{
		return(valid_dates[n])
	}
}


#-----------------------------------------------
# works on a "result" object from DoPhenology, and retrieves
# the ith threshold phenometric of the nth vegetation cycle
# with i_key phenometric between start_date and end_date
ExtractPhenoMetric <- function(x, i, i_key, start_date, end_date, n=1){
	if(is.na(x)){
		return(NA)
	}

	# something funky with certain row blocks and the number of nested lists
	# this corrects for those problems during the extraction process
	tmp <- try(x[[1]], silent=T)
	if(inherits(tmp, "try-error")){
		return(NA)
	}else if(typeof(tmp) == "list"){
		# tmp <- tmp[[i]]
		# x <- unlist(unlist(x, rec=F), rec=F)
		x <- unlist(x, rec=F, use.names=F)
	}else{
		# tmp <- x[[i]]
		# x <- unlist(x, rec=F)
		x <- x
	}

	# convert to dates
	key_dates <- try(as.Date(x[[i_key]], origin="1970-1-1"), silent=T)
	if(inherits(key_dates, "try-error")){
		return(NA)
	}

	# find indices for i_key phenometrics in the valid window
	valid_indices <- try(which(key_dates >= start_date & key_dates <= end_date), silent=T)
	if(inherits(valid_indices, "try-error")){
		return(NA)
	}
	# subset_dates <- my_dates[my_dates >= start_date & my_dates <= end_date]

	# get the ith phenometric(s)
	valid_dates <- try(x[[i]][valid_indices], silent=T)
	if(inherits(valid_dates, "try-error")){
		return(NA)
	}

	# get the specified cycle, checking if it exists first
	if(length(valid_dates) < n){
		return(NA)
	}else{
		return(valid_dates[n])
	}
}

#-----------------------------------------------
# works on a "result" object from DoPhenology, returns
# between [start.date, end.date)
ExtractNumCyclesWindow <- function(x, i, start_date, end_date){
	if(is.na(x)){
		return(NA)
	}

	# try and get the specified phenometric list

	# something funky with certain row blocks and the number of nested lists
	# this corrects for those problems during the extraction process
	tmp <- try(x[[1]], silent=T)
	if(inherits(tmp, "try-error")){
		return(NA)
	}
	if(typeof(tmp) == "list"){
		tmp <- tmp[[i]]
	}else{
		tmp <- tmp
	}


	# tmp <- try(x[[1]][[i]], silent=T)
	# if(inherits(tmp, "try-error")){
	# 	return(NA)
	# }

	# try to convert the phenometric list to dates
	my_dates <- try(
		as.Date(
			tmp,
			origin=as.Date("1970-1-1")
		),
		silent=T
	)
	if(inherits(my_dates, "try-error")){
		return(NA)
	}

	# select only those dates between start and end date
	win_dates <- try(my_dates[my_dates >= start_date & my_dates < end_date])
	if(inherits(win_dates, "try-error")){
		return(NA)
	}

	# calculate the number of cycles
	win_num_cycles <- try(length(win_dates))
	if(inherits(win_num_cycles, "try-error")){
		return(NA)
	}else{
		return(win_num_cycles)
	}
}

#-----------------------------------------------
ExtractDates <- function(x, i){
	if(is.na(x)){
		return(NA)
	}

	tmp <- try(x[[1]][[i]], silent=T)
	if(inherits(tmp, "try-error")){
		return(NA)
	}

	my_dates <- try(
		as.Date(
			tmp,
			origin=as.Date("1970-1-1")
		),
		silent=T
	)
	if(inherits(my_dates, "try-error")){
		return(NA)
	}else{
		return(my_dates)
	}

}

#-----------------------------------------------
ExtractGSL <- function(x, i_gup=1, i_gdown=2){
	if(is.na(x)){
		return(NA)
	}

	tmp_gup <- try(x[[1]][[i_gup]], silent=T)
	if(inherits(tmp_gup, "try-error")){
		return(NA)
	}

	tmp_gdown <- try(x[[1]][[i_gdown]], silent=T)
	if(inherits(tmp_gdown, "try-error")){
		return(NA)
	}

	gup_dates <- try(
		as.Date(
			tmp_gup,
			origin=as.Date("1970-1-1")
		),
		silent=T
	)
	if(inherits(gup_dates, "try-error")){
		return(NA)
	}

	gdown_dates <- try(
		as.Date(
			tmp_gdown,
			origin=as.Date("1970-1-1")
		),
		silent=T
	)
	if(inherits(gdown_dates, "try-error")){
		return(NA)
	}

	gsl <- try(gdown_dates - gup_dates)
	if(inherits(gsl, "try-error")){
		return(NA)
	}else{
		return(gsl)
	}
}

#-----------------------------------------------
# works on a "result" object from DoPhenology, returns
# between [start.date, end.date)
ExtractGSLWindow <- function(x, i_gup, i_gdown, start.date, end.date, cycle=1){
	if(is.na(x)){
		return(NA)
	}

	# try and get the g_up/g_down phenometrics
	tmp_gup <- try(x[[1]][[i_gup]], silent=T)
	if(inherits(tmp_gup, "try-error")){
		return(NA)
	}

	tmp_gdown <- try(x[[1]][[i_gdown]], silent=T)
	if(inherits(tmp_gdown, "try-error")){
		return(NA)
	}


	# try to convert the phenometric lists to dates
	my_dates_gup <- try(
		as.Date(
			tmp_gup,
			origin=as.Date("1970-1-1")
		),
		silent=T
	)
	if(inherits(my_dates_gup, "try-error")){
		return(NA)
	}

	my_dates_gdown <- try(
		as.Date(
			tmp_gdown,
			origin=as.Date("1970-1-1")
		),
		silent=T
	)
	if(inherits(my_dates_gdown, "try-error")){
		return(NA)
	}


	# select only those gup dates between start and end date
	gup_dates_inds <- try(which(my_dates_gup >= start.date & my_dates_gup < end.date))
	if(inherits(gup_dates_inds, "try-error")){
		return(NA)
	}

	# get the g_uyp/g_down dates for the specified cycle in the window
	gup_date <- my_dates_gup[gup_dates_inds[cycle]]
	gdown_date <- my_dates_gdown[gup_dates_inds[cycle]]

	# calculate GSL
	tmp_gsl <- gdown_date - gup_date

	return(tmp_gsl)
}



#-----------------------------------------------
ExtractNumCycles <- function(x, i=1){
	if(is.na(x)){
		return(NA)
	}

	tmp <- try(x[[1]][[i]], silent=T)
	if(inherits(tmp, "try-error")){
		return(NA)
	}

	num_cycles <- try(length(tmp))
	if(inherits(num_cycles, "try-error")){
		return(NA)
	}else{
		return(num_cycles)
	}
}

#==============================================================================
#==============================================================================
# plotting functions
#==============================================================================
#==============================================================================

# plot with a linear % cutoff stretch and proper legend
PlotStretch <- function(r, thresh_percent=c(2, 2), thresh_q=NA, na.rm=T, legend_ticks=5, dates=T, legend_round=3, ...){
	if(is.na(thresh_q)){
		q <- quantile(r, c(0, thresh_percent[1] / 100, 1 - (thresh_percent[2] / 100), 1), na.rm=na.rm)
	}else{
		if(length(thresh_q) != 4){
			print("If specified, thresh_q must have length 4: (min, low_cut, high_cut, max)")
			return(NA)
		}else{
			q <- thresh_q
		}
	}

	if(dates){
		breaks <- unique(round(c(q[1], seq(q[2], q[3], len=254), q[4])))
	}else{
		breaks <- unique(c(q[1], seq(q[2], q[3], len=254), q[4]))
	}

	pal <-  colorRampPalette(brewer.pal(11, "Spectral"))
	plot(r, breaks=breaks, col=pal(length(breaks) - 1), legend=F, ...)

	# make the legend
	# tmp_r <- raster(matrix(q[2]:q[3]))
	tmp_r <- raster(matrix(seq(q[2], q[3], len=255)))
	# legend_at <- round(seq(q[2], q[3], len=legend_ticks))
	legend_at <- seq(q[2], q[3], len=legend_ticks)
	# legend_at <- round(seq(q[2], q[3], len=legend_ticks), legend_round)

	if(dates){
		legend_at <- round(legend_at)
		legend_labels <- c(
							paste("<", as.Date(legend_at[1], origin="1970-1-1")),
							as.character(as.Date(legend_at[2:(length(legend_at) - 1)], origin="1970-1-1")),
							paste("<", as.Date(legend_at[length(legend_at)], origin="1970-1-1"))
						)
	}else{
		legend_labels_tmp <- round(legend_at, legend_round)
		legend_labels <- c(paste("<", legend_labels_tmp[1]), legend_labels_tmp[2:(length(legend_labels_tmp) - 1)], paste(">", legend_labels_tmp[length(legend_labels_tmp)]))
	}

	plot(
		tmp_r,
		legend.only=T,
		col=pal(length(breaks) - 1),
		axis.args=list(at=legend_at, labels=legend_labels)
	)
}

#---------------------------------------------------------------------
# function for plotting an individual time series row w/ or w/o seg metrics
PlotSeries <- function(x, dates, NA_value=32767, legend=T, seg_metrics=NA){
  mycols <- c("#636363", "#DE2D26", "#3182BD", "#31A354", "#756BB1", "#E6550D")
  x[x==NA_value] <- NA
  data_length <- length(x) / 3
  filtered <- x[1:data_length]
  resids <- x[(data_length + 1):(2 * data_length)]
  snowflags <- x[(2 * data_length + 1):(3 * data_length)]
  cols <- mycols[snowflags + 1]
  pchs <- snowflags + 1

  # start plotting
  plot(c(dates, dates), c(filtered, filtered + resids), type="n", xlab="", ylab="EVI2")

  # more extensive plotting if seg_metrics are available
  if(!is.na(seg_metrics)){
    for(seg_metric in seg_metrics){
      # plot the segment data
      gup_x_tmp <- as.Date(seg_metric$peak - seg_metric$length_gup, origin="1970-1-1"):as.Date(seg_metric$peak, origin="1970-1-1")
      gup_x <- c(gup_x_tmp, rev(gup_x_tmp))
      gup_y_tmp <- filtered[which(dates == as.Date(seg_metric$peak - seg_metric$length_gup, origin="1970-1-1")):which(dates == as.Date(seg_metric$peak, origin="1970-1-1"))]
      gup_y <- c(gup_y_tmp, rep(par()$usr[3], length(gup_x_tmp)))
      polygon(x=gup_x, y=gup_y, border=NA, col=rgb(0.87, 0.87, 0.87))

      gdown_x_tmp <- as.Date(seg_metric$peak, origin="1970-1-1"):as.Date(seg_metric$peak + seg_metric$length_gdown, origin="1970-1-1")
      gdown_x <- c(gdown_x_tmp, rev(gdown_x_tmp))
      gdown_y_tmp <- filtered[which(dates == as.Date(seg_metric$peak, origin="1970-1-1")):which(dates == as.Date(seg_metric$peak + seg_metric$length_gdown, origin="1970-1-1"))]
      gdown_y <- c(gdown_y_tmp, rep(par()$usr[3], length(gdown_x_tmp)))
      polygon(x=gdown_x, y=gdown_y, border=NA, col=rgb(0.82, 0.82, 0.82))

      # this creates simple rectangle fills:
      # polygon(x=c(rep(as.Date(seg_metric$peak - seg_metric$length_gup, origin="1970-1-1"), 2), rep(as.Date(seg_metric$peak, origin="1970-1-1"), 2)), y=c(par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3]), border=NA, col=rgb(0.87, 0.87, 0.87))
      # polygon(x=c(rep(as.Date(seg_metric$peak, origin="1970-1-1"), 2), rep(as.Date(seg_metric$peak + seg_metric$length_gdown, origin="1970-1-1"), 2)), y=c(par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3]), border=NA, col=rgb(0.82, 0.82, 0.82))
    } # end seg loop

    # add the data/spline with fill/miss/snow info
    points(dates, filtered + resids, type="p", pch=pchs, cex=0.75, col=cols)
    points(dates, filtered, type="l", lwd=2, col="darkgrey")

    for(seg_metric in seg_metrics){
      # add vertical lines for phenometrics
      abline(v=as.Date(seg_metric$ogi, origin="1970-1-1"), lty=2, col=mycols[1])
      abline(v=as.Date(seg_metric$midgup, origin="1970-1-1"), lty=2, col=mycols[1])
      abline(v=as.Date(seg_metric$mat, origin="1970-1-1"), lty=2, col=mycols[1])
      abline(v=as.Date(seg_metric$peak, origin="1970-1-1"), lty=2, col=mycols[1])
      abline(v=as.Date(seg_metric$sen, origin="1970-1-1"), lty=2, col=mycols[1])
      abline(v=as.Date(seg_metric$midgdown, origin="1970-1-1"), lty=2, col=mycols[1])
      abline(v=as.Date(seg_metric$dor, origin="1970-1-1"), lty=2, col=mycols[1])

      # annotate with QA information

    } # end seg loop

    # create a legend
    legend("topleft", legend=c("Spline-smoothed", "Flag=0", "Flag=1", "Flag=2", "Flag=3", "Flag=4"), pch=c(NA, 1:5), col=c("darkgrey", mycols[1:5]), lwd=c(1, rep(NA, 5)), bg="white")
  }else{
    # plot(c(dates, dates), c(filtered, filtered + resids), type="n", xlab="", ylab="EVI2")
    points(dates, filtered + resids, type="p", pch=pchs, cex=0.75, col=cols)
    points(dates, filtered, type="l", lwd=2, col="darkgrey")
    legend("topleft", legend=c("Spline-smoothed", "Flag=0", "Flag=1", "Flag=2", "Flag=3", "Flag=4"), pch=c(NA, 1:5), col=c("darkgrey", mycols[1:5]), lwd=c(1, rep(NA, 5)), bg="white")
  }
}

#---------------------------------------------------------------------
# PlotSeries <- function(x, dates, out_quant=0.99, spline_spar=0.3, out_iterations=1, min_peak_to_peak_distance=180, min_seg_amplitude=0.15, max_seg_length=180, gup_threshes=c(0.1, 0.95), gdown_threshes=c(0.5, 0.05)){
PlotSeriesPhenocam <- function(x, dates, pheno_pars, pred_dates=NA){
	# PLOT
	par(mar=c(3, 4, 3, 1))

	# if prediction dates aren't given, assume a daily time series
	if(is.na(pred_dates)) pred_dates <- min(dates, na.rm=T):max(dates, na.rm=T)

	# get a time series of snow-filtered EVI2
	filt_evi2 <- x
	# filt_evi2 <- try(SnowFilterEVI2(x, pheno_pars$evi2_snow_quant, pheno_pars$ndsi_thresh), silent=T)
	# if(inherits(filt_evi2, 'try-error'))
	# 	return(NA)

	# PLOT
	# plot(dates, filt_evi2, type="n", col=rgb(0.2, 0.2, 0.2), pch=4, cex=1, xlab="", ylab="EVI2")
	plot(dates, filt_evi2, type="n", col=rgb(0.2, 0.2, 0.2), pch=4, cex=1, xlab="", ylab="EVI2", ylim=c(0, 0.9))

	# smooth and eliminate outliers
	# smooth_evi2 <- try(SplineAndOutlierRemoval(filt_evi2, dates, out_sigma=out_sigma, spline_spar=spline_spar, out_iterations=out_iterations, pred_dates=pred_dates), silent=T)
	smooth_evi2 <- try(SplineAndOutlierRemoval(filt_evi2, dates, out_quant=pheno_pars$out_quant, spline_spar=pheno_pars$spline_spar, out_iterations=pheno_pars$out_iterations, pred_dates=pred_dates), silent=T)
	if(inherits(smooth_evi2, 'try-error'))
		return(NA)

	# PLOT
	# points(pred_dates, smooth_evi2, type="l", col=brewer.pal(5, "Blues")[4], lwd=1.5)

	# find valid peaks in the time series
	valid_peaks <- try(FilterPeaks(smooth_evi2, FindPotentialPeaks_minmax(smooth_evi2), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance), silent=T)

	if(inherits(valid_peaks, 'try-error'))
		return(NA)

	# find full segments
	full_segs <- try(GetFullSegs(smooth_evi2, valid_peaks, max_seg_length=pheno_pars$max_seg_length, min_seg_amplitude=pheno_pars$min_seg_amplitude, pheno_pars$agg_amp_frac))
	if(inherits(full_segs, 'try-error')){
		return(NA)
	}else if(!is.na(full_segs)){
		# PLOT
		parcoords <- par()$usr #get plot coordinates
		# greenup.color <- "lightgreen"
		# greendown.color <- "pink"
		greenup.color <- brewer.pal(5, "Greens")[3]
		greendown.color <- brewer.pal(5, "Reds")[3]

		transp=0.5
		greenup.color <- rgb(t(col2rgb(greenup.color))/255, alpha=transp)
		greendown.color <- rgb(t(col2rgb(greendown.color))/255, alpha=transp)

		for(seg in full_segs){
			gup.poly.x <- c(rep(pred_dates[seg[1]], 2), rep(pred_dates[seg[2]], 2))
			gup.poly.y <- c(parcoords[3], parcoords[4], parcoords[4], parcoords[3])
			gdown.poly.x <- c(rep(pred_dates[seg[2]], 2), rep(pred_dates[seg[3]], 2))
			gdown.poly.y <- c(parcoords[3], parcoords[4], parcoords[4], parcoords[3])

			polygon(gup.poly.x, gup.poly.y, col=greenup.color, border=NA)
			polygon(gdown.poly.x, gdown.poly.y, col=greendown.color, border=NA)
		}
	}

	# get PhenoDates
	pheno_dates <- try(GetPhenoDates(full_segs, smooth_evi2, pred_dates, pheno_pars$gup_threshes, pheno_pars$gdown_threshes))
	if(!inherits(pheno_dates, "try-error")){
		coords <- par()$usr
		for(i in 1:length(pheno_pars$gup_threshes)){
			abline(v=pheno_dates[[1]][[i]], lty=2, col="darkgreen")
			text(pheno_dates[[1]][[i]] - 13, y=coords[4] - (0.1 * (coords[4] - coords[3])), labels=as.Date(pheno_dates[[1]][[i]], origin=as.Date("1970-1-1")), srt=90, col="darkgreen", cex=0.7)
		}
		for(i in 1:length(pheno_pars$gdown_threshes)){
			abline(v=pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]], lty=2, col="darkred")
			text(pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]] + 13, y=coords[4] - (0.1 * (coords[4] - coords[3])), labels=as.Date(pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]], origin=as.Date("1970-1-1")), srt=90, col="darkred", cex=0.7)
		}
	}

	points(dates, filt_evi2, type="p", col=rgb(0.2, 0.2, 0.2), pch=16, cex=0.75)
	points(pred_dates, smooth_evi2, type="l", col=brewer.pal(5, "Blues")[4], lwd=1.5)
}

PlotSeriesLandsat <- function(x, dates, pheno_pars, pred_dates=NA){
	# PLOT
	par(mar=c(3, 3, 3, 1))

	# if prediction dates aren't given, assume a daily time series
	# if(is.na(pred_dates)) pred_dates <- min(dates, na.rm=T):max(dates, na.rm=T)

	# get a time series of snow-filtered EVI2
	filt_evi2 <- x
	# filt_evi2 <- try(SnowFilterEVI2(x, pheno_pars$evi2_snow_quant, pheno_pars$ndsi_thresh), silent=T)
	# if(inherits(filt_evi2, 'try-error'))
	# 	return(NA)
	# tmp <- LandsatScreenSmooth(x, dates)
	tmp <- LandsatSmoother(x, dates, fill_quant=pheno_pars$LandsatFillQuant, spike_thresh=pheno_pars$LandsatSpikeThresh, min_resid=pheno_pars$LandsatMinResid, fill_doy=pheno_pars$LandsatFillDOY)
	pred_dates <- tmp[[1]]
	smooth_evi2 <- tmp[[2]]

	# PLOT
	plot(dates, filt_evi2, type="n", col=rgb(0.2, 0.2, 0.2), pch=4, cex=1, xlab="", ylab="EVI2")

	# smooth and eliminate outliers
	# smooth_evi2 <- try(SplineAndOutlierRemoval(filt_evi2, dates, out_sigma=out_sigma, spline_spar=spline_spar, out_iterations=out_iterations, pred_dates=pred_dates), silent=T)
	# smooth_evi2 <- try(SplineAndOutlierRemoval(filt_evi2, dates, out_quant=pheno_pars$out_quant, spline_spar=pheno_pars$spline_spar, out_iterations=pheno_pars$out_iterations, pred_dates=pred_dates), silent=T)
	# if(inherits(smooth_evi2, 'try-error'))
	# 	return(NA)

	# PLOT
	# points(pred_dates, smooth_evi2, type="l", col=brewer.pal(5, "Blues")[4], lwd=1.5)

	# find valid peaks in the time series
	valid_peaks <- try(FilterPeaks(smooth_evi2, FindPotentialPeaks_minmax(smooth_evi2), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance), silent=T)

	if(inherits(valid_peaks, 'try-error'))
		return(NA)

	# find full segments
	full_segs <- try(GetFullSegs(smooth_evi2, valid_peaks, max_seg_length=pheno_pars$max_seg_length, min_seg_amplitude=pheno_pars$min_seg_amplitude, pheno_pars$agg_amp_frac))
	if(inherits(full_segs, 'try-error')){
		return(NA)
	}else if(!is.na(full_segs)){
		# PLOT
		parcoords <- par()$usr #get plot coordinates
		# greenup.color <- "lightgreen"
		# greendown.color <- "pink"
		greenup.color <- brewer.pal(5, "Greens")[3]
		greendown.color <- brewer.pal(5, "Reds")[3]

		transp=0.5
		greenup.color <- rgb(t(col2rgb(greenup.color))/255, alpha=transp)
		greendown.color <- rgb(t(col2rgb(greendown.color))/255, alpha=transp)

		for(seg in full_segs){
			gup.poly.x <- c(rep(pred_dates[seg[1]], 2), rep(pred_dates[seg[2]], 2))
			gup.poly.y <- c(parcoords[3], parcoords[4], parcoords[4], parcoords[3])
			gdown.poly.x <- c(rep(pred_dates[seg[2]], 2), rep(pred_dates[seg[3]], 2))
			gdown.poly.y <- c(parcoords[3], parcoords[4], parcoords[4], parcoords[3])

			polygon(gup.poly.x, gup.poly.y, col=greenup.color, border=NA)
			polygon(gdown.poly.x, gdown.poly.y, col=greendown.color, border=NA)
		}
	}

	# get PhenoDates
	pheno_dates <- try(GetPhenoDates(full_segs, smooth_evi2, pred_dates, pheno_pars$gup_threshes, pheno_pars$gdown_threshes))
	if(!inherits(pheno_dates, "try-error")){
		coords <- par()$usr
		for(i in 1:length(pheno_pars$gup_threshes)){
			abline(v=pheno_dates[[1]][[i]], lty=2, col="darkgreen")
			text(pheno_dates[[1]][[i]] - 13, y=coords[4] - (0.1 * (coords[4] - coords[3])), labels=as.Date(pheno_dates[[1]][[i]], origin=as.Date("1970-1-1")), srt=90, col="darkgreen", cex=0.7)
		}
		for(i in 1:length(pheno_pars$gdown_threshes)){
			abline(v=pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]], lty=2, col="darkred")
			text(pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]] + 13, y=coords[4] - (0.1 * (coords[4] - coords[3])), labels=as.Date(pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]], origin=as.Date("1970-1-1")), srt=90, col="darkred", cex=0.7)
		}
	}

	points(dates, filt_evi2, type="p", col=rgb(0.2, 0.2, 0.2), pch=16, cex=0.75)
	points(pred_dates, smooth_evi2, type="l", col=brewer.pal(5, "Blues")[4], lwd=1.5)
}

#---------------------------------------------------------------------
# PlotSeries <- function(x, dates, out_quant=0.99, spline_spar=0.3, out_iterations=1, min_peak_to_peak_distance=180, min_seg_amplitude=0.15, max_seg_length=180, gup_threshes=c(0.1, 0.95), gdown_threshes=c(0.5, 0.05)){
PlotSeries_old <- function(x, dates, pheno_pars){
	# PLOT
	par(mar=c(3, 3, 3, 1))

	pred_dates <- min(dates, na.rm=T):max(dates, na.rm=T)

	# get a time series of snow-filtered EVI2
	# filt_evi2 <- try(SnowFilterEVI2_mv(x, pheno_pars$evi2_snow_quant, pheno_pars$ndsi_thresh), silent=T)
	# if(inherits(filt_evi2, 'try-error'))
	# 	return(NA)

	filt_return <- try(SnowFilterEVI2_mv(x, dates, evi2_snow_quant=pheno_pars$evi2_snow_quant, ndsi_thresh=pheno_pars$ndsi_thresh, max_snow_fill_ratio=pheno_pars$max_snow_fill_ratio), silent=T)
	if(inherits(filt_return, 'try-error'))
		return(NA)

	filt_evi2 <- filt_return[[1]] # The EVI2 time series
	snow_fills <- filt_return[[2]] # Boolean vector indicating snow fills

	# PLOT
	plot(dates, filt_evi2, type="n", col=rgb(0.2, 0.2, 0.2), pch=4, cex=1, xlab="", ylab="EVI2")

	# smooth and eliminate outliers
	# smooth_evi2 <- try(SplineAndOutlierRemoval(filt_evi2, dates, out_sigma=out_sigma, spline_spar=spline_spar, out_iterations=out_iterations, pred_dates=pred_dates), silent=T)
	smooth_evi2 <- try(SplineAndOutlierRemoval(filt_evi2, dates, out_quant=pheno_pars$out_quant, spline_spar=pheno_pars$spline_spar, out_iterations=pheno_pars$out_iterations, pred_dates=pred_dates), silent=T)
	if(inherits(smooth_evi2, 'try-error'))
		return(NA)

	# PLOT
	# points(pred_dates, smooth_evi2, type="l", col=brewer.pal(5, "Blues")[4], lwd=1.5)

	# find valid peaks in the time series
	valid_peaks <- try(FilterPeaks(smooth_evi2, FindPotentialPeaks_minmax(smooth_evi2), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance), silent=T)

	if(inherits(valid_peaks, 'try-error'))
		return(NA)

	# find full segments
	full_segs <- try(GetFullSegs(smooth_evi2, valid_peaks, max_seg_length=pheno_pars$max_seg_length, min_seg_amplitude=pheno_pars$min_seg_amplitude, pheno_pars$agg_amp_frac))
	if(inherits(full_segs, 'try-error')){
		return(NA)
	}else if(!is.na(full_segs)){
		# PLOT
		parcoords <- par()$usr #get plot coordinates
		# greenup.color <- "lightgreen"
		# greendown.color <- "pink"
		greenup.color <- brewer.pal(5, "Greens")[3]
		greendown.color <- brewer.pal(5, "Reds")[3]

		transp=0.5
		greenup.color <- rgb(t(col2rgb(greenup.color))/255, alpha=transp)
		greendown.color <- rgb(t(col2rgb(greendown.color))/255, alpha=transp)

		for(seg in full_segs){
			gup.poly.x <- c(rep(pred_dates[seg[1]], 2), rep(pred_dates[seg[2]], 2))
			gup.poly.y <- c(parcoords[3], parcoords[4], parcoords[4], parcoords[3])
			gdown.poly.x <- c(rep(pred_dates[seg[2]], 2), rep(pred_dates[seg[3]], 2))
			gdown.poly.y <- c(parcoords[3], parcoords[4], parcoords[4], parcoords[3])

			polygon(gup.poly.x, gup.poly.y, col=greenup.color, border=NA)
			polygon(gdown.poly.x, gdown.poly.y, col=greendown.color, border=NA)
		}
	}

	# get PhenoDates
	pheno_dates <- try(GetPhenoDates(full_segs, smooth_evi2, pred_dates, pheno_pars$gup_threshes, pheno_pars$gdown_threshes))
	if(!inherits(pheno_dates, "try-error")){
		coords <- par()$usr
		for(i in 1:length(pheno_pars$gup_threshes)){
			abline(v=pheno_dates[[1]][[i]], lty=2, col="darkgreen")
			text(pheno_dates[[1]][[i]] - 13, y=coords[4] - (0.1 * (coords[4] - coords[3])), labels=as.Date(pheno_dates[[1]][[i]], origin=as.Date("1970-1-1")), srt=90, col="darkgreen", cex=0.7)
		}
		for(i in 1:length(pheno_pars$gdown_threshes)){
			abline(v=pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]], lty=2, col="darkred")
			text(pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]] + 13, y=coords[4] - (0.1 * (coords[4] - coords[3])), labels=as.Date(pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]], origin=as.Date("1970-1-1")), srt=90, col="darkred", cex=0.7)
		}
	}

	points(dates, filt_evi2, type="p", col=rgb(0.2, 0.2, 0.2), pch=16, cex=0.75)
	points(pred_dates, smooth_evi2, type="l", col=brewer.pal(5, "Blues")[4], lwd=1.5)
}

#---------------------------------------------------------------------
# PlotSeries <- function(x, dates, out_quant=0.99, spline_spar=0.3, out_iterations=1, min_peak_to_peak_distance=180, min_seg_amplitude=0.15, max_seg_length=180, gup_threshes=c(0.1, 0.95), gdown_threshes=c(0.5, 0.05)){
PlotSeriesSIF <- function(x, dates, pheno_pars, throttle_zero=T){
	# PLOT
	par(mar=c(3, 3, 3, 1))

	pred_dates <- min(dates, na.rm=T):max(dates, na.rm=T)

	if(throttle_zero){
		# x[x < 0] <- 0
		x[x < 0] <- NA
	}
	filt_evi2 <- x
	# # get a time series of snow-filtered EVI2
	# filt_evi2 <- try(SnowFilterEVI2(x, pheno_pars$evi2_snow_quant, pheno_pars$ndsi_thresh), silent=T)
	# if(inherits(filt_evi2, 'try-error'))
	# 	return(NA)

	# PLOT
	plot(dates, filt_evi2, type="n", col=rgb(0.2, 0.2, 0.2), pch=4, cex=1, xlab="", ylab="EVI2")

	# smooth and eliminate outliers
	# smooth_evi2 <- try(SplineAndOutlierRemoval(filt_evi2, dates, out_sigma=out_sigma, spline_spar=spline_spar, out_iterations=out_iterations, pred_dates=pred_dates), silent=T)
	smooth_evi2 <- try(SplineAndOutlierRemoval(filt_evi2, dates, out_quant=pheno_pars$out_quant, spline_spar=pheno_pars$spline_spar, out_iterations=pheno_pars$out_iterations, pred_dates=pred_dates), silent=T)
	if(inherits(smooth_evi2, 'try-error'))
		return(NA)

	# PLOT
	# points(pred_dates, smooth_evi2, type="l", col=brewer.pal(5, "Blues")[4], lwd=1.5)

	# find valid peaks in the time series
	valid_peaks <- try(FilterPeaks(smooth_evi2, FindPotentialPeaks_minmax(smooth_evi2), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance), silent=T)

	if(inherits(valid_peaks, 'try-error'))
		return(NA)

	# find full segments
	full_segs <- try(GetFullSegs(smooth_evi2, valid_peaks, max_seg_length=pheno_pars$max_seg_length, min_seg_amplitude=pheno_pars$min_seg_amplitude, pheno_pars$agg_amp_frac))
	if(inherits(full_segs, 'try-error')){
		return(NA)
	}else if(!is.na(full_segs)){
		# PLOT
		parcoords <- par()$usr #get plot coordinates
		# greenup.color <- "lightgreen"
		# greendown.color <- "pink"
		greenup.color <- brewer.pal(5, "Greens")[3]
		greendown.color <- brewer.pal(5, "Reds")[3]

		transp=0.5
		greenup.color <- rgb(t(col2rgb(greenup.color))/255, alpha=transp)
		greendown.color <- rgb(t(col2rgb(greendown.color))/255, alpha=transp)

		for(seg in full_segs){
			gup.poly.x <- c(rep(pred_dates[seg[1]], 2), rep(pred_dates[seg[2]], 2))
			gup.poly.y <- c(parcoords[3], parcoords[4], parcoords[4], parcoords[3])
			gdown.poly.x <- c(rep(pred_dates[seg[2]], 2), rep(pred_dates[seg[3]], 2))
			gdown.poly.y <- c(parcoords[3], parcoords[4], parcoords[4], parcoords[3])

			polygon(gup.poly.x, gup.poly.y, col=greenup.color, border=NA)
			polygon(gdown.poly.x, gdown.poly.y, col=greendown.color, border=NA)
		}
	}

	# get PhenoDates
	pheno_dates <- try(GetPhenoDates(full_segs, smooth_evi2, pred_dates, pheno_pars$gup_threshes, pheno_pars$gdown_threshes))
	if(!inherits(pheno_dates, "try-error")){
		coords <- par()$usr
		for(i in 1:length(pheno_pars$gup_threshes)){
			abline(v=pheno_dates[[1]][[i]], lty=2, col="darkgreen")
			text(pheno_dates[[1]][[i]] - 13, y=coords[4] - (0.1 * (coords[4] - coords[3])), labels=as.Date(pheno_dates[[1]][[i]], origin=as.Date("1970-1-1")), srt=90, col="darkgreen", cex=0.7)
		}
		for(i in 1:length(pheno_pars$gdown_threshes)){
			abline(v=pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]], lty=2, col="darkred")
			text(pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]] + 13, y=coords[4] - (0.1 * (coords[4] - coords[3])), labels=as.Date(pheno_dates[[1]][[length(pheno_pars$gup_threshes) + i]], origin=as.Date("1970-1-1")), srt=90, col="darkred", cex=0.7)
		}
	}

	points(dates, filt_evi2, type="p", col=rgb(0.2, 0.2, 0.2), pch=16, cex=0.75)
	points(pred_dates, smooth_evi2, type="l", col=brewer.pal(5, "Blues")[4], lwd=1.5)
}


# # testing the segmentation function
# x <- v[1505,]
# PlotSeries(x,dates,max_seg_length=300)
# filt_evi2 <- try(SnowFilterEVI2(x, evi2_snow_quant=0.02), silent=T)
# smooth_evi2 <- try(SplineAndOutlierRemoval(filt_evi2, dates, out_sigma=3, spline_spar=0.3, out_iterations=1, pred_dates=pred_dates), silent=T)
# valid_peaks <- try(FilterPeaks(smooth_evi2, FindPotentialPeaks_minmax(smooth_evi2), min_peak_to_peak_distance=180), silent=T)
# full_segs <- GetFullSegs(smooth_evi2, valid_peaks, max_seg_length=180, min_seg_amplitude=0.10)
# plot(dates, filt_evi2)
# points(pred_dates, smooth_evi2, type="l", col=4)
# lty <- 1
# for(seg in full_segs){
# 	abline(v=pred_dates[seg], lty=lty, lwd=2)
# 	lty <- lty + 1
# }
#


#--------------------------------------------------------
GetSnowFreeQuants <- function(start_row, datasets, quants=c(0.02, 0.05, 0.5, 0.95, 0.98), rows_to_do, ndsi_thresh=-0.2, scale_factor=1e4){
	# load a data chunk
	system.time(x <- getValuesGDAL(datasets, start_row=start_row, n=rows_to_do))

	# calculate indices of particular SDS
	band_length <- dim(x)[2] / 5
	red_inds <- 1:band_length
	nir_inds <- (band_length + 1):(2 * band_length)
	band4_inds <- ((3 * band_length) + 1):(4 * band_length)
	band6_inds <- ((4 * band_length) + 1):(5 * band_length)

	# calculate EVI2 and NDSI
	evi2 <- 2.5 * ((x[, nir_inds] - x[, red_inds]) / (x[, nir_inds] + (2.4 * x[, red_inds]) + 1))
	ndsi <- (x[, band4_inds] - x[, band6_inds]) / (x[, band4_inds] + x[, band6_inds])
	rm(x) # clearup memory, hopefully

	# eliminate snow values
	evi2[ndsi > ndsi_thresh] <- NA
	rm(ndsi) # clearup memory, hopefully

	# get the requested quantiles in each time series
	tmp <- t(apply(evi2, 1, quantile, probs=quants, na.rm=T))
	rm(evi2) # clearup memory, hopefully

	# convert to integer
	tmp_int <- matrix(as.integer(scale_factor * tmp), nrow=dim(tmp)[1])

	return(tmp_int)
}

#---------------------------------------------------------------------
WriteQuants <- function(out_file, out_values, nbands, tile){
	# Appends a chunk of phenology data to a file, or creates and writes a chunk if the file doesn't exist.
	# Output is BIP Int16. Thanks to DSM for ENVI header code

	# write the data to a file as binary integers (4 bytes!)
	ff <- file(out_file, 'wb') # create the file and open for writing
	writeBin(out_values, ff)
	close(ff)

	# this code courtesy of Damien Sulla Menashe
	out_hdr <- paste(out_file, ".hdr", sep="")
	if(!file.exists(out_hdr)){
		# variables specifying the upper left corner of the modis grid and the tile/pixel lengths
		uly_map = 10007554.677
		ulx_map = -20015109.354
		lry_map = -10007554.677
		lrx_map = 20015109.354
		pix = 463.312716525
		dims = 2400

		# calculate the upper left corner of the current tile
		tile_h <- as.integer(substr(tile, 2, 3))
		tile_v <- as.integer(substr(tile, 5, 6))
		ulx = ulx_map + (tile_h * pix * dims)
		uly = uly_map - (tile_v * pix * dims)

		# create the header text
		temp_txt = paste("ENVI description = { C5 MCD43A4 NBAR-EVI2 Quantiles }\nlines = ", dims, "\nsamples = ", dims, "\nbands = ", nbands, "\nheader offset = 0\nfile type = ENVI Standard\ndata type = 3\ninterleave = bsq\nbyte order = 0\nmap info = {Sinusoidal, 1, 1,", ulx, ", ", uly, ", ", pix, ", ", pix, "}", "\ncoordinate system string = {PROJCS[\"Sinusoidal\",GEOGCS[\"GCS_unnamed ellipse\",DATUM[\"D_unknown\",SPHEROID[\"Unknown\",6371007.181,0]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Sinusoidal\"],PARAMETER[\"central_meridian\",0],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"Meter\",1]]}", sep="")

		# write the header to a file
		sink(out_hdr)
		cat(temp_txt)
		sink()
	}
}
