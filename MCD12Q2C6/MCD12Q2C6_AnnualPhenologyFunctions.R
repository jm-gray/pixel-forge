#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Josh Gray, Boston University, 2016
#
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Functions for reading/setting algorithm parameters, creating/scaling return values, etc.
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#---------------------------------------------------------------------
# DefaultPhenoParameters <- function(year_of_interest){
DefaultPhenoParameters <- function(){
	# Default algorithm parameters for AnnualPhenologyC6

	pheno_pars <- list(
		## valid segment definition ##
		min_seg_amplitude=0.15, # absolute amplitude threshold, not used if rel_amp_frac and rel_peak_frac are provided
		min_increase_length=30,
		max_increase_length=250,
		min_decrease_length=30,
		max_decrease_length=250,
		rel_amp_frac= 0.35, # segment amplitude must be >= rel_amp_frac * global_amp (max - min)
		rel_peak_frac=0.7, # EVI value at peak must be >= rel_peak_frac * global max

		# min_peak_to_peak_distance=120, # not used
		# agg_amp_frac=0.15, # not used
		# max_seg_length=200, # not used
		## date thresholds ##
	  ogi_thresh=0.15,
	  midgup_thresh=0.5,
	  mat_thresh=0.9, #???
	  sen_thresh=0.9, #???
	  midgdown_thresh=0.5,
	  dor_thresh=0.15,
		## QA parameters ##
		qual_buffer_days=14, # buffer around phenometric for QA calc
		qual_r2_weight=1, # rel weight of r2 in composite QA
		qual_fill_weight=4, # rel weight of fill fraction in composite QA
		## data in/out parameters ##
	  nbar_scale_factor=1e4,
	  nbar_NA_value=32767,
		out_float_scale=1e4,
		out_NA_value=32767
	)
	return(pheno_pars)
}

#---------------------------------------------------------------------
ReadPhenoParameters <- function(parameter_file, year_of_interest){
	# NOTE: just stubbed here, but for future LC-specific version

	return(DefaultPhenoParameters(year_of_interest))
}

#---------------------------------------------------------------------
PhenoReturnValue <- function(default_value=NA){
	# returns a blank AnnualPhenologyC6 return value

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
ScaleToIntegerAndSetNA <- function(annual_pheno_metrics, scale=1e4, outNA=32767){
	# this function takes a version of annual_pheno_metrics with floating point actual data values and scales & casts for output

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


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Functions for reading splined NBAR-EVI2, spline residual, and snow flag
# binary data files (signed, 2-byte integers, 365 layers, BIP)
# and writing chunks of algorithm output to binary files
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#---------------------------------------------------------------------
ReadSplinedNBAR <- function(in_file, lines_to_read=2400, start_line=1, samples=2400, dsize=2, nbands=365, s_flag=T, bsq_flag=F){
	# For reading from the binary, splined output files which contain 365 bands (days) of spline-smoothed EVI2 data/resids/QA

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
Get3YearDataChunk <- function(evi2_files, resid_files, snow_files, year_of_interest, start_line, lines_to_read){
	# retrieves 3 full years of splined C6 NBAR-EVI2 data, residuals, and snow flags as a data.frame where each row
	# is a time series of evi, resid, and snow: evi001, ..., evi365, resid001, ..., resid365, flag001, ..., flag365

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
		ReadSplinedNBAR(evi2_files[evi2_file_inds[1]], lines_to_read=lines_to_read, start_line=start_line),
		ReadSplinedNBAR(evi2_files[evi2_file_inds[2]], lines_to_read=lines_to_read, start_line=start_line),
		ReadSplinedNBAR(evi2_files[evi2_file_inds[3]], lines_to_read=lines_to_read, start_line=start_line),
		ReadSplinedNBAR(resid_files[resid_file_inds[1]], lines_to_read=lines_to_read, start_line=start_line),
		ReadSplinedNBAR(resid_files[resid_file_inds[2]], lines_to_read=lines_to_read, start_line=start_line),
		ReadSplinedNBAR(resid_files[resid_file_inds[3]], lines_to_read=lines_to_read, start_line=start_line),
		ReadSplinedNBAR(snow_files[snow_file_inds[1]], lines_to_read=lines_to_read, start_line=start_line),
		ReadSplinedNBAR(snow_files[snow_file_inds[2]], lines_to_read=lines_to_read, start_line=start_line),
		ReadSplinedNBAR(snow_files[snow_file_inds[3]], lines_to_read=lines_to_read, start_line=start_line)
	) # end data read/cbind data

	return(tmp)
}

#---------------------------------------------------------------------
WritePhenologyData <- function(out_file, pheno_values_bip, tile, nbands=44){
	# Appends a chunk of phenology data to a file, or creates and writes a chunk if the file doesn't exist.
	# Output is BIP Int16. Thanks to DSM for ENVI header code

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
	# write OGI output, converting from Int32 to Int16 and stacking as 2 bands:
	# gdal_translate -ot Int16 -of ENVI -b 10 -b 31 MCD12Q2C6_h12v04_2003 MCD12Q2C6_h12v04_2003_OGI
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Function for time series segmentation
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#------------------------------------------------------------------
FindPeaks <- function(x, mag_order=T){
	# Function to identify peaks in time series x (or troughs if x=-x), supports "flat top" peaks
  # if mag_order is TRUE, peaks are returned in order of increasing magnitude (of x)
  d <- diff(x)
  d_code <- (d > 0) + (2 * (d < 0)) # 0=no change, 1=inc, 2=dec
  peaks <- unlist(gregexpr("12", paste(d_code, collapse=""))) # no match is -1
  if(peaks[1] == -1) peaks <- NULL
  flat_peaks <- unlist(gregexpr("10+2", paste(d_code, collapse=""))) # no match is -1
  if(flat_peaks[1] == -1) flat_peaks <- NULL
  d_code_rle <- rle(d_code)
  flat_peaks <- flat_peaks + round(d_code_rle$l[match(flat_peaks, cumsum(d_code_rle$l)) + 1] / 2)
  # all_peaks <- c(ifelse(peaks[1] == -1, NULL, peaks + 1), ifelse(flat_peaks[1] == -1, NULL, flat_peaks + 1))
  peaks <- sort(c(peaks + 1, flat_peaks + 1))
  if(mag_order) return(peaks[order(x[peaks])])
  return(peaks)
}

#------------------------------------------------------------------
GetSegs <- function(peaks, x, pars, peak=NA){
  # identifies valid increasing-decreasing segments in x subject to the parameters in pars
  # returns a list of segments: c(start, peak, end). DON'T call directly w/ peak!=NA
  # NOTE: returned segments will not necessarily be in order, and may not completely partition x

  # ensure that peaks are in increasing order of x's magnitude
  tmp_peaks <- peaks[order(x[peaks])] # so we only have to sort once if they're in the wrong order
  if(!identical(tmp_peaks, peaks)) peaks <- tmp_peaks

  # if no peak is specified, we start at the beginnning
  if(is.na(peak)) peak <- peaks[1]

  # get the next largest peak; will be NA if this peak is the highest (last one to do)
  next_highest_peak <- peaks[which(peaks == peak) + 1]

  # check if we're doing C5-style relative amplitude and peak identification
  if(!is.na(pars$rel_amp_frac) & !is.na(pars$rel_peak_frac)){
    global_max <- max(x, na.rm=T)
    seg_thresh <- (global_max - min(x, na.rm=T)) * pars$rel_amp_frac
    peak_thresh <- global_max * pars$rel_peak_frac
  }else{
    seg_thresh <- pars$min_seg_amplitude
    peak_thresh <- 0
  }

	# checks if the period preceding the peak covers enough amplitude
	# search before the peak up to the maximum of: previous peak, the head of x, or the peak - max_increase_length
  previous_peaks <- peaks[peaks - peak < 0]
  previous_peak <- NA
  if(length(previous_peaks) > 0) previous_peak <- max(previous_peaks)
  search_start <- max(1, peak - pars$max_increase_length, previous_peak, na.rm=T)
  search_end <- peak
	# get the index of the closest minimum value within the search window
  # NOTE: should maybe retrieve the troughs here with FindPeaks(-x) instead
  # in the event of repeated minimum values, we take the closest one here
  inc_min_ind <- max(which(x[search_start:search_end] == min(x[search_start:search_end], na.rm=T)) + search_start - 1, na.rm=T)
  seg_amp <- x[peak] - x[inc_min_ind] # get the increasing segment amplitude
  # if(seg_amp > pars$min_seg_amplitude){
  if((seg_amp >= seg_thresh) & (x[peak] >= peak_thresh)){
    # check for a valid decreasing segment
    next_peaks <- peaks[peaks - peak > 0]
    next_peak <- NA
    if(length(next_peaks) > 0) next_peak <- min(next_peaks)
		# search after the peak up to the minimum of: next peak, the tail of x, or the max_decrease_length
    search_start <- peak
    search_end <- min(length(x), peak + pars$max_decrease_length, next_peak, na.rm=T)
		# get the index of the closest minimum value within the search window
    # NOTE: see above note about finding troughs instead
    dec_min_ind <- min(which(x[search_start:search_end] == min(x[search_start:search_end], na.rm=T)) + search_start - 1, na.rm=T)
    seg_amp <- x[peak] - x[dec_min_ind] # get the decreasing segment amplitude
    # if(seg_amp > pars$min_seg_amplitude){
    if(seg_amp >= seg_thresh){
      # we found a valid segment, store it as a list with a single vector: c(start, peak, end)
      tmp_seg <- list(c(inc_min_ind, peak, dec_min_ind))
      # if this isn't the last peak, then call CheckSegRec again w/ next highest peak
      if(!is.na(next_highest_peak)){
        return(c(tmp_seg, GetSegs(peaks, x, pars, peak=next_highest_peak)))
      }else{
				# that was the last peak, and it was valid
        return(tmp_seg) # covers the case where there's only one valid peak
      }
    }else{
      # increase was valid, but decrease was not
      peaks <- peaks[-which(peaks == peak)] # remove peak from peaks list
      # if this isn't the last peak, then call CheckSegRec again w/ next highest peak
      if(!is.na(next_highest_peak)){
        return(GetSegs(peaks, x, pars, peak=next_highest_peak))
      }else{
				# that was the last peak, and it was invalid
        return(NULL)
      }
    }
  }else{
    # increase segment not valid
    peaks <- peaks[-which(peaks == peak)] # remove peak from peaks list
    # if this isn't the last peak, then call CheckSegRec again w/ next highest peak
    if(!is.na(next_highest_peak)){
      return(GetSegs(peaks, x, pars, peak=next_highest_peak))
    }else{
			# that was the last peak, and it was invalid
      return(NULL)
    }
  }
}

#------------------------------------------------------------------
# GetSegs <- function(peaks, x, pars, peak=NA){
#   # identifies valid increasing-decreasing segments in x subject to the parameters in pars
#   # returns a list of segments: c(start, peak, end). DON'T call directly w/ peak!=NA
#   # NOTE: returned segments will not necessarily be in order, and may not completely partition x
#
#   # ensure that peaks are in increasing order of x's magnitude
#   tmp_peaks <- peaks[order(x[peaks])] # so we only have to sort once if they're in the wrong order
#   if(!identical(tmp_peaks, peaks)) peaks <- tmp_peaks
#
#   # if no peak is specified, we start at the beginnning
#   if(is.na(peak)) peak <- peaks[1]
#
#   # get the next largest peak; will be NA if this peak is the highest (last one to do)
#   next_highest_peak <- peaks[which(peaks == peak) + 1]
#
# 	# checks if the period preceding the peak covers enough amplitude
# 	# search before the peak up to the maximum of: previous peak, the head of x, or the max_increase_length
#   previous_peaks <- peaks[peaks - peak < 0]
#   previous_peak <- NA
#   if(length(previous_peaks) > 0) previous_peak <- max(previous_peaks)
#   search_start <- max(1, peak - pars$max_increase_length, previous_peak, na.rm=T)
#   search_end <- peak
# 	# get the index of the closest minimum value within the search window
#   # NOTE: should maybe retrieve the troughs here with FindPeaks(-x) instead
#   # in the event of repeated minimum values, we take the closest one here
#   inc_min_ind <- max(which(x[search_start:search_end] == min(x[search_start:search_end], na.rm=T)) + search_start - 1, na.rm=T)
#   seg_amp <- x[peak] - x[inc_min_ind] # get the increasing segment amplitude
#   if(seg_amp > pars$min_seg_amplitude){
#     # check for a valid decreasing segment
#     next_peaks <- peaks[peaks - peak > 0]
#     next_peak <- NA
#     if(length(next_peaks) > 0) next_peak <- min(next_peaks)
# 		# search after the peak up to the minimum of: next peak, the tail of x, or the max_decrease_length
#     search_start <- peak
#     search_end <- min(length(x), peak + pars$max_decrease_length, next_peak, na.rm=T)
# 		# get the index of the closest minimum value within the search window
#     # NOTE: see above note about finding troughs instead
#     dec_min_ind <- min(which(x[search_start:search_end] == min(x[search_start:search_end], na.rm=T)) + search_start - 1, na.rm=T)
#     seg_amp <- x[peak] - x[dec_min_ind] # get the decreasing segment amplitude
#     if(seg_amp > pars$min_seg_amplitude){
#       # we found a valid segment, store it as a list with a single vector: c(start, peak, end)
#       tmp_seg <- list(c(inc_min_ind, peak, dec_min_ind))
#       # if this isn't the last peak, then call CheckSegRec again w/ next highest peak
#       if(!is.na(next_highest_peak)){
#         return(c(tmp_seg, GetSegs(peaks, x, pars, peak=next_highest_peak)))
#       }else{
# 				# that was the last peak, and it was valid
#         return(tmp_seg) # covers the case where there's only one valid peak
#       }
#     }else{
#       # increase was valid, but decrease was not
#       peaks <- peaks[-which(peaks == peak)] # remove peak from peaks list
#       # if this isn't the last peak, then call CheckSegRec again w/ next highest peak
#       if(!is.na(next_highest_peak)){
#         return(GetSegs(peaks, x, pars, peak=next_highest_peak))
#       }else{
# 				# that was the last peak, and it was invalid
#         return(NULL)
#       }
#     }
#   }else{
#     # increase segment not valid
#     peaks <- peaks[-which(peaks == peak)] # remove peak from peaks list
#     # if this isn't the last peak, then call CheckSegRec again w/ next highest peak
#     if(!is.na(next_highest_peak)){
#       return(GetSegs(peaks, x, pars, peak=next_highest_peak))
#     }else{
# 			# that was the last peak, and it was invalid
#       return(NULL)
#     }
#   }
# }


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Function for retrieving vegetation cycle metrics (phenodates, QA, etc.)
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#---------------------------------------------------------------------
SegOverlapsPeriod <- function(seg, dates, period_start, period_end){
	# Checks if a single seg overlaps the period interest at all, if not then we don't want to bother processing
	if(dates[seg[3]] < period_start | dates[seg[1]] > period_end){
		return(FALSE)
	}else{
		return(TRUE)
	}
}

#---------------------------------------------------------------------
PeakInPeriod <- function(seg_met, start_date, end_date){
	# Checks if the peak of the segment with metrics "segmet" are between "start_date" and "end_date"
	if(seg_met$peak >= start_date & seg_met$peak <= end_date){
		return(TRUE)
	}else{
		return(FALSE)
	}
}

#------------------------------------------------------------------
GetThresh <- function(thresh_value, x, first_greater=T, gup=T){
	# returns the index of the first/last value  of x that is greater/less than the value of thresh.
	# If gup is False (greendown) then it returns the first/last value of x that is less/greater than
	# the value of thresh. first/last and greater/less determined by first_greater
	# NOTE: if thresh is 1 or 0, rounding error can be a problem. Now we round the threshold and each
	# of the evi values to 6 decimal places to compensate

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
R2 <- function(x, resids){
	# calculate coefficient of determination from x (evi) and residuals

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

	#--------------------------------------
	# new snowflag values:
	# 0 - observed value, no snow
	# 1 - interpolated value, no snow. Is used to initialize array so values between stride are going to be 1 even if they are snow.
	# 2 - truncated to minimum value, no snow
	# 3 - observed value, snow
	# 4 - interpolated value, snow
	# 5 - interpolated value off-stride, snow

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

#---------------------------------------------------------------------
SetReturnValues <- function(annual_pheno_metrics, seg_met, cycle=1, num_cycles=NA, fill_code=NA){
	# sets return object annual_pheno_metrics values for cycle 1 or 2 using the SegMet return object seg_met
	# and num_cycles and fill_code arguments

	######################
	# NOTE: here is where we collapse quality scores from continuous to 0-3, calculate
	# overall quality, and do the bitpacking
	# quals order: gup, midgup, mat, peak, sen, midsen, dor, overall:
	# quals <- c(1, 0, 3, 2, 2, 0, 1, 2) # example
	# QualBit <- function(x) as.integer(intToBits(x))[1:2]
	# QualPack <- function(quals) sum(unlist(lapply(quals, QualBit)) * (2^(0:15)))


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


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Highest level processing function
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#---------------------------------------------------------------------
AnnualPhenologyC6 <- function(x, dates, pheno_pars, pheno_period_start, pheno_period_end, plot=F){
	# processes MCD12Q2 phenology for a single pixel's time series
	# Up to two cycles are returned (but total # of peaks in period is too)
	# if there were >2 segs in the period, then the two with the highest magnitude
	# peaks are returned, but in chronological order

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
	# valid_peaks <- try(FilterPeaks(evi, FindPotentialPeaks_minmax(evi), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance), silent=T)
	valid_peaks <- try(FindPeaks(evi), silent=T)
	# if(inherits(valid_peaks, 'try-error') | is.na(valid_peaks)){
	if(inherits(valid_peaks, 'try-error') | (length(valid_peaks) == 0)){
		annual_pheno_metrics$fill_code <- 1 # no valid peaks
		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
		return(annual_pheno_metrics)
	}

	# find full segments
	full_segs <- try(GetSegs(valid_peaks, evi, pheno_pars), silent=T)
	if(inherits(full_segs, 'try-error') | is.null(full_segs)){
		annual_pheno_metrics$fill_code <- 2 # no valid segments
		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
		return(annual_pheno_metrics)
	}

	# eliminate segs that end before, or start after the period of interest
	seg_overlaps <- unlist(lapply(full_segs, SegOverlapsPeriod, dates, pheno_period_start, pheno_period_end))
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
	is_in_period <- unlist(lapply(seg_metrics, PeakInPeriod, start_date=pheno_period_start, end_date=pheno_period_end))
	segs_in_period <- seg_metrics[is_in_period]

	# assign the seg metric values to the output return value
	num_cycles <- length(segs_in_period) # total number of valid peaks within the period of interest
	if(num_cycles == 0){
		# no segments in period of interest, return default annual pheno metrics
		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
		return(annual_pheno_metrics)
	}else{
		# so, there is at least one valid segment.
		segs_in_period <- rev(segs_in_period) # reverse order so highest magnitude segs are first
		segs_in_period <- segs_in_period[order(unlist(segs_in_period)[which(names(unlist(segs_in_period)) == "peak")])] # put into chronological order

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

# #---------------------------------------------------------------------
# AnnualPhenologyC6_old <- function(x, dates, pheno_pars, pheno_period_start, pheno_period_end, plot=F){
# 	# processes MCD12Q2 phenology for a single pixel's time series
#
# 	# get a template return value
# 	annual_pheno_metrics <- PhenoReturnValue()
#
# 	# split the data up: first third are smoothed nbar-evi2, second third are residuals, and last third are snowflags
# 	# also scale, and fill for NA
# 	x[x==pheno_pars$nbar_NA_value] <- NA
# 	data_length <- length(x) / 3
# 	evi <- x[1:data_length]
# 	evi <- evi / pheno_pars$nbar_scale_factor
# 	resids <- x[(data_length + 1):(2 * data_length)]
# 	resids <- resids / pheno_pars$nbar_scale_factor
# 	snowflags <- x[(2 * data_length + 1):(3 * data_length)]
#
# 	# find valid peaks in the time series
# 	valid_peaks <- try(FilterPeaks(evi, FindPotentialPeaks_minmax(evi), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance), silent=T)
# 	if(inherits(valid_peaks, 'try-error') | is.na(valid_peaks)){
# 		annual_pheno_metrics$fill_code <- 1 # no valid peaks
# 		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
# 		return(annual_pheno_metrics)
# 	}
#
# 	# find full segments
# 	full_segs <- try(GetFullSegs(evi, valid_peaks, max_seg_length=pheno_pars$max_seg_length, min_seg_amplitude=pheno_pars$min_seg_amplitude, agg_amp_frac=pheno_pars$agg_amp_frac), silent=T)
# 	if(inherits(full_segs, 'try-error') | is.na(full_segs)){
# 		annual_pheno_metrics$fill_code <- 2 # no valid segments
# 		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
# 		return(annual_pheno_metrics)
# 	}
#
# 	# eliminate segs that end before, or start after the period of interest
# 	seg_overlaps <- unlist(lapply(full_segs, SegOverlapsPeriod, dates, pheno_period_start, pheno_period_end))
# 	if(all(!seg_overlaps)){
# 		# no segs within period of interest
# 		annual_pheno_metrics$fill_code <- 2 # no valid segments
# 		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
# 		return(annual_pheno_metrics)
# 	}
# 	full_segs <- full_segs[seg_overlaps]
#
# 	# get the segment metrics
# 	seg_metrics <- lapply(full_segs, SegMet, x=x, dates=dates, pheno_pars=pheno_pars)
#
# 	# limit to segments where the peak is in the period of interest
# 	is_in_period <- unlist(lapply(seg_metrics, PeakInPeriod, start_date=pheno_period_start, end_date=pheno_period_end))
# 	segs_in_period <- seg_metrics[is_in_period]
#
# 	# assign the seg metric values to the output return value
# 	num_cycles <- length(segs_in_period) # total number of valid peaks within the period of interest
# 	if(num_cycles == 0){
# 		# no segments in period of interest, return default annual pheno metrics
# 		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
# 		return(annual_pheno_metrics)
# 	}else{
# 		# set first cycle return metrics
# 		cycle1_seg_met <- segs_in_period[[1]]
# 		annual_pheno_metrics <- SetReturnValues(annual_pheno_metrics, cycle1_seg_met, cycle=1, num_cycles=num_cycles, fill_code=0)
# 		if(num_cycles > 1){
# 			# set the second cycle
# 			cycle2_seg_met <- segs_in_period[[2]]
# 			annual_pheno_metrics <- SetReturnValues(annual_pheno_metrics, cycle2_seg_met, cycle=2)
# 		}
# 		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
# 		return(annual_pheno_metrics)
# 	}
# }

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Functions for visualizing the data
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#---------------------------------------------------------------------
PlotSeries <- function(x, dates, NA_value=32767, scale_value=1e4, plot_legend=T, seg_metrics=NA, ylim=NULL, grid=T){
	# function for plotting an individual time series row w/ or w/o seg metrics
  mycols <- c("#636363", "#DE2D26", "#3182BD", "#31A354", "#756BB1", "#E6550D")
  x[x==NA_value] <- NA
  data_length <- length(x) / 3
  filtered <- x[1:data_length] / scale_value
  resids <- x[(data_length + 1):(2 * data_length)] / scale_value
  snowflags <- x[(2 * data_length + 1):(3 * data_length)]
  cols <- mycols[snowflags + 1]
  pchs <- snowflags + 1

  # start plotting
	if(!is.null(ylim)) ylim <- ylim
  plot(c(dates, dates), c(filtered, filtered + resids), type="n", xlab="", ylab="EVI2", ylim=ylim)
	if(grid){
		grid(ny=NULL, nx=NA)
		abline(v=c(as.Date(paste(unique(strftime(dates, format="%Y")), "-1-1", sep="")), as.Date(paste(unique(strftime(dates, format="%Y")), "-7-1", sep=""))), lty="dotted", col="lightgray")
	}

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
		# points(dates[is.na(resids) & (snowflags == 1 | snowflags == 5)], filtered[is.na(resids) & (snowflags == 1 | snowflags == 5)], type="p", pch=pchs[is.na(resids) & (snowflags == 1 | snowflags == 5)], cex=0.75, col=cols[is.na(resids) & (snowflags == 1 | snowflags == 5)])
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
		if(plot_legend) legend("topleft", legend=c("Spline-smoothed", "Flag=0", "Flag=1", "Flag=2", "Flag=3", "Flag=4"), pch=c(NA, 1:5), col=c("darkgrey", mycols[1:5]), lwd=c(1, rep(NA, 5)), bg="white")
  }else{
    # plot(c(dates, dates), c(filtered, filtered + resids), type="n", xlab="", ylab="EVI2")
    points(dates, filtered + resids, type="p", pch=pchs, cex=0.75, col=cols)
    points(dates, filtered, type="l", lwd=2, col="darkgrey")
    if(plot_legend) legend("topleft", legend=c("Spline-smoothed", "Flag=0", "Flag=1", "Flag=2", "Flag=3", "Flag=4"), pch=c(NA, 1:5), col=c("darkgrey", mycols[1:5]), lwd=c(1, rep(NA, 5)), bg="white")
  }
}

#---------------------------------------------------------------------
CheckForTile <- function(tile, year, data_dir, prefix="MCD12Q2C6"){
  # check for the existence of a MCC12Q2C6 output tile year in data_dir and return it

  in_files <- dir(data_dir, pattern=paste(prefix, ".*", tile, ".*", year, "$", sep=""), full=T)
  return(in_files[1])
}

#---------------------------------------------------------------------
BuildVRT <- function(file_list, out_file, band=NULL, vrtnodata=0){
  # Builds a VRT mosaic of files in file_list, possibly w/ band specification

  if(!is.null(band)){
    sys_cmd <- paste("gdalbuildvrt", paste("-vrtnodata", vrtnodata), paste("-b", band), out_file, paste(file_list, collapse=" "))
  }else{
    sys_cmd <- paste("gdalbuildvrt", paste("-vrtnodata", vrtnodata), out_file, paste(file_list, collapse=" "))
  }
  system(sys_cmd)
}

#---------------------------------------------------------------------
ContinentWGS <- function(in_file, out_file, subset=NA){
	# gdalwarp -dstnodata 0 -s_srs "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m" -t_srs "+proj=latlong +datum=WGS84" -tr 0.025 0.025 -te -170 10 -50 75 NAmerica_midgup_cycle1.vrt NAmerica_midgup_cycle1_wgs.tif
}

#---------------------------------------------------------------------
PlotDOY <- function(r, year, cutoffs=c(0, 365), pdf_out=NULL, plot_height=12, plotNEW=T, garish=F, continents=NA,...){
  # plots an MCD12Q2 date layer as doy; linear stretch between cutoff values

  doy_offset <- as.numeric(as.Date(paste(year, "-1-1", sep="")))
	if(garish){
		pal <- colorRampPalette(rev(rainbow(12)))
	}else{
		pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
	}
  breaks <- c(-1e6, seq(cutoffs[1], cutoffs[2]), 1e6)
  breaks <- unique(breaks)

  if(!is.null(pdf_out)){
    # create PDF output
    pdf(file=pdf_out, h=plot_height, w=1.2 * (plot_height * (dim(r)[2] / dim(r)[1])))
  }else{
    # print to screen w/ x11 device
		plot_height
		if(plotNEW) x11(h=plot_height, w=1.2 * (plot_height * (dim(r)[2] / dim(r)[1])))
  }

  par(col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white", bg="black", mar=c(3.5, 3.5, 1, 1), oma=rep(0, 4))
  # plot(r - doy_offset, breaks=breaks, col=pal(length(breaks) - 1), legend=F, colNA="black", ...)
	if(!is.na(continents)){
		plot(continents, border="NA", col=rgb(0.3, 0.3, 0.3), xlim=extent(midgup)[1:2], ylim=extent(midgup)[3:4])
		plot(r - doy_offset, breaks=breaks, col=pal(length(breaks) - 1), legend=F, xaxt="n", yaxt="n", box=FALSE, bty="n", add=T, ...)
	}else{
		plot(r - doy_offset, breaks=breaks, col=pal(length(breaks) - 1), legend=F, colNA="black", xaxt="n", yaxt="n", box=FALSE, bty="n", ...)
	}


  tmp_r <- raster(matrix(seq(breaks[2], breaks[length(breaks) - 1], len=length(breaks) - 1)))
  legend_at <- seq(breaks[2], breaks[length(breaks) - 1], len=9)
  legend_labels_tmp <- round(legend_at, 0)
  legend_labels <- c(paste("<", legend_labels_tmp[1]), legend_labels_tmp[2:(length(legend_labels_tmp) - 1)], paste(">", legend_labels_tmp[length(legend_labels_tmp)]))

  plot(
    tmp_r,
    legend.only=T,
    col=pal(length(breaks) - 1),
    axis.args=list(at=legend_at, labels=legend_labels)
  )

  if(!is.null(pdf_out)) dev.off()
}

#---------------------------------------------------------------------
PlotNumCycles <- function(r, pdf_out=NULL, plot_height=12, plotNEW=T, garish=F, continents=NA,...){
  # plots an MCD12Q2 date layer as doy; linear stretch between cutoff values

  # doy_offset <- as.numeric(as.Date(paste(year, "-1-1", sep="")))
	if(garish){
		pal <- colorRampPalette(rev(rainbow(12)))
	}else{
		# pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
		cols <- brewer.pal(4, "Greens")[2:4]
	}
	breaks <- c(0, 1.5, 2.5, 1e3)
  # breaks <- c(-1e6, seq(cutoffs[1], cutoffs[2]), 1e6)
  # breaks <- unique(breaks)

  if(!is.null(pdf_out)){
    # create PDF output
    pdf(file=pdf_out, h=plot_height, w=1.2 * (plot_height * (dim(r)[2] / dim(r)[1])))
  }else{
    # print to screen w/ x11 device
		# plot_height
		if(plotNEW) x11(h=plot_height, w=1.2 * (plot_height * (dim(r)[2] / dim(r)[1])))
  }

  par(col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white", bg="black", mar=c(3.5, 3.5, 1, 1), oma=rep(0, 4))
  # plot(r - doy_offset, breaks=breaks, col=pal(length(breaks) - 1), legend=F, colNA="black", ...)
	if(!is.na(continents)){
		plot(continents, border="NA", col=rgb(0.3, 0.3, 0.3), xlim=extent(midgup)[1:2], ylim=extent(midgup)[3:4])
		# plot(r - doy_offset, breaks=breaks, col=pal(length(breaks) - 1), legend=F, xaxt="n", yaxt="n", box=FALSE, bty="n", add=T, ...)
		plot(r, breaks=breaks, col=cols, legend=F, xaxt="n", yaxt="n", box=FALSE, bty="n", add=T, ...)
	}else{
		# plot(r - doy_offset, breaks=breaks, col=pal(length(breaks) - 1), legend=F, colNA="black", xaxt="n", yaxt="n", box=FALSE, bty="n", ...)
		plot(r, breaks=breaks, col=cols, legend=F, xaxt="n", yaxt="n", box=FALSE, bty="n", add=F, ...)
	}


  # tmp_r <- raster(matrix(seq(breaks[2], breaks[length(breaks) - 1], len=length(breaks) - 1)))
  # legend_at <- seq(breaks[2], breaks[length(breaks) - 1], len=9)
  # legend_labels_tmp <- round(legend_at, 0)
  # legend_labels <- c(paste("<", legend_labels_tmp[1]), legend_labels_tmp[2:(length(legend_labels_tmp) - 1)], paste(">", legend_labels_tmp[length(legend_labels_tmp)]))
	#
  # plot(
  #   tmp_r,
  #   legend.only=T,
  #   col=pal(length(breaks) - 1),
  #   axis.args=list(at=legend_at, labels=legend_labels)
  # )

  if(!is.null(pdf_out)) dev.off()
}

#---------------------------------------------------------------------
PlotStretch <- function(r, year, cutoffs=c(0, 1), pdf_out=NULL, plot_height=12, plotNEW=T, garish=F, ...){
  # plots an MCD12Q2 date layer as doy; linear stretch between cutoff values

  # doy_offset <- as.numeric(as.Date(paste(year, "-1-1", sep="")))
	if(garish){
		pal <- colorRampPalette(rev(rainbow(12)))
	}else{
		pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
	}
  breaks <- c(-1e6, seq(cutoffs[1], cutoffs[2]), 1e6)
  breaks <- unique(breaks)

  if(!is.null(pdf_out)){
    # create PDF output
    pdf(file=pdf_out, h=plot_height, w=1.2 * (plot_height * (dim(r)[2] / dim(r)[1])))
  }else{
    # print to screen w/ x11 device
		plot_height
		if(plotNEW) x11(h=plot_height, w=1.2 * (plot_height * (dim(r)[2] / dim(r)[1])))
  }

  par(col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white", bg="black", mar=c(3.5, 3.5, 1, 1), oma=rep(0, 4))
  # plot(r - doy_offset, breaks=breaks, col=pal(length(breaks) - 1), legend=F, colNA="black", ...)
  plot(r, breaks=breaks, col=pal(length(breaks) - 1), legend=F, colNA="black", xaxt="n", yaxt="n", box=FALSE, bty="n", ...)

  tmp_r <- raster(matrix(seq(breaks[2], breaks[length(breaks) - 1], len=length(breaks) - 1)))
  legend_at <- seq(breaks[2], breaks[length(breaks) - 1], len=9)
  legend_labels_tmp <- round(legend_at, 0)
  legend_labels <- c(paste("<", legend_labels_tmp[1]), legend_labels_tmp[2:(length(legend_labels_tmp) - 1)], paste(">", legend_labels_tmp[length(legend_labels_tmp)]))

  plot(
    tmp_r,
    legend.only=T,
    col=pal(length(breaks) - 1),
    axis.args=list(at=legend_at, labels=legend_labels)
  )

  if(!is.null(pdf_out)) dev.off()
}

#---------------------------------------------------------------------
# PlotStretch <- function(r, thresh_percent=c(2, 2), thresh_q=NA, na.rm=T, legend_ticks=5, dates=T, legend_round=3, ...){
# 	# plot with a linear % cutoff stretch and proper legend
#
# 	library(RColorBrewer)
# 	if(is.na(thresh_q)){
# 		q <- quantile(r, c(0, thresh_percent[1] / 100, 1 - (thresh_percent[2] / 100), 1), na.rm=na.rm)
# 	}else{
# 		if(length(thresh_q) != 4){
# 			print("If specified, thresh_q must have length 4: (min, low_cut, high_cut, max)")
# 			return(NA)
# 		}else{
# 			q <- thresh_q
# 		}
# 	}
#
# 	if(dates){
# 		breaks <- unique(round(c(q[1], seq(q[2], q[3], len=254), q[4])))
# 	}else{
# 		breaks <- unique(c(q[1], seq(q[2], q[3], len=254), q[4]))
# 	}
#
# 	pal <-  colorRampPalette(brewer.pal(11, "Spectral"))
# 	plot(r, breaks=breaks, col=pal(length(breaks) - 1), legend=F, ...)
#
# 	# make the legend
# 	# tmp_r <- raster(matrix(q[2]:q[3]))
# 	tmp_r <- raster(matrix(seq(q[2], q[3], len=255)))
# 	# legend_at <- round(seq(q[2], q[3], len=legend_ticks))
# 	legend_at <- seq(q[2], q[3], len=legend_ticks)
# 	# legend_at <- round(seq(q[2], q[3], len=legend_ticks), legend_round)
#
# 	if(dates){
# 		legend_at <- round(legend_at)
# 		legend_labels <- c(
# 							paste("<", as.Date(legend_at[1], origin="1970-1-1")),
# 							as.character(as.Date(legend_at[2:(length(legend_at) - 1)], origin="1970-1-1")),
# 							paste("<", as.Date(legend_at[length(legend_at)], origin="1970-1-1"))
# 						)
# 	}else{
# 		legend_labels_tmp <- round(legend_at, legend_round)
# 		legend_labels <- c(paste("<", legend_labels_tmp[1]), legend_labels_tmp[2:(length(legend_labels_tmp) - 1)], paste(">", legend_labels_tmp[length(legend_labels_tmp)]))
# 	}
#
# 	plot(
# 		tmp_r,
# 		legend.only=T,
# 		col=pal(length(breaks) - 1),
# 		axis.args=list(at=legend_at, labels=legend_labels)
# 	)
# }


#---------------------------------------------------------------------
PlotPhenology <- function(x, dates, pheno_pars=DefaultPhenoParameters(), NA_value=32767, scale_value=1e4, ylim=NULL, plot_legend=T, grid=T, plot_dates=T, blackBG=F){
# PlotSeries <- function(x, dates, NA_value=32767, scale_value=1e4, plot_legend=T, seg_metrics=NA, ylim=NULL, grid=T){
	# function for plotting an individual time series row w/ or w/o seg metrics
  mycols <- c("#636363", "#DE2D26", "#3182BD", "#31A354", "#756BB1", "#E6550D")
  # gup_poly_color <- "lightgreen"
	gup_poly_color <- rgb(144,238,144,155,max=255)
  # gdown_poly_color <- "pink"
	gdown_poly_color <- rgb(255,192,203,155,max=255)
  poly_density <- 25
  x[x==NA_value] <- NA
  data_length <- length(x) / 3
  filtered <- x[1:data_length] / scale_value
  resids <- x[(data_length + 1):(2 * data_length)] / scale_value
  snowflags <- x[(2 * data_length + 1):(3 * data_length)]
  cols <- mycols[snowflags + 1]
  pchs <- snowflags + 1

  # # check if pheno_period_start/end are provided, if not use the whole time series
  # if(is.na(pheno_period_start) | is.na(pheno_period_end)){
  #   pheno_period_start <- min(dates, na.rm=T)
  #   pheno_period_end <- max(dates, na.rm=T)
  # }

  # retrieve phenology
  # pheno_values <- AnnualPhenologyC6(x, dates, pheno_pars, pheno_period_start, pheno_period_end)
  valid_peaks <- try(FindPeaks(filtered), silent=T)
  full_segs <- try(GetSegs(valid_peaks, filtered, pheno_pars), silent=T)
  seg_metrics <- lapply(full_segs, SegMet, x=x, dates=dates, pheno_pars=pheno_pars)

  # start plotting
	if(!is.null(ylim)) ylim <- ylim
	if(blackBG){
		par(col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white", bg="black", mar=c(4, 4, 2, 2), oma=rep(1, 4))
	}
  plot(c(dates, dates), c(filtered, filtered + resids), type="n", xlab="", ylab="EVI2", ylim=ylim)
	if(grid){
		grid(nx=NA, ny=NULL) # plot horizontal grid lines at tick marks
		abline(v=as.Date(paste(unique(as.numeric(strftime(dates, format="%Y"))), "-1-1", sep="")), col="lightgray", lty="dotted")
		abline(v=as.Date(paste(unique(as.numeric(strftime(dates, format="%Y"))), "-7-1", sep="")), col="lightgray", lty="dotted")
	}

  # more extensive plotting if seg_metrics are available
  if(length(seg_metrics) > 0){
    for(seg_metric in seg_metrics){
      # plot the segment data
      # gup_x_tmp <- as.Date(seg_metric$peak - seg_metric$length_gup, origin="1970-1-1"):as.Date(seg_metric$peak, origin="1970-1-1")
			gup_x_tmp <- dates[which(dates == as.Date(seg_metric$peak - seg_metric$length_gup, origin="1970-1-1")):which(dates == as.Date(seg_metric$peak, origin="1970-1-1"))]
      gup_x <- c(gup_x_tmp, rev(gup_x_tmp))
      gup_y_tmp <- filtered[which(dates == as.Date(seg_metric$peak - seg_metric$length_gup, origin="1970-1-1")):which(dates == as.Date(seg_metric$peak, origin="1970-1-1"))]
      gup_y <- c(gup_y_tmp, rep(par()$usr[3], length(gup_x_tmp)))
      # polygon(x=gup_x, y=gup_y, border=NA, col=gup_poly_color, density=poly_density, angle=45)
			polygon(x=gup_x, y=gup_y, border=NA, col=gup_poly_color)

      # gdown_x_tmp <- as.Date(seg_metric$peak, origin="1970-1-1"):as.Date(seg_metric$peak + seg_metric$length_gdown, origin="1970-1-1")
			gdown_x_tmp <- dates[which(dates == as.Date(seg_metric$peak, origin="1970-1-1")):which(dates == as.Date(seg_metric$peak + seg_metric$length_gdown, origin="1970-1-1"))]
      gdown_x <- c(gdown_x_tmp, rev(gdown_x_tmp))
      gdown_y_tmp <- filtered[which(dates == as.Date(seg_metric$peak, origin="1970-1-1")):which(dates == as.Date(seg_metric$peak + seg_metric$length_gdown, origin="1970-1-1"))]
      gdown_y <- c(gdown_y_tmp, rep(par()$usr[3], length(gdown_x_tmp)))
      # polygon(x=gdown_x, y=gdown_y, border=NA, col=gdown_poly_color, density=poly_density, angle=-45)
			polygon(x=gdown_x, y=gdown_y, border=NA, col=gdown_poly_color)

      # this creates simple rectangle fills:
      # polygon(x=c(rep(as.Date(seg_metric$peak - seg_metric$length_gup, origin="1970-1-1"), 2), rep(as.Date(seg_metric$peak, origin="1970-1-1"), 2)), y=c(par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3]), border=NA, col=rgb(0.87, 0.87, 0.87))
      # polygon(x=c(rep(as.Date(seg_metric$peak, origin="1970-1-1"), 2), rep(as.Date(seg_metric$peak + seg_metric$length_gdown, origin="1970-1-1"), 2)), y=c(par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3]), border=NA, col=rgb(0.82, 0.82, 0.82))
    } # end seg loop

    # add the data/spline with fill/miss/snow info
    points(dates, filtered + resids, type="p", pch=pchs, cex=0.75, col=cols)
    points(dates, filtered, type="l", lwd=2, col="darkgrey")

		if(plot_dates){
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
		}


    # create a legend
		if(plot_legend) legend("topleft", legend=c("Spline-smoothed", "Flag=0", "Flag=1", "Flag=2", "Flag=3", "Flag=4"), pch=c(NA, 1:5), col=c("darkgrey", mycols[1:5]), lwd=c(1, rep(NA, 5)), bg="white", horiz=F)
  }else{
    # plot(c(dates, dates), c(filtered, filtered + resids), type="n", xlab="", ylab="EVI2")
    points(dates, filtered + resids, type="p", pch=pchs, cex=0.75, col=cols)
    points(dates, filtered, type="l", lwd=2, col="darkgrey")
    if(plot_legend) legend("topleft", legend=c("Spline-smoothed", "Flag=0", "Flag=1", "Flag=2", "Flag=3", "Flag=4"), pch=c(NA, 1:5), col=c("darkgrey", mycols[1:5]), lwd=c(1, rep(NA, 5)), bg="white", horiz=F)
  }
}
