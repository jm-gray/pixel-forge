library(raster)
library(rgdal)
library(RColorBrewer)

source("~/Documents/pixel-forge/MCD12Q2C6/MCD12Q2C6_AnnualPhenologyFunctions.R")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
GetCellTimeSeries <- function(cell, r, evi2_files, resid_files, snow_files){
  the_line <- rowFromCell(r, cell)
  the_col <- colFromCell(r, cell)
  tmp <- Get3YearDataChunk_tmp(evi2_files, resid_files, snow_files, year_to_do, the_line, 1)
  return(tmp[the_col,])
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
lc_dir <- "~/Desktop/MCD12Q2_analysis/LC/"
data_dir="~/Desktop/MCD12Q2_analysis/"
splined_nbar_data_dir <- file.path(data_dir, "MCD12I2")
splined_qa_data_dir <- file.path(data_dir, "MCD12I3")
splined_resids_data_dir <- file.path(data_dir, "MCD12I4")
phen_dir_v0 <- file.path(data_dir, "MCD12I6")
phen_dir_v1 <- file.path(data_dir, "MCD12I6_v1")

tile_to_do <- "h08v05"
year_to_do <- 2009
metric <- "Dormancy"
c6_dates <- as.Date(paste(rep((year_to_do - 1):(year_to_do + 1), each=365), rep(1:365, 3), sep="-"), format="%Y-%j")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Select some places where dormancy is delayed in v1
r_v0 <- raster(dir(phen_dir_v0, pattern=paste(metric, ".*", tile_to_do, ".*", year_to_do, "$", sep=""), full=T))
r_v1 <- raster(dir(phen_dir_v1, pattern=paste(metric, ".*", tile_to_do, ".*", year_to_do, "$", sep=""), full=T))
NAvalue(r_v0) <- 32767
NAvalue(r_v1) <- 32767

r_diff <- r_v0 - r_v1
v_diff <- values(r_diff)
not_zero <- which(v_diff != 0)
cells_to_do <- not_zero[c(1, 3, 9, 7, 12, 14, 15, 17, 26)] # hand selected
cell_lines <- rowFromCell(r_v0, cells_to_do)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# get the phenology data for v0
pheno_in_files <- dir(phen_dir_v0, pattern=paste("(Greenup|MidGreenup|Maturity|Peak|Senescence|MidGreendown|Dormancy).*", tile_to_do, "_[0-9]{4}$", sep=""), full=T)
pheno_years <- as.integer(gsub(".*([0-9]{4})$", "\\1", basename(pheno_in_files)))
pheno_in_files <- pheno_in_files[pheno_years %in% c(year_to_do - 1, year_to_do, year_to_do + 1)]
reorder_inds <- c(4:6, 13:15, 7:9, 16:18, 19:21, 10:12, 1:3)
pheno_s_v0 <- stack(pheno_in_files[reorder_inds])
NAvalue(pheno_s_v0) <- 32767
pheno_v_v0 <- values(pheno_s_v0)

# get the phenology data for v1
pheno_in_files <- dir(phen_dir_v1, pattern=paste("(Greenup|MidGreenup|Maturity|Peak|Senescence|MidGreendown|Dormancy).*", tile_to_do, "_[0-9]{4}$", sep=""), full=T)
pheno_years <- as.integer(gsub(".*([0-9]{4})$", "\\1", basename(pheno_in_files)))
pheno_in_files <- pheno_in_files[pheno_years %in% c(year_to_do - 1, year_to_do, year_to_do + 1)]
reorder_inds <- c(4:6, 13:15, 7:9, 16:18, 19:21, 10:12, 1:3)
pheno_s_v1 <- stack(pheno_in_files[reorder_inds])
NAvalue(pheno_s_v1) <- 32767
pheno_v_v1 <- values(pheno_s_v1)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# retrieve the splined input data
patt <- paste("MCD12", ".*", tile_to_do, ".*.bip$", sep="")
evi2_files <- dir(splined_nbar_data_dir, pattern=patt, full=T)
resid_files <- dir(splined_resids_data_dir, pattern=patt, full=T)
snow_files <- dir(splined_qa_data_dir, pattern=patt, full=T)

tile_ts <- t(sapply(cells_to_do, GetCellTimeSeries, r=r_v0, evi2_files=evi2_files, resid_files=resid_files, snow_files=snow_files))

mycols <- c(brewer.pal(5, "Reds")[4], brewer.pal(5, "Blues")[4], brewer.pal(5, "Greens")[4])
layout(matrix(1:nrow(tile_ts), ncol=1, byrow=T))
par(mar=c(1, 4, 1, 1))
for(i in 1:nrow(tile_ts)){
  # PlotSeries(tile_ts[i,], c6_dates, igbp_names)
  PlotSeries(tile_ts[i,], c6_dates)
  abline(v=as.Date(pheno_v_v0[cells_to_do[i],], origin=as.Date("1970-1-1")), lty=rep(1:2, 21), col=mycols[1], lwd=1.5)
  abline(v=as.Date(pheno_v_v1[cells_to_do[i],], origin=as.Date("1970-1-1")), lty=rep(1:2, 21), col=mycols[2])
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
which_series <- 9
pars <- DefaultPhenoParameters()
pars$min_seg_amplitude <- 0.15
pars$rel_amp_frac <- 0.35
pars$rel_peak_frac <- 0.67
# pars$rel_peak_frac <- 0.7
pheno_results_v0 <- AnnualPhenologyC6(as.numeric(tile_ts[which_series,]), c6_dates, pars, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)
pars <- DefaultPhenoParameters()
pars$min_seg_amplitude <- 0.1
pheno_results_v1 <- AnnualPhenologyC6(as.numeric(tile_ts[which_series,]), c6_dates, pars, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)
layer_names <- c("num_cycles", "evi_area_cycle1", "evi_amp_cycle1", "evi_min_cycle1", "ogi_cycle1", "midgup_cycle1", "mat_cycle1", "peak_cycle1", "sen_cycle1", "midgdown_cycle1", "dor_cycle1", "overall_qa_cycle1", "detailed_qa_cycle1",  "evi_area_cycle2", "evi_amp_cycle2", "evi_min_cycle2", "ogi_cycle2", "midgup_cycle2", "mat_cycle2", "peak_cycle2", "sen_cycle2", "midgdown_cycle2", "dor_cycle2", "overall_qa_cycle2", "detailed_qa_cycle2")

ind_pattern <- c(rep(NA, 4), 3, 9, 15, 21, 27, 33, 39, rep(NA, 5), 4, 10, 16, 22, 28, 34, 40, NA, NA)
data.frame(name=layer_names, v0_mine=pheno_results_v0, v0_giss=pheno_v_v0[cells_to_do[which_series],][ind_pattern], v1_mine=pheno_results_v1, v1_giss=pheno_v_v1[cells_to_do[which_series],][ind_pattern])


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Testing the peak-to-peak distance problem

x[x==pheno_pars$nbar_NA_value] <- NA
data_length <- length(x) / 3
evi <- x[1:data_length]
evi <- evi / pheno_pars$nbar_scale_factor
resids <- x[(data_length + 1):(2 * data_length)]
resids <- resids / pheno_pars$nbar_scale_factor
snowflags <- x[(2 * data_length + 1):(3 * data_length)]

which_series <- 3
pars <- DefaultPhenoParameters()
pars$min_seg_amplitude <- 0.1
pheno_no_filter <- AnnualPhenologyC6(as.numeric(tile_ts[which_series,]), c6_dates, pars, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)
pars$min_peak_to_peak_distance=90
pheno_filtered <- AnnualPhenologyC6_filter(as.numeric(tile_ts[which_series,]), c6_dates, pars, as.Date("2009-1-1"), as.Date("2009-12-31"), plot=T)
# data.frame(layer_names, no_filt=pheno_no_filter, filt=pheno_filtered)
layout(matrix(1:2, nrow=2))
par(mar=c(0, 4, 0, 2), oma=c(4, 0, 1, 1))
PlotSeries(tile_ts[which_series,], c6_dates)
abline(v=as.Date(pheno_no_filter[c(5:11, 17:23)], origin=as.Date("1970-1-1")), lty=c(rep(1,7), rep(2,7)), col=mycols[1], lwd=1.5)
PlotSeries(tile_ts[which_series,], c6_dates)
abline(v=as.Date(pheno_filtered[c(5:11, 17:23)], origin=as.Date("1970-1-1")), lty=c(rep(1,7), rep(2,7)), col=mycols[2])

FilterPeaks <- function(peaks, min_peak_to_peak_distance=NA){
  # this version takes advantage of the fact that the peaks are in magnitude order, so no need to check EVI
  peaks <- rev(peaks)
  i <- 1
  while(i < length(peaks)){
    if(!is.na(peaks[i])){
      peaks_too_close <- which(abs(peaks[i] - peaks) < min_peak_to_peak_distance)
      peaks_too_close <- peaks_too_close[peaks_too_close != i]
      peaks[peaks_too_close] <- NA
      print(peaks)
    }
    i <- i + 1
  }
  peaks <- peaks[!is.na(peaks)]
  if(length(peaks) == 0) return(NA)
  return(peaks)
}


# FilterPeaks <- function(x, potential_peaks, min_peak_to_peak_distance=NA){
#
#   potential_peaks <- sort(potential_peaks)
# 	# check that potential_peaks is not NA
# 	if(all(is.na(potential_peaks))) return(NA)
#
#   # return unmodified peaks if the min_peak_to_peak_distance is NA or 0
#   if(is.na(min_peak_to_peak_distance) | min_peak_to_peak_distance == 0) return(potential_peaks)
#
# 	# we loop through all the potential peaks, only moving on when we have satisfied
# 	# the minimum distance between peaks requirement
# 	i <- 1
#
# 	while(i < length(potential_peaks)){
#     # print(paste(i, potential_peaks[i]))
# 		# find the distance to next peak
# 		peak_to_peak_dist <- potential_peaks[i + 1] - potential_peaks[i]
#
# 		if(peak_to_peak_dist < min_peak_to_peak_distance){
# 			# eliminate the smaller of the two potential peaks
# 			if(x[potential_peaks[i + 1]] >= x[potential_peaks[i]]){
# 				potential_peaks <- potential_peaks[-1 * i]
# 			}else{
# 				potential_peaks <- potential_peaks[-1 * (i + 1)]
# 			}
# 		}else{
# 			# distance is fine, move on
# 			i <- i + 1
# 		}
# 		# if we've eliminated the last peak, return NA
# 		if(length(potential_peaks) == 0) return(NA)
# 	}
# 	return(potential_peaks)
# }

AnnualPhenologyC6_filter <- function(x, dates, pheno_pars, pheno_period_start, pheno_period_end, plot=F){
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
  valid_peaks <- FilterPeaks(valid_peaks, min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance)
	# if(inherits(valid_peaks, 'try-error') | is.na(valid_peaks)){
	if(inherits(valid_peaks, 'try-error') | (length(valid_peaks) == 0)){
		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
		return(annual_pheno_metrics)
	}

	# find full segments
	full_segs <- try(GetSegs(valid_peaks, evi, pheno_pars), silent=T)
	if(inherits(full_segs, 'try-error') | is.null(full_segs)){
		annual_pheno_metrics <- ScaleToIntegerAndSetNA(annual_pheno_metrics, pheno_pars$out_float_scale, pheno_pars$out_NA_value)
		return(annual_pheno_metrics)
	}

	# eliminate segs that end before, or start after the period of interest
	seg_overlaps <- unlist(lapply(full_segs, SegOverlapsPeriod, dates, pheno_period_start, pheno_period_end))
	if(all(!seg_overlaps)){
		# no segs within period of interest
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
