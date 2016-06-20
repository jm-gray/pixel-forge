#---------------------------------------------------------------------
# Prelims
rm(list=ls())
library(plyr)
library(raster)
library(rgdal)
library(doParallel)
library(foreach)

#---------------------------------------------------------------------
# Functions

#---------------------------------------------------------------------
# these functions return an SDS filename string so the raster package can handle the HDF inputs
Get_SDS_name_band1 <- function(mcd43a4_file_path, fill_tile){
	if(is.na(mcd43a4_file_path)){
		# Fill tile
		data_set_name <- fill_tile
	}else{
		data_set_name <- paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":MOD_Grid_BRDF:Nadir_Reflectance_Band1", sep = "")
	}

	return(data_set_name)
}

Get_SDS_name_band2 <- function(mcd43a4_file_path, fill_tile){
	if(is.na(mcd43a4_file_path)){
		# Fill tile
		data_set_name <- fill_tile
	}else{
		data_set_name <- paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":MOD_Grid_BRDF:Nadir_Reflectance_Band2", sep = "")
	}

	return(data_set_name)
}

Get_SDS_name_band3 <- function(mcd43a4_file_path, fill_tile){
	if(is.na(mcd43a4_file_path)){
		# Fill tile
		data_set_name <- fill_tile
	}else{
		data_set_name <- paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":MOD_Grid_BRDF:Nadir_Reflectance_Band3", sep = "")
	}

	return(data_set_name)
}

Get_SDS_name_band4 <- function(mcd43a4_file_path, fill_tile){
	if(is.na(mcd43a4_file_path)){
		# Fill tile
		data_set_name <- fill_tile
	}else{
		data_set_name <- paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":MOD_Grid_BRDF:Nadir_Reflectance_Band4", sep = "")
	}

	return(data_set_name)
}

Get_SDS_name_band6 <- function(mcd43a4_file_path, fill_tile){
	if(is.na(mcd43a4_file_path)){
		# Fill tile
		data_set_name <- fill_tile
	}else{
		data_set_name <- paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":MOD_Grid_BRDF:Nadir_Reflectance_Band6", sep = "")
	}

	return(data_set_name)
}

Get_SDS_name_NBAR_QA <- function(mcd43a2_file_path, fill_tile){
	if(is.na(mcd43a2_file_path)){
		# Fill tile
		data_set_name <- fill_tile
	}else{
		data_set_name <- paste("HDF4_EOS:EOS_GRID:\"", mcd43a2_file_path, "\":MOD_Grid_BRDF:BRDF_Albedo_Quality", sep = "")
	}

	return(data_set_name)
}

Get_SDS_name_NBAR_snow <- function(mcd43a2_file_path, fill_tile){
	if(is.na(mcd43a2_file_path)){
		# Fill tile
		data_set_name <- fill_tile
	}else{
		data_set_name <- paste("HDF4_EOS:EOS_GRID:\"", mcd43a2_file_path, "\":MOD_Grid_BRDF:Snow_BRDF_Albedo", sep = "")
	}

	return(data_set_name)
}

Get_SDS_name_NBAR_bandquality <- function(mcd43a2_file_path, fill_tile){
	if(is.na(mcd43a2_file_path)){
		# Fill tile
		data_set_name <- fill_tile
	}else{
		data_set_name <- paste("HDF4_EOS:EOS_GRID:\"", mcd43a2_file_path, "\":MOD_Grid_BRDF:BRDF_Albedo_Band_Quality", sep = "")
	}

	return(data_set_name)
}

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
# function to return the date from a file path, works as long as the file paths are from our GEO MCD43A4 data directory
DateFunc <- function(x){
	as.Date(strsplit(basename(x), split="\\.")[[1]][2], format="A%Y%j")
}


#---------------------------------------------------------------------
# this function gets all time series within a particular tile:
# dsets is a list of datasets (filenames), sites is a two column matrix
# with "row", and "col" attributes describing the sample locations in dsets
GetTimeSeriesInTile <- function(dsets, sites, max_open_datasets=2.75e3){
	# need to open files in blocks b/c of max file open limits
	num_blocks <- ceiling(length(dsets) / max_open_datasets)
	out_vals <- matrix(NA, nrow=dim(sites)[1], ncol=length(dsets))

	for(dset_block in 1:num_blocks){
		dset_start <- ((dset_block - 1) * max_open_datasets) + 1
		dset_end <- min((dset_block * max_open_datasets), length(dsets))
		gds <- list()

		for (i in dset_start:dset_end) { gds <- c(gds, GDAL.open(dsets[i])) }

		# loop through each site
		siteids <- unique(sites$siteid)
		for (siteid in siteids){
			# DEBUG
			print(paste("Site:", which(siteids == siteid), "of", length(siteids)))

			# Sites can include multiple pixels within an area, but not all of the pixels
			# are in the site because of the sinusoidal projection. We extract the area
			# enclosing all of the pixels, then extract just the site's pixels
			#
			# First, determine the offset and the region.dim for enclosing area
			offset <- c(min(sites$row[sites$siteid == siteid], na.rm=T), min(sites$col[sites$siteid == siteid], na.rm=T))
			region.dim <- c(max(sites$row[sites$siteid == siteid], na.rm=T), max(sites$col[sites$siteid == siteid], na.rm=T)) - offset + c(1, 1)

			# next we create a mask to subset the extracted values to just the pixels
			# for the site.
			# tmp_df <- subset(sites, siteid == siteid)
			tmp_df <- sites[sites$siteid == siteid, ]
			# tmp_df$row <- factor(tmp_df$row, levels=(offset[1] + 1):(offset[1] + region.dim[1]))
			# tmp_df$col <- factor(tmp_df$col, levels=(offset[2] + 1):(offset[2] + region.dim[2]))
			tmp_df$row <- factor(tmp_df$row, levels=(offset[1]):(offset[1] + region.dim[1] - 1))
			tmp_df$col <- factor(tmp_df$col, levels=(offset[2]):(offset[2] + region.dim[2] - 1))
			mask <- table(tmp_df$row, tmp_df$col)
			mask <- t(mask == 1)

			# loop through each dataset and get values, appending to the output matrix
			for (j in 1:length(gds)) {
				val <- getRasterData(gds[[j]], offset=offset, region.dim=region.dim)
				out_vals[sites$siteid == siteid, (dset_start + j - 1)] <- c(val[mask])
			}
		}

		# for (k in 1:dim(sites)[1]){
		# 	# # DEBUG
		# 	# print(paste("Site: ", k))
		#
		# 	data_row <- sites$row[k]
		# 	data_col <- sites$col[k]
		#
		# 	# loop through each dataset and get values, appending to the output matrix
		# 	for (j in 1:length(gds)) {
		# 		val <- getRasterData(gds[[j]], offset=c(data_row, data_col), region.dim=c(1,1))
		# 		out_vals[k, (dset_start + j - 1)] <- c(val)
		# 	}
		# }

		# close all open files
		for (i in 1:length(gds)) { GDAL.close(gds[[i]]) }
	}

	return(out_vals)
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
# top-level function to retrieve all step sites for single tile
GetTileTimeSeries <- function(this_tile, step, start_date, end_date, save_int=F, int_dir="/projectnb/modislc/users/joshgray/MCD12Q2C6/STEP_ts"){
	# DEBUG
	print(paste("Working on tile: ", this_tile, sep=""))
	sites <- subset(step, tile==this_tile, select=c(3, 4, 6))
	names(sites) <- c("col", "row", "siteid")
	step_sub <- subset(step, tile==this_tile)

	# create the dataset lists and extract the time series
	tile_data <- GetTileDataSets(this_tile, start_date, end_date)
	data_sets <- tile_data[[1]]
	dates <- tile_data[[2]]
	time_series <- GetTimeSeriesInTile(data_sets, sites)

	# combine with step site data
	tmp <- cbind(step_sub, time_series)

	# assign variable names
	var_names <- c("b1", "b2", "b3", "b4", "b6", "qa", "snow", "bandquality")
	names(tmp)[(dim(step)[2] + 1):dim(tmp)[2]] <- paste(rep(strftime(dates, format="%Y%j"), length(var_names)), rep(var_names, each=length(dates)), sep="_")

	# save an intermediate file if requested
	if(save_int){
		# save an intermediate file, in case things go south...
		tmp_out_file <- file.path(int_dir, paste("step_ts_", this_tile, ".Rdata", sep=""))
		save(tmp, file=tmp_out_file)
		return(tmp)
	}else{
		return(tmp)
	}
}

#------------------------------
# for Eli
load("/projectnb/modislc/users/joshgray/MCD12Q2C6/panel_info.Rdata")
site_df <- out_dat[, 1:5]
# names(site_df) <- c("row", "col")
site_df$siteid <- paste("panel", 1:dim(site_df)[1], sep="_")
tile <- "h11v03"
start_date <- as.Date("2001-11-1")
end_date <- as.Date("2003-2-28")
# start_date <- as.Date("2008-11-1")
# end_date <- as.Date("2010-2-28")
panel_MODIS_data <- GetTileTimeSeries(tile, site_df, start_date, end_date)
save(panel_MODIS_data, file="/projectnb/modislc/users/joshgray/MCD12Q2C6/eli_h11v03_2002_output.Rdata")

#------------------------------
# get phenocam lat/long
phenocams <- read.table("http://phenocam.sr.unh.edu/static/txt/sitelocs_active.txt", sep="|", head=T)[, c(3,2,1)]
sin_proj <- projection(raster("/projectnb/modislc/users/joshgray/MCD12Q2C6/filltiles/h12v04fill.tif"))
sin_coords <- project(as.matrix(phenocams[, c(3, 2)]), sin_proj)
phenocams$sin_x <- sin_coords[, 1]
phenocams$sin_y <- sin_coords[, 2]

GetTile <- function(x){
	sin_x <- as.numeric(x[4])
	sin_y <- as.numeric(x[5])
	horiz_bins <- c(-18:18) * 2400 * 463.3127
	horiz_labels <- 0:35
	vert_bins <- c(-9:9) * 2400 * 463.3127
	vert_labels <- 17:0
	v_tile <- formatC(vert_labels[findInterval(sin_y, vert_bins, all.inside=T)], width=2, flag="0")
	h_tile <- formatC(horiz_labels[findInterval(sin_x, horiz_bins, all.inside=T)], width=2, flag="0")
	tile <- paste("h", h_tile, "v", v_tile, sep="")
	return(tile)
}

phenocams$tile <- apply(phenocams, 1, GetTile)

phenocams$row <- phenocams$col <- rep(NA, dim(phenocams)[1])
for(this_tile in unique(phenocams$tile)){
	print(paste("Doing tile:", this_tile))
	# tmp <- phenocams[phenocams$tile == this_tile, ]
	tmp_r <- Sys.glob(paste("/projectnb/modislc/users/joshgray/MCD12Q2C6/filltiles/*", this_tile, "*", sep=""))
	if(!length(tmp_r) == 0){
		tmp_r <- raster(tmp_r)
		phenocams$col[phenocams$tile == this_tile] <- colFromX(tmp_r, phenocams$sin_x[phenocams$tile == this_tile]) - 1

		phenocams$row[phenocams$tile == this_tile] <- rowFromY(tmp_r, phenocams$sin_y[phenocams$tile == this_tile]) - 1
	}else{
		print("No tile found for ", this_tile)
	}
}

# write.csv(phenocams, file="/projectnb/modislc/users/joshgray/MCD12Q2C6/phenocam_locations.csv", row.names=F, quote=F)

# check for successfully completed tiles and subset to only missing ones
output_missing <- function(x) { !length(Sys.glob(paste("/projectnb/modislc/users/joshgray/MCD12Q2C6/PHENOCAM_ts/*", x[1], "*", sep="")) > 0) }
tiles_to_do <- unique(phenocams$tile)[unlist(lapply(unique(phenocams$tile), output_missing))]

start_date <- as.Date("2001-1-1")
end_date <- as.Date("2014-12-31")

# now do the extraction, parallelized over tiles using dopar and foreach
num_cores <- 16
registerDoParallel(num_cores)
all_phenocam_ts <- foreach(x=tiles_to_do, .combine=rbind) %dopar% GetTileTimeSeries_PhenoCam(x, phenocams, start_date, end_date, save_int=T)
# save(all_step_ts, file="/projectnb/modislc/users/joshgray/MCD12Q2C6/phenocam_all_ts.Rdata")



#---------------------------------------------------------------------
# Main operations: now with parallel!
step <- read.csv("/projectnb/modislc/users/joshgray/MCD12Q2C6/step_export.csv",head=F)
names(step) <- c("tile_h", "tile_v", "pix_x", "pix_y", "gid", "siteid", "nveg", "veg", "cov", "lf_type", "lf_phen", "cult", "irig", "ag_type", "wet", "igbb", "struct", "tcov", "domht", "year_start", "year_end")
step$tile <- paste("h", formatC(step$tile_h, width=2, flag="0"), "v", formatC(step$tile_v, width=2, flag="0"), sep="")

# get number of sites and number of pixels per tile
step <- ddply(step, .(tile), transform, num_sites=length(unique(siteid)), num_pixels=length(siteid))

# data frame of sites per tile
step_tile <- data.frame(tile=unique(step$tile), num_sites=step$num_sites[match(unique(step$tile), step$tile)], num_pixels=step$num_pixels[match(unique(step$tile), step$tile)])
step_tile <- step_tile[rev(order(step_tile$num_sites)), ]

start_date <- as.Date("2001-1-1")
# end_date <- as.Date("2014-12-31")
end_date <- as.Date("2015-12-31")

# check for successfully completed tiles and subset to only missing ones
output_missing <- function(x) { !length(Sys.glob(paste("/projectnb/modislc/users/joshgray/MCD12Q2C6/STEP_ts/*", x[1], "*", sep="")) > 0) }
step_tile <- subset(step_tile, apply(step_tile, 1, output_missing))

# now do the extraction, parallelized over tiles using dopar and foreach
num_cores <- 16
registerDoParallel(num_cores)
all_step_ts <- foreach(x=step_tile$tile, .combine=rbind) %dopar% GetTileTimeSeries(x, step, start_date, end_date, save_int=T)
save(all_step_ts, file="/projectnb/modislc/users/joshgray/MCD12Q2C6/step_all_ts.Rdata")

#---------------------------------------------------------------------
# Get h19v03 problem pixels for spline testing
# 1433,981 - funky in winter 2004-2005
# 419,1936 - looks good both winters 2004-2005 and 2005-2006
# 1913,2254 - funky in winter 2004-2005
tile <- "h19v03"
start_date <- as.Date("2004-1-1")
end_date <- as.Date("2007-12-31")
step[1:3,1] <- 19; step[1:3,2] <- 3
step[1:3,22] <- "h19v03"
step[1,3] <- 1433; step[1,4] <- 981; step[1,6] <- 1
step[2,3] <- 419; step[2,4] <- 1936; step[2,6] <- 2
step[2,3] <- 1913; step[2,4] <- 2254; step[3,6] <- 3

step <- step[1:3,]
tmp <- GetTileTimeSeries("h19v03", step, start_date, end_date, save_int=F)


#---------------------------------------------------------------------
# Analyses

b1 <- grep("b1", names(tmp))
b2 <- grep("b2", names(tmp))
b3 <- grep("b3", names(tmp))
b4 <- grep("b4", names(tmp))
b6 <- grep("b6", names(tmp))
snow <- grep("snow", names(tmp))
# ndvi <- (tmp[, b2] - tmp[, b1]) / (tmp[, b2] + tmp[, b1])
evi2 <- 2.5 * ((tmp[, b2] - tmp[, b1]) / (tmp[, b2] + (2.4 * tmp[, b1]) + 1))
ndsi <- (tmp[, b4] - tmp[, b6]) / (tmp[, b4] + tmp[, b6])
# ndii <- (tmp[, b2] - tmp[, b6]) / (tmp[, b2] + tmp[, b6])
# pheno_index <- rep(NA, length(dates))
# pheno_index[ndvi < 0 | ndii < 0] <- 0
# pheno_index <- (ndvi^2) - (ndii^2)
# pheno_index[pheno_index < 0] <- 0

dates <- as.Date(names(tmp)[b1],format="%Y%j")


#---------------------
# for spline testing:
par(mfrow=c(3,1), mar=rep(2,4))
pred_dates <- seq(dates[1], dates[length(dates)], by=1)
for(i in 1:3){
	tmp_evi2 <- as.numeric(evi2[i,])
	tmp_ndsi <- as.numeric(ndsi[i,])
	tmp_snow <- as.numeric(tmp[i, snow])
	snow_fills <- rep(F, length(tmp_ndsi))
	snow_fills[tmp_ndsi > 0 | tmp_ndsi > 0] <- T

	plot(dates, tmp_evi2, type="n")
	points(dates[tmp_snow == 0], tmp_evi2[tmp_snow == 0], pch=16, col=1)
	points(dates[tmp_snow == 1], tmp_evi2[tmp_snow == 1], pch=4, col=4)
	points(dates[tmp_ndsi > -0.2], tmp_evi2[tmp_ndsi > -0.2], pch=1, col=2, ex=1.5)
	points(pred_dates, SmoothSnowFill(tmp_evi2, snow_fills, as.numeric(dates), as.numeric(pred_dates), pheno_pars, plot=F), type="l", col=2)

}
legend("topleft", legend=c("snow flagged", "no snow flag", "ndsi > 0"), col=c(4,1,2), pch=c(4,16,1), cex=2)

SmoothSnowFill <- function(evi2, snow_fills, dates, pred_dates, pheno_pars, plot=F){

i <- 3
tmp_evi2 <- as.numeric(evi2[i,])
tmp_ndsi <- as.numeric(ndsi[i,])
tmp_snow <- as.numeric(tmp[i, snow])
snow_fills <- rep(F, length(tmp_ndsi))
snow_fills[tmp_ndsi > 0 | tmp_ndsi > 0] <- T
pred_dates <- seq(dates[1], dates[length(dates)], by=1)
trash <- SmoothSnowFill(tmp_evi2, snow_fills, as.numeric(dates), as.numeric(pred_dates), pheno_pars, plot=T)

# calculate site means
evi2_mean <- aggregate(evi2, by=list(gid=tmp$gid), FUN=mean)

############################
# Make a data.frame that looks like the thing that DoPhenology expects:
# b1, b2, b3, b4, b6
tmp_df <- data.frame(tmp[, c(b1, b2, b3, b4, b6)])

# Do some phenology
library(zoo)
library(RColorBrewer)
source("/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_Functions.R")
pheno_pars <- list(
	out_quant=0.99,
	# spline_spar=0.15,
	spline_spar=0,
	out_iterations=0,
	min_peak_to_peak_distance=180,
	# min_seg_amplitude=0.2,
	min_seg_amplitude=0.15, # what was used for CLMv3
	agg_amp_frac=0.15,
	max_seg_length=200,
	evi2_snow_quant=0.02,
	ndsi_thresh=0,
	gup_threshes=c(0.2,0.5,1),
	gdown_threshes=c(0.5,0.2),
	max_snow_fill_ratio=1.25
)

DoPhenology(as.numeric(tmp_df[1,]), dates, pheno_pars)


############################
# plot all cropland site pixel EVI2 time series and mean
EvenOrNextHighest <- function(x) {ifelse(x %% 2 == 0, x, x + 1) }
croplands_gids <- unique(tmp$gid[which(tmp$igbb == 12)])
all_gids <- unique(tmp$gid)
layout(matrix(1:EvenOrNextHighest(length(croplands_gids)), ncol=2))
par(mar=rep(1,4), oma=rep(0.5, 4))
for(i in 1:length(croplands_gids)){
	plot(dates, as.numeric(evi2_mean[which(evi2_mean$gid == croplands_gids[i]), 2:dim(evi2_mean)[2]]), type="n", xlab="", ylab="NBAR-EVI2")
	for(k in which(tmp$gid == croplands_gids[i])) { points(dates, as.numeric(evi2[k, ]), col="grey", type="l") }
	points(dates, as.numeric(evi2_mean[which(evi2_mean$gid == croplands_gids[i]), 2:dim(evi2_mean)[2]]), type="l", lwd=2)
}

############################
# Checking that the projection information works
h <- subset(step,tile_h == 12 & tile_v == 4)
point_i <- 1
project(matrix(c(xFromCol(r, h$pix_x[point_i] - 1), yFromRow(r, h$pix_y[point_i]) - 1), nrow=1), projection(r), inv=T)
