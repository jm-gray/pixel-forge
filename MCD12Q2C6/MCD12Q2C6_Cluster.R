# example call
# R --args -tile h12v04 -gup 0.1 0.5 1 -gdown 0.1 0.5 -min_seg_amplitude 0.1 -out_prefix new_

#------------------------------------------
#------------------------------------------
# Prelims
rm(list=ls())
library(raster)
library(rgdal)
library(parallel)
library(doParallel)
library(foreach)
library(RColorBrewer)
library(zoo)
library(argparse)
library(pracma)
source("/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_Functions.R")

#------------------------------------------
#------------------------------------------
# parse the command line arguments
arg_parser <- ArgumentParser()

# general arguments
arg_parser$add_argument("-tile", type="character")
arg_parser$add_argument("-start_date", type="character", default="2000-1-1")
arg_parser$add_argument("-end_date", type="character", default="2014-12-31")
arg_parser$add_argument("-out_dir", type="character", default="/projectnb/modislc/users/joshgray/MCD12Q2C6/output")
arg_parser$add_argument("-out_prefix", type="character", default="pheno_c6_")
arg_parser$add_argument("-data_dir", type="character", default="/projectnb/modislc/data/mcd12_in/c5/mcd43a4")
arg_parser$add_argument("-qa_dir", type="character", default="/projectnb/modislc/data/mcd12_in/c5/mcd43a2")
arg_parser$add_argument("-fill_dir", type="character", default="/projectnb/modislc/users/joshgray/MCD12Q2C6/filltiles")

# pheno parameters
arg_parser$add_argument("-out_quant", type="double", default=0.99)
arg_parser$add_argument("-spline_spar", type="double", default=0.2)
arg_parser$add_argument("-out_iterations", type="integer", default=0)
arg_parser$add_argument("-min_peak_to_peak_distance", type="integer", default=90)
arg_parser$add_argument("-min_seg_amplitude", type="double", default=0.15)
arg_parser$add_argument("-agg_amp_frac", type="double", default=0.15)
arg_parser$add_argument("-max_seg_length", type="integer", default=200)
arg_parser$add_argument("-evi2_snow_quant", type="double", default=0.05)
arg_parser$add_argument("-ndsi_thresh", type="double", default=0)
arg_parser$add_argument("-max_snow_fill_ratio", type="double", default=1.25)
arg_parser$add_argument("-gup_threshes", type="double", nargs="*", default=0.5)
arg_parser$add_argument("-gdown_threshes", type="double", nargs="*", default=0.5)

# temporary testing parameters
arg_parser$add_argument("-GDALopen", action="store_true") # flag to use faster GDAL.open
arg_parser$add_argument("-C6input", action="store_true") # flag to use faster GDAL.open

# parse the args, change dates to Date type, and check if out_file is missing
args <- arg_parser$parse_args()
args$start_date <- as.Date(args$start_date)
args$end_date <- as.Date(args$end_date)
args$out_file <- file.path(args$out_dir, paste(args$out_prefix, args$tile, ".Rdata", sep=""))

# special case: if spline_spar is 0, then we set to NULL (GCV based smoothing)
if(args$spline_spar == 0) args$spline_spar <- NULL

# DEBUG: for parsing the *.sh.o* file to find what tile is associated with a particular job
print(paste("Doing tile:", args$tile))

#------------------------------------------
#------------------------------------------
# create a RasterStack of NBAR data from which to calculate EVI, EVI2, NDSI, etc

# C6 vs C5 data directories are setup differently
if(!args$C6input){
	# get all datasets within date range
	datasets <- GetTileDataSets(args$tile, start_date=args$start_date, end_date=args$end_date, nbar_data_dir=args$data_dir, nbarqa_data_dir=args$qa_dir, tmp_out_dir=args$fill_dir)
	data_sets <- datasets[[1]]
	my_dates <- datasets[[2]]
	# GetTileDataSets also returns 3 QA datasets, but they aren't currently used, so we eliminate it
	data_sets <- data_sets[1:(length(my_dates) * 5)]
}else{
	# this still uses the older way of doing thins...
	nbar_files <- Sys.glob(file.path(args$data_dir, "*", paste("*MCD43A4*", args$tile, "*hdf", sep=""))) # for c6 NBARS
	dates <- do.call("c", lapply(nbar_files, DateFunc)) # get the dates for each file
	nbar_files <- nbar_files[(dates >= args$start_date) & (dates <= args$end_date)]
	my_dates <- dates[(dates >= args$start_date) & (dates <= args$end_date)]
	# make a vector of the SDS names
	data_sets <- c(unlist(lapply(nbar_files, Get_SDS_name_b1)), unlist(lapply(nbar_files, Get_SDS_name_b2)), unlist(lapply(nbar_files, Get_SDS_name_b3)), unlist(lapply(nbar_files, Get_SDS_name_band4)), unlist(lapply(nbar_files, Get_SDS_name_band6)))
}

# create a stack if we're not using the GDAL.open method
if(!args$GDALopen){
	# make a stack: red, nir, blue, b4, b6
	nbar_s <- stack(data_sets)
}else{
	# get a temporary raster
	nbar_s <- raster(data_sets[1])
}

#------------------------------------------
#------------------------------------------
# Do phenology
registerDoParallel(16)
rows_to_do <- 30
i_seq <- 1:(nrow(nbar_s) / rows_to_do)

# getValues() and stack() based function to apply to each row block
foreach_phenology_function_stack <- function(i, nbar_s, dates, rows_to_do, args){
	v <- getValues(nbar_s, ((i - 1) * rows_to_do) + 1, rows_to_do)
	tmp_result <- apply(v, 1, DoPhenology, dates, args)
}

# GDAL.open based function to apply to each row block
foreach_phenology_function_GDAL <- function(i, data_sets, dates, rows_to_do, args){
	v <- getValuesGDAL(data_sets, i, rows_to_do)
	tmp_result <- apply(v, 1, DoPhenology, dates, args)
}

# apply the phenology function with foreach/dopar
if(!args$GDALopen){
	pheno_elapsed_time <- system.time(
		result <- foreach(i=i_seq, .combine=c) %dopar% foreach_phenology_function_stack(i, nbar_s, my_dates, rows_to_do, args)
	)
}else{
	pheno_elapsed_time <- system.time(
		result <- foreach(i=i_seq, .combine=c) %dopar% foreach_phenology_function_GDAL(i, data_sets, my_dates, rows_to_do, args)
	)
}


# save the result object as an Rdata file
save(pheno_elapsed_time, result, file=args$out_file)

# # extract phenophase transition DOY for a particular year and make a raster
# # NOTE: this way would apparently take 43 hours...
# system.time(
# 	doy_2011 <- foreach(pheno_date=result, .combine=c) %dopar% extract_doy(pheno_date, 1, 2011)
# )


# #------------------------------------------
# #------------------------------------------
# # Extract a threshold date for each year and write to a multi-layer raster file
# cl <- makeCluster(16) # make cluster
#
# # extract year in years
# thresh_to_extract <- 2
# # thresh_to_extract <- 5
# tmp_r <- raster(nbar_s, 1)
# # for(year_to_extract in 2008:2012){
# for(year_to_extract in 2000:2013){
# 	print(paste("Doing year: ", year_to_extract))
# 	pheno_doy_v <- unlist(parLapply(cl, result, ExtractDOYYear, i=thresh_to_extract, year=year_to_extract, cycle=1))
# 	tmp_r <- setValues(tmp_r, pheno_doy_v)
# 	# if(year_to_extract == 2008){
# 	if(year_to_extract == 2000){
# 		mid_gup_doy_s <- stack(tmp_r)
# 	}else{
# 		mid_gup_doy_s <- stack(mid_gup_doy_s, tmp_r)
# 	}
# }
#
# # write to a raster
# # out_raster_prefix <- "pheno_c6_"
# # out_raster_suffix <- "_2000_2013_gup_stack.tif"
# # out_raster_dir <- "/projectnb/modislc/users/joshgray/MCD12Q2C6/output_tiff"
# # out_raster_file <- file.path(out_raster_dir, paste(out_raster_prefix, args$tile, out_raster_suffix, sep=""))
# # writeRaster(file=out_raster_file, mid_gup_doy_s)
# writeRaster(file=out_file, mid_gup_doy_s)
# # writeRaster(file="/projectnb/modislc/users/joshgray/MCD12Q2C6/pheno_c6_h11v03_outstack.tif", mid_gup_doy_s)


# #------------------------------------------
# #------------------------------------------
# # Extract number of cycles within window
# thresh_to_extract <- 2
# tmp_r <- raster(nbar_s, 1)
# for(year in 2000:2011){
# 	win_start <- as.Date(paste(year, "-7-1", sep=""))
# 	win_end <- as.Date(paste(year + 3, "-6-30", sep=""))
# 	print(paste("Start:", win_start, "End:", win_end))
#
# 	pheno_num_cycles_v <- unlist(parLapply(cl, result, ExtractNumCyclesWindow, i=thresh_to_extract, start.date=win_start, end.date=win_end))
# 	tmp_r <- setValues(tmp_r, pheno_num_cycles_v)
# 	if(year == 2000){
# 		num_cycles_s <- stack(tmp_r)
# 	}else{
# 		num_cycles_s <- stack(num_cycles_s, tmp_r)
# 	}
# }
# # write to a raster
# out_raster_prefix <- "pheno_c6_"
# out_raster_suffix <- "_num_cycles_3yr_v2.tif"
# out_raster_dir <- "/projectnb/modislc/users/joshgray/MCD12Q2C6"
# out_raster_file <- file.path(out_raster_dir, paste(out_raster_prefix, tile, out_raster_suffix, sep=""))
# writeRaster(file=out_raster_file, num_cycles_s)
#
#
# tmp_s <- stack("pheno_c6_h12v11_num_cycles_3yr.tif")
# for(i in 1:nlayers(tmp_s)){
# 	print(paste("doing:", i))
# 	tmp_v <- values(raster(tmp_s, i))
# 	if(!exists("num_greater_than_three")){
# 		num_greater_than_three <- sum(tmp_v > 3, na.rm=T)
# 	}else{
# 		num_greater_than_three <- c(num_greater_than_three, sum(tmp_v > 3, na.rm=T))
# 	}
# }
#
# for(i in 1:nlayers(tmp_s_crop)){
# 	print(paste("doing:", i))
# 	tmp_v <- values(raster(tmp_s_crop, i))
# 	if(!exists("num_greater_than_three_crop")){
# 		num_greater_than_three_crop <- sum(tmp_v > 3, na.rm=T)
# 	}else{
# 		num_greater_than_three_crop <- c(num_greater_than_three_crop, sum(tmp_v > 3, na.rm=T))
# 	}
# }
#
# par(mar=c(2, 4, 1, 1))
# plot(2002:2013, num_greater_than_three_crop / ncell(tmp_s_crop), xlab="", ylab="Fraction of pixels > 1 cycle/yr", type="l", col=1, lwd=2)
# points(2002:2013, num_greater_than_three / ncell(tmp_s), type="l", col=2, lwd=2, lty=2)
# legend("bottomleft", legend=c("Tile h12v11", "Taboada Belgrano"), lty=c(2, 1), col=c(2, 1))
#
# my_breaks <- c(0, 2, 3, 10)
# my_colors <- c("burlywood", "lightgreen", "darkgreen")
# plot(raster(tmp_s, 1), breaks=my_breaks, col=my_colors, legend=F)
# plot(my_x, add=T, col=2)
# plot(raster(matrix(1:3)), breaks=c(0,1,2,3), col=my_colors, legend.only=T, axis.args=list(at=c(0.5, 1.5, 2.5), labels=c("<3","3",">3")))
# title("tile h12v11 number of veg cycles 2000/7/1-2003/6/30")
#
# # for the subset
# cur_dir <- getwd()
# setwd("/projectnb/modislc/users/joshgray/MCD12Q2C6/TaboadaBelgrano_shapefile")
# tmp_shp <- readOGR("departabandera_sin.shp", "departabandera_sin")
# setwd(cur_dir)
#
# my_breaks <- c(0, 2, 3, 10)
# my_colors <- c("burlywood", "lightgreen", "darkgreen")
# layout(matrix(1:12, nrow=3))
# par(mar=c(0.5, 0.5, 2, 0.5))
# for(i in 1:12){
# 	plot(raster(tmp_s_crop, i), breaks=my_breaks, col=my_colors, legend=F, xlab="", ylab="", axes=F, maxpixels=253000)
# 	plot(tmp_shp, add=T, border=2)
# 	title(paste((1999+i), "-7-1 to ", (1999 + i + 3), "-6-30", sep=""))
# }
#
# # plot(tmp_s_crop, breaks=my_breaks, col=my_colors, legend=F, xlab="", ylab="")



# #------------------------------------------
# #------------------------------------------
# # Extract annual GSL for first cycle
# # lapply(result[1:10],ExtractGSLWindow, i_gup=2, i_gdown=5, start.date="2000-7-1", end.date="2001-6-30")
#
# i_gup <- 2
# i_gdown <- 5
# tmp_r <- raster(nbar_s, 1)
# for(year in 2000:2014){
# 	win_start <- as.Date(paste(year, "-7-1", sep=""))
# 	win_end <- as.Date(paste(year + 1, "-6-30", sep=""))
# 	print(paste("Start:", win_start, "End:", win_end))
#
# 	gsl_v <- unlist(parLapply(cl, result, ExtractGSLWindow, i_gup=i_gup, i_gdown=i_gdown, start.date=win_start, end.date=win_end, cycle=1))
# 	tmp_r <- setValues(tmp_r, gsl_v)
# 	if(year == 2000){
# 		gsl_s <- stack(tmp_r)
# 	}else{
# 		gsl_s <- stack(gsl_s, tmp_r)
# 	}
# }
# # write to a raster
# out_raster_prefix <- "pheno_c6_"
# out_raster_suffix <- "_gsl_annual_v2.tif"
# out_raster_dir <- "/projectnb/modislc/users/joshgray/MCD12Q2C6"
# out_raster_file <- file.path(out_raster_dir, paste(out_raster_prefix, tile, out_raster_suffix, sep=""))
# writeRaster(file=out_raster_file, gsl_s)




# #------------------------------------------
# #------------------------------------------
# # Extract other phenometrics
# # extract number of cycles
# system.time(
# 	pheno_num_cycles_v <- unlist(parLapply(cl, result, ExtractNumCycles))
# )
#
# # extract GUP Dates
# system.time(
# 	pheno_gup_v <- parLapply(cl, result, ExtractDates, i=1)
# )
#
# # extract GDOWN Dates
# system.time(
# 	pheno_gdown_v <- parLapply(cl, result, ExtractDates, i=2)
# )
#
# # extract GSL Dates
# system.time(
# 	pheno_gsl_v <- parLapply(cl, result, ExtractGSL, i_gup=1, i_gdown=2)
# )


# #------------------------------------------
# #------------------------------------------
# # write number of cycles to raster
# pheno_num_cycles_r <- raster(nbar_s, 1)
# pheno_num_cycles_r <- setValues(pheno_num_cycles_r, pheno_num_cycles_v)
#
# pheno_gsl_r <- raster(nbar_s, 1)
# pheno_gsl_r <- setValues(pheno_gsl_r, pheno_gsl_v)

# ------------------------------------------
# ------------------------------------------
# merge rasters
# i <- 1 # i=1 is 2008
# r_merge <- merge(
# 	raster("pheno_c6_h08v05_outstack.tif", i),
# 	raster("pheno_c6_h08v06_outstack.tif", i),
# 	raster("pheno_c6_h09v04_outstack.tif", i),
# 	raster("pheno_c6_h09v05_outstack.tif", i),
# 	raster("pheno_c6_h09v06_outstack.tif", i),
# 	raster("pheno_c6_h10v03_outstack.tif", i),
# 	raster("pheno_c6_h10v04_outstack.tif", i),
# 	raster("pheno_c6_h10v05_outstack.tif", i),
# 	raster("pheno_c6_h10v06_outstack.tif", i),
# 	raster("pheno_c6_h11v03_outstack.tif", i),
# 	raster("pheno_c6_h11v04_outstack.tif", i),
# 	raster("pheno_c6_h11v05_outstack.tif", i),
# 	raster("pheno_c6_h12v02_outstack.tif", i),
# 	raster("pheno_c6_h12v03_outstack.tif", i),
# 	raster("pheno_c6_h12v04_outstack.tif", i),
# 	raster("pheno_c6_h12v05_outstack.tif", i),
# 	raster("pheno_c6_h13v03_outstack.tif", i),
# 	raster("pheno_c6_h13v04_outstack.tif", i),
# 	raster("pheno_c6_h14v03_outstack.tif", i)
# )
# r_merge_v <- values(r_merge) # get values for statistics


#------------------------------------------
#------------------------------------------
# plot pheno date

# library(RColorBrewer)
# lin_cutoff <- 0.02
# # my_q <- quantile(pheno_doy_v, c(lin_cutoff, (1 - lin_cutoff)), na.rm=T)
# my_q <- quantile(r_merge_v, c(lin_cutoff, (1 - lin_cutoff)), na.rm=T)
# if((my_q[2] - my_q[1]) >= 254){
# 	break_len <- 256
# }else{
# 	break_len <- (my_q[2] - my_q[1]) + 2
# }
#
# gup_breaks <- c(1, seq(my_q[1], my_q[2], len=(break_len - 2)), 366)
# my_pal <- colorRampPalette(brewer.pal(11, "Spectral"))
# # plot(pheno_doy_r, breaks=gup_breaks, col=my_pal(break_len))
#
#
# pdf(file="NA_2008.pdf", h=17, w=25)
# plot(r_merge, breaks=gup_breaks, col=my_pal(break_len), maxpixels=45e6, legend=F)
# tmp_r <- raster(matrix(my_q[1]:my_q[2]))
# plot(tmp_r, col=my_pal(break_len), legend.only=T)
# title("2008 50% Greenup Amplitude")
# dev.off()

#------------------------------------------
#------------------------------------------












#-------------------------------------------------------------
#-------------------------------------------------------------
#-------------------------------------------------------------
# parApply version
#-------------------------------------------------------------
#-------------------------------------------------------------
#-------------------------------------------------------------

# oldstuff{
# #-------------------------------------------------------------
# # Running the function
# # cl <- makeCluster(16)
# # clusterEvalQ(cl, {library(raster); source("/projectnb/modislc/users/joshgray/LA_MC/JoshSeg.R")})
# # clusterEvalQ(cl, {library(raster); source("/projectnb/modislc/users/joshgray/MCD12Q2C6_Functions.R")})
# beginCluster()
# cl <- getCluster()
# clusterExport(cl, c("my_dates"))
# clusterEvalQ(cl, {source("/projectnb/modislc/users/joshgray/MCD12Q2C6_Functions.R")})
#
# start_time <- Sys.time()
# tmp <- clusfun(nbar_s, dates=my_dates, cl)
# print(Sys.time() - start_time)
# elapsed_time <- Sys.time() - start_time
# save(tmp, elapsed_time, file="/projectnb/modislc/users/joshgray/ValidPeaksTmp_new.Rdata")
#
# #-------------------------------------------------------------
# # this example function takes a RasterStack of DSM's composites and returns the maximum NDVI across all years
# clusfun <- function(x, dates, cl, ...){
# 	# node_mult <- 4 # should probably be determined by available memory and size of x...
# 	# out <- raster(x)
# 	# out <- array(NA, nrow(x), ncol(x))
# 	out <- array(NA, nrow(x) * ncol(x))
# 	# cl <- getCluster()
# 	on.exit(returnCluster())
# 	nodes <- length(cl)
# 	# bs <- blockSize(x, minblocks=nodes*node_mult)
# 	bs <- blockSize(x)
# 	pb <- pbCreate(bs$n, style=3, progress="text", label="Progress")
#
# 	# the function to be used (simple example)
# 	# clFun <- function(i, dates){
# 	clFun <- function(i){
# 		#-----------------------------------------------------
# 		# the actual phenology function
# 		#-----------------------------------------------------
#
# 		# # source the required functions
# 		# source("/projectnb/modislc/users/joshgray/MCD12Q2C6_Functions.R")
# 		#
# 		# # get the values as a data.frame
# 		v <- getValues(x, bs$row[i], bs$nrows[i])
# 		#
# 		# # run the testing function to gather valid peaks over the block
# 		# # filt_evi2 <- t(apply(v, 1, SnowFilterEVI2, evi2_snow_quant=evi2_snow_quant))
# 		peaks <- apply(v, 1, TestingFunction, dates)
# 		return(peaks)
#
# 		#-----------------------------------------------------
# 		# return 1-4 random values from x
# 		#-----------------------------------------------------
# 		# randFunc <- function(x){
# 		# 	# returns samples of x of random length in 1:4
# 		# 	val <- sample(1:length(x), round(runif(1, 1, 4)))
# 		# 	return(val)
# 		# }
# 		#
# 		# v <- getValues(x, bs$row[i], bs$nrows[i])
# 		# randSamples <- apply(v, 1, randFunc)
# 		#
# 		# # this returns a list of lists, one per pixel
# 		# return(randSamples)
#
# 		#-----------------------------------------------------
# 		# just a simple copy for testing and timing
# 		#-----------------------------------------------------
# 		# v <- getValues(x, bs$row[i], bs$nrows[i])
# 		# return(v)
# 	}
#
# 	# get all nodes going
# 	for (i in 1:nodes) {
# 		sendCall(cl[[i]], clFun, i, tag=i)
# 		# sendCall(cl[[i]], clFun, list(i=i, dates=dates), tag=i)
# 		# sendCall(cl[[i]], clFun, i=i, dates=dates, tag=i)
# 		# sendCall(cl[[i]], clFun, args=list(i=i, dates=dates), tag=i)
# 		# sendCall(cl[[i]], clFun, args=list(i=i, x=x, dates=dates), tag=i)
# 		# sendCall(cl[[i]], clFun, args=list(i=i, dates=dates), tag=i)
# 	}
#
# 	for (i in 1:bs$n){
# 		# receive results from a node
# 		d <- recvOneData(cl)
#
# 		# error?
# 		if(!d$value$success){
# 			# DEBUG
# 			return(d)
#
# 			stop('cluster error')
# 		}
#
# 		# which block is this?
# 		b <- d$value$tag
#
# 		# DEBUG
# 		cat('received block: ',b,'\n'); flush.console();
# 		return(d)
#
# 		# write block results to the output array
# 		out_start <- (nrow(x) * (bs$row[b] - 1)) + 1
# 		out_end <- (nrow(x) * bs$nrows[b]) + (nrow(x) * (bs$row[b] - 1))
# 		out[out_start:out_end] <- d$value$value
#
# 		# need to send more data?
# 		ni <- nodes + i
# 		if(ni <= bs$n){
# 			sendCall(cl[[d$node]], clFun, ni, tag=ni)
# 			# sendCall(cl[[d$node]], clFun, list(i=ni, dates=dates), tag=ni)
# 			# sendCall(cl[[d$node]], clFun, i=ni, dates=dates, tag=ni)
# 			# sendCall(cl[[d$node]], clFun, args=list(i=ni, dates=dates), tag=ni)
# 			# sendCall(cl[[d$node]], clFun, args=list(i=ni, x=x, dates=dates), tag=ni)
# 			# sendCall(cl[[d$node]], clFun, args=list(i=ni, dates=dates), tag=ni)
# 		}
# 		pbStep(pb)
# 	}
#
# 	pbClose(pb)
#
# 	return(out)
# }
#
#
# }
