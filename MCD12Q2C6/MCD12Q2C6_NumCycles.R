# example call
# R --args -tile h12v04 -i 2 -prefix numcycles
# R --args -tile h12v04 -i 3 -prefix numcycles

#------------------------------------------
#------------------------------------------
# Prelims
rm(list=ls())
library(raster)
library(parallel)
library(doParallel)
library(foreach)
library(RColorBrewer)
library(zoo)
library(argparse)
source("/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_Functions.R")

#------------------------------------------
#------------------------------------------
# parse the command line arguments
arg_parser <- ArgumentParser()
arg_parser$add_argument("-tile", type="character")
arg_parser$add_argument("-in_dir", type="character", default="/projectnb/modislc/users/joshgray/MCD12Q2C6/output_new")
arg_parser$add_argument("-out_dir", type="character", default="/projectnb/modislc/users/joshgray/MCD12Q2C6/output_tiff_new")
arg_parser$add_argument("-data_dir", type="character", default="/projectnb/modislc/data/mcd12_in/c5/mcd43a4")
arg_parser$add_argument("-yearrange", type="integer", nargs=2, default=c(2000, 2011))
arg_parser$add_argument("-win_years", type="integer", nargs=1, default=3)
arg_parser$add_argument("-i", type="integer", nargs=1, default=1)
arg_parser$add_argument("-prefix", type="character", nargs=1, default="num_cycles")
arg_parser$add_argument("-num_nodes", type="integer", default=8)
args <- arg_parser$parse_args()

#------------------------------------------
#------------------------------------------
# get the Rdata file and load
in_file <- dir(args$in_dir, pattern=paste(".*", args$tile, ".*.Rdata", sep=""), full=T)
load(in_file)

#------------------------------------------
#------------------------------------------
# get a temporary RasterLayer for this tile in order to set extent, geo, etc
nbar_files <- Sys.glob(file.path(args$data_dir, "*", "*", paste("*", args$tile, "*", sep="")))
tmp_r <- raster(Get_SDS_name_red(nbar_files[1]))

#------------------------------------------
#------------------------------------------
# Extract a threshold date for each year and write to a multi-layer raster file
cl <- makeCluster(args$num_nodes) # make cluster

# extract for each threshold
for(start_year in args$yearrange[1]:args$yearrange[2]){
	# get the vertical tile number, and choose start/end dates based on N/S Hemisphere
	if(as.integer(substr(args$tile, 5, 6)) < 9){
		# Northern Hemisphere dates
		start_date <- as.Date(paste(start_year, "-1-1", sep=""))
		end_date <- as.Date(paste((start_year + (args$win_years - 1)), "-12-31", sep=""))
	}else{
		# Southern Hemisphere
		start_date <- as.Date(paste(start_year, "-7-1", sep=""))
		end_date <- as.Date(paste((start_year + args$win_years), "-6-30", sep=""))
	}
	
	# DEBUG
	print(paste("Start:", start_date, "End:", end_date))
	
	# do the extraction
	num_cycles <- unlist(parLapply(cl, result, ExtractNumCyclesWindow, i=args$i, start_date=start_date, end_date=end_date))
	
	# set Raster values
	tmp_r <- setValues(tmp_r, num_cycles)
	
	# check for first time through, and create, or add to RasterStack
	if(start_year == args$yearrange[1]){
		num_cycles_stack <- stack(tmp_r)
	}else{
		num_cycles_stack <- stack(num_cycles_stack, tmp_r)
	}
}

# write output raster
out_file <- file.path(args$out_dir, paste(args$prefix, "_", args$tile, ".tif", sep=""))
writeRaster(file=out_file, num_cycles_stack)
