# example call
# R --args -tile h12v04 -i 2 -prefix midgup
# R --args -tile h12v04 -i 2 4 -i_key 2 -prefix new_midgup new_midgdown -in_dir /projectnb/modislc/users/joshgray/MCD12Q2C6/output_new -out_dir /projectnb/modislc/users/joshgray/MCD12Q2C6/output_tiff_new

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
arg_parser$add_argument("-in_dir", type="character", default="/projectnb/modislc/users/joshgray/MCD12Q2C6/output")
arg_parser$add_argument("-out_dir", type="character", default="/projectnb/modislc/users/joshgray/MCD12Q2C6/output_tiff")
arg_parser$add_argument("-data_dir", type="character", default="/projectnb/modislc/data/mcd12_in/c5/mcd43a4")
arg_parser$add_argument("-data_prefix", type="character", default="")
arg_parser$add_argument("-yearrange", type="integer", nargs=2, default=c(2000, 2014))
arg_parser$add_argument("-i", type="integer", nargs="*", default=1)
arg_parser$add_argument("-i_key", type="integer", default=2)
arg_parser$add_argument("-prefix", type="character", nargs="*", default="pheno")
arg_parser$add_argument("-n", type="integer", nargs="*", default=1)
arg_parser$add_argument("-num_nodes", type="integer", default=8)
arg_parser$add_argument("-annual", type="logical", default=T)
arg_parser$add_argument("-overwrite", type="logical", default=T)
args <- arg_parser$parse_args()

#------------------------------------------
#------------------------------------------
# get the Rdata file and load
# in_file <- dir(args$in_dir, pattern=paste(".*", args$tile, ".*.Rdata", sep=""), full=T)
in_file <- dir(args$in_dir, pattern=paste(".*", args$data_prefix, ".*", args$tile, ".*.Rdata", sep=""), full=T)
load(in_file)

# DEBUG
print(paste("Loading file:", in_file))

#------------------------------------------
#------------------------------------------
# get a temporary RasterLayer for this tile in order to set extent, geo, etc
nbar_files <- Sys.glob(file.path(args$data_dir, "*", "*", paste("*", args$tile, "*", sep="")))
tmp_r <- raster(Get_SDS_name_band1(nbar_files[1]))

#------------------------------------------
#------------------------------------------
# Extract each requested threshold date for each year and write to a multi-layer raster file
# year range
cl <- makeCluster(args$num_nodes) # make cluster

# extract for each threshold
n <- 1 # counter for args$prefix
for(i in args$i){
	# extract year in years
	for(year_to_extract in args$yearrange[1]:args$yearrange[2]){
		print(paste("Metric:", i, "Year: ", year_to_extract))

		# get the vertical tile number, and choose start/end dates based on N/S Hemisphere
		if(as.integer(substr(args$tile, 5, 6)) < 9){
			# Northern Hemisphere dates
			start_date <- as.Date(paste(year_to_extract, "-1-1", sep=""))
			end_date <- as.Date(paste(year_to_extract, "-12-31", sep=""))
		}else{
			# Southern Hemisphere
			start_date <- as.Date(paste(year_to_extract, "-7-1", sep=""))
			end_date <- as.Date(paste(year_to_extract + 1, "-6-30", sep=""))
		}

		# do the extraction
		# pheno_dates <- unlist(parLapply(cl, result, ExtractSingleDate, i=i, start_date=start_date, end_date=end_date, cycle=args$cycle))
		pheno_dates <- unlist(parLapply(cl, result, ExtractSingleDate, i=i, i_key=args$i_key, start_date=start_date, end_date=end_date, n=args$n))

		# set Raster values
		tmp_r <- setValues(tmp_r, pheno_dates)

		# check if annual or total output is requested, if annual, then write the raster and move on, otherwise add it to the stack
		if(args$annual){
			out_file <- file.path(args$out_dir, paste(args$prefix[n], "_", args$tile, "_", year_to_extract,".tif", sep=""))
			writeRaster(tmp_r, file=out_file, overwrite=args$overwrite)
		}else{
			# check for first time through, and create, or add to RasterStack
			if(year_to_extract == args$yearrange[1]){
				pheno_stack <- stack(tmp_r)
			}else{
				pheno_stack <- stack(pheno_stack, tmp_r)
			}
		}
	}

	# write output raster as a big stack if annual values were note requested
	if(!args$annual){
		out_file <- file.path(args$out_dir, paste(args$prefix[n], "_", args$tile, ".tif", sep=""))
		writeRaster(file=out_file, pheno_stack, overwrite=args$overwrite)
	}

	n <- n + 1
}
