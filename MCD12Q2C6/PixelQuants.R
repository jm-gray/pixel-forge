#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Josh Gray, June 2016
# Here, the computation is low but the number of files are large, so we parallelize
# across data chunks: each "node" loads an equal number of rows and processes the
# snowfree EVI2 quantiles. Results are bound together and written as binary integers
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Submitted to the GridEngine job handler w/ shell script (subPixelQuants.sh):
# #!/bin/bash
# echo Submitting tile: $1
# R --vanilla < /projectnb/modislc/users/joshgray/MCD12Q2C6/PixelQuants.R --args -tile $1

# example submission command using default parameters:
# qsub -V -pe omp 16 -l h_rt=02:00:00 -l mem_total=98G subPixelQuants.sh h11v04

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # batch submit many tiles:
# namerica_tiles <- c('h10v02', 'h10v03', 'h10v04', 'h10v05', 'h10v06', 'h10v07', 'h10v08', 'h11v02', 'h11v03', 'h11v04', 'h11v05', 'h11v06', 'h11v07', 'h12v01', 'h12v02', 'h12v03', 'h12v04', 'h12v05', 'h12v07', 'h13v01', 'h13v02', 'h13v03', 'h13v04', 'h14v01', 'h14v02', 'h14v03', 'h14v04', 'h15v01', 'h15v02', 'h15v03', 'h16v00', 'h16v01', 'h16v02', 'h17v00', 'h17v01', 'h17v02', 'h28v03', 'h29v03', 'h06v03', 'h07v03', 'h07v05', 'h07v06', 'h07v07', 'h08v03', 'h08v04', 'h08v05', 'h08v06', 'h08v07', 'h09v02', 'h09v03', 'h09v04', 'h09v05', 'h09v06', 'h09v07', 'h09v08')
# tiles <- namerica_tiles
# f <- function(x){
#   existing_file <- dir("/projectnb/modislc/users/joshgray/MCD12Q2C6/PixelQuants", pattern=paste(x, ".*365$", sep=""), full=T)
#   good_size <- 115200000
#   if(length(existing_file) < 1){
#     sys_cmd <- paste("qsub -V -pe omp 16 -l h_rt=02:00:00 -l mem_total=98G /projectnb/modislc/users/joshgray/MCD12Q2C6/subPixelQuants.sh", x)
#     print(paste("Resubmitting tile:", x))
#     system(sys_cmd)
#   }else if(!file.info(existing_file)$size == good_size){
#     # resub
#     sys_cmd <- paste("qsub -V -pe omp 16 -l h_rt=02:00:00 -l mem_total=98G /projectnb/modislc/users/joshgray/MCD12Q2C6/subPixelQuants.sh", x)
#     print(paste("Resubmitting tile:", x))
#     system(sys_cmd)
#   }
# }
# trash <- lapply(tiles, f)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# make continental mosaic of output
# tiles <- namerica_tiles
# in_files <- unlist(lapply(tiles, function(x) dir("/projectnb/modislc/users/joshgray/MCD12Q2C6/PixelQuants", pattern=paste(x, ".*365$", sep=""), full=T)))
# mosaic_file <- "/projectnb/modislc/users/joshgray/MCD12Q2C6/PixelQuants/NAmerica_mosaic.vrt"
# sys_cmd <- paste("gdalbuildvrt", mosaic_file, paste(in_files, collapse=" "))
# system(sys_cmd)
# gdalwarp -dstnodata 0 -s_srs "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m" -t_srs "+proj=latlong +datum=WGS84" -tr 0.025 0.025 -te -170 10 -50 75 /projectnb/modislc/users/joshgray/MCD12Q2C6/PixelQuants/NAmerica_mosaic.vrt /projectnb/modislc/users/joshgray/MCD12Q2C6/PixelQuants/NAmericaQuants_wgs.tif
#
# s <- stack("/projectnb/modislc/users/joshgray/MCD12Q2C6/PixelQuants/NAmericaQuants_wgs.tif")
# s[s<0]<-0
# amp <- raster(s,5)-raster(s,1)
# plotRGB(s,r=1,g=3,b=5,stretch="lin")


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# prelims
library(raster)
library(rgdal)
library(doParallel)
library(foreach)
library(argparse)
source("/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_Functions.R")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# parse the command line arguments
arg_parser <- ArgumentParser()
arg_parser$add_argument("-tile", type="character") # tile to process
arg_parser$add_argument("-start_date", type="character", default="2011-1-1")
arg_parser$add_argument("-end_date", type="character", default="2015-12-31")
arg_parser$add_argument("-out_dir", type="character", default="/projectnb/modislc/users/joshgray/MCD12Q2C6/PixelQuants")
arg_parser$add_argument("-num_cores", type="integer", default=16) # number of CPU cores to use for processing
args <- arg_parser$parse_args()

args$start_date <- as.Date(args$start_date)
args$end_date <- as.Date(args$end_date)

# DEBUG (makes it possible to grep the *.sh.o* files to identify a particular job)
print(paste("Processing tile", args$tile))

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Get input data_sets
tile_dsets <- GetTileDataSets(args$tile, start_date=args$start_date, end_date=args$end_date)
dates <- tile_dsets[[2]]
in_files <- tile_dsets[[1]][1:(length(dates) * 5)] # restrict to non-QA values
tmp_r <- raster(tile_dsets[[1]][1]) # temporary raster for sizing...as if it's not MODIS ;)
rows_to_do <- nrow(tmp_r) / args$num_cores
i_seq <- 1:(nrow(tmp_r) / rows_to_do)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# split row chunks across processors with library(foreach) and library(doParallel)
registerDoParallel(args$num_cores)
system.time(all_quants <- foreach(i=i_seq, .combine=rbind) %dopar% GetSnowFreeQuants(start_row=i, datasets=in_files, rows_to_do=rows_to_do))

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# check for proper size of all_quants
if(dim(all_quants)[1] == dim(tmp_r)[1]^2){
	# write to disk
	out_file <- file.path(args$out_dir, paste(args$tile, "SNOWFREE_EVI2_QUANTS", strftime(args$start_date, format="%Y%j"), strftime(args$end_date, format="%Y%j"), sep="_"))
	WriteQuants(out_file, c(all_quants), nbands=dim(all_quants)[2], tile=args$tile)
}else{
  # NOTE: could also resubmit here
	print("Failure!")
}
