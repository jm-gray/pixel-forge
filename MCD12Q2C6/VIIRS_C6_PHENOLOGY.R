#============================================================
#============================================================
# Read in xyz's VIIRS NBAR EVI2 data
# The data are 3-day composites in BIP format. The band 1 is EVI2, band2 is QA, and band3 is NDVI.
# QA=0   good quality
# QA=1 snow flag
# QA=9 fill value
# QA=other values   other quality.
# scaled by 1e4

#============================================================
#============================================================
# preliminaries and constants
library(raster)
library(rgdal)
library(parallel)
library(zoo)
library(doParallel)
library(foreach)
library(RColorBrewer)

nbands = 3 #
dsize = 2
s_flag = T
NA_VALUE <- 32767
NUM_PIXELS <- 5760000
YEARS_TO_DO <- 3
TILE <- "H12V04"
# DATA_DIR <- "/Users/joshuagray/Desktop/new_xyz_viirs/VIIRS_3daycomposite"
DATA_DIR <- "/projectnb/modislc/users/joshgray/xyz_pheno_comparison/VIIRS_NBAR/VIIRS_3daycomposite"
SCALE_FACTOR <- 1e4
OUT_FILE <- "/projectnb/modislc/users/joshgray/xyz_pheno_comparison/VIIRS_NBAR/C6_RESULTS/h12v04_second_run_output.Rdata"

#============================================================
#============================================================
# Functions
# function from DSM to read a binary file
read_int_file <- function(in_file, dsize, nbands, s_flag=T, bsq_flag=F){
  f_size <- file.info(in_file)$size # get the total dimensions of file = nrow*ncol*dsize*nbands
  f <- file(in_file,"rb") # open the file
  tot_size <- (f_size / dsize) # dsize is typically 2 for integers
  temp <- readBin(f, integer(), n=tot_size, endian= "little", size=dsize, signed=s_flag) # read all of the data into an array
  close(f) # close the file
  # re-order the temp array into a matrix; if bsq_flag then read by row, otherwise by column
  ifelse(bsq_flag, byr_flag<-FALSE, byr_flag<-TRUE)
  temp <- matrix(temp, ncol=nbands, byrow=byr_flag)
  return(temp) # return the data
}

#============================================================
#============================================================
# get input files
evi_in_files <- dir(DATA_DIR, pattern=paste("VIIRS.*", TILE, ".*BIP", sep=""), full=T)
file_years <- as.integer(unlist(lapply(evi_in_files, function(x) as.numeric(substr(unlist(strsplit(basename(x), split="\\."))[2], 1, 4)))))

start_year <- min(file_years)
end_year <- start_year + (YEARS_TO_DO - 1)

evi_dates <- as.Date(unlist(lapply(evi_in_files, function(x) as.Date(paste(as.numeric(substr(unlist(strsplit(basename(x), split="\\."))[2], 1, 4)), as.integer(substr(unlist(strsplit(basename(x), split="_"))[3], 1, 3)), sep="-"), format="%Y-%j"))), origin="1970-1-1")

#============================================================
#============================================================
# set up the cluster
source("/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_Functions.R")
# source("/Users/joshuagray/Dropbox/RCode/c6code/MCD12Q2C6_Functions.R")
cl <- makeCluster(16)
clusterEvalQ(cl, {source("/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_Functions.R")})
clusterExport(cl, c("evi_dates"))

#============================================================
#============================================================
# set up the data arrays
evi_mat <- matrix(NA, nrow=NUM_PIXELS, ncol=length(evi_in_files))
qa_mat <- matrix(NA, nrow=NUM_PIXELS, ncol=length(evi_in_files))

# loop through files and add to evi_mat array
n <- 1 # column counter
for(in_file in evi_in_files){
  print(paste("Working on file:", in_file))
  tmp <- read_int_file(in_file, dsize=2, nbands=3, s_flag=T, bsq_flag=F)
  tmp[tmp == NA_VALUE] <- NA
  # tmp <- tmp / SCALE_FACTOR
  evi_mat[, n] <- tmp[, 1] / SCALE_FACTOR
  qa_mat[, n] <- tmp[, 2]
  n <- n + 1 # increment column counter
}

# append evi and qa information
data_mat <- cbind(evi_mat, qa_mat)
rm(evi_mat, qa_mat) # clean-up

#============================================================
#============================================================
# do phenology for the tile
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
	gup_threshes=c(0.1, 0.2 ,0.5, 1),
	gdown_threshes=c(0.5, 0.2, 0.1),
	max_snow_fill_ratio=1.25,
  snow_free_min_quant=0.1,
  snow_max_quant=0.98
)

#============================================================
#============================================================
# timing single vs multithreaded:
# system.time(mean_pheno <- parApply(cl, data_mat[1:10000,], 1, DoPhenologyVIIRS, dates=evi_dates, pheno_pars=pheno_pars))
# system.time(mean_pheno <- apply(data_mat[1:10000,], 1, DoPhenologyVIIRS, dates=evi_dates, pheno_pars=pheno_pars))
# object.size(x=lapply(ls(), get)) and print(object.size(x=lapply(ls(), get)), units="Mb")

# do the entire thing
# system.time(pheno_results <- parApply(cl, data_mat, 1, DoPhenologyVIIRS, dates=evi_dates, pheno_pars=pheno_pars))
# save the output
# save(pheno_results, file=OUT_FILE)

# single-thread, chunked up
# system.time(pheno_results1 <- apply(data_mat[1:1e6,], 1, DoPhenologyVIIRS, dates=evi_dates, pheno_pars=pheno_pars))
# gc()
# system.time(pheno_results2 <- apply(data_mat[1000001:2e6,], 1, DoPhenologyVIIRS, dates=evi_dates, pheno_pars=pheno_pars))
# gc()
# system.time(pheno_results3 <- apply(data_mat[2000001:3e6,], 1, DoPhenologyVIIRS, dates=evi_dates, pheno_pars=pheno_pars))
# gc()
# system.time(pheno_results4 <- apply(data_mat[3000001:4e6,], 1, DoPhenologyVIIRS, dates=evi_dates, pheno_pars=pheno_pars))
# gc()
# system.time(pheno_results5 <- apply(data_mat[4000001:5e6,], 1, DoPhenologyVIIRS, dates=evi_dates, pheno_pars=pheno_pars))
# gc()
# system.time(pheno_results6 <- apply(data_mat[5000001:NUM_PIXELS,], 1, DoPhenologyVIIRS, dates=evi_dates, pheno_pars=pheno_pars))
# # save
# save(pheno_results1, pheno_results2, pheno_results3, pheno_results4, pheno_results5, pheno_results6, file=OUT_FILE)
# save(pheno_results1, pheno_results2, pheno_results3, file="/projectnb/modislc/users/joshgray/xyz_pheno_comparison/VIIRS_NBAR/C6_RESULTS/h12v04_first_three_chunks.Rdata")

#===========================================
# multithreaded, chunked up
system.time(pheno_results1 <- parApply(cl, data_mat[1:1e6,], 1, DoPhenologyVIIRS, dates=evi_dates, pheno_pars=pheno_pars))
gc()
system.time(pheno_results2 <- parApply(cl, data_mat[1000001:2e6,], 1, DoPhenologyVIIRS, dates=evi_dates, pheno_pars=pheno_pars))
gc()
system.time(pheno_results3 <- parApply(cl, data_mat[2000001:3e6,], 1, DoPhenologyVIIRS, dates=evi_dates, pheno_pars=pheno_pars))
gc()
system.time(pheno_results4 <- parApply(cl, data_mat[3000001:4e6,], 1, DoPhenologyVIIRS, dates=evi_dates, pheno_pars=pheno_pars))
gc()
system.time(pheno_results5 <- parApply(cl, data_mat[4000001:5e6,], 1, DoPhenologyVIIRS, dates=evi_dates, pheno_pars=pheno_pars))
gc()
system.time(pheno_results6 <- parApply(cl, data_mat[5000001:NUM_PIXELS,], 1, DoPhenologyVIIRS, dates=evi_dates, pheno_pars=pheno_pars))

pheno_results <- append(pheno_results1, pheno_results2)
pheno_results <- append(pheno_results, pheno_results3)
pheno_results <- append(pheno_results, pheno_results4)
pheno_results <- append(pheno_results, pheno_results5)
pheno_results <- append(pheno_results, pheno_results6)
# save(pheno_results1, pheno_results2, pheno_results3, pheno_results4, pheno_results5, pheno_results6, file=OUT_FILE)
save(pheno_results, file=OUT_FILE)


# put the output into a raster:
# r_evi2 <- raster(nrow=2400, ncol=2400)
# r_qa <- raster(nrow=2400, ncol=2400)
# r_ndvi <- raster(nrow=2400, ncol=2400)
# tmp <- read_int_file(evi_in_files[170], dsize=2, nbands=3, s_flag=T, bsq_flag=F)
# tmp[tmp==32767] <- NA
# values(r_evi2) <- tmp[, 1]
# values(r_qa) <- tmp[, 2]
# values(r_ndvi) <- tmp[, 3]

#============================================================
#============================================================
# Extract
args <- list(
  tile=TILE,
  in_dir="/projectnb/modislc/users/joshgray/xyz_pheno_comparison/VIIRS_NBAR/C6_RESULTS",
  out_dir="/projectnb/modislc/users/joshgray/xyz_pheno_comparison/VIIRS_NBAR/C6_RESULTS",
  data_dir="/projectnb/modislc/data/mcd12_in/c5/mcd43a4",
  yearrange=c(2012, 2014),
  data_prefix="",
  i=1:11,
  i_key=3,
  prefix=c("tengup", "twentygup", "fiftygup", "peak", "fiftygdown", "twentygdown", "tengdown", "evi2min", "evi2max", "evi2int", "guprsq"),
  n=1,
  num_nodes=8,
  annual=TRUE,
  overwrite=TRUE
)

#------------------------------------------
#------------------------------------------
# get a temporary RasterLayer for this tile in order to set extent, geo, etc
# NOTE: this should search for a SINGLE random file, and repeat until it finds one...
nbar_files <- Sys.glob(file.path(args$data_dir, "*", "*", paste("*", tolower(args$tile), "*", sep="")))
tmp_r <- raster(Get_SDS_name_band1(nbar_files[1]))

#------------------------------------------
#------------------------------------------
# Extract each requested threshold date for each year and write to a multi-layer raster file
# year range
# cl <- makeCluster(args$num_nodes) # make cluster

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
		pheno_dates <- unlist(parLapply(cl, pheno_results, ExtractSingleDate, i=i, i_key=args$i_key, start_date=start_date, end_date=end_date, n=args$n))

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
