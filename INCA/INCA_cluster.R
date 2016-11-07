# the sub_INCA_cluster.R script looks like this:
#
# #!/bin/csh
#
# source /usr/local/apps/R/R-312.csh
# R --vanilla < /share/jmgray2/INCA_cluster.R --args -tile $1


# To submit all tiles, we use this script:
#
# #!/bin/bash
#
# # submit all CONUS tiles as INCA tile summary jobs on henry2
# SCRATCHDIR="/share/jmgray2/"
# SCRIPTNAME="sub_INCA_cluster.sh"
# SCRATCHROOT="INCA"
#
# declare -a tiles=("h08v04" "h09v04" "h10v04" "h11v04" "h12v04" "h13v04" "h08v05" "h09v05" "h10v05" "h11v05" "h12v05" "h09v06" "h10v06")
# for i in "${tiles[@]}"
# do
#   bsub_cmd="bsub -q cnr -W 4:00 -n 16 -R \"oc span[ptile=16]\" -o ${SCRATCHDIR}${SCRATCHROOT}.$i.out.%J -e $SCRATCHDIR$SCRATCHROOT.$i.err.%J \"csh ${SCRATCHDIR}${SCRIPTNAME} $i\""
#   # echo $bsub_cmd
#   eval $bsub_cmd
# done

.libPaths("/home/jmgray2/R/x86_64-unknown-linux-gnu-library/3.1")
library(raster)
library(rgdal)
library(parallel)
library(argparse)
library(mblm)

#===============================================================================
# Functions
#===============================================================================

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
PhenoNormals <- function(x, years=NULL){
  # function takes a vector of DOY's (x) and input years (years) and returns:
  # x's median, MAD, average, std dev, Theil-Sen slope, T-S p-value
  # then the rounded median anomalies of x for all years, then that anomaly divided
  # by the MAD (how many median absolute deviations from the median is the year?)

  if(is.null(years)) years <- 1:length(x) # if no years are given, assume they are in the right order and have no gaps
  med <- median(x, na.rm=T)
  MAD <- mad(x, na.rm=T)
  avg <- mean(x, na.rm=T)
  stddev <- sd(x, na.rm=T)

  # calculate annual median anomalies: raw and as multiples of MAD
  ann.median.anoms <- round(x - med)
  ann.median.anoms.mad <- (x - med) / MAD

  # mblm can't handle missing data
  x.na <- x[!is.na(x)]
  years.na <- years[!is.na(x)]

  # fit the median based linear model (Theil-Sen slope estimator)
  library(mblm) # required for Theil-Sen
  mblm.fit <- try(mblm(x.na ~ years.na), silent=TRUE)
  if(inherits(mblm.fit, 'try-error')){
    ts.slope <- NA
    ts.pvalue <- NA
  }else{
    ts.vals <- summary(mblm.fit)$coeff[2, c(1, 4)]
    ts.slope <- ts.vals[1]
    ts.pvalue <- ts.vals[2]
  }

  return(c(round(med), MAD, round(avg), stddev, ts.slope, ts.pvalue, ann.median.anoms, ann.median.anoms.mad))
}

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
PhenoNormalsTile <- function(tile, cl, data_dir, output_dir){
  # soup-to-nuts processing of INCA variables for a single tile

  # first, export HDF as TIFF and calculate: ogi, half_spring, half_fall, dormancy, and gsl
  # NOTE: we don't need to do this, b/c we have the HDF4 driver on henry2
  # this is the leftover workaround that will work on the PC's
  mcd12q2_names <- dir(data_dir, patt=paste("MCD12.*", tile, ".*hdf$", sep=""), full=T)
  for(mcd12q2_name in mcd12q2_names){
    tmp_out_name <- paste("temp_out_", tile, ".tif", sep="")
    # gdal_translate(mcd12q2_name, dst_dataset = file.path(output_dir, tmp_out_name), sds=T)
    gdal_cmd <- paste("gdal_translate -sds", mcd12q2_name, file.path(output_dir, tmp_out_name))
    system(gdal_cmd)
    out_files <- dir(output_dir, patt=paste("temp_out_", tile, ".*tif$", sep=""), full=T)
    # out_files <- dir(output_dir, patt="temp_out.*.tif$", full=T)
    s <- stack(out_files[1:4])
    s <- subset(s, seq(1, nlayers(s), by=2)) # get rid of 2nd layer of each year
    gsl <- raster(s, 4) - raster(s, 1)
    half_spring <- (raster(s, 1) + raster(s, 2)) / 2
    half_fall <- (raster(s, 3) + raster(s, 4)) / 2
    annual_out_s <- stack(raster(s, 1), half_spring, half_fall, raster(s, 4), gsl)
    out_file <- file.path(output_dir, paste(paste(unlist(strsplit(basename(mcd12q2_name), split="\\."))[1:3], collapse="_"), "_INCA.tif", sep=""))
    writeRaster(annual_out_s, filename=out_file, overwrite=TRUE)
  }

  # now calculate the pheno normals
  metrics <- c("ogi", "halfspring", "halffall", "dormancy", "gsl")
  in_files <- dir(output_dir, patt=paste("MCD12.*", tile, ".*INCA.tif$", sep=""), full=T)
  in_years <- as.integer(gsub(".*A([0-9]{4})001.*$", "\\1", basename(in_files))) # get data years from file name
  doy_offset <- as.integer(as.Date(paste(in_years, "-1-1", sep="")) - as.Date("2000-1-1")) # Jan 1 of each data year as days since 2000-1-1
  pheno_s <- stack(in_files)

  # DEBUG: this is just b/c of weird job failures, do full processing for h11v05
  # but only gsl for any others
  mets_to_do <- 5
  if(tile == "h11v05") mets_to_do <- 1:5
  for(i in mets_to_do){
    tmp_s <- subset(pheno_s, seq(i, nlayers(pheno_s), by=5))
    # only need to subtract DOY offset from first four metrics!
    if(i != 5) tmp_s <- tmp_s - doy_offset # convert days since 2000-1-1 to DOY
    pheno_stack_v <- values(tmp_s)
    # system.time(pheno_output <- parApply(cl, pheno_stack_v[1:1e4,], 1, PhenoNormals)) # do just first 1e4 pixels, for testing
    system.time(pheno_output <- parApply(cl, pheno_stack_v, 1, PhenoNormals))
    output_s <- stack(raster(tmp_s, 1))
    # Definitely NOT the best way to do this, but sort of stumped and moving on...
    for(n in 1:(dim(pheno_output)[1] - 1)){
      output_s <- addLayer(output_s, raster(tmp_s, 1))
    }
    values(output_s) <- t(pheno_output)
    out_file_name <- file.path(output_dir, paste(paste("INCA_summary", tile, metrics[i], sep="_"), ".tiff", sep=""))
    writeRaster(file=out_file_name, output_s)
  }
}

#===============================================================================
# Do the actual processing...
#===============================================================================
# # Print start time
# start_time <- Sys.time()
# print(start_time)

# parse the command line arguments
arg_parser <- ArgumentParser()
arg_parser$add_argument("-tile", type="character") # tile to process
arg_parser$add_argument("-output_dir", type="character", default="/share/jmgray2/INCA/output") # output directory
arg_parser$add_argument("-data_dir", type="character", default="/share/jmgray2/MODIS/MCD12Q2") # input binary splined evi data directory
arg_parser$add_argument("-cluster_size", type="integer", default=16) # number of CPU cores to use for processing
args <- arg_parser$parse_args()
# args <- arg_parser$parse_args(c("-tile","h12v04")) # example to test parser

# make a cluster
cl <- makeCluster(args$cluster_size)

# compute INCA summaries
system.time(PhenoNormalsTile(tile=args$tile, cl=cl, data_dir=args$data_dir, output_dir=args$output_dir))

# # Print end time
# end_time <- Sys.time()
# print(end_time)
# print(end_time - start_time)
