# use INCA_suball.sh to submit all jobs

# library(raster)
library(raster, lib.loc="/home/jmgray2/libs/R_libs")
# library(rgdal)
library(rgdal, lib.loc="/home/jmgray2/libs/R_libs")
library(parallel)
# library(argparse)
library(argparse, lib.loc="/home/jmgray2/libs/R_libs")
# library(mblm)
library(mblm, lib.loc="/home/jmgray2/libs/R_libs")


#===============================================================================
# Functions
#===============================================================================

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
PhenoNormals_full <- function(x, years=NULL){
  # NOTE: NOT UPDATED LIKE PHENONORMALS() TO HANDLE ROUND_VALS!!!
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
PhenoNormals <- function(x, years=NULL, round_vals=T){
  # function takes a vector of DOY's (x) and input years (years) and returns:
  # x's median, MAD, average, std dev, Theil-Sen slope, T-S p-value, number of years

  if(is.null(years)) years <- 1:length(x) # if no years are given, assume they are in the right order and have no gaps
  med <- median(x, na.rm=T)
  MAD <- mad(x, na.rm=T)
  avg <- mean(x, na.rm=T)
  stddev <- sd(x, na.rm=T)
  num_vals <- sum(!is.na(x))

  # calculate annual median anomalies: raw and as multiples of MAD
  # ann.median.anoms <- round(x - med)
  # ann.median.anoms.mad <- (x - med) / MAD

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

  # return(c(round(med), MAD, round(avg), stddev, ts.slope, ts.pvalue, ann.median.anoms, ann.median.anoms.mad))
  if(round_vals){
    return(c(round(med), MAD, round(avg), stddev, ts.slope, ts.pvalue, num_vals))
  }else{
    return(c(med, MAD, avg, stddev, ts.slope, ts.pvalue, num_vals))
  }
  
}

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
GetSDS <- function(file_path, sds=NULL){
  # gets the SDS names for EOS HDF data access
  # valid SDS: Greenup, MidGreenup, Peak, Senescence, MidGreendown, Dormancy, EVI_Minimum, EVI_Amplitude, NumCycles, QA_Detailed, QA_Overall
  all_sds <- c("NumCycles", "Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy", "EVI_Minimum", "EVI_Amplitude", "EVI_Area", "QA_Overall", "QA_Detailed")
  if(is.null(sds)){
    return(paste("HDF4_EOS:EOS_GRID:\"", file_path, "\":MCD12Q2:", all_sds, sep = ""))
  }else{
    return(paste("HDF4_EOS:EOS_GRID:\"", file_path, "\":MCD12Q2:", sds, sep = ""))
  }
}

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
PhenoNormalsTile <- function(tile, cl, data_dir, output_dir, doGSL=FALSE){
  # soup-to-nuts processing of INCA variables for a single tile
  # added doGSL flag: if this is TRUE, then ONLY growing season length summaries are calculated
  # these are calculated as MidGreenup to MidGreendown.

  # the collection of SDS to calculate normals and trends over
  metrics_to_do <- c("Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy", "EVI_Minimum", "EVI_Amplitude", "EVI_Area")
  doy_metric <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) # is this a DOY value?

  # gather and sort input files
  in_files <- dir(data_dir, patt=paste("MCD12.*", tile, ".*hdf$", sep=""), full=T, rec=T)
  in_years <- as.integer(gsub(".*A([0-9]{4})001.*$", "\\1", basename(in_files))) # get data years from file name
  in_files <- in_files[order(in_years)]
  in_years <- sort(in_years)

  # calculate offset to convert from days since 1970-1-1 to DOY
  doy_offset <- as.integer(as.Date(paste(in_years, "-1-1", sep="")) - as.Date("1970-1-1")) # Jan 1 of each data year as days since 2000-1-1

  # if doing GSL
  if(doGSL){
    # get MidGreenup and MidGreendown SDS, calculate their difference in days
    midgup_s <- stack(GetSDS(in_files, "MidGreenup"))
    midgup_s <- subset(midgup_s, seq(1, nlayers(midgup_s), by=2))
    NAvalue(midgup_s) <- 32767
    midgup_s[midgup_s == -32768] <- NA
    midgdown_s <- stack(GetSDS(in_files, "MidGreendown"))
    midgdown_s <- subset(midgdown_s, seq(1, nlayers(midgdown_s), by=2))
    NAvalue(midgdown_s) <- 32767
    midgdown_s[midgdown_s == -32768] <- NA
    gsl_s <- midgdown_s - midgup_s
    gsl_v <- values(gsl_s)
    
    # calculate pheno normals, trends, and anomalies
    pheno_output_v <- parApply(cl, gsl_v, 1, PhenoNormals, round_vals=TRUE)

    # create a blank output raster, assign the values, and write to a file
    pheno_output_s <- do.call(stack, replicate(dim(pheno_output_v)[1], raster(midgup_s, 1)))
    values(pheno_output_s) <- t(pheno_output_v)
    out_file_name <- file.path(output_dir, paste("INCA", tile, "GSL", "tiff", sep="."))
    # writeRaster(file=out_file_name, pheno_output_s)
    writeRaster(file=out_file_name, pheno_output_s, NAflag=32767, overwrite=T)

    # testing...
    # data_dir <- "/rsstu/users/j/jmgray2/SEAL/INCA/MCD12Q2C6/MCD12Q2"
    # tile <- "h12v04"
    # in_files <- dir(data_dir, patt=paste("MCD12.*", tile, ".*hdf$", sep=""), full=T, rec=T)
    # in_years <- as.integer(gsub(".*A([0-9]{4})001.*$", "\\1", basename(in_files))) # get data years from file name
    # in_files <- in_files[order(in_years)]
    # in_years <- sort(in_years)
    # set.seed(42)
    # x <- sample(1:nrow(gsl_v), 10)
    # round_flag=TRUE
    # pheno_output_v <- apply(gsl_v[x,], 1, PhenoNormals, round_vals=round_flag)
    # system.time(pheno_output_v <- parApply(cl, gsl_v, 1, PhenoNormals, round_vals=round_flag))
    
  }else{
    for(this_metric in metrics_to_do){
      # get a RasterStack for all years of this metric
      pheno_s <- stack(GetSDS(in_files, this_metric))
      # subset to first layer in each year
      pheno_s <- subset(pheno_s, seq(1, nlayers(pheno_s), by=2))
      # set NA value
      NAvalue(pheno_s) <- 32767
      # some tiles also have -32768 as an out-of-range value
      pheno_s[pheno_s == -32768] <- NA

      # convert to DOY if necessary, and set rounding flag
      if(doy_metric[which(metrics_to_do == this_metric)]){ 
        pheno_s <- pheno_s - doy_offset
        round_flag <- TRUE
      }else{
        round_flag <- FALSE
      }

      # get the raw values in a matrix
      pheno_v <- values(pheno_s)

      # calculate pheno normals, trends, and anomalies
      pheno_output_v <- parApply(cl, pheno_v, 1, PhenoNormals, round_vals=round_flag)

      # create a blank output raster, assign the values, and write to a file
      pheno_output_s <- do.call(stack, replicate(dim(pheno_output_v)[1], raster(pheno_s, 1)))
      values(pheno_output_s) <- t(pheno_output_v)
      out_file_name <- file.path(output_dir, paste("INCA", tile, this_metric, "tiff", sep="."))
      # writeRaster(file=out_file_name, pheno_output_s)
      writeRaster(file=out_file_name, pheno_output_s, NAflag=32767, overwrite=T)
    }
  }
}

#===============================================================================
# Do the actual processing...
#===============================================================================
# parse the command line arguments
arg_parser <- ArgumentParser()
arg_parser$add_argument("-tile", type="character") # tile to process
# arg_parser$add_argument("-output_dir", type="character", default="/rsstu/users/j/jmgray2/SEAL/INCA/INCAoutput") # output directory
arg_parser$add_argument("-output_dir", type="character", default="/rsstu/users/j/jmgray2/SEAL/INCA/INCAglobaloutput") # output directory
arg_parser$add_argument("-data_dir", type="character", default="/rsstu/users/j/jmgray2/SEAL/INCA/MCD12Q2C6/MCD12Q2") # input binary splined evi data directory
# arg_parser$add_argument("-cluster_size", type="integer", default=16) # number of CPU cores to use for processing
arg_parser$add_argument("-cluster_size", type="integer", default=8) # number of CPU cores to use for processing
args <- arg_parser$parse_args()
# args <- arg_parser$parse_args(c("-tile","h12v04")) # example to test parser

# make a cluster
cl <- makeCluster(args$cluster_size)
clusterExport(cl, c("GetSDS", "PhenoNormals"))
# clusterEvalQ(cl, {library(raster); library(rgdal); library(mblm)})
clusterEvalQ(cl, {library(raster, lib.loc="/home/jmgray2/libs/R_libs"); library(rgdal, lib.loc="/home/jmgray2/libs/R_libs"); library(mblm, lib.loc="/home/jmgray2/libs/R_libs")})

# compute INCA summaries
# system.time(PhenoNormalsTile(tile=args$tile, cl=cl, data_dir=args$data_dir, output_dir=args$output_dir))

# do it just for GSL
system.time(PhenoNormalsTile(tile=args$tile, cl=cl, data_dir=args$data_dir, output_dir=args$output_dir, doGSL=TRUE))

# # use INCA_suball.sh to submit all jobs

# library(raster)
# library(rgdal)
# library(parallel)
# library(argparse)
# library(mblm)

# #===============================================================================
# # Functions
# #===============================================================================

# #=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
# PhenoNormals <- function(x, years=NULL){
#   # function takes a vector of DOY's (x) and input years (years) and returns:
#   # x's median, MAD, average, std dev, Theil-Sen slope, T-S p-value
#   # then the rounded median anomalies of x for all years, then that anomaly divided
#   # by the MAD (how many median absolute deviations from the median is the year?)

#   if(is.null(years)) years <- 1:length(x) # if no years are given, assume they are in the right order and have no gaps
#   med <- median(x, na.rm=T)
#   MAD <- mad(x, na.rm=T)
#   avg <- mean(x, na.rm=T)
#   stddev <- sd(x, na.rm=T)

#   # calculate annual median anomalies: raw and as multiples of MAD
#   ann.median.anoms <- round(x - med)
#   ann.median.anoms.mad <- (x - med) / MAD

#   # mblm can't handle missing data
#   x.na <- x[!is.na(x)]
#   years.na <- years[!is.na(x)]

#   # fit the median based linear model (Theil-Sen slope estimator)
#   library(mblm) # required for Theil-Sen
#   mblm.fit <- try(mblm(x.na ~ years.na), silent=TRUE)
#   if(inherits(mblm.fit, 'try-error')){
#     ts.slope <- NA
#     ts.pvalue <- NA
#   }else{
#     ts.vals <- summary(mblm.fit)$coeff[2, c(1, 4)]
#     ts.slope <- ts.vals[1]
#     ts.pvalue <- ts.vals[2]
#   }

#   return(c(round(med), MAD, round(avg), stddev, ts.slope, ts.pvalue, ann.median.anoms, ann.median.anoms.mad))
# }

# #=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
# GetSDS <- function(file_path, sds=NULL){
#   # gets the SDS names for EOS HDF data access
#   # valid SDS: Greenup, MidGreenup, Peak, Senescence, MidGreendown, Dormancy, EVI_Minimum, EVI_Amplitude, NumCycles, QA_Detailed, QA_Overall
#   all_sds <- c("NumCycles", "Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy", "EVI_Minimum", "EVI_Amplitude", "EVI_Area", "QA_Overall", "QA_Detailed")
#   if(is.null(sds)){
#     return(paste("HDF4_EOS:EOS_GRID:\"", file_path, "\":MCD12Q2:", all_sds, sep = ""))
#   }else{
#     return(paste("HDF4_EOS:EOS_GRID:\"", file_path, "\":MCD12Q2:", sds, sep = ""))
#   }
# }

# #=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
# PhenoNormalsTile <- function(tile, cl, data_dir, output_dir){
#   # soup-to-nuts processing of INCA variables for a single tile

#   # the collection of SDS to calculate normals and trends over
#   metrics_to_do <- c("Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy", "EVI_Minimum", "EVI_Amplitude", "EVI_Area")
#   doy_metric <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) # is this a DOY value?

#   # gather and sort input files
#   in_files <- dir(data_dir, patt=paste("MCD12.*", tile, ".*hdf$", sep=""), full=T, rec=T)
#   in_years <- as.integer(gsub(".*A([0-9]{4})001.*$", "\\1", basename(in_files))) # get data years from file name
#   in_files <- in_files[order(in_years)]
#   in_years <- sort(in_years)

#   # calculate offset to convert from days since 1970-1-1 to DOY
#   doy_offset <- as.integer(as.Date(paste(in_years, "-1-1", sep="")) - as.Date("1970-1-1")) # Jan 1 of each data year as days since 2000-1-1

#   for(this_metric in metrics_to_do){
#     # get a RasterStack for all years of this metric
#     pheno_s <- stack(GetSDS(in_files, this_metric))
#     # subset to first layer in each year
#     pheno_s <- subset(pheno_s, seq(1, nlayers(pheno_s), by=2))
#     # set NA value
#     NAvalue(pheno_s) <- 32767
#     # some tiles also have -32768 as an out-of-range value
#     pheno_s[pheno_s == -32768] <- NA

#     # convert to DOY if necessary
#     if(doy_metric[which(metrics_to_do == this_metric)]){
#       pheno_s <- pheno_s - doy_offset
#     }

#     # get the raw values in a matrix
#     pheno_v <- values(pheno_s)

#     # calculate pheno normals, trends, and anomalies
#     pheno_output_v <- parApply(cl, pheno_v, 1, PhenoNormals)

#     # create a blank output raster, assign the values, and write to a file
#     pheno_output_s <- do.call(stack, replicate(dim(pheno_output_v)[1], raster(pheno_s, 1)))
#     values(pheno_output_s) <- t(pheno_output_v)
#     out_file_name <- file.path(output_dir, paste("INCA", tile, this_metric, "tiff", sep="."))
#     # writeRaster(file=out_file_name, pheno_output_s)
#     writeRaster(file=out_file_name, pheno_output_s, NAflag=32767)
#   }
# }

# #===============================================================================
# # Do the actual processing...
# #===============================================================================
# # parse the command line arguments
# arg_parser <- ArgumentParser()
# arg_parser$add_argument("-tile", type="character") # tile to process
# arg_parser$add_argument("-output_dir", type="character", default="/rsstu/users/j/jmgray2/SEAL/INCA/INCAoutput") # output directory
# arg_parser$add_argument("-data_dir", type="character", default="/rsstu/users/j/jmgray2/SEAL/INCA/MCD12Q2C6/MCD12Q2") # input binary splined evi data directory
# arg_parser$add_argument("-cluster_size", type="integer", default=16) # number of CPU cores to use for processing
# args <- arg_parser$parse_args()
# # args <- arg_parser$parse_args(c("-tile","h12v04")) # example to test parser

# # make a cluster
# cl <- makeCluster(args$cluster_size)
# clusterExport(cl, c("GetSDS", "PhenoNormals"))
# clusterEvalQ(cl, {library(raster); library(rgdal); library(mblm)})

# # compute INCA summaries
# system.time(PhenoNormalsTile(tile=args$tile, cl=cl, data_dir=args$data_dir, output_dir=args$output_dir))
