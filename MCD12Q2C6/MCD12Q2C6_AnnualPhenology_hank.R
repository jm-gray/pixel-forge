# MCD12Q2C6 Code adapted for use on henry2 cluster
#
# Example bsub submission:
# bsub -q cnr -W 6:00 -n 16 -R "span[ptile=16]" -o /rsstu/users/j/jmgray2/SEAL/hank_code/Pheno_par1_h08v05_2009.out.%J -e /rsstu/users/j/jmgray2/SEAL/hank_code/Pheno_par1_h08v05_2009.err.%J R CMD BATCH --vanilla '--args -tile h12v01 -year 2003 -out_prefix par1 -params /rsstu/users/j/jmgray2/SEAL/hank_code/PhenoPars1.txt' MCD12Q2C6_AnnualPhenology_hank.R Pheno_par1_h08v05_2009.log
#
# bsub -q cnr -W 0:05 -o /rsstu/users/j/jmgray2/SEAL/hank_code/Test_o.out.%J -e /rsstu/users/j/jmgray2/SEAL/hank_code/Test_o.err.%J R CMD BATCH --vanilla '--args -tile h08v05 -year 2009 -params /rsstu/users/j/jmgray2/SEAL/hank_code/PhenoPars1.txt' Test.R Test_o.log

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# prelims
library(parallel)
library(argparse)
library(rgdal)
source("/rsstu/users/j/jmgray2/SEAL/hank_code/MCD12Q2C6_AnnualPhenologyFunctions_hank.R")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# parse the command line arguments
arg_parser <- ArgumentParser()
arg_parser$add_argument("-tile", type="character") # tile to process
arg_parser$add_argument("-year", type="integer") # year of interest
arg_parser$add_argument("-pheno_period_start", type="character") # limit analysis to cycles between pheno_period_start/end
arg_parser$add_argument("-pheno_period_end", type="character")
arg_parser$add_argument("-out_dir", type="character", default="/rsstu/users/j/jmgray2/SEAL/MCD12Q2C6_output") # output directory
arg_parser$add_argument("-out_prefix", type="character", default="MCD12Q2C6") # prefix for output files
arg_parser$add_argument("-data_dir", type="character", default="/rsstu/users/j/jmgray2/SEAL/") # input binary splined evi data directory
arg_parser$add_argument("-data_prefix", type="character", default="MCD12I") # prefix of input binary files
arg_parser$add_argument("-params", type="character") # phenology processing parameters file
arg_parser$add_argument("-chunk_line_size", type="integer", default=480) # lines to process at once; consider memory; default does tile in 5 equal chunks; NOTE: limited to factors of 2400
arg_parser$add_argument("-cluster_size", type="integer", default=16) # number of CPU cores to use for processing

args <- arg_parser$parse_args()
# args <- arg_parser$parse_args(c("-tile","h12v04","-year",2003))

# if a parameter file is specified, retrieve the parameters from there, use defaults otherwise
if(is.null(args$params)){
  pheno_pars <- DefaultPhenoParameters()
}else{
  pheno_pars <- ReadPhenoParameters(args$params)
}

# if the pheno_period_start and pheno_period_end are not given as char strings (e.g. "2007-1-1")
# then we assume it is the calendar year of interest
if(is.null(args$pheno_period_start) | is.null(args$pheno_period_end)){
  args$pheno_period_start <- as.numeric(as.Date(paste(args$year, "-1-1", sep="")))
  args$pheno_period_end <- as.numeric(as.Date(paste(args$year, "-12-31", sep="")))
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# DEBUG: allows us to grep the *.sh.o* files to find a particular tile-year job
print(paste("Processing tile", args$tile, "for year", args$year))
print(paste(names(pheno_pars), pheno_pars, sep="=", collapse=", ")) # to record the pheno parameters used for processing

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# set up the cluster
cl <- makeCluster(args$cluster_size)
clusterEvalQ(cl, {source("/rsstu/users/j/jmgray2/SEAL/hank_code/MCD12Q2C6_AnnualPhenologyFunctions_hank.R"); library(rgdal); library(raster)})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# get all data files for tile
patt <- paste(args$data_prefix, ".*", args$tile, ".*bip$", sep="")
evi2_files <- dir(file.path(args$data_dir, "MCD12I2"), pattern=patt, full=T)
resid_files <- dir(file.path(args$data_dir, "MCD12I4"), pattern=patt, full=T)
snow_files <- dir(file.path(args$data_dir, "MCD12I3"), pattern=patt, full=T)

# create 3 years worth of daily dates, IGNORING LEAP YEARS!
c6_dates <- as.Date(paste(rep((args$year - 1):(args$year + 1), each=365), rep(1:365, 3), sep="-"), format="%Y-%j")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# do the phenology and write to output in chunks of an equal number of lines
for(i in 1:(2400 / args$chunk_line_size)){
	# load the data
	print(paste("Loading chunk", i))
	start_time <- Sys.time()
	start_line <- ((i - 1) * args$chunk_line_size) + 1 # start reading data chunk at this line
	tmp <- Get3YearDataChunk(evi2_files, resid_files, snow_files, args$year, start_line, args$chunk_line_size)
	print(Sys.time() - start_time)

  # apply the phenology algorithm
	print(paste("Retrieving phenology of chunk", i))
  system.time(pheno_values <- parApply(cl, tmp, 1, AnnualPhenologyC6, dates=c6_dates, pheno_pars=pheno_pars, pheno_period_start=args$pheno_period_start, pheno_period_end=args$pheno_period_end))
	print(Sys.time() - start_time)

  # convert to vector and write to file
	print(paste("Writing chunk", i))
	out_file <- file.path(args$out_dir, paste(args$out_prefix, args$tile, args$year, sep="_"))
	WritePhenologyData(out_file=out_file, pheno_values_bip=c(pheno_values), tile=args$tile, nbands=dim(pheno_values)[1])
	print(Sys.time() - start_time)

	# clean up to free memory; NOTE: maybe not useful...
	rm(tmp)
	rm(pheno_values)
	rm(out_file)
}
