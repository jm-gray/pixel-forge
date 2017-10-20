#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Josh Gray, Boston University, 2016
#
# environment includes the R_earth/3.1 module
# NOTE: scale factor for evi_min, evi_amp, frac_filled_gup, and frac_filled_gdown and all qual scores are 1e4, but evi_amp is 10!
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Submitted to the GridEngine job handler w/ shell script (subAnnualPheno.sh):
# #!/bin/bash
# echo Submitting tile: $1
# echo Processing year: $2
# R --vanilla < /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_AnnualPhenology.R --args -tile $1 -year $2

# example submission command using default parameters:
# qsub -V -pe omp 16 -l h_rt=04:00:00 -l mem_total=98G subAnnualPheno.sh h11v04 2003

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Batch submission codes:
# asia_tiles <- c('h10v02', 'h11v02', 'h12v01', 'h19v00', 'h19v01', 'h19v04', 'h20v01', 'h20v02', 'h20v03', 'h20v04', 'h20v05', 'h21v01', 'h21v02', 'h21v03', 'h21v04', 'h21v05', 'h21v06', 'h21v07', 'h22v01', 'h22v02', 'h22v03', 'h22v04', 'h22v05', 'h22v06', 'h22v07', 'h23v01', 'h23v02', 'h23v03', 'h23v04', 'h23v05', 'h23v06', 'h23v07', 'h24v02', 'h24v03', 'h24v04', 'h24v05', 'h24v06', 'h24v07', 'h25v02', 'h25v03', 'h25v04', 'h25v05', 'h25v06', 'h25v07', 'h25v08', 'h25v09', 'h26v02', 'h26v03', 'h26v04', 'h26v05', 'h26v06', 'h26v07', 'h26v08', 'h27v03', 'h27v04', 'h27v05', 'h27v06', 'h27v07', 'h27v08', 'h27v09', 'h28v10', 'h28v04', 'h28v05', 'h28v06', 'h28v07', 'h28v08', 'h28v09', 'h29v10', 'h29v05', 'h29v06', 'h29v07', 'h29v08', 'h29v09', 'h30v10', 'h30v07', 'h30v08', 'h30v09', 'h31v09', 'h32v10', 'h32v09')
# namerica_tiles <- c('h10v02', 'h10v03', 'h10v04', 'h10v05', 'h10v06', 'h10v07', 'h10v08', 'h11v02', 'h11v03', 'h11v04', 'h11v05', 'h11v06', 'h11v07', 'h12v01', 'h12v02', 'h12v03', 'h12v04', 'h12v05', 'h12v07', 'h13v01', 'h13v02', 'h13v03', 'h13v04', 'h14v01', 'h14v02', 'h14v03', 'h14v04', 'h15v01', 'h15v02', 'h15v03', 'h16v00', 'h16v01', 'h16v02', 'h17v00', 'h17v01', 'h17v02', 'h28v03', 'h29v03', 'h06v03', 'h07v03', 'h07v05', 'h07v06', 'h07v07', 'h08v03', 'h08v04', 'h08v05', 'h08v06', 'h08v07', 'h09v02', 'h09v03', 'h09v04', 'h09v05', 'h09v06', 'h09v07', 'h09v08')
# europe_tiles <- c('h15v05', 'h16v02', 'h16v05', 'h17v01', 'h17v02', 'h17v03', 'h17v04', 'h17v05', 'h18v00', 'h18v01', 'h18v02', 'h18v03', 'h18v04', 'h18v05', 'h19v00', 'h19v01', 'h19v02', 'h19v03', 'h19v04', 'h19v05', 'h20v01', 'h20v02', 'h20v03', 'h20v04', 'h20v05', 'h21v03', 'h21v04')
# africa_tiles <- c('h15v07', 'h16v05', 'h16v06', 'h16v07', 'h16v08', 'h17v10', 'h17v05', 'h17v06', 'h17v07', 'h17v08', 'h18v05', 'h18v06', 'h18v07', 'h18v08', 'h18v09', 'h19v10', 'h19v11', 'h19v12', 'h19v05', 'h19v06', 'h19v07', 'h19v08', 'h19v09', 'h20v10', 'h20v11', 'h20v12', 'h20v05', 'h20v06', 'h20v07', 'h20v08', 'h20v09', 'h21v10', 'h21v11', 'h21v05', 'h21v06', 'h21v07', 'h21v08', 'h21v09', 'h22v10', 'h22v11', 'h22v07', 'h22v08', 'h22v09', 'h23v10', 'h23v11', 'h23v07', 'h23v08', 'h23v09')
# samerica_tiles <- c('h10v10', 'h10v07', 'h10v08', 'h10v09', 'h11v10', 'h11v11', 'h11v12', 'h11v07', 'h11v08', 'h11v09', 'h12v10', 'h12v11', 'h12v12', 'h12v13', 'h12v08', 'h12v09', 'h13v10', 'h13v11', 'h13v12', 'h13v13', 'h13v14', 'h13v08', 'h13v09', 'h14v10', 'h14v11', 'h14v14', 'h14v09', 'h08v08', 'h08v09', 'h09v08', 'h09v09')
# oceania_tiles <- c('h00v10', 'h00v08', 'h01v10', 'h01v11', 'h01v07', 'h01v09', 'h02v10', 'h02v06', 'h02v08', 'h27v10', 'h28v14', 'h29v13', 'h03v10', 'h03v11', 'h03v06', 'h03v07', 'h30v13', 'h31v12', 'h31v13', 'h31v08', 'h32v11', 'h32v12', 'h32v07', 'h32v09', 'h33v10', 'h33v11', 'h33v07', 'h33v08', 'h33v09', 'h34v10', 'h34v07', 'h34v08', 'h34v09', 'h35v10', 'h35v08', 'h35v09', 'h04v09', 'h05v13', 'h06v11', 'h08v11')
# australia_tiles <- c('h27v11', 'h27v12', 'h27v14', 'h28v11', 'h28v12', 'h28v13', 'h29v10', 'h29v11', 'h29v12', 'h29v13', 'h30v10', 'h30v11', 'h30v12', 'h31v10', 'h31v11', 'h31v12', 'h32v10')
# all_tiles <- scan("/projectnb/modislc/users/joshgray/MCD12Q2C6/gltiles.txt", what=character(), quiet=T)
#
# years <- 2003
# tiles <- africa_tiles
# submitTileYear <- function(x){
#   sys_cmd <- paste("qsub -V -pe omp 16 -l h_rt=08:00:00 -l mem_total=98G /projectnb/modislc/users/joshgray/MCD12Q2C6/subAnnualPheno.sh", x[1], x[2])
#   print(paste("Submitting", x[1], "for year", x[2]))
#   system(sys_cmd)
# }
#
# output_dir <- "/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT"
# checkExists <- function(tile, year, output_dir){
#   out_files <- dir(output_dir, paste(".*", tile, ".*", year, ".hdr$", sep=""))
#   if(length(out_files) == 0){
#     return(FALSE)
#   }else{
#     return(TRUE)
#   }
# }
#
# # check for existence
# df <- expand.grid(tiles, years)
# df$Exists <- apply(df, 1, checkExists, year=years, output_dir=output_dir)
# df_missing <- subset(df, !Exists)
#
# # submit all missing tile-years
# apply(df_missing, 1, submitTileYear)

# this is for convenience...
# layer_names <- c("num_cycles", "fill_code", "evi_area_cycle1", "evi_amp_cycle1", "evi_min_cycle1", "frac_filled_gup_cycle1", "frac_filled_gdown_cycle1", "length_gup_cycle1", "length_gdown_cycle1", "ogi_cycle1", "midgup_cycle1", "mat_cycle1", "peak_cycle1", "sen_cycle1", "midgdown_cycle1", "dor_cycle1", "ogi_qual_cycle1", "midgup_qual_cycle1", "mat_qual_cycle1", "peak_qual_cycle1", "sen_qual_cycle1", "midgdown_qual_cycle1", "dor_qual_cycle1", "evi_area_cycle2", "evi_amp_cycle2", "evi_min_cycle2", "frac_filled_gup_cycle2", "frac_filled_gdown_cycle2", "length_gup_cycle2", "length_gdown_cycle2", "ogi_cycle2", "midgup_cycle2", "mat_cycle2", "peak_cycle2", "sen_cycle2", "midgdown_cycle2", "dor_cycle2", "ogi_qual_cycle2", "midgup_qual_cycle2", "mat_qual_cycle2", "peak_qual_cycle2", "sen_qual_cycle2", "midgdown_qual_cycle2", "dor_qual_cycle2")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# prelims
library(parallel)
library(argparse)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# parse the command line arguments
arg_parser <- ArgumentParser()
arg_parser$add_argument("-tile", type="character") # tile to process
arg_parser$add_argument("-year", type="integer") # year of interest
arg_parser$add_argument("-pheno_period_start", type="character") # limit analysis to cycles between pheno_period_start/end
arg_parser$add_argument("-pheno_period_end", type="character")
arg_parser$add_argument("-out_dir", type="character", default="/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT") # output directory
arg_parser$add_argument("-out_prefix", type="character", default="MCD12Q2C6") # prefix for output files
arg_parser$add_argument("-data_dir", type="character", default="/projectnb/modislc/users/dsm/spline_codes/spline_one_v3/outputs_global") # input binary splined evi data directory
arg_parser$add_argument("-data_prefix", type="character", default="c6_str5") # prefix of input binary files
arg_parser$add_argument("-params", type="character") # phenology processing parameters file
arg_parser$add_argument("-chunk_line_size", type="integer", default=480) # lines to process at once; consider memory; default does tile in 5 equal chunks; NOTE: limited to factors of 2400
arg_parser$add_argument("-cluster_size", type="integer", default=16) # number of CPU cores to use for processing
# arg_parser$add_argument("-code_dir", type="character", default="/usr2/faculty/dsm/codes/pixel-forge/MCD12Q2C6") # input data where codes are stored
arg_parser$add_argument("-code_dir", type="character", default="/projectnb/modislc/users/joshgray/MCD12Q2C6") # input data where codes are stored
arg_parser$add_argument("-lw_mask_dir", type="character", default="/projectnb/modislc/data/mcd12_in/c6/ancillary_layers/C6_LW_Mask/lw_mask_500m") # input data where codes are stored

args <- arg_parser$parse_args()
# args <- arg_parser$parse_args(c("-tile","h12v04","-year",2002))

## load the functions we need
if(!is.null(args$code_dir)) {
  path_to_fxns = paste(args$code_dir,"/MCD12Q2C6_AnnualPhenologyFunctions.R",sep="")
  source(path_to_fxns)
} else {
  print("No code dir specified. Abort!")
  quit(save="no",status=0)
}

# if a parameter file is specified, retrieve the parameters from there, use defaults otherwise
if(is.null(args$params)){
  pheno_pars <- DefaultPhenoParameters()
}else{
  pheno_pars <- ReadPhenoParameters(args$params, args$year)
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

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# set up the cluster
cl <- makeCluster(args$cluster_size)
clusterExport(cl, "path_to_fxns")
# NOTE: not sure if this next line is necessary...
clusterEvalQ(cl, { source(path_to_fxns) })

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# get all data files for tile
### sanjeeb - TODO where you need to edit
patt <- paste(args$data_prefix, "\\..*", args$tile, ".*evi2.bip", sep="")
#patt <- paste(args$data_prefix, "\\..*", args$tile, ".*MCD12I2*.bip", sep="")
evi2_files <- dir(args$data_dir, pattern=patt, full=T)
patt <- paste(args$data_prefix, "\\..*", args$tile, ".*evi2.residual.bip", sep="")
#patt <- paste(args$data_prefix, "\\..*", args$tile, ".*MCD12I4*.bip", sep="")
resid_files <- dir(args$data_dir, pattern=patt, full=T)
patt <- paste(args$data_prefix, "\\..*", args$tile, ".*evi2.flag.bip", sep="")
#patt <- paste(args$data_prefix, "\\..*", args$tile, ".*MCD12I3*.bip", sep="")
snow_files <- dir(args$data_dir, pattern=patt, full=T)

# create 3 years worth of daily dates, IGNORING LEAP YEARS!
c6_dates <- as.Date(paste(rep((args$year - 1):(args$year + 1), each=365), rep(1:365, 3), sep="-"), format="%Y-%j")

## create the output names and dimensions - these will be ~13 in number
pheno_metric_names <- c("NumCycles", "EVI_Area", "EVI_Amplitude", "EVI_Minimum", "Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy", "QA_Overall", "QA_Detailed")
num_out_mets <- length(pheno_metric_names)
## the first one has 1 band and the rest have 2
nbands_out <- c(1, rep(2,(num_out_mets-1)))

## make output directories for each of these layers
for(b in 1:num_out_mets) {
	dir.create(file.path(args$out_dir, pheno_metric_names[b]), showWarnings = FALSE)
}

## get the LW mask for this tile
# NOTE: this search pattern may need to be updated for GISS implementation!
lwmask_file <- dir(args$lw_mask_dir, pattern=paste("map.*", args$tile, "$", sep=""), full=T)
lwmask_v <- ReadSplinedNBAR(lwmask_file, nbands=1, s_flag=F, dsize=1)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# do the phenology and write to output in chunks of an equal number of lines
for(i in 1:(2400 / args$chunk_line_size)){
	# load the data
	print(paste("Loading chunk", i))
	start_time <- Sys.time()
	start_line <- ((i - 1) * args$chunk_line_size) + 1 # start reading data chunk at this line
	tmp <- Get3YearDataChunk(evi2_files, resid_files, snow_files, args$year, start_line, args$chunk_line_size)
	# if(i ==1) {
	# print(evi2_files)
	# print(dim(tmp)) }
	print(Sys.time() - start_time)

  # apply the phenology algorithm
	print(paste("Retrieving phenology of chunk", i))
  	system.time(pheno_values <- parApply(cl, tmp, 1, AnnualPhenologyC6, dates=c6_dates, pheno_pars=pheno_pars, pheno_period_start=args$pheno_period_start, pheno_period_end=args$pheno_period_end))
	# if(i ==1 ) {
	# print(dim(pheno_values))
	# }
	print(Sys.time() - start_time)

  # convert to vector and write to file - for each layer
  lwmask_inds <- ((start_line - 1) * 2400 + 1):((start_line + args$chunk_line_size - 1) * 2400) # indices for this chunk within the LW mask file
  # DEBUG
  print(paste("i", i, "start_line", start_line, "range lwmask inds", paste(range(lwmask_inds), collapse=", "), sep=": "))

	print(paste("Writing chunk", i))
	for(b in 1:num_out_mets) {
		out_file <- file.path(paste(args$out_dir,"/",pheno_metric_names[b],sep=''), paste(pheno_metric_names[b], args$tile, args$year, sep="_"))
		if(b==1) {
      layer_1_vals <- c(pheno_values[b,])
      layer_1_vals[lwmask_v[lwmask_inds] != 2] <- 32767 # fill NA where the LW mask is not LAND
			# WritePhenologyData(out_file=out_file, pheno_values_bip=c(pheno_values[b,]), tile=args$tile, nbands=nbands_out[b])
      WritePhenologyData(out_file=out_file, pheno_values_bip=layer_1_vals, tile=args$tile, nbands=nbands_out[b])
		} else {
      layer_1_vals <- pheno_values[b,]
      layer_1_vals[lwmask_v[lwmask_inds] != 2] <- 32767 # fill NA where the LW mask is not LAND
      layer_2_vals <- pheno_values[(b+12),]
      layer_2_vals[lwmask_v[lwmask_inds] != 2] <- 32767 # fill NA where the LW mask is not LAND
			# temp_in = rbind(pheno_values[b,],pheno_values[(b+12),])
      temp_in <- rbind(layer_1_vals, layer_2_vals)
			# WritePhenologyData(out_file=out_file, pheno_values_bip=as.vector(temp_in), tile=args$tile, nbands=nbands_out[b])
      WritePhenologyData(out_file=out_file, pheno_values_bip=as.vector(temp_in), tile=args$tile, nbands=nbands_out[b])
		}
	}
	print(Sys.time() - start_time)

	# clean up to free memory; NOTE: maybe not useful...
	rm(tmp)
	rm(pheno_values)
	rm(out_file)
}
