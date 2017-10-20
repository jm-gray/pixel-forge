library(argparse)
library(raster)
library(RColorBrewer)
library(rgdal)

source("/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_AnnualPhenologyFunctions.R")

#---------------------------------------------------------------------
ReadSplinedNBAR <- function(in_file, lines_to_read=2400, start_line=1, samples=2400, dsize=2, nbands=365, s_flag=T, bsq_flag=F){
	# For reading from the binary, splined output files which contain 365 bands (days) of spline-smoothed EVI2 data/resids/QA

  # start_line begins at 1!!!
  where_to_start <- (start_line - 1) * samples * dsize * nbands
  size_to_read <- lines_to_read * nbands * samples
  f <- file(in_file, "rb") # open the file
  seek(f, where_to_start) # position the file connection
  temp <- readBin(f, integer(), n=size_to_read, endian= "little", size=dsize, signed=s_flag) # read all of the data into an array
  close(f) # close the file
  # re-order the temp array into a matrix; if bsq_flag then read by row, otherwise by column
  ifelse(bsq_flag, byr_flag<-FALSE, byr_flag<-TRUE)
  temp <- matrix(temp, ncol=nbands, byrow=byr_flag)
  return(temp) # return the data
}

#---------------------------------------------------------------------
Get3YearDataChunk_tmp <- function(evi2_files, resid_files, snow_files, year_of_interest, start_line, lines_to_read){
	# retrieves 3 full years of splined C6 NBAR-EVI2 data, residuals, and snow flags as a data.frame where each row
	# is a time series of evi, resid, and snow: evi001, ..., evi365, resid001, ..., resid365, flag001, ..., flag365

	# small function to retrieve the year from the nbar file names
	# NOTE: this is FRAGILE! A more robust solution would be better
	year_function <- function(x) gsub(".*\\.A([0-9]{4}).*", "\\1", basename(x))

	evi2_years <- as.integer(unlist(lapply(evi2_files, year_function)))
	resid_years <- as.integer(unlist(lapply(resid_files, year_function)))
	snow_years <- as.integer(unlist(lapply(snow_files, year_function)))
	years_of_interest <- c(year_of_interest - 1, year_of_interest, year_of_interest + 1)

	# get indices for year before, year of interest, and year after in the files
	evi2_file_inds <- match(years_of_interest, evi2_years)
	resid_file_inds <- match(years_of_interest, resid_years)
	snow_file_inds <- match(years_of_interest, snow_years)

	# check that all needed data is available
	if(any(is.na(c(evi2_file_inds, resid_file_inds, snow_file_inds)))){
		print("Some data are missing, aborting")
		return(NA)
	}

	tmp <- cbind(
		ReadSplinedNBAR(evi2_files[evi2_file_inds[1]], lines_to_read=lines_to_read, start_line=start_line),
		ReadSplinedNBAR(evi2_files[evi2_file_inds[2]], lines_to_read=lines_to_read, start_line=start_line),
		ReadSplinedNBAR(evi2_files[evi2_file_inds[3]], lines_to_read=lines_to_read, start_line=start_line),
		ReadSplinedNBAR(resid_files[resid_file_inds[1]], lines_to_read=lines_to_read, start_line=start_line),
		ReadSplinedNBAR(resid_files[resid_file_inds[2]], lines_to_read=lines_to_read, start_line=start_line),
		ReadSplinedNBAR(resid_files[resid_file_inds[3]], lines_to_read=lines_to_read, start_line=start_line),
		ReadSplinedNBAR(snow_files[snow_file_inds[1]], lines_to_read=lines_to_read, start_line=start_line),
		ReadSplinedNBAR(snow_files[snow_file_inds[2]], lines_to_read=lines_to_read, start_line=start_line),
		ReadSplinedNBAR(snow_files[snow_file_inds[3]], lines_to_read=lines_to_read, start_line=start_line)
	) # end data read/cbind data

	return(tmp)
}
#---------------------------------------------------------------------

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

# args <- arg_parser$parse_args()
# args <- arg_parser$parse_args(c("-tile","h09v07", "-year", 2006))
args <- arg_parser$parse_args(c("-tile","h12v04", "-year", 2004))

# New Data Directory for run 2.0 processing:
# data_dir="/Users/jmgray2/Desktop/NewSpline"
data_dir="/projectnb/modislc/users/dsm/eval_modis_lc_061917"
splined_nbar_data_dir <- file.path(data_dir, "MCD12I2")
splined_qa_data_dir <- file.path(data_dir, "MCD12I3")
splined_resids_data_dir <- file.path(data_dir, "MCD12I4")
# splined_nbar_data_dir <- "/projectnb/modislc/users/joshgray/C6_Diagnostics/SplinedOutput/splined_nbar"
# splined_qa_data_dir <- "/projectnb/modislc/users/joshgray/C6_Diagnostics/SplinedOutput/splined_qa"
# splined_resids_data_dir <- "/projectnb/modislc/users/joshgray/C6_Diagnostics/SplinedOutput/splined_resids"

patt <- paste("MCD12", ".*", args$tile, ".*.bip$", sep="")
evi2_files <- dir(splined_nbar_data_dir, pattern=patt, full=T)
resid_files <- dir(splined_resids_data_dir, pattern=patt, full=T)
snow_files <- dir(splined_qa_data_dir, pattern=patt, full=T)

c6_dates <- as.Date(paste(rep((args$year - 1):(args$year + 1), each=365), rep(1:365, 3), sep="-"), format="%Y-%j")
# v <- ReadSplinedNBAR(evi2_files[2])
# i <- 1199
# start_line <- ((i - 1) * args$chunk_line_size) + 1 # start reading data chunk at this line
start_line <- 1200
tmp <- Get3YearDataChunk_tmp(evi2_files, resid_files, snow_files, args$year, start_line, args$chunk_line_size)
pheno_pars <- DefaultPhenoParameters()
pheno_period_start <- as.numeric(as.Date(paste(args$year, "-1-1", sep="")))
pheno_period_end <- as.numeric(as.Date(paste(args$year, "-12-31", sep="")))
AnnualPhenologyC6(tmp[1,], c6_dates, pheno_pars, pheno_period_start, pheno_period_end, plot=T)

##########################################
mod_inds <- seq(1, ncol(tmp), by=length(c6_dates))
cell_number <- (2400*400) + 1600
tmp_evi2 <- tmp[cell_number, mod_inds[1]:(mod_inds[1] + length(c6_dates) - 1)]
tmp_resids <- tmp[cell_number, mod_inds[2]:(mod_inds[2] + length(c6_dates) - 1)]
tmp_evi2[tmp_evi2 == 32767] <- NA
tmp_resids[tmp_resids == 32767] <- NA
plot(c6_dates, tmp_evi2 + tmp_resids, xlab="", ylab="EVI2", col=2)
points(c6_dates, tmp_evi2, type="l", lty=2)
legend("topleft", legend=c("Splined EVI2", "NBAR Valuea"), lty=c(NA, 2), pch=c(1, NA), col=c(2, 1))

abline(v=seq.Date(as.Date("2003-6-1-"), as.Date("2005-6-1"),by="year"))






# scratch plotting
# tile <- "h09v07"
# in_file <- "/projectnb/modislc/users/joshgray/C6_Diagnostics/SplinedOutput/splined_nbar/MCD12I2.A2006001.h09v07.1467.20161203132004.bip"
tmp_r <- raster("/projectnb/modislc/users/joshgray/MCD12Q2C6/CLM_output/CLMextractv2/twentygup_h09v07_2005.tif")
lwmask_dir <- "/projectnb/modislc/data/mcd12_in/c6/ancillary_layers/C6_LW_Mask/lw_mask_500m"
tile_lwmask <- raster(dir(lwmask_dir, pattern=paste("map.*", args$tile, "$", sep=""), full=T))

which_file <- 1
v <- ReadSplinedNBAR(evi2_files[which_file])
v_resid <- ReadSplinedNBAR(resid_files[which_file])
v_qa <- ReadSplinedNBAR(snow_files[which_file])
v[v == 32767] <- NA
v_resid[v_resid == 32767] <- NA

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# Plot whole tile splined NBAR for a particular date
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
r <- raster(t(matrix(v[,151], nrow=2400, ncol=2400, byrow=F)))
extent(tile_lwmask) <- extent(r)
r[tile_lwmask == 1] <- NA
qs <- quantile(r, c(0,0.02,0.98,1), na.rm=T)
pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
breaks <- unique(c(qs[1], seq(qs[2], qs[3], len=254), qs[4]))
plot(r, breaks=breaks, col=pal(length(breaks) - 1), legend=F)



#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# Plot time series for select pixels
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
layout(matrix(1:2, nrow=2))
par(mar=c(2,4,1,1))
plot_dates <- as.Date(paste(2006, 1:365, sep="-"), format="%Y-%j")

the_row1 <- cellFromXY(r, myc[3,1:2])
the_row2 <- cellFromXY(r, myc[4,1:2])
ylim <- c(min(c(v[the_row1, ] + v_resid[the_row1, ], v[the_row2, ] + v_resid[the_row2, ]), na.rm=T), max(c(v[the_row1, ] + v_resid[the_row1, ], v[the_row2, ] + v_resid[the_row2, ]), na.rm=T))
plot(plot_dates, rep(NA, length(plot_dates)), xlab="", ylab="EVI2", ylim=ylim, type="n")
points(plot_dates, v[the_row1, ], type="l", lwd=2, col="tomato", lty=2)
points(plot_dates, v[the_row1, ] + v_resid[the_row1, ], type="p", col="tomato", pch=16, cex=1.5)
points(plot_dates, v[the_row2, ], type="l", lwd=2, col="seagreen", lty=2)
points(plot_dates, v[the_row2, ] + v_resid[the_row2, ], type="p", col="seagreen", pch=16, cex=1.5)
legend("topleft", legend=c(paste("line:", rowFromY(r, myc[3,2])), paste("line:", rowFromY(r, myc[4,2]))), pch=16, lty=2, lwd=2, col=c("tomato", "seagreen"))

the_row1 <- cellFromXY(r, myc[1,1:2])
the_row2 <- cellFromXY(r, myc[2,1:2])
ylim <- c(min(c(v[the_row1, ] + v_resid[the_row1, ], v[the_row2, ] + v_resid[the_row2, ]), na.rm=T), max(c(v[the_row1, ] + v_resid[the_row1, ], v[the_row2, ] + v_resid[the_row2, ]), na.rm=T))
plot(plot_dates, rep(NA, length(plot_dates)), xlab="", ylab="EVI2", ylim=ylim, type="n")
points(plot_dates, v[the_row1, ], type="l", lwd=2, col="tomato", lty=2)
points(plot_dates, v[the_row1, ] + v_resid[the_row1, ], type="p", col="tomato", pch=16, cex=1.5)
points(plot_dates, v[the_row2, ], type="l", lwd=2, col="seagreen", lty=2)
points(plot_dates, v[the_row2, ] + v_resid[the_row2, ], type="p", col="seagreen", pch=16, cex=1.5)
legend("topleft", legend=c(paste("line:", rowFromY(r, myc[1,2])), paste("line:", rowFromY(r, myc[2,2]))), pch=16, lty=2, lwd=2, col=c("tomato", "seagreen"))


#------------------------------------------------------------------------
GetSDSName <- function(mcd43a4_file_path, band){
  return(paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":", paste("MOD_Grid_BRDF:Nadir_Reflectance_Band", band, sep=""), sep = ""))
}

#------------------------------------------------------------------------
CalcEVI2 <- function(nir_r, red_r, scale_factor=1){
  # s <- stack(in_file)
  nir_v <- values(nir_r)
  red_v <- values(red_r)
  evi2 <- 2.5 * (((nir_v * scale_factor) - (red_v * scale_factor)) / ((nir_v * scale_factor) + (2.4 * (red_v * scale_factor)) + 1))
	tmp_r <- nir_r
	values(tmp_r) <- evi2
  return(tmp_r)
}


# spline_v <- ReadSplinedNBAR("~/Desktop/MCD43A4/SplinedOutput/MCD12I2.A2005001.h09v07.1467.20161203074847.bip")
# resid_v <- ReadSplinedNBAR("~/Desktop/MCD43A4/SplinedOutput/MCD12I4.A2005001.h09v07.1467.20161203074847.bip")
lw_mask <- raster("LW.map.h09v07")

# spline_v <- ReadSplinedNBAR("~/Desktop/MCD43A4/SplinedOutput/MCD12I2.A2006001.h09v07.1467.20161203132004.bip")
# resid_v <- ReadSplinedNBAR("~/Desktop/MCD43A4/SplinedOutput/MCD12I4.A2006001.h09v07.1467.20161203132004.bip")
# spline_v <- ReadSplinedNBAR("~/Desktop/MCD43A4/SplinedOutput/c5_str8.h09v07.2006.evi2.bip")
# resid_v <- ReadSplinedNBAR("~/Desktop/MCD43A4/SplinedOutput/c5_str8.h09v07.2006.evi2.residual.bip")
# spline_v <- ReadSplinedNBAR("~/Desktop/MCD43A4/SplinedOutput/MCD12I2.A2002001.h09v07.1467.20161201201234.bip")
# resid_v <- ReadSplinedNBAR("~/Desktop/MCD43A4/SplinedOutput/MCD12I4.A2002001.h09v07.1467.20161201201234.bip")
spline_v <- ReadSplinedNBAR("~/Desktop/MCD43A4/SplinedOutput/BugFixed/c6_str5.h09v07.2006.evi2.bip")
resid_v <- ReadSplinedNBAR("~/Desktop/MCD43A4/SplinedOutput/BugFixed/c6_str5.h09v07.2006.evi2.residual.bip")

geo_spline_v <- ReadSplinedNBAR("~/Desktop/MCD43A4/SplinedOutput/c6_str5.h09v07.2002.evi2.bip")
geo_resid_v <- ReadSplinedNBAR("~/Desktop/MCD43A4/SplinedOutput/c6_str5.h09v07.2002.evi2.residual.bip")

doy <- 151
geo_evi_spline <- tmp_r
values(geo_evi_spline) <- geo_spline_v[,doy]
geo_evi_spline[geo_evi_spline == 32767] <- NA
geo_evi_spline[lw_mask != 2] <- NA
geo_evi_spline <- geo_evi_spline / 1e4


doy <- 151
evi_spline <- tmp_r
values(evi_spline) <- spline_v[,doy]
evi_spline[evi_spline == 32767] <- NA
evi_spline[lw_mask != 2] <- NA
evi_spline <- evi_spline / 1e4

spline_tmp <- spline_v[,doy]
resid_tmp <- resid_v[,doy]
spline_tmp[spline_tmp == 32767] <- NA
resid_tmp[resid_tmp == 32767] <- NA
# evi_tmp <- spline_tmp + resid_tmp
evi_tmp <- spline_tmp - resid_tmp
# evi_151[evi_151 == 32767] <- NA
evi_spline_resids <- tmp_r
values(evi_spline_resids) <- evi_tmp
evi_spline_resids[lw_mask != 2] <- NA
evi_spline_resids <- evi_spline_resids / 1e4

qs <- quantile(evi_spline, c(0, 0.02, 0.98, 1), na.rm=T)
breaks <- c(qs[1], seq(qs[2], qs[3], len=254), qs[4])

layout(matrix(1:3, nrow=1))
par(mar=rep(1, 4))
plot(evi_spline, maxpixels=6e6, col=pal(255), breaks=breaks, legend=F)
title("Splined Out 2006-151")
plot(evi_spline_resids, maxpixels=6e6, col=pal(255), breaks=breaks, legend=F)
title("Spline Out 2006-151 w/ Resids")
plot(evi2_r, maxpixels=6e6, col=pal(255), breaks=breaks, legend=F)
title("MCD43A4 C6 2006-151")

set.seed(42)
mysamp <- sample(1:ncell(tmp_r), 1e6)
plot(evi2_r[mysamp], evi_spline_resids[mysamp], xlab="EVI2 from NBAR", ylab="EVI2 Values Recovered from Spline", xlim=c(0,1), ylim=c(0,1), pch=16, cex=0.7, col=rgb(0.2,0.2,0.2,0.5))
abline(a=0,b=1,lty=2,col=2,lwd=2)

in_files <- normalizePath(dir(".",pattern="MCD43",full=T))
pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

file_ind <- 10
red_r <- raster(GetSDSName(in_files[file_ind], band=1))
red_r[lw_mask != 2] <- NA
nir_r <- raster(GetSDSName(in_files[file_ind], band=2))
nir_r[lw_mask != 2] <- NA

evi2_r <- CalcEVI2(nir_r, red_r)
# evi2_r_sub <- crop(evi2_r, extent(myp))

layout(matrix(1:2, nrow=1))
par(mar=rep(0.5,4))
qs <- quantile(evi2_r, c(0, 0.02, 0.98, 1), na.rm=T)
breaks <- c(qs[1], seq(qs[2], qs[3], len=254), qs[4])
plot(evi2_r, breaks=breaks, col=pal(255), maxpixels=6e6, legend=F)
plot(extent(myp), lwd=3, add=T, col="white")
plot(evi2_r_sub, breaks=breaks, col=pal(255), maxpixels=6e6, legend=F)
