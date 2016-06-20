#============================================================
# VIIRS Phenology comparison with C6 Phenology on VIIRS data
#============================================================
# Josh Gray Feb 1, 2016
#============================================================

library(raster)
library(RColorBrewer)
source("/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_Visualize.R")

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

# finds the LW mask for 500 m C5 MODIS tiles
FindLWMask <- function(tile){
  mask_dir <- "/projectnb/modislc/data/mcd12_in/c5/ancillary_layers/C5_LW_Mask/lw_mask_500m"
  # mask_dir <- "/projectnb/modislc/data/mcd12_in/c6/ancillary_layers/C6_LW_Mask/lw_mask_500m"
  mask_file <- dir(mask_dir, pattern=paste(".*", tile, ".*", "bin$", sep=""), full=T)
  return(mask_file)
}


#============================================================
YEAR_TO_DO <- 2013
TILE <- "H12V04"
# TILE <- "H11V04"
# TILE <- "H08V05"
# C6_RESULTS_DIR <- "/projectnb/modislc/users/joshgray/xyz_pheno_comparison/VIIRS_NBAR/C6_RESULTS/"
C6_RESULTS_DIR <- "/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT"

# HPLM_RESULTS_DIR <- "/projectnb/modislc/users/joshgray/xyz_pheno_comparison/VIIRS_NBAR/RESULTS/H12V04/HPLM"
HPLM_RESULTS_DIR <- file.path("/projectnb/modislc/users/joshgray/xyz_pheno_comparison/VIIRS_NBAR/RESULTS", toupper(TILE), "HPLM")
# C6_PHENO_METRIC <- "ogi_cycle1"
# HPLM_PHENO_METRIC <- "greenup"
C6_PHENO_METRIC <- "midgup_cycle1"
HPLM_PHENO_METRIC <- "midtime_grow"
# C6_PHENO_METRIC <- "mat_cycle1"
# HPLM_PHENO_METRIC <- "maturity"
# C6_PHENO_METRIC <- "sen_cycle1"
# HPLM_PHENO_METRIC <- "senescence"
# C6_PHENO_METRIC <- "midgdown_cycle1"
# HPLM_PHENO_METRIC <- "mid_sene"
# C6_PHENO_METRIC <- "dor_cycle1"
# HPLM_PHENO_METRIC <- "Dormant"

LIN_STRETCH_PERCENTILES <- c(0.05, 0.95)
LAYER_NAMES <- c("num_cycles", "fill_code", "evi_area_cycle1", "evi_amp_cycle1", "evi_min_cycle1", "frac_filled_gup_cycle1", "frac_filled_gdown_cycle1", "length_gup_cycle1", "length_gdown_cycle1", "ogi_cycle1", "midgup_cycle1", "mat_cycle1", "peak_cycle1", "sen_cycle1", "midgdown_cycle1", "dor_cycle1", "ogi_qual_cycle1", "midgup_qual_cycle1", "mat_qual_cycle1", "peak_qual_cycle1", "sen_qual_cycle1", "midgdown_qual_cycle1", "dor_qual_cycle1", "evi_area_cycle2", "evi_amp_cycle2", "evi_min_cycle2", "frac_filled_gup_cycle2", "frac_filled_gdown_cycle2", "length_gup_cycle2", "length_gdown_cycle2", "ogi_cycle2", "midgup_cycle2", "mat_cycle2", "peak_cycle2", "sen_cycle2", "midgdown_cycle2", "dor_cycle2", "ogi_qual_cycle2", "midgup_qual_cycle2", "mat_qual_cycle2", "peak_qual_cycle2", "sen_qual_cycle2", "midgdown_qual_cycle2", "dor_qual_cycle2")
#============================================================
# read in the C6 VIIRS data, convert to DOY
# c6_raster <- raster(dir(C6_RESULTS_DIR, pattern=paste(C6_PHENO_METRIC, ".*", TILE, ".*", YEAR_TO_DO, ".tif", sep=""), full=T))
# c6_raster <- raster(dir(C6_RESULTS_DIR, pattern=paste(C6_PHENO_METRIC, ".*(", TILE, "|", tolower(TILE), ").*", YEAR_TO_DO, ".tif", sep=""), full=T))
c6_file <- dir(C6_RESULTS_DIR, pattern=paste(".*", tolower(TILE), ".*", YEAR_TO_DO, "$", sep=""), full=T)
s <- stack(c6_file)
NAvalue(s) <- 32767
# c6_offset <- as.numeric(as.Date(paste(year, "-1-1", sep="")))
c6_raster <- raster(s, which(LAYER_NAMES == C6_PHENO_METRIC)) - as.numeric(as.Date("2013-1-1"))
lw <- raster(FindLWMask(tolower(TILE)))
c6_raster[lw != 1] <- NA

# get the C6 quality
C6_QUAL_METRIC <- paste(c(unlist(strsplit(C6_PHENO_METRIC, split="_")), "qual")[c(1,3,2)], collapse="_")
c6_qual_raster <- raster(s, which(LAYER_NAMES == C6_QUAL_METRIC))
c6_qual_raster[lw != 1] <- NA

# Read in the VIIRS data, take care of NoData, make a raster
hplm_in_file <- dir(HPLM_RESULTS_DIR, pattern=paste(".*", HPLM_PHENO_METRIC, ".*", YEAR_TO_DO, sep=""), full=T)
tmp <- read_int_file(hplm_in_file, dsize=2, nbands=1, s_flag=F, bsq_flag=F)
tmp[tmp == 32767] <- NA # take care of NoData values
hplm_raster <- c6_raster
values(hplm_raster) <- tmp

#============================================================
# Make a plot with three rasters: HPLM, C6, difference and histogram of differences
diff_raster <- hplm_raster - c6_raster
doy_pal <- colorRampPalette(brewer.pal(11, "Spectral"))
diff_pal <- colorRampPalette(brewer.pal(9, "RdBu"))

# setup breaks for a linear contrast stretch
doy_breaks <- c(0, seq(min(quantile(values(c6_raster), LIN_STRETCH_PERCENTILES[1], na.rm=T), quantile(values(hplm_raster), LIN_STRETCH_PERCENTILES[1], na.rm=T)), max(quantile(values(c6_raster), LIN_STRETCH_PERCENTILES[2], na.rm=T), quantile(values(hplm_raster), LIN_STRETCH_PERCENTILES[2], na.rm=T)), by=1), 365)
doy_breaks <- sort(unique(doy_breaks))

max_abs_diff_break <- max(abs(min(values(diff_raster), na.rm=T)), max(values(diff_raster), na.rm=T))
diff_breaks <- c(-1 * max_abs_diff_break, seq(quantile(values(diff_raster), LIN_STRETCH_PERCENTILES[1], na.rm=T), quantile(values(diff_raster), LIN_STRETCH_PERCENTILES[2], na.rm=T), by=2), max_abs_diff_break)

# plot to a PDF
pdf_out <- file.path("/projectnb/modislc/users/joshgray/MCD12Q2C6/", paste("HPLM_C6_comparison_", HPLM_PHENO_METRIC, ".pdf", sep=""))
pdf(h=8, w=25, file=pdf_out)
# setup the plot device
par(mar=c(2,1,2,4), oma=rep(1,4))
nf <- layout(matrix(c(1, 2, 3), nrow=1, byrow=T))
PlotPheno(hplm_raster, breaks=doy_breaks, pal=doy_pal, maxpixels=3e6)
title(paste("HPLM", HPLM_PHENO_METRIC))
PlotPheno(c6_raster, breaks=doy_breaks, pal=doy_pal, maxpixels=3e6)
title(paste("C6", C6_PHENO_METRIC))
PlotDiff(diff_raster, maxpixels=3e6)
title("HPLM - C6")
dev.off()

# x11();
pdf_out <- file.path("/projectnb/modislc/users/joshgray/MCD12Q2C6/", paste("HPLM_C6_comparison_", HPLM_PHENO_METRIC, "_diff_hist.pdf", sep=""))
pdf(h=10, w=10, file=pdf_out)
diff_v <- values(diff_raster)
hist(diff_v, breaks=diff_breaks, xlim=c(-45, 60), main="HPLM - C6")
abline(v=mean(diff_v, na.rm=T), col=2, lwd=2, lty=2)
print(mean(diff_v, na.rm=T))
dev.off()


# PlotDoyReadableLegend(hplm_raster, doy_breaks)
# PlotDOYLegend(hplm_raster, breaks=doy_breaks, pal=doy_pal)

#---------------------------------
# get residuals vs quality
NUM_SAMPLE <- 1e5
lw_v <- values(lw)
cell_nums <- 1:ncell(lw)
sample_cells <- sample(cell_nums[lw_v == 1], NUM_SAMPLE) # get a sample of land pixels
hplm_sample <- hplm_raster[sample_cells]
c6_sample <- c6_raster[sample_cells]
c6_sample_qual <- c6_qual_raster[sample_cells]
# plot(c6_sample_qual/1e4, hplm_sample - c6_sample, pch=16, cex=0.75, col=rgb(0.3, 0.3, 0.3, 0.5), xlim=c(0.8,1), ylim=c(-30, 30))
# abline(h=0, col=2, lty=1)

par(col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white", bg="black", mar=c(4, 4, 2, 2), oma=rep(1, 4))
smoothScatter(c6_sample_qual/1e4, hplm_sample - c6_sample, xlab="C6 Quality", ylab="HPLM-C6", colramp=colorRampPalette(c("black", rev(brewer.pal(11, "Spectral")))), pch=1, nrpoints=0, xlim=c(0.7,1), ylim=c(-45, 45))
abline(h=0, lty=2, col="white", lwd=1.5)

#=============================================================================
# Create a C5 EVI2 stack
start_date <- as.Date("2012-1-1")
end_date <- as.Date("2014-12-31")
c5_data_dir <- "/projectnb/modislc/data/mcd12_in/c5/mcd43a4"
nbar_files <- Sys.glob(file.path(c5_data_dir, "*", "*", paste("*MCD43A4*", tolower(TILE), "*hdf", sep="")))
dates <- do.call("c", lapply(nbar_files, DateFunc)) # get the dates for each file
nbar_files <- nbar_files[(dates >= start_date) & (dates <= end_date)]
my_dates <- dates[(dates >= start_date) & (dates <= end_date)]
data_sets <- c(unlist(lapply(nbar_files, Get_SDS_name_band1)), unlist(lapply(nbar_files, Get_SDS_name_band2)), unlist(lapply(nbar_files, Get_SDS_name_band3)), unlist(lapply(nbar_files, Get_SDS_name_band4)), unlist(lapply(nbar_files, Get_SDS_name_band6)))

nbar_s <- stack(data_sets)
v <- getValues(nbar_s)

red_inds <- 1:length(my_dates)
nir_inds <- (length(my_dates) + 1):(2 * length(my_dates))
evi2 <- 2.5 * ((v[, nir_inds] - v[, red_inds]) / (v[, nir_inds] + (2.4 * v[, red_inds]) + 1))

tmp_r <- c6_raster
samp_ind <- 72
values(tmp_r) <- evi2[, samp_ind]
PlotStretch(tmp_r,dates=F)
myc <- click(tmp_r, cell=T)
points(xyFromCell(r, myc$cell), pch=1)
text(xyFromCell(r, myc$cell), labels=1:dim(myc)[1], pos=4, cex=1.5)

nf <- layout(matrix(1:8, nrow=4, byrow=F))
par(mar=c(3,4,1,1), oma=rep(1,4))
for(i in 1:8){
  plot(my_dates, evi2[myc$cell[i],], xlab="", ylab="NBAR-EVI2", type="n")
  xs <- c(as.Date("2013-1-1"), as.Date("2013-1-1"), as.Date("2013-12-31"), as.Date("2013-12-31"))
  ys <- c(par()$usr[3], par()$usr[4], par()$usr[4], par()$usr[3])
  polygon(xs,ys, col="lightgrey", border=NA)
  points(my_dates, evi2[myc$cell[i],])
  text(par()$usr[1], par()$usr[4], labels=i, cex=2.5, adj=c(0, 1))
}
