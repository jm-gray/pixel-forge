# /projectnb/modislc/users/dsm/spline_codes/spline_one_v3/outputs2
#
# You want to read in the c5.h19v03 files, then the files starting with c6_str# where the # is how often files were read in either every file for 1 up to every 16.  So you want to read in the c5, str1, str2, str4, str8, and str16 files and see how these are different.

library(raster)
library(parallel)
library(RColorBrewer)
source("/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_Functions.R")

#-----------------------------------------------------------------------
data_dir <- "/projectnb/modislc/users/dsm/spline_codes/spline_one_v3/outputs"
# tile <- "h12v04"
# tile <- "h25v06"
# tile <- "h20v08"
# tile <- "h23v02"
tile <- "h11v04"

# c5_evi2_files <- dir(data_dir, pattern="c5.*evi2.bip", full=T)
# patt <- paste("c6_str1\\..*", tile, ".*evi2.bip", sep="")
# c6_1day_evi2_files <- dir(data_dir, pattern=patt, full=T)

# patt <- paste("c6_str2\\..*", tile, ".*evi2.bip", sep="")
# c6_2day_evi2_files <- dir(data_dir, pattern=patt, full=T)
# patt <- paste("c6_str2\\..*", tile, ".*evi2.residual.bip", sep="")
# c6_2day_resid_files <- dir(data_dir, pattern=patt, full=T)
#
# patt <- paste("c6_str3\\..*", tile, ".*evi2.bip", sep="")
# c6_3day_evi2_files <- dir(data_dir, pattern=patt, full=T)
# patt <- paste("c6_str3\\..*", tile, ".*evi2.residual.bip", sep="")
# c6_3day_resid_files <- dir(data_dir, pattern=patt, full=T)
#
# patt <- paste("c6_str4\\..*", tile, ".*evi2.bip", sep="")
# c6_4day_evi2_files <- dir(data_dir, pattern=patt, full=T)
# patt <- paste("c6_str4\\..*", tile, ".*evi2.residual.bip", sep="")
# c6_4day_resid_files <- dir(data_dir, pattern=patt, full=T)

# patt <- paste("c6_str8\\..*", tile, ".*evi2.bip", sep="")
# c6_8day_evi2_files <- dir(data_dir, pattern=patt, full=T)
# patt <- paste("c6_str16\\..*", tile, ".*evi2.bip", sep="")
# c6_16day_evi2_files <- dir(data_dir, pattern=patt, full=T)

patt <- paste("c6_str4t\\..*", tile, ".*evi2.bip", sep="")
c6_4day_evi2_files <- dir(data_dir, pattern=patt, full=T)
patt <- paste("c6_str4t\\..*", tile, ".*evi2.residual.bip", sep="")
c6_4day_resid_files <- dir(data_dir, pattern=patt, full=T)
patt <- paste("c6_str4t\\..*", tile, ".*evi2.flag.bip", sep="")
c6_4day_snow_files <- dir(data_dir, pattern=patt, full=T)


patt <- paste("c6_str5\\..*", tile, ".*evi2.bip", sep="")
c6_5day_evi2_files <- dir(data_dir, pattern=patt, full=T)
patt <- paste("c6_str5\\..*", tile, ".*evi2.residual.bip", sep="")
c6_5day_resid_files <- dir(data_dir, pattern=patt, full=T)
patt <- paste("c6_str5\\..*", tile, ".*evi2.flag.bip", sep="")
c6_5day_snow_files <- dir(data_dir, pattern=patt, full=T)


# patt <- paste("c6_str4t\\..*", tile, ".*evi2.bip", sep="")
# c6_4day_evi2_files_new <- dir(data_dir, pattern=patt, full=T)
# patt <- paste("c6_str4t\\..*", tile, ".*evi2.residual.bip", sep="")
# c6_4day_resid_files_new <- dir(data_dir, pattern=patt, full=T)

# read in the central year's data
# c5_v <- read_int_file(c5_evi2_files[3], 2, 365)
# c6_1day_v <- read_int_file(c6_1day_evi2_files[3], 2, 365)
# c6_2day_v <- read_int_file(c6_2day_evi2_files[3])
# c6_3day_v <- read_int_file(c6_3day_evi2_files[3])
# c6_4day_v <- read_int_file(c6_4day_evi2_files[3])
# c6_4day_v_new <- read_int_file(c6_4day_evi2_files[3])
# c6_2day_resid_v <- read_int_file(c6_2day_resid_files[3])
# c6_3day_resid_v <- read_int_file(c6_3day_resid_files[3])
# c6_4day_resid_v <- read_int_file(c6_4day_resid_files[3])
# c6_4day_resid_v_new <- read_int_file(c6_4day_resid_files[3])
# c6_8day_v <- read_int_file(c6_8day_evi2_files[3], 2, 365)
# c6_16day_v <- read_int_file(c6_16day_evi2_files[3], 2, 365)

c6_4day_v <- read_int_file(c6_4day_evi2_files[3])
c6_4day_resid_v <- read_int_file(c6_4day_resid_files[3])
c6_5day_v <- read_int_file(c6_5day_evi2_files[3])
c6_5day_resid_v <- read_int_file(c6_5day_resid_files[3])

# c6_4day_resid_v <- read_int_file(c6_4day_resid_files[3])


#-----------------------------------------------------------------------
# plotting
mycols <- c("#636363", "#DE2D26", "#3182BD", "#31A354", "#756BB1", "#E6550D")
mypal <- colorRampPalette(brewer.pal(11, "Spectral"))
evi2_breaks <- c(-1e6, seq(0.1e4, 0.55e4, by=50), 1e6)
r <- raster(c6_4day_evi2_files[3], 240)
# par(mfrow=c(1,2),mar=c(3,1,1,1))
# hist(r, breaks=evi2_breaks,, xlim=c(0, 1e4))
plot(r, breaks=evi2_breaks, col=mypal(length(evi2_breaks) - 1), maxpixels=1e6)
# myc <- click(r, cell=T)
points(xyFromCell(r, myc$cell), pch=1)
text(xyFromCell(r, myc$cell), labels=1:dim(myc)[1], pos=4, cex=1.5)

# plot the t.s.
layout(matrix(1:6, nrow=2, byrow=T))
par(mar=c(1, 4, 1, 1), oma=rep(0.5, 4))
for(i in 1:6){
  # vs <- c(c5_v[myc$cell[i], ], c6_1day_v[myc$cell[i], ], c6_4day_v[myc$cell[i], ], c6_4day_v[myc$cell[i], ], c6_8day_v[myc$cell[i], ], c6_16day_v[myc$cell[i], ])
  # vs <- c(c6_1day_v[myc$cell[i], ], c6_4day_v[myc$cell[i], ], c6_4day_v[myc$cell[i], ], c6_8day_v[myc$cell[i], ], c6_16day_v[myc$cell[i], ])
  # vs <- c(c6_4day_v[myc$cell[i], ], c6_3day_v[myc$cell[i], ], c6_4day_v[myc$cell[i], ])
  vs <- c(c6_4day_v[myc$cell[i], ], c6_5day_v[myc$cell[i], ])
  plot(vs, type="n", xlab="", ylab="EVI2", xlim=c(0, 365))
  # points(c5_v[myc$cell[i], ], type="l", col=mycols[1], lwd=1.5)
  # points(c6_1day_v[myc$cell[i], ], type="l", col=mycols[1], lwd=1.5)

  # filt_data_2day <- c6_2day_v[myc$cell[i], ]
  # filt_data_2day[filt_data_2day == 32767] <- NA
  # resid_data_2day <- c6_2day_resid_v[myc$cell[i], ]
  # resid_data_2day[resid_data_2day == 32767] <- NA
  # og_data_2day <- filt_data_2day - resid_data_2day
  # points(filt_data_2day, type="l", col=mycols[1], lwd=1.5)
  # points(og_data_2day, type="p", col=mycols[1], pch=1, cex=0.5)
  #
  # filt_data_3day <- c6_3day_v[myc$cell[i], ]
  # filt_data_3day[filt_data_3day == 32767] <- NA
  # resid_data_3day <- c6_3day_resid_v[myc$cell[i], ]
  # resid_data_3day[resid_data_3day == 32767] <- NA
  # og_data_3day <- filt_data_3day - resid_data_3day
  # points(filt_data_3day, type="l", col=mycols[2], lwd=1.5)
  # points(og_data_3day, type="p", col=mycols[2], pch=2, cex=0.5)

  filt_data_4day <- c6_4day_v[myc$cell[i], ]
  filt_data_4day[filt_data_4day == 32767] <- NA
  resid_data_4day <- c6_4day_resid_v[myc$cell[i], ]
  resid_data_4day[resid_data_4day == 32767] <- NA
  og_data_4day <- filt_data_4day - resid_data_4day
  points(filt_data_4day, type="l", col=mycols[1], lwd=1.5)
  points(og_data_4day, type="p", col=mycols[1], pch=3, cex=0.5)

  filt_data_5day <- c6_5day_v[myc$cell[i], ]
  filt_data_5day[filt_data_5day == 32767] <- NA
  resid_data_5day <- c6_5day_resid_v[myc$cell[i], ]
  resid_data_5day[resid_data_5day == 32767] <- NA
  og_data_5day <- filt_data_5day - resid_data_5day
  points(filt_data_5day, type="l", col=mycols[2], lwd=1.5)
  points(og_data_5day, type="p", col=mycols[2], pch=3, cex=0.5)


  # filt_data_4day_new <- c6_4day_v_new[myc$cell[i], ]
  # filt_data_4day_new[filt_data_4day_new == 32767] <- NA
  # resid_data_4day_new <- c6_4day_resid_v_new[myc$cell[i], ]
  # resid_data_4day_new[resid_data_4day_new == 32767] <- NA
  # og_data_4day_new <- filt_data_4day_new - resid_data_4day_new
  # points(filt_data_4day_new, type="l", col=mycols[4], lwd=1.5)
  # points(og_data_4day_new, type="p", col=mycols[4], pch=3, cex=0.5)

  # points(c6_4day_v[myc$cell[i], ], type="l", col=mycols[1], lwd=1.5)
  # points(c6_3day_v[myc$cell[i], ], type="l", col=mycols[2], lwd=1.5)
  # points(c6_4day_v[myc$cell[i], ], type="l", col=mycols[3], lwd=1.5)
  # points(c6_8day_v[myc$cell[i], ], type="l", col=mycols[4], lwd=1.5)
  # points(c6_16day_v[myc$cell[i], ], type="l", col=mycols[5], lwd=1.5)
  # if(i==4) legend("topleft", legend=c("C5", "C6 1-day", "C6 2-day", "C6 4-day", "C6 8-day", "C6 16-day"), lwd=1.5, lty=1, pch=NA, col=mycols)
  # if(i==4) legend("topleft", legend=c("C6 1-day", "C6 2-day", "C6 4-day", "C6 8-day", "C6 16-day"), lwd=1.5, lty=1, pch=NA, col=mycols[1:5])
  # if(i==4) legend("topleft", legend=c("C6 2-day", "C6 3-day", "C6 4-day"), lwd=1.5, lty=1, pch=NA, col=mycols[1:3])
  # if(i==4) legend("topleft", legend=c("C6 2-day", "C6 3-day", "C6 4-day", "New 4-day"), lwd=1.5, lty=1, pch=NA, col=mycols[1:4])
  if(i==4) legend("topleft", legend=c("4-day", "5-day"), lwd=1.5, lty=1, pch=NA, col=mycols[1:2])
}

#-----------------------------------------------------------------------
# check on snow filtering vs no filtering
tile <- "h12v04"
patt <- paste("c6_str4\\..*", tile, ".*evi2.bip", sep="")
snow_filt_files <- dir(data_dir, pattern=patt, full=T)
patt <- paste("c6_str4\\..*", tile, ".*evi2.residual.bip", sep="")
snow_filt_resid_files <- dir(data_dir, pattern=patt, full=T)

patt <- paste("c6_str4_nf\\..*", tile, ".*evi2.bip", sep="")
nofilt_files <- dir(data_dir, pattern=patt, full=T)

# snow_filt_v <- read_int_file(snow_filt_files[3], 2, 365)
# snow_filt_resid_v <- read_int_file(snow_filt_resid_files[3], 2, 365)
# snow_filt_resid_v[snow_filt_resid_v == 32767] <- NA
# nofilt_v <- read_int_file(nofilt_files[3], 2, 365)

snow_filt_v <- read_int_file(snow_filt_files[4], 2, 365)
snow_filt_resid_v <- read_int_file(snow_filt_resid_files[4], 2, 365)
snow_filt_resid_v[snow_filt_resid_v == 32767] <- NA
nofilt_v <- read_int_file(nofilt_files[4], 2, 365)

evi2_breaks <- c(-1e6, seq(0.1e4, 0.4e4, by=50), 1e6)
r <- raster(snow_filt_files[3], 240)
plot(r, breaks=evi2_breaks, col=mypal(length(evi2_breaks) - 1), maxpixels=1e6)
myc <- click(r, cell=T)
points(xyFromCell(r, myc$cell), pch=1)
text(xyFromCell(r, myc$cell), labels=1:dim(myc)[1], pos=4, cex=1.5)

layout(matrix(1:6, nrow=2, byrow=T))
par(mar=c(1, 4, 1, 1), oma=rep(0.5, 4))
for(i in 1:6){
  vs <- c(snow_filt_v[myc$cell[i], ], nofilt_v[myc$cell[i], ], (snow_filt_v[myc$cell[i], ] + snow_filt_resid_v[myc$cell[i], ]))
  filt_data <- snow_filt_v[myc$cell[i], ]
  resid_data <- snow_filt_resid_v[myc$cell[i], ]
  og_data <- snow_filt_v[myc$cell[i], ] + snow_filt_resid_v[myc$cell[i], ]
  og_data[is.na(snow_filt_resid_v[myc$cell[i], ])] <- NA
  plot(vs, type="n", xlab="", ylab="EVI2", xlim=c(0, 365))
  points(snow_filt_v[myc$cell[i], ], type="l", col=mycols[1], lwd=1.5)
  points(nofilt_v[myc$cell[i], ], type="l", col=mycols[2], lwd=1.5)
  points(og_data, type="p", col=1)
  if(i==4) legend("topleft", legend=c("snow filtered", "not filtered", "data"), lwd=c(1.5, 1.5, NA), lty=c(1, 1, NA), pch=c(NA,NA,1), col=c(mycols[1:2], 1))
}

#================================================================================
source("/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_Functions.R")
data_dir <- "/projectnb/modislc/users/dsm/spline_codes/spline_one_v3/outputs"
tile <- "h11v04"
year_of_interest <- 2003

cl <- makeCluster(16)
clusterEvalQ(cl, {source("/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_Functions.R")})

pheno_pars <- list(
	min_peak_to_peak_distance=120,
	min_seg_amplitude=0.15,
	agg_amp_frac=0.15,
	max_seg_length=200,
  ogi_thresh=0.15,
  midgup_thresh=0.5,
  mat_thresh=0.95, #???
  sen_thresh=0.8, #???
  midgdown_thresh=0.5,
  dor_thresh=0.15,
  nbar_scale_factor=1e4,
  nbar_NA_value=32767,
  qual_buffer_days=14,
  qual_r2_weight=1,
  qual_fill_weight=1,
  pheno_period_start=as.numeric(as.Date(paste(year_of_interest, "-1-1", sep=""))),
  pheno_period_end=as.numeric(as.Date(paste(year_of_interest, "-12-31", sep="")))

  # unused for C6 methodology
  # out_quant=0.99,
  # spline_spar=0.15,
  # spline_spar=0,
  # out_iterations=0,
  # evi2_snow_quant=0.02,
  # ndsi_thresh=0,
  # max_snow_fill_ratio=1.25,
  # snow_free_min_quant=0.1,
  # snow_max_quant=0.98,
)

patt <- paste("c6_str4t\\..*", tile, ".*evi2.bip", sep="")
evi2_files <- dir(data_dir, pattern=patt, full=T)
patt <- paste("c6_str4t\\..*", tile, ".*evi2.residual.bip", sep="")
resid_files <- dir(data_dir, pattern=patt, full=T)
patt <- paste("c6_str4t\\..*", tile, ".*evi2.flag.bip", sep="")
snow_files <- dir(data_dir, pattern=patt, full=T)

# get 3 full years of data
c6_dates <- as.Date(paste(rep((year_of_interest - 1):(year_of_interest + 1), each=365), rep(1:365, 3), sep="-"), format="%Y-%j")

# loop through all of the tile data in 5 chunks
# NOTE: the size of the chunk has not been optimized in any way!
# NOTE: what is faster: parallelizing the pheno processing for each chunk, or have each thread read in a max chunk and single-thread process it?
# NOTE: these years are hard coded, need to get year from filename and match to year of interest
lines_to_read <- 480
# output_list <- list()
for(i in 1:5){
  start_line <- ((i - 1) * lines_to_read) + 1 # start reading data chunk at this line
  # load a chunk of data
  tmp <- cbind(
    read_int_file(evi2_files[2], lines_to_read=lines_to_read, start_line=start_line),
    read_int_file(evi2_files[3], lines_to_read=lines_to_read, start_line=start_line),
    read_int_file(evi2_files[4], lines_to_read=lines_to_read, start_line=start_line),
    read_int_file(resid_files[2], lines_to_read=lines_to_read, start_line=start_line),
    read_int_file(resid_files[3], lines_to_read=lines_to_read, start_line=start_line),
    read_int_file(resid_files[4], lines_to_read=lines_to_read, start_line=start_line),
    read_int_file(snow_files[2], lines_to_read=lines_to_read, start_line=start_line),
    read_int_file(snow_files[3], lines_to_read=lines_to_read, start_line=start_line),
    read_int_file(snow_files[4], lines_to_read=lines_to_read, start_line=start_line)
  ) # end data read/cbind data

  # apply the phenology algorithm
  system.time(pheno_values <- parApply(cl, tmp, 1, AnnualPhenologyC6, dates=c6_dates, pheno_pars=pheno_pars))

  # convert to vector and write to file
  pheno_values_bip <- unlist(pheno_values, use.names=F)
  # WritePhenologyData(out_file, pheno_values_bip, overwrite=F)

  # create and write a header

}

system.time(pheno_values <- parApply(cl, tmp[1:1e5,], 1, AnnualPhenologyC6, dates=c6_dates, pheno_pars=pheno_pars))

# investigate a particular pixel's time series
i <- 1200
x <- tmp[i, ]
x[x==pheno_pars$nbar_NA_value] <- NA
data_length <- length(x) / 3
evi <- x[1:data_length]
evi <- evi / pheno_pars$nbar_scale_factor
resids <- x[(data_length + 1):(2 * data_length)]
resids <- resids / pheno_pars$nbar_scale_factor
snowflags <- x[(2 * data_length + 1):(3 * data_length)]
valid_peaks <- try(FilterPeaks(evi, FindPotentialPeaks_minmax(evi), min_peak_to_peak_distance=pheno_pars$min_peak_to_peak_distance), silent=T)
full_segs <- try(GetFullSegs(evi, valid_peaks, max_seg_length=pheno_pars$max_seg_length, min_seg_amplitude=pheno_pars$min_seg_amplitude, agg_amp_frac=pheno_pars$agg_amp_frac), silent=T)
segmets <- lapply(full_segs, SegMet, x=x, dates=c6_dates, pheno_pars=pheno_pars)
PlotSeries(tmp[i,], c6_dates, seg_metrics=segmets)

# for(i in 1:1e3){
#   print(paste("Doing", i))
#   trash <- DoPhenologyC6_3year(tmp[i,], dates=c6_dates, pheno_pars=pheno_pars)
# }

#================================================================================
year_of_interest <- 2003
c6_dates <- as.Date(paste(rep((year_of_interest - 1):(year_of_interest + 1), each=365), rep(1:365, 3), sep="-"), format="%Y-%j")
c64day <- Get3YearDataChunk(c6_4day_evi2_files, c6_4day_resid_files, c6_4day_snow_files, year_of_interest, 1, 20)
c65day <- Get3YearDataChunk(c6_5day_evi2_files, c6_5day_resid_files, c6_5day_snow_files, year_of_interest, 1, 20)
layout(matrix(1:2, nrow=2))
par(mar=rep(1,4))
PlotSeries(c64day[1,], c6_dates)
PlotSeries(c65day[1,], c6_dates)
