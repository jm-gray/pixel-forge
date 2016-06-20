source("/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_AnnualPhenologyFunctions.R")
library(rgdal)
library(raster)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# PhenoCam Processing
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
GetTile <- function(x){
	sin_x <- as.numeric(x[4])
	sin_y <- as.numeric(x[5])
	horiz_bins <- c(-18:18) * 2400 * 463.3127
	horiz_labels <- 0:35
	vert_bins <- c(-9:9) * 2400 * 463.3127
	vert_labels <- 17:0
	v_tile <- formatC(vert_labels[findInterval(sin_y, vert_bins, all.inside=T)], width=2, flag="0")
	h_tile <- formatC(horiz_labels[findInterval(sin_x, horiz_bins, all.inside=T)], width=2, flag="0")
	tile <- paste("h", h_tile, "v", v_tile, sep="")
	return(tile)
}

# get phenocam sinusoidal coordinates
phenocams <- read.table("/projectnb/modislc/users/joshgray/MCD12Q2C6/siteloc_active.txt", sep="|", head=T)[, c(3,2,1)]
sin_proj <- projection(raster("/projectnb/modislc/users/joshgray/MCD12Q2C6/filltiles/h12v04fill.tif"))
sin_coords <- project(as.matrix(phenocams[, c(3, 2)]), sin_proj)
phenocams$sin_x <- sin_coords[, 1]
phenocams$sin_y <- sin_coords[, 2]
phenocams$tile <- apply(phenocams, 1, GetTile)
phenocams$row <- phenocams$col <- rep(NA, dim(phenocams)[1])
# get line/sample information
for(this_tile in unique(phenocams$tile)){
	print(paste("Doing tile:", this_tile))
	# tmp <- phenocams[phenocams$tile == this_tile, ]
	tmp_r <- Sys.glob(paste("/projectnb/modislc/users/joshgray/MCD12Q2C6/filltiles/*", this_tile, "*", sep=""))
	if(!length(tmp_r) == 0){
		tmp_r <- raster(tmp_r)
		phenocams$col[phenocams$tile == this_tile] <- colFromX(tmp_r, phenocams$sin_x[phenocams$tile == this_tile]) - 1

		phenocams$row[phenocams$tile == this_tile] <- rowFromY(tmp_r, phenocams$sin_y[phenocams$tile == this_tile]) - 1
	}else{
		print("No tile found for ", this_tile)
	}
}


# grab some PC t.s.
# umichbiological
# merbleue
# hubbardbrook (DB2)
# groundhog (XX)
tile_of_interest <- "h03v06" # for kamuela
# tile_of_interest <- "h12v04" # for harvard, umichbiological, merbleue, hubbardbrook, groundhog
year_of_interest <- 2013

# extract time series
# get all data files for tile
data_prefix <- "c6_str5"
spline_data_dir <- "/projectnb/modislc/users/dsm/spline_codes/spline_one_v3/outputs_global"
patt <- paste(data_prefix, "\\..*", tile_of_interest, ".*evi2.bip", sep="")
evi2_files <- dir(spline_data_dir, pattern=patt, full=T)
patt <- paste(data_prefix, "\\..*", tile_of_interest, ".*evi2.residual.bip", sep="")
resid_files <- dir(spline_data_dir, pattern=patt, full=T)
patt <- paste(data_prefix, "\\..*", tile_of_interest, ".*evi2.flag.bip", sep="")
snow_files <- dir(spline_data_dir, pattern=patt, full=T)

# create 3 years worth of daily dates, IGNORING LEAP YEARS!
dates <- as.Date(paste(rep((year_of_interest - 1):(year_of_interest + 1), each=365), rep(1:365, 3), sep="-"), format="%Y-%j")

tmp_df <-  phenocams[phenocams$title=="kamuela",]
# tmp_df <-  phenocams[phenocams$title=="harvard",]
# tmp_df <-  phenocams[phenocams$title=="groundhog",]
# tmp_df <-  phenocams[phenocams$title=="merbleue",]
# tmp_df <-  phenocams[phenocams$title=="umichbiological",]
# tmp_df <-  phenocams[phenocams$title=="hubbardbrook",]
pixel_time_series <- Get3YearDataChunk(evi2_files, resid_files, snow_files, year_of_interest, start_line=tmp_df$row, lines_to_read=1)[tmp_df$col, ]

# get the PC t.s.
phenocam_ts <- read.csv("/projectnb/modislc/users/joshgray/MCD12Q2C6/kamuela_GR_0001_3day_v4.csv", comment.char="#", head=T)
# phenocam_ts <- read.csv("/projectnb/modislc/users/joshgray/MCD12Q2C6/harvard_DB_0001_3day_v4.csv", comment.char="#", head=T)
# phenocam_ts <- read.csv("/projectnb/modislc/users/joshgray/MCD12Q2C6/groundhog_XX_0001_3day_v4.csv", comment.char="#", head=T)
# phenocam_ts <- read.csv("/projectnb/modislc/users/joshgray/MCD12Q2C6/merbleue_SH_0001_3day_v4.csv", comment.char="#", head=T)
# phenocam_ts <- read.csv("/projectnb/modislc/users/joshgray/MCD12Q2C6/umichbiological_DB_0001_3day_v4.csv", comment.char="#", head=T)
# phenocam_ts <- read.csv("/projectnb/modislc/users/joshgray/MCD12Q2C6/hubbardbrook_DB_0002_3day_v4.csv", comment.char="#", head=T)
phenocam_ts$date <- as.Date(phenocam_ts$date)

# Plot the overlapping time period in the time series
PlotPhenology(pixel_time_series, dates, plot_legend=F, plot_dates=F, ylim=c(0.15, 0.4))
par(new=T)
plot(phenocam_ts$date[phenocam_ts$date >= "2012-1-1" & phenocam_ts$date < "2015-1-1"], phenocam_ts$gcc_90[phenocam_ts$date >= "2012-1-1" & phenocam_ts$date < "2015-1-1"], type="l", lwd=1.5, col="seagreen", yaxt="n", xaxt="n", xlab="", ylab="", xlim=c(min(dates), max(dates)))
legend(
	"topleft",
	legend=c("C6 NBAR-EVI2", "Splined NBAR-EVI2", "PhenoCam GCC"),
	lty=c(NA, 1, 1),
	lwd=c(NA, 1, 1),
	col=c("#636363", "#636363", "seagreen"),
	bg="white",
	pch=c(1, NA, NA)
)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Plot Mosaics
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
source("/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_AnnualPhenologyFunctions.R")

FindLWMask <- function(tile){
  mask_dir <- "/projectnb/modislc/data/mcd12_in/c5/ancillary_layers/C5_LW_Mask/lw_mask_500m"
  # mask_dir <- "/projectnb/modislc/data/mcd12_in/c6/ancillary_layers/C6_LW_Mask/lw_mask_500m"
  mask_file <- dir(mask_dir, pattern=paste(".*", tile, ".*", "bin$", sep=""), full=T)
  return(mask_file)
}

#-----------------------------------------------------------
data_dir <- "/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT"
out_dir <- "/projectnb/modislc/users/joshgray/MCD12Q2C6/LW_MASKS"

asia_tiles <- c('h10v02', 'h11v02', 'h12v01', 'h19v00', 'h19v01', 'h19v04', 'h20v01', 'h20v02', 'h20v03', 'h20v04', 'h20v05', 'h21v01', 'h21v02', 'h21v03', 'h21v04', 'h21v05', 'h21v06', 'h21v07', 'h22v01', 'h22v02', 'h22v03', 'h22v04', 'h22v05', 'h22v06', 'h22v07', 'h23v01', 'h23v02', 'h23v03', 'h23v04', 'h23v05', 'h23v06', 'h23v07', 'h24v02', 'h24v03', 'h24v04', 'h24v05', 'h24v06', 'h24v07', 'h25v02', 'h25v03', 'h25v04', 'h25v05', 'h25v06', 'h25v07', 'h25v08', 'h25v09', 'h26v02', 'h26v03', 'h26v04', 'h26v05', 'h26v06', 'h26v07', 'h26v08', 'h27v03', 'h27v04', 'h27v05', 'h27v06', 'h27v07', 'h27v08', 'h27v09', 'h28v10', 'h28v04', 'h28v05', 'h28v06', 'h28v07', 'h28v08', 'h28v09', 'h29v10', 'h29v05', 'h29v06', 'h29v07', 'h29v08', 'h29v09', 'h30v10', 'h30v07', 'h30v08', 'h30v09', 'h31v09', 'h32v10', 'h32v09')
namerica_tiles <- c('h10v02', 'h10v03', 'h10v04', 'h10v05', 'h10v06', 'h10v07', 'h10v08', 'h11v02', 'h11v03', 'h11v04', 'h11v05', 'h11v06', 'h11v07', 'h12v01', 'h12v02', 'h12v03', 'h12v04', 'h12v05', 'h12v07', 'h13v01', 'h13v02', 'h13v03', 'h13v04', 'h14v01', 'h14v02', 'h14v03', 'h14v04', 'h15v01', 'h15v02', 'h15v03', 'h16v00', 'h16v01', 'h16v02', 'h17v00', 'h17v01', 'h17v02', 'h28v03', 'h29v03', 'h06v03', 'h07v03', 'h07v05', 'h07v06', 'h07v07', 'h08v03', 'h08v04', 'h08v05', 'h08v06', 'h08v07', 'h09v02', 'h09v03', 'h09v04', 'h09v05', 'h09v06', 'h09v07', 'h09v08')
europe_tiles <- c('h15v05', 'h16v02', 'h16v05', 'h17v01', 'h17v02', 'h17v03', 'h17v04', 'h17v05', 'h18v00', 'h18v01', 'h18v02', 'h18v03', 'h18v04', 'h18v05', 'h19v00', 'h19v01', 'h19v02', 'h19v03', 'h19v04', 'h19v05', 'h20v01', 'h20v02', 'h20v03', 'h20v04', 'h20v05', 'h21v03', 'h21v04')
africa_tiles <- c('h15v07', 'h16v05', 'h16v06', 'h16v07', 'h16v08', 'h17v10', 'h17v05', 'h17v06', 'h17v07', 'h17v08', 'h18v05', 'h18v06', 'h18v07', 'h18v08', 'h18v09', 'h19v10', 'h19v11', 'h19v12', 'h19v05', 'h19v06', 'h19v07', 'h19v08', 'h19v09', 'h20v10', 'h20v11', 'h20v12', 'h20v05', 'h20v06', 'h20v07', 'h20v08', 'h20v09', 'h21v10', 'h21v11', 'h21v05', 'h21v06', 'h21v07', 'h21v08', 'h21v09', 'h22v10', 'h22v11', 'h22v07', 'h22v08', 'h22v09', 'h23v10', 'h23v11', 'h23v07', 'h23v08', 'h23v09')
samerica_tiles <- c('h10v10', 'h10v07', 'h10v08', 'h10v09', 'h11v10', 'h11v11', 'h11v12', 'h11v07', 'h11v08', 'h11v09', 'h12v10', 'h12v11', 'h12v12', 'h12v13', 'h12v08', 'h12v09', 'h13v10', 'h13v11', 'h13v12', 'h13v13', 'h13v14', 'h13v08', 'h13v09', 'h14v10', 'h14v11', 'h14v14', 'h14v09', 'h08v08', 'h08v09', 'h09v08', 'h09v09')
oceania_tiles <- c('h00v10', 'h00v08', 'h01v10', 'h01v11', 'h01v07', 'h01v09', 'h02v10', 'h02v06', 'h02v08', 'h27v10', 'h28v14', 'h29v13', 'h03v10', 'h03v11', 'h03v06', 'h03v07', 'h30v13', 'h31v12', 'h31v13', 'h31v08', 'h32v11', 'h32v12', 'h32v07', 'h32v09', 'h33v10', 'h33v11', 'h33v07', 'h33v08', 'h33v09', 'h34v10', 'h34v07', 'h34v08', 'h34v09', 'h35v10', 'h35v08', 'h35v09', 'h04v09', 'h05v13', 'h06v11', 'h08v11')
australia_tiles <- c('h27v11', 'h27v12', 'h27v14', 'h28v11', 'h28v12', 'h28v13', 'h29v10', 'h29v11', 'h29v12', 'h29v13', 'h30v10', 'h30v11', 'h30v12', 'h31v10', 'h31v11', 'h31v12', 'h32v10')

#-----------------------------------------------------------
# create headers for all the LW MASK files
mask_dir <- "/projectnb/modislc/data/mcd12_in/c5/ancillary_layers/C5_LW_Mask/lw_mask_500m"
lw_mask_files <- dir(mask_dir, pattern="*.bin", full=T)

for(in_file in lw_mask_files){
  print(paste("Writing header for:", basename(in_file)))
  WriteHDR(in_file)
}

#-----------------------------------------------------------
# create a continental LW mosaic
year_to_mosaic <- 2003
out_prefix <- "Europe"
# out_prefix <- "NAmerica"
# out_prefix <- "Africa"
out_sep <- "_"
layer_name <- "LWMask"
# band_to_mosaic <- which(layer_names == layer_name)
out_file <- file.path(out_dir, paste(paste(out_prefix, layer_name, sep=out_sep), ".vrt", sep=""))
# tiles_to_mosaic <- namerica_tiles
# tiles_to_mosaic <- africa_tiles
tiles_to_mosaic <- europe_tiles
mosaic_files <- unlist(lapply(tiles_to_mosaic, FindLWMask))
BuildVRT(mosaic_files, out_file=out_file, vrtnodata=32767)

# reproject
gdalwarp -r mode -s_srs "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m" -t_srs "+proj=latlong +datum=WGS84" -tr 0.025 0.025 -te -170 10 -50 75 /projectnb/modislc/users/joshgray/MCD12Q2C6/LW_MASKS/NAmerica_LWMask.vrt /projectnb/modislc/users/joshgray/MCD12Q2C6/LW_MASKS/NAmerica_LWMask_wgs.tif

gdalwarp -r mode -s_srs "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m" -t_srs "+proj=latlong +datum=WGS84" -tr 0.025 0.025 /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT/Africa_LWMask.vrt /projectnb/modislc/users/joshgray/MCD12Q2C6/LW_MASKS/Africa_LWMask_wgs.tif

gdalwarp -s_srs "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m" -t_srs "+proj=latlong +datum=WGS84" -tr 0.025 0.025 /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT/Africa_evi_amp_cycle1.vrt /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT/Africa_evi_amp_wgs.tif

# Maps:
# - N America midgup, number of cycles
# - Africa midgup, number of cycles
# - Europe midgup, number of cycles
#
# Time series:
# - desert SW
# - boreal N America
# - Florida live oaks
# - scandanavia evergreen
gdalwarp -s_srs "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m" -t_srs "+proj=latlong +datum=WGS84" -tr 0.025 0.025 -te -170 10 -50 75 /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT/NAmerica_evi_amp_cycle1.vrt /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT/NAmerica_evi_amp_wgs.tif

gdalwarp -dstnodata 0 -s_srs "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m" -t_srs "+proj=latlong +datum=WGS84" -tr 0.025 0.025 Europe_midgup_cycle1.vrt Europe_midgup_cycle1_wgs.tif
gdalwarp -dstnodata 0 -s_srs "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m" -t_srs "+proj=latlong +datum=WGS84" -tr 0.025 0.025 Africa_midgup_cycle1.vrt Africa_midgup_cycle1_wgs.tif

in_dir <- "/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT"
in_files <- dir(in_dir, pattern="MCD12Q2C6.*2003$", full=T)
good_size <- 1013760000
for(in_file in in_files){
  if(file.info(in_file)$size != good_size){
    print(paste(basename(in_file), file.info(in_file)$size))
    del_cmd <- paste("rm -f", in_file)
    system(del_cmd)
    tile <- unlist(strsplit(basename(in_file),split="_"))[2]
    sub_cmd <- paste("qsub -V -pe omp 16 -l h_rt=14:00:00 -l mem_total=98G /projectnb/modislc/users/joshgray/MCD12Q2C6/subAnnualPheno.sh", tile, "2003")
    system(sub_cmd)
  }
}
