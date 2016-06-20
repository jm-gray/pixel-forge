#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Josh Gray, Boston University, 2016
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Prelims
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
library(raster)
library(RColorBrewer)

# source plotting functions
source("/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_AnnualPhenologyFunctions.R")

# tile lists
asia_tiles <- c('h10v02', 'h11v02', 'h12v01', 'h19v00', 'h19v01', 'h19v04', 'h20v01', 'h20v02', 'h20v03', 'h20v04', 'h20v05', 'h21v01', 'h21v02', 'h21v03', 'h21v04', 'h21v05', 'h21v06', 'h21v07', 'h22v01', 'h22v02', 'h22v03', 'h22v04', 'h22v05', 'h22v06', 'h22v07', 'h23v01', 'h23v02', 'h23v03', 'h23v04', 'h23v05', 'h23v06', 'h23v07', 'h24v02', 'h24v03', 'h24v04', 'h24v05', 'h24v06', 'h24v07', 'h25v02', 'h25v03', 'h25v04', 'h25v05', 'h25v06', 'h25v07', 'h25v08', 'h25v09', 'h26v02', 'h26v03', 'h26v04', 'h26v05', 'h26v06', 'h26v07', 'h26v08', 'h27v03', 'h27v04', 'h27v05', 'h27v06', 'h27v07', 'h27v08', 'h27v09', 'h28v10', 'h28v04', 'h28v05', 'h28v06', 'h28v07', 'h28v08', 'h28v09', 'h29v10', 'h29v05', 'h29v06', 'h29v07', 'h29v08', 'h29v09', 'h30v10', 'h30v07', 'h30v08', 'h30v09', 'h31v09', 'h32v10', 'h32v09')
namerica_tiles <- c('h10v02', 'h10v03', 'h10v04', 'h10v05', 'h10v06', 'h10v07', 'h10v08', 'h11v02', 'h11v03', 'h11v04', 'h11v05', 'h11v06', 'h11v07', 'h12v01', 'h12v02', 'h12v03', 'h12v04', 'h12v05', 'h12v07', 'h13v01', 'h13v02', 'h13v03', 'h13v04', 'h14v01', 'h14v02', 'h14v03', 'h14v04', 'h15v01', 'h15v02', 'h15v03', 'h16v00', 'h16v01', 'h16v02', 'h17v00', 'h17v01', 'h17v02', 'h28v03', 'h29v03', 'h06v03', 'h07v03', 'h07v05', 'h07v06', 'h07v07', 'h08v03', 'h08v04', 'h08v05', 'h08v06', 'h08v07', 'h09v02', 'h09v03', 'h09v04', 'h09v05', 'h09v06', 'h09v07', 'h09v08')
europe_tiles <- c('h15v05', 'h16v02', 'h16v05', 'h17v01', 'h17v02', 'h17v03', 'h17v04', 'h17v05', 'h18v00', 'h18v01', 'h18v02', 'h18v03', 'h18v04', 'h18v05', 'h19v00', 'h19v01', 'h19v02', 'h19v03', 'h19v04', 'h19v05', 'h20v01', 'h20v02', 'h20v03', 'h20v04', 'h20v05', 'h21v03', 'h21v04')
africa_tiles <- c('h15v07', 'h16v05', 'h16v06', 'h16v07', 'h16v08', 'h17v10', 'h17v05', 'h17v06', 'h17v07', 'h17v08', 'h18v05', 'h18v06', 'h18v07', 'h18v08', 'h18v09', 'h19v10', 'h19v11', 'h19v12', 'h19v05', 'h19v06', 'h19v07', 'h19v08', 'h19v09', 'h20v10', 'h20v11', 'h20v12', 'h20v05', 'h20v06', 'h20v07', 'h20v08', 'h20v09', 'h21v10', 'h21v11', 'h21v05', 'h21v06', 'h21v07', 'h21v08', 'h21v09', 'h22v10', 'h22v11', 'h22v07', 'h22v08', 'h22v09', 'h23v10', 'h23v11', 'h23v07', 'h23v08', 'h23v09')
samerica_tiles <- c('h10v10', 'h10v07', 'h10v08', 'h10v09', 'h11v10', 'h11v11', 'h11v12', 'h11v07', 'h11v08', 'h11v09', 'h12v10', 'h12v11', 'h12v12', 'h12v13', 'h12v08', 'h12v09', 'h13v10', 'h13v11', 'h13v12', 'h13v13', 'h13v14', 'h13v08', 'h13v09', 'h14v10', 'h14v11', 'h14v14', 'h14v09', 'h08v08', 'h08v09', 'h09v08', 'h09v09')
oceania_tiles <- c('h00v10', 'h00v08', 'h01v10', 'h01v11', 'h01v07', 'h01v09', 'h02v10', 'h02v06', 'h02v08', 'h27v10', 'h28v14', 'h29v13', 'h03v10', 'h03v11', 'h03v06', 'h03v07', 'h30v13', 'h31v12', 'h31v13', 'h31v08', 'h32v11', 'h32v12', 'h32v07', 'h32v09', 'h33v10', 'h33v11', 'h33v07', 'h33v08', 'h33v09', 'h34v10', 'h34v07', 'h34v08', 'h34v09', 'h35v10', 'h35v08', 'h35v09', 'h04v09', 'h05v13', 'h06v11', 'h08v11')
australia_tiles <- c('h27v11', 'h27v12', 'h27v14', 'h28v11', 'h28v12', 'h28v13', 'h29v10', 'h29v11', 'h29v12', 'h29v13', 'h30v10', 'h30v11', 'h30v12', 'h31v10', 'h31v11', 'h31v12', 'h32v10')

# MCD12Q2C6 band names
layer_names <- c("num_cycles", "fill_code", "evi_area_cycle1", "evi_amp_cycle1", "evi_min_cycle1", "frac_filled_gup_cycle1", "frac_filled_gdown_cycle1", "length_gup_cycle1", "length_gdown_cycle1", "ogi_cycle1", "midgup_cycle1", "mat_cycle1", "peak_cycle1", "sen_cycle1", "midgdown_cycle1", "dor_cycle1", "ogi_qual_cycle1", "midgup_qual_cycle1", "mat_qual_cycle1", "peak_qual_cycle1", "sen_qual_cycle1", "midgdown_qual_cycle1", "dor_qual_cycle1", "evi_area_cycle2", "evi_amp_cycle2", "evi_min_cycle2", "frac_filled_gup_cycle2", "frac_filled_gdown_cycle2", "length_gup_cycle2", "length_gdown_cycle2", "ogi_cycle2", "midgup_cycle2", "mat_cycle2", "peak_cycle2", "sen_cycle2", "midgdown_cycle2", "dor_cycle2", "ogi_qual_cycle2", "midgup_qual_cycle2", "mat_qual_cycle2", "peak_qual_cycle2", "sen_qual_cycle2", "midgdown_qual_cycle2", "dor_qual_cycle2")

# input/output data directories
data_dir <- "/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT"
out_dir <- "/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT"


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Create mosaics and plot
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# create mosaic of europe
year_to_mosaic <- 2003
# out_prefix <- "Africa"
out_prefix <- "Europe"
# out_prefix <- "NAmerica"
out_sep <- "_"
layer_name <- "midgup_cycle1"
# layer_name <- "evi_amp_cycle1"
band_to_mosaic <- which(layer_names == layer_name)
out_file <- file.path(out_dir, paste(paste(out_prefix, layer_name, sep=out_sep), ".vrt", sep=""))
# tiles_to_mosaic <- namerica_tiles
# tiles_to_mosaic <- africa_tiles
tiles_to_mosaic <- europe_tiles
mosaic_files <- unlist(lapply(tiles_to_mosaic, CheckForTile, year=year_to_mosaic, data_dir=data_dir))
BuildVRT(mosaic_files, out_file=out_file, band=band_to_mosaic, vrtnodata=32767)

# resample/reproject
# gdalwarp -s_srs "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m" -t_srs "+proj=latlong +datum=WGS84" -tr 0.025 0.025 -te -170 10 -50 75 /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT/NAmerica_evi_amp_cycle1.vrt /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT/NAmerica_evi_amp_wgs.tif

# gdalwarp -s_srs "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m" -t_srs "+proj=latlong +datum=WGS84" -tr 0.025 0.025 /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT/Africa_midgup_cycle1.vrt /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT/Africa_midgup_cycle1_wgs.tif

gdalwarp -s_srs "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m" -t_srs "+proj=latlong +datum=WGS84" -tr 0.025 0.025 -te -15 40 35 72 /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT/Europe_evi_amp_cycle1.vrt /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT/Europe_evi_amp_cycle1_wgs.tif

gdalwarp -s_srs "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m" -t_srs "+proj=latlong +datum=WGS84" -tr 0.025 0.025 -te -15 40 35 72 /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT/Europe_midgup_cycle1.vrt /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT/Europe_midgup_cycle1_wgs.tif

gdalwarp -r mode -s_srs "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m" -t_srs "+proj=latlong +datum=WGS84" -tr 0.025 0.025 -te -15 40 35 72 /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT/Europe_LWMask.vrt /projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT/Europe_LWMask_wgs.tif


# plot Europe
r_mosaic <- raster(out_file)
NAvalue(r_mosaic) <- 32767
# PlotDOY(europe_mosaic, 2003, cutoffs=c(30, 200), maxpixels=1e6)
out_pdf <- file.path(out_dir, paste(paste(out_prefix, layer_name, year_to_mosaic, sep=out_sep), ".pdf", sep=""))
PlotDOY(r_mosaic, 2003, cutoffs=c(30, 200), pdf_out=out_pdf, maxpixels=4e6)

# plot NAmerica midgup
# midgup <- raster("~/Desktop/NAmerica_midgup_cycle1_wgs.tif")
# lw <- raster("~/Desktop/NAmerica_LWMask_wgs.tif")
# amp <- raster("~/Desktop/NAmerica_evi_amp_wgs.tif")
# midgup <- raster("~/Desktop/Africa_midgup_cycle1_wgs.tif")
# lw <- raster("~/Desktop/Africa_LWMask_wgs.tif")
# amp <- raster("~/Desktop/Africa_evi_amp_wgs.tif")
midgup <- raster("~/Desktop/Europe_midgup_cycle1_wgs.tif")
lw <- raster("~/Desktop/Europe_LWMask_wgs.tif")
amp <- raster("~/Desktop/Europe_evi_amp_cycle1_wgs.tif")


midgup[lw != 1] <- NA
midgup[amp <= 350] <- NA
# out_pdf <- "~/Desktop/NAmerica_midgup.pdf"
# out_pdf <- "~/Desktop/Africa_midgup.pdf"
out_pdf <- "~/Desktop/Europe_midgup.pdf"
continents <- shapefile("~/Downloads/continents/continent.shp")
# plot(subset(continents, CONTINENT=="North America"), border="NA", col="gray", xlim=extent(midgup)[1:2], ylim=extent(midgup)[3:4])
# PlotDOY(midgup, 2003, cutoffs=c(60,190), pdf_out=out_pdf, maxpixels=6e6, continents=subset(continents, CONTINENT=="North America"))
# PlotDOY(midgup, 2003, cutoffs=c(-89,288), pdf_out=out_pdf, maxpixels=6e6, continents=subset(continents, CONTINENT=="Africa"))
PlotDOY(midgup, 2003, cutoffs=c(67,183), pdf_out=out_pdf, maxpixels=6e6, continents=subset(continents, CONTINENT=="Europe"))


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Plotting for single tiles and time series
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
year_of_interest <- 2003
tile_of_interest <- "h12v04"
# layer_name <- "ogi_qual_cycle1"
layer_name <- "num_cycles"
out_sep <- "_"
out_pdf <- file.path(out_dir, paste(paste(tile_of_interest, layer_name, year_of_interest, sep=out_sep), ".pdf", sep=""))
input_file <- CheckForTile(tile=tile_of_interest, year=year_of_interest, data_dir=data_dir)
band_to_plot <- which(layer_names == layer_name)
s <- stack(input_file)
r <- raster(s, band_to_plot)
NAvalue(r) <- 32767

mask_dir <- "/projectnb/modislc/data/mcd12_in/c5/ancillary_layers/C5_LW_Mask/lw_mask_500m"
lw_mask_file <- dir(mask_dir, pattern=paste(".*", tile_of_interest, ".*", ".bin$", sep=""), full=T)
lw <- raster(lw_mask_file)
r[lw != 1] <- NA

# create a 2% linear stretch for plotting
plot_cutoffs <- quantile(r, c(0.02, 0.98), na.rm=T) - as.numeric(as.Date(paste(year_of_interest, "-1-1", sep="")))

PlotDOY(r, year_of_interest, cutoffs=plot_cutoffs, pdf_out=out_pdf, plot_height=20, maxpixels=6e6)

# PlotDOY(r, year_of_interest, cutoffs=plot_cutoffs, plot_height=8, maxpixels=2.5e6)
# for non-DOY:
# plot_cutoffs <- quantile(r, c(0.02, 0.98), na.rm=T)
# PlotStretch(r, cutoffs=plot_cutoffs, plot_height=8, maxpixels=2.5e6)

# Plot the QA
qa_breaks <- c(0, 0.75, 0.85, 0.95, 1)
# qa_cols <- c("green", "yellow", "orange", "red")
qa_cols <- rev(brewer.pal(4, "PuRd"))
pdf(h=10, w=12, file="/projectnb/modislc/users/joshgray/MCD12Q2C6/MCD12Q2C6_OUTPUT/h13v03_ogi_qual_plot.pdf")
par(col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white", bg="black", mar=c(3.5, 3.5, 1, 1), oma=rep(0, 4))
plot(r / 1e4, breaks=qa_breaks, col=qa_cols, maxpixels=5e6)
dev.off()

#------------------------------------
# click and plot some time series

myc <- click(r, cell=T)
points(xyFromCell(r, myc$cell), pch=1, col="black")
text(xyFromCell(r, myc$cell), labels=1:dim(myc)[1], pos=4, cex=1.5, col="black")

# get all data files for tile
data_prefix <- "c6_str5"
# spline_data_dir <- "/projectnb/modislc/users/dsm/spline_codes/spline_one_v3/outputs_global"
spline_data_dir <- "/projectnb/modislc/users/dsm/spline_codes/spline_one_v3/outputs"
patt <- paste(data_prefix, "\\..*", tile_of_interest, ".*evi2.bip", sep="")
evi2_files <- dir(spline_data_dir, pattern=patt, full=T)
patt <- paste(data_prefix, "\\..*", tile_of_interest, ".*evi2.residual.bip", sep="")
resid_files <- dir(spline_data_dir, pattern=patt, full=T)
patt <- paste(data_prefix, "\\..*", tile_of_interest, ".*evi2.flag.bip", sep="")
snow_files <- dir(spline_data_dir, pattern=patt, full=T)

# create 3 years worth of daily dates, IGNORING LEAP YEARS!
c6_dates <- as.Date(paste(rep((year_of_interest - 1):(year_of_interest + 1), each=365), rep(1:365, 3), sep="-"), format="%Y-%j")

pixel_time_series <- list()
for(i in 1:length(myc$cell)){
  tmp_row <- rowFromCell(r, myc$cell[i])
  tmp_col <- colFromCell(r, myc$cell[i])

  print(paste("Reading pixel", i, "of", length(myc$cell)))
  pixel_time_series[[i]] <- Get3YearDataChunk(evi2_files, resid_files, snow_files, year_of_interest, start_line=tmp_row, lines_to_read=1)[tmp_col, ]
}

x11(h=8,w=14)
layout(matrix(1:length(pixel_time_series), ncol=2, byrow=T))
par(mar=c(2,4,0,1), oma=rep(1,4))
for(i in 1:length(pixel_time_series)){
  PlotSeries(pixel_time_series[[i]], c6_dates, ylim=c(0,0.7), plot_legend=F)
  text(par()$usr[2], par()$usr[4], label=i, cex=3, adj=c(1.2, 1.2))
}
