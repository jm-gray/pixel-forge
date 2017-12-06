library(raster)
library(RColorBrewer)
library(tools)
library(argparse)
#---------------------------------------------------------------------
# run all continents:
# qsub -V -l mem_total=24G -l h_rt=24:00:00 ./run_continent_diagnostics.sh namerica
# qsub -V -l mem_total=24G -l h_rt=24:00:00 ./run_continent_diagnostics.sh asia
# qsub -V -l h_rt=24:00:00 ./run_continent_diagnostics.sh europe
# qsub -V -l h_rt=24:00:00 ./run_continent_diagnostics.sh samerica
# qsub -V -l h_rt=24:00:00 ./run_continent_diagnostics.sh africa
# qsub -V -l h_rt=24:00:00 ./run_continent_diagnostics.sh australia
# qsub -V -l h_rt=24:00:00 ./run_continent_diagnostics.sh oceania

# R --vanilla < /projectnb/modislc/users/joshgray/C6_Diagnostics/Mosaics/MCD12Q2C6_DiagnosticsContinents.R --args -continent $1

# #!/bin/bash
#
# echo Submitting continent $1
# R --vanilla < /projectnb/modislc/users/joshgray/C6_Diagnostics/Mosaics/MCD12Q2C6_DiagnosticsContinents.R --args -continent $1

#---------------------------------------------------------------------
BuildVRT <- function(file_list, out_file, band=NULL, vrtnodata=0){
  # Builds a VRT mosaic of files in file_list, possibly w/ band specification
  if(!is.null(band)){
    sys_cmd <- paste("gdalbuildvrt", paste("-vrtnodata", vrtnodata), paste("-b", band), out_file, paste(file_list, collapse=" "))
  }else{
    sys_cmd <- paste("gdalbuildvrt", paste("-vrtnodata", vrtnodata), out_file, paste(file_list, collapse=" "))
  }
  system(sys_cmd)
}

#---------------------------------------------------------------------
GetTile <- function(tile, metric, year, data_dir) list.files(file.path(data_dir, metric), pattern=paste(metric, "_", tile, "_", year, "$", sep=""), full=T)

#---------------------------------------------------------------------
GetLWTile <- function(tile) list.files("/projectnb/modislc/data/mcd12_in/c6/ancillary_layers/C6_LW_Mask/lw_mask_500m", pattern=paste("LW.map.", tile, "$", sep=""), full=T)

#---------------------------------------------------------------------
WriteHDR <- function(in_file){
  tile <- gsub("LW.map.(h[0-9]{2}v[0-9]{2})$", "\\1", basename(in_file))
  tile_h <- as.numeric(substr(tile, 2, 3))
  tile_v <- as.numeric(substr(tile, 5, 6))
  out_hdr <- file.path(dirname(in_file), paste(basename(in_file), ".hdr", sep=""))

  uly_map = 10007554.677
  ulx_map = -20015109.354
  lry_map = -10007554.677
  lrx_map = 20015109.354
  pix = 463.312716525
  dims = 2400

  # calculate the upper left corner of the current tile
  ulx = ulx_map + (tile_h * pix * dims)
  uly = uly_map - (tile_v * pix * dims)

  nbands=1

  temp_txt = paste("ENVI\nENVI description = { MODIS 500 m LW Mask }\nlines = ", dims, "\nsamples = ", dims, "\nbands = ", nbands, "\nheader offset = 0\nfile type = ENVI Standard\ndata type = 1\ninterleave = bip\nbyte order = 0\nmap info = {Sinusoidal, 1, 1,", ulx, ", ", uly, ", ", pix, ", ", pix, "}", "\ncoordinate system string = {PROJCS[\"Sinusoidal\",GEOGCS[\"GCS_unnamed ellipse\",DATUM[\"D_unknown\",SPHEROID[\"Unknown\",6371007.181,0]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Sinusoidal\"],PARAMETER[\"central_meridian\",0],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"Meter\",1]]}", sep="")


  # write the header to a file
  f <- file(out_hdr)
  writeLines(c(temp_txt), f)
  close(f)
}

PlotTile <- function(r, lwmask, cutoffs=c(0, 365), breaks=NULL, round_digs=0, pal=colorRampPalette(rev(brewer.pal(11, "Spectral"))), scale=1, title=NA, plot_legend=T, legend_only=F, legend_title="", pdf_out=NULL, plot_height=12, plotNEW=F, garish=F, continents=NA, MAXPIXELS=1.44e6,...){
  # MAXPIXELS <- 5.7e6
  WATERCOLOR <- rgb(0.7, 0.7, 0.7)
  LANDCOLOR <- rgb(0.2, 0.2, 0.2)
  LEGENDAXISCEX <- 1
  LEGENDMAINCEX <- 1
  LEGENDWIDTH <- 2.5

  if(is.null(breaks)){
    breaks <- round(c(minValue(r), seq(cutoffs[1], cutoffs[2]), maxValue(r)), round_digs)
    breaks <- unique(breaks)
  }

  if(legend_only){
    legend_at <- round(seq(breaks[2], breaks[length(breaks) - 1], len=7))
    legend_at_scaled <- round(legend_at / scale, round_digs)
    legend_labels <- c(paste("<", legend_at_scaled[1]), as.character(legend_at_scaled[2:(length(legend_at_scaled) - 1)]), paste(">", legend_at_scaled[length(legend_at_scaled)]))
    plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=pal(length(breaks)-1), legend.width=LEGENDWIDTH, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=LEGENDAXISCEX), legend.args=list(text=legend_title, side=3, font=2, line=0.5, cex=LEGENDMAINCEX), smallplot=c(0.2,0.5,0.2,0.9))
    # plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=pal(length(breaks)-1), legend.width=5, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=2), legend.args=list(text="", side=3, font=2, line=0.5, cex=2),  smallplot=c(0.1,0.5,0.2,0.9))
  }else{
    plot(lwmask, breaks=c(-1, 0.5, 1.5, 1e6), col=c(WATERCOLOR, LANDCOLOR, WATERCOLOR), xaxt="n", yaxt="n", legend=F, bty="n", box=FALSE, maxpixels=MAXPIXELS)
    plot(r, breaks=breaks, col=pal(length(breaks) - 1), maxpixels=MAXPIXELS, legend=F, xaxt="n", yaxt="n", bty="n", box=FALSE, add=T, ...)
    if(plot_legend){
      legend_at <- round(seq(breaks[2], breaks[length(breaks) - 1], len=7))
      legend_at_scaled <- round(legend_at / scale, round_digs)
      legend_labels <- c(paste("<", legend_at_scaled[1]), as.character(legend_at_scaled[2:(length(legend_at_scaled) - 1)]), paste(">", legend_at_scaled[length(legend_at_scaled)]))
      plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=pal(length(breaks)-1), legend.width=LEGENDWIDTH, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=LEGENDAXISCEX), legend.args=list(text="", side=3, font=2, line=0.5, cex=LEGENDMAINCEX))
    }
    title(title)
  }
}[]

#------------------------------------------------
# First create proper headers for all LW mask tiles
# all_in_files <- list.files("/projectnb/modislc/data/mcd12_in/c6/ancillary_layers/C6_LW_Mask/lw_mask_500m", pattern="LW.map.h[0-9]{2}v[0-9]{2}$", full=T)
# for(in_file in all_in_files){
#   # in_file <- GetLWTile(tile)
#   WriteHDR(in_file)
# }

#------------------------------------------------
# tile lists
# asia_tiles <- c('h10v02', 'h11v02', 'h12v01', 'h19v00', 'h19v01', 'h19v04', 'h20v01', 'h20v02', 'h20v03', 'h20v04', 'h20v05', 'h21v01', 'h21v02', 'h21v03', 'h21v04', 'h21v05', 'h21v06', 'h21v07', 'h22v01', 'h22v02', 'h22v03', 'h22v04', 'h22v05', 'h22v06', 'h22v07', 'h23v01', 'h23v02', 'h23v03', 'h23v04', 'h23v05', 'h23v06', 'h23v07', 'h24v02', 'h24v03', 'h24v04', 'h24v05', 'h24v06', 'h24v07', 'h25v02', 'h25v03', 'h25v04', 'h25v05', 'h25v06', 'h25v07', 'h25v08', 'h25v09', 'h26v02', 'h26v03', 'h26v04', 'h26v05', 'h26v06', 'h26v07', 'h26v08', 'h27v03', 'h27v04', 'h27v05', 'h27v06', 'h27v07', 'h27v08', 'h27v09', 'h28v10', 'h28v04', 'h28v05', 'h28v06', 'h28v07', 'h28v08', 'h28v09', 'h29v10', 'h29v05', 'h29v06', 'h29v07', 'h29v08', 'h29v09', 'h30v10', 'h30v07', 'h30v08', 'h30v09', 'h31v09', 'h32v10', 'h32v09')
# this asia tile leaves out ones that cross the meridian for plotting purposes
asia_tiles <- c('h19v00', 'h19v01', 'h19v04', 'h20v01', 'h20v02', 'h20v03', 'h20v04', 'h20v05', 'h21v01', 'h21v02', 'h21v03', 'h21v04', 'h21v05', 'h21v06', 'h21v07', 'h22v01', 'h22v02', 'h22v03', 'h22v04', 'h22v05', 'h22v06', 'h22v07', 'h23v01', 'h23v02', 'h23v03', 'h23v04', 'h23v05', 'h23v06', 'h23v07', 'h24v02', 'h24v03', 'h24v04', 'h24v05', 'h24v06', 'h24v07', 'h25v02', 'h25v03', 'h25v04', 'h25v05', 'h25v06', 'h25v07', 'h25v08', 'h25v09', 'h26v02', 'h26v03', 'h26v04', 'h26v05', 'h26v06', 'h26v07', 'h26v08', 'h27v03', 'h27v04', 'h27v05', 'h27v06', 'h27v07', 'h27v08', 'h27v09', 'h28v10', 'h28v04', 'h28v05', 'h28v06', 'h28v07', 'h28v08', 'h28v09', 'h29v10', 'h29v05', 'h29v06', 'h29v07', 'h29v08', 'h29v09', 'h30v10', 'h30v07', 'h30v08', 'h30v09', 'h31v09', 'h32v10', 'h32v09')
namerica_tiles <- c('h10v02', 'h10v03', 'h10v04', 'h10v05', 'h10v06', 'h10v07', 'h10v08', 'h11v02', 'h11v03', 'h11v04', 'h11v05', 'h11v06', 'h11v07', 'h12v01', 'h12v02', 'h12v03', 'h12v04', 'h12v05', 'h12v07', 'h13v01', 'h13v02', 'h13v03', 'h13v04', 'h14v01', 'h14v02', 'h14v03', 'h14v04', 'h15v01', 'h15v02', 'h15v03', 'h16v00', 'h16v01', 'h16v02', 'h17v00', 'h17v01', 'h17v02', 'h28v03', 'h29v03', 'h06v03', 'h07v03', 'h07v05', 'h07v06', 'h07v07', 'h08v03', 'h08v04', 'h08v05', 'h08v06', 'h08v07', 'h09v02', 'h09v03', 'h09v04', 'h09v05', 'h09v06', 'h09v07', 'h09v08')
europe_tiles <- c('h15v05', 'h16v02', 'h16v05', 'h17v01', 'h17v02', 'h17v03', 'h17v04', 'h17v05', 'h18v00', 'h18v01', 'h18v02', 'h18v03', 'h18v04', 'h18v05', 'h19v00', 'h19v01', 'h19v02', 'h19v03', 'h19v04', 'h19v05', 'h20v01', 'h20v02', 'h20v03', 'h20v04', 'h20v05', 'h21v03', 'h21v04')
africa_tiles <- c('h15v07', 'h16v05', 'h16v06', 'h16v07', 'h16v08', 'h17v10', 'h17v05', 'h17v06', 'h17v07', 'h17v08', 'h18v05', 'h18v06', 'h18v07', 'h18v08', 'h18v09', 'h19v10', 'h19v11', 'h19v12', 'h19v05', 'h19v06', 'h19v07', 'h19v08', 'h19v09', 'h20v10', 'h20v11', 'h20v12', 'h20v05', 'h20v06', 'h20v07', 'h20v08', 'h20v09', 'h21v10', 'h21v11', 'h21v05', 'h21v06', 'h21v07', 'h21v08', 'h21v09', 'h22v10', 'h22v11', 'h22v07', 'h22v08', 'h22v09', 'h23v10', 'h23v11', 'h23v07', 'h23v08', 'h23v09')
samerica_tiles <- c('h10v10', 'h10v07', 'h10v08', 'h10v09', 'h11v10', 'h11v11', 'h11v12', 'h11v07', 'h11v08', 'h11v09', 'h12v10', 'h12v11', 'h12v12', 'h12v13', 'h12v08', 'h12v09', 'h13v10', 'h13v11', 'h13v12', 'h13v13', 'h13v14', 'h13v08', 'h13v09', 'h14v10', 'h14v11', 'h14v14', 'h14v09', 'h08v08', 'h08v09', 'h09v08', 'h09v09')
oceania_tiles <- c('h00v10', 'h00v08', 'h01v10', 'h01v11', 'h01v07', 'h01v09', 'h02v10', 'h02v06', 'h02v08', 'h27v10', 'h28v14', 'h29v13', 'h03v10', 'h03v11', 'h03v06', 'h03v07', 'h30v13', 'h31v12', 'h31v13', 'h31v08', 'h32v11', 'h32v12', 'h32v07', 'h32v09', 'h33v10', 'h33v11', 'h33v07', 'h33v08', 'h33v09', 'h34v10', 'h34v07', 'h34v08', 'h34v09', 'h35v10', 'h35v08', 'h35v09', 'h04v09', 'h05v13', 'h06v11', 'h08v11')
australia_tiles <- c('h27v11', 'h27v12', 'h27v14', 'h28v11', 'h28v12', 'h28v13', 'h29v10', 'h29v11', 'h29v12', 'h29v13', 'h30v10', 'h30v11', 'h30v12', 'h31v10', 'h31v11', 'h31v12', 'h32v10')

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Create mosaics and plot
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

arg_parser <- ArgumentParser()
arg_parser$add_argument("-continent", type="character") # tile to process
arg_parser$add_argument("-metrics", type="character", nargs="*", default="all") # tile to process
arg_parser$add_argument("-years", type="integer", nargs="*", default=0) # tile to process
args <- arg_parser$parse_args()
# args <- arg_parser$parse_args(c("-continent","europe", "-metrics","Greenup", "MidGreenup", "Dormancy","-years","2011", "2014"))
all_metrics <- c("Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy")

# check for special value of years "0" which means "all years"
if(args$years == 0){
  years <- 2001:2015
}else{
  years <- args$years
}

# check for special value of metrics "all" which means "all metrics"
if(tolower(args$metrics) == "all"){
  metrics <- all_metrics
}else{
  metrics <- args$metrics
}

if(tolower(args$continent) == "asia"){
  continent_list <- asia_tiles
}else if(tolower(args$continent) == "namerica"){
  continent_list <- namerica_tiles
}else if(tolower(args$continent) == "europe"){
  continent_list <- europe_tiles
}else if(tolower(args$continent) == "africa"){
  continent_list <- africa_tiles
}else if(tolower(args$continent) == "samerica"){
  continent_list <- samerica_tiles
}else if(tolower(args$continent) == "oceania"){
  continent_list <- oceania_tiles
}else if(tolower(args$continent) == "australia"){
  continent_list <- australia_tiles
}else{
  print("not a recognized continent")
  quit()
}

print(paste("Doing:", tolower(args$continent)))

# continent_lists <- list(asia_tiles, namerica_tiles, europe_tiles, africa_tiles, samerica_tiles, oceania_tiles, australia_tiles)
# continent_names <- c("Asia", "N America", "Europe", "Africa", "S America", "Oceania", "Australia")
# doy_metrics <- c("Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy")
# data_dir <- "/projectnb/modislc/data/mcd12_out/phen_out/c6"
data_dir="/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6"
out_dir <- "/projectnb/modislc/users/joshgray/C6_Diagnostics/Mosaics"
band_to_mosaic <- 1
# years <- 2001:2014
i <- 1
tiles_to_mosaic <- continent_list
out_prefix <- tolower(args$continent)
# pdf(file.path(out_dir, paste("C6_Diagnostics_", out_prefix, ".pdf", sep="")), height=15, width=15)
for(metric_name in metrics){
  for(year in years){
    doy_offset <- as.Date(paste(year, "-1-1", sep="")) - as.Date("1970-1-1")
    out_file <- file.path(out_dir, paste(paste(out_prefix, metric_name, year, sep="_"), ".vrt", sep=""))
    mosaic_files <- unlist(lapply(tiles_to_mosaic, GetTile, metric=metric_name, year=year, data_dir=data_dir))
    BuildVRT(mosaic_files, out_file=out_file, band=band_to_mosaic, vrtnodata=32767)
    out_lwmask_file <- file.path(out_dir, paste(paste(out_prefix, "LWMASK", sep="_"), ".vrt", sep=""))
    if(!file.exists(out_lwmask_file)){
      lwmask_mosaic_files <- unlist(lapply(tiles_to_mosaic, GetLWTile))
      BuildVRT(lwmask_mosaic_files, out_file=out_lwmask_file, band=1, vrtnodata=32767)
    }
    # make the plot
    pdf(file.path(out_dir, paste("C6_Diagnostics_", out_prefix, "_", metric_name, "_", year, ".pdf", sep="")), height=15, width=15)
    r <- raster(out_file) - as.integer(doy_offset)
    lw_mask <- raster(out_lwmask_file)
    qs <- quantile(r, c(0.02, 0.98), na.rm=T)
    plot_title <- paste(out_prefix, metric_name, year)
    PlotTile(r, lw_mask, cutoffs=c(qs[1], qs[2]), MAXPIXELS=5e6, title=plot_title)
    dev.off()
  }
}
# dev.off() # close the continent file
