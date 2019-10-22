library(raster)
library(RColorBrewer)
library(tools)
library(argparse)
library(viridis)
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
  return(out_file)
}

#---------------------------------------------------------------------
ReprojectContinent <- function(vrt_file, out_dir, cut_continent_path, out_res=5000){
    cut_continent <- shapefile(cut_continent_path)

    # first gdalwarp: go to resampled resolution
    resample_out_file <- file.path(out_dir, paste(file_path_sans_ext(basename(vrt_file)), "_resampled.tif", sep=""))
    gdal_cmd <- paste("gdalwarp -tr", out_res, out_res, vrt_file, resample_out_file)
    system(gdal_cmd)
    # second gdalwarp: reproject and crop to cutline
    final_out_file <- file.path(out_dir, paste(file_path_sans_ext(basename(vrt_file)), "_reproj_cut.tif", sep=""))
    gdal_cmd <- paste("gdalwarp -t_srs '", projection(cut_continent), "' -crop_to_cutline -cutline ", cut_continent_path, " ", resample_out_file, " ", final_out_file, sep="")
    system(gdal_cmd)

    system(paste("rm", resample_out_file))
    return(final_out_file)
}

#---------------------------------------------------------------------
GetTile <- function(tile, metric, data_dir) dir(data_dir, pattern=paste(tile, "\\.", metric, "\\.tif$", sep=""), full=T)

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
eurasia_tiles <- c(asia_tiles, europe_tiles)
oceania_tiles <- c(australia_tiles, oceania_tiles)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Create mosaics and plot
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

arg_parser <- ArgumentParser()
arg_parser$add_argument("-continent", type="character") # tile to process
arg_parser$add_argument("-metrics", type="character", nargs="*", default="all") # tile to process
arg_parser$add_argument("-years", type="integer", nargs="*", default=0) # tile to process
args <- arg_parser$parse_args()
# args <- arg_parser$parse_args(c("-continent","europe", "-metrics","Greenup", "MidGreenup", "Dormancy","-years","2011", "2014"))
# all_metrics <- c("Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy")
all_metrics <- c("Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy", "EVI_Area", "EVI_Amplitude", "EVI_Minimum")
doy_metrics <- c(rep(TRUE, 7), rep(FALSE, 3))
layer_names <- c("Median", "MAD", "Average", "Std. Dev.", "TS-slope", "TS p-value", "NumVals")

# check for special value of metrics "all" which means "all metrics"
if(tolower(args$metrics) == "all"){
  metrics <- all_metrics
}else{
  metrics <- args$metrics
}

cut_continent_dir <- "/rsstu/users/j/jmgray2/SEAL/INCA/CutContinents"
# cut_continent_dir <- "/Volumes/users/j/jmgray2/SEAL/INCA/CutContinents"
dir(cut_continent_dir)
if(tolower(args$continent) == "asia"){
  continent_list <- asia_tiles
  the_continent_shp_file <- ""
}else if(tolower(args$continent) == "namerica"){
  continent_list <- namerica_tiles
  the_continent_shp_file <- file.path(cut_continent_dir, "INCA_NAmerica_aea_cutline.shp")
}else if(tolower(args$continent) == "europe"){
  continent_list <- europe_tiles
  the_continent_shp_file <- ""
}else if(tolower(args$continent) == "africa"){
  continent_list <- africa_tiles
  the_continent_shp_file <- ""
}else if(tolower(args$continent) == "samerica"){
  continent_list <- samerica_tiles
  the_continent_shp_file <- ""
}else if(tolower(args$continent) == "oceania"){
  continent_list <- oceania_tiles
  the_continent_shp_file <- ""
}else if(tolower(args$continent) == "australia"){
  continent_list <- australia_tiles
  the_continent_shp_file <- ""
}else{
  print("not a recognized continent")
  quit()
}

print(paste("Doing:", tolower(args$continent)))

# continent_lists <- list(asia_tiles, namerica_tiles, europe_tiles, africa_tiles, samerica_tiles, oceania_tiles, australia_tiles)
# continent_names <- c("Asia", "N America", "Europe", "Africa", "S America", "Oceania", "Australia")
# doy_metrics <- c("Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy")
# data_dir <- "/projectnb/modislc/data/mcd12_out/phen_out/c6"
# data_dir="/projectnb/modislc/users/dsm/eval_modis_lc_061917/MCD12I6"
data_dir <- "/rsstu/users/j/jmgray2/SEAL/INCA/INCAglobaloutput"
# data_dir <- "/Volumes/users/j/jmgray2/SEAL/INCA/INCAglobaloutput"
out_dir <- "/rsstu/users/j/jmgray2/SEAL/INCA/INCAglobaloutput/continent_mosaics"
# out_dir <- "/Volumes/users/j/jmgray2/SEAL/INCA/INCAglobaloutput/continent_mosaics"

tiles_to_mosaic <- continent_list
out_prefix <- tolower(args$continent)

MAXPIXELS <- 2.5e6
PVALUETHRESH <- 0.05
MEDIANCUTOFFS <- NULL
MEDIANCUTOFFSQUANTS <- c(0.02, 0.98)
SLOPECUTOFFS <- NULL
SLOPECUTOFFSQUANTS <- c(0.05, 0.95)
MADCUTOFFS <- c(0, 21)
MADCUTOFFSQUANTS <- c(0, 0.98)
LANDCOLOR <- rgb(0.5, 0.5, 0.5)
WATERCOLOR <- rgb(0.2, 0.2, 0.2)
LEGENDAXISCEX <- 1
LEGENDMAINCEX <- 1
LEGENDWIDTH <- 2.5


# for(metric_name in metrics){
for(metric_name in metrics[!doy_metrics]){
    # check if we're doing a DOY metric or not
    if(metric_name %in% metrics[doy_metrics]){
        doing_doy <- TRUE
    }else{
        doing_doy <- FALSE
    }
    # mosaic, resample, and reproject
    out_res <- 2500
    out_vrt_file <- file.path(out_dir, paste(paste(out_prefix, metric_name, sep="_"), ".vrt", sep=""))
    mosaic_files <- unlist(lapply(tiles_to_mosaic, GetTile, metric=metric_name, data_dir=data_dir))
    out_vrt_file <- BuildVRT(mosaic_files, out_file=out_vrt_file, vrtnodata=32767)
    out_file <- ReprojectContinent(vrt_file=out_vrt_file, out_dir=out_dir, cut_continent_path=the_continent_shp_file, out_res=out_res)
    
    # make the plots
    s <- stack(out_file)
    continent_shp <- shapefile(the_continent_shp_file)
    pdf_height <- nrow(s) / 120
    pdf_width <- ncol(s) / 120
    pdf(file.path(out_dir, paste("INCA_", out_prefix, "_", metric_name, ".pdf", sep="")), height=pdf_height, width=pdf_width)
    
    if(doing_doy){
        #----------------------------------------
        # plot median
        med_qs <- quantile(raster(s, 1), c(0, 0.02, 0.98, 1))
        med_breaks <- unique(round(c(med_qs[1], seq(med_qs[2], med_qs[3], len=254), med_qs[4])))
        par(mar=rep(1, 4), oma=c(0, 0, 2, 0), bg=WATERCOLOR, col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white", cex.main=1.5)
        plot(continent_shp, border=NA, col=LANDCOLOR, box=F)
        plot(raster(s, 1), breaks=med_breaks, col=plasma(length(med_breaks) - 1), maxpixels=ncell(s), legend=F, xaxt="n", yaxt="n", bty="n", box=FALSE, add=T)
        legend_at <- round(seq(med_breaks[2], med_breaks[length(med_breaks) - 1], len=7))
        legend_labels <- c(paste("<", legend_at[1], sep=""), as.character(legend_at[2:(length(legend_at) - 1)]), paste(">", legend_at[length(legend_at)], sep=""))
        plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=plasma(length(med_breaks) - 1), legend.width=LEGENDWIDTH, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=LEGENDAXISCEX), legend.args=list(text="DOY", side=3, font=2, line=0.5, cex=LEGENDMAINCEX))
        title(paste("MCD12Q2", metric_name, "Median"), cex=2)

        #----------------------------------------
        # plot MAD
        mad_qs <- quantile(raster(s, 2), c(0, 0.02, 0.98, 1))
        mad_breaks <- c(seq(0, 42), mad_qs[4])
        par(mar=rep(1, 4), oma=c(0, 0, 2, 0), bg=WATERCOLOR, col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white", cex.main=1.5)
        plot(continent_shp, border=NA, col=LANDCOLOR, box=F)
        plot(raster(s, 2), breaks=mad_breaks, col=viridis(length(mad_breaks) - 1), maxpixels=ncell(s), legend=F, xaxt="n", yaxt="n", bty="n", box=FALSE, add=T)
        legend_at <- round(seq(0, mad_breaks[length(mad_breaks) - 1], len=7))
        legend_labels <- c(legend_at[1], as.character(legend_at[2:(length(legend_at) - 1)]), paste(">", legend_at[length(legend_at)], sep=""))
        plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=viridis(length(med_breaks) - 1), legend.width=LEGENDWIDTH, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=LEGENDAXISCEX), legend.args=list(text="days", side=3, font=2, line=0.5, cex=LEGENDMAINCEX))
        title(paste("MCD12Q2", metric_name, "MAD"), cex=2)

        #----------------------------------------
        # plot TS slope
        slope_r <- raster(s, 5)
        pvalue_r <- raster(s, 6)
        slope_r[pvalue_r >= 0.05] <- NA
        slope_qs <- quantile(slope_r, c(0, 0.02, 0.98, 1))
        neg_breaks <- c(-max(abs(slope_qs[c(1, 4)])), seq(-3, 0, len=128))
        pos_breaks <- c(seq(0, 3, len=128), max(abs(slope_qs[c(1, 4)])))
        slope_breaks <- c(neg_breaks[1:(length(neg_breaks) - 1)], pos_breaks[2:length(pos_breaks)])
        slope_pal <- colorRampPalette(rev(brewer.pal(9, "RdBu")))
        par(mar=rep(1, 4), oma=c(0, 0, 2, 0), bg=WATERCOLOR, col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white", cex.main=1.5)
        plot(continent_shp, border=NA, col=LANDCOLOR, box=F)
        plot(slope_r, breaks=slope_breaks, col=slope_pal(length(slope_breaks) - 1), maxpixels=ncell(s), legend=F, xaxt="n", yaxt="n", bty="n", box=FALSE, add=T)
        legend_at <- round(seq(slope_breaks[2], slope_breaks[length(slope_breaks) - 1], len=7))
        legend_labels <- c(paste("<", legend_at[1], sep=""), as.character(legend_at[2:(length(legend_at) - 1)]), paste(">", legend_at[length(legend_at)], sep=""))
        plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=slope_pal(length(slope_breaks) - 1), legend.width=LEGENDWIDTH, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=LEGENDAXISCEX), legend.args=list(text="days/yr", side=3, font=2, line=0.5, cex=LEGENDMAINCEX))
        title(paste("MCD12Q2", metric_name, "Sig. TS-slopes"), cex=2)

        dev.off()
    }else{
        if(metric_name == "EVI_Area"){
            legend_text <- "sum(EVI2)"
        }else{
            legend_text <- "EVI2"
        }
        #----------------------------------------
        # plot median
        med_qs <- quantile(raster(s, 1), c(0, 0.02, 0.98, 1))
        med_breaks <- c(med_qs[1], seq(med_qs[2], med_qs[3], len=254), med_qs[4])
        par(mar=rep(1, 4), oma=c(0, 0, 2, 0), bg=WATERCOLOR, col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white", cex.main=1.5)
        plot(continent_shp, border=NA, col=LANDCOLOR, box=F)
        plot(raster(s, 1), breaks=med_breaks, col=plasma(length(med_breaks) - 1), maxpixels=ncell(s), legend=F, xaxt="n", yaxt="n", bty="n", box=FALSE, add=T)
        # legend_at <- round(seq(med_breaks[2], med_breaks[length(med_breaks) - 1], len=7))
        legend_at <- round(seq(med_breaks[2], med_breaks[length(med_breaks) - 1], len=7), digits=2)
        legend_labels <- c(paste("<", legend_at[1], sep=""), as.character(legend_at[2:(length(legend_at) - 1)]), paste(">", legend_at[length(legend_at)], sep=""))
        plot(raster(matrix(c(legend_at[1], legend_at[length(legend_at)]))), legend.only=T, col=plasma(length(med_breaks) - 1), legend.width=LEGENDWIDTH, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=LEGENDAXISCEX), legend.args=list(text=legend_text, side=3, font=2, line=0.5, cex=LEGENDMAINCEX))
        title(paste("MCD12Q2", metric_name, "Median"), cex=2)

        #----------------------------------------
        # plot MAD
        mad_qs <- quantile(raster(s, 2), c(0, 0.02, 0.98, 1))
        # mad_breaks <- c(seq(0, 42), mad_qs[4])
        mad_breaks <- c(seq(0, mad_qs[3], len=254), mad_qs[4])
        par(mar=rep(1, 4), oma=c(0, 0, 2, 0), bg=WATERCOLOR, col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white", cex.main=1.5)
        plot(continent_shp, border=NA, col=LANDCOLOR, box=F)
        plot(raster(s, 2), breaks=mad_breaks, col=viridis(length(mad_breaks) - 1), maxpixels=ncell(s), legend=F, xaxt="n", yaxt="n", bty="n", box=FALSE, add=T)
        legend_at <- signif(seq(0, mad_breaks[length(mad_breaks) - 1], len=7), digits=2)
        legend_labels <- c(legend_at[1], as.character(legend_at[2:(length(legend_at) - 1)]), paste(">", legend_at[length(legend_at)], sep=""))
        plot(raster(matrix(c(legend_at[1], legend_at[length(legend_at)]))), legend.only=T, col=viridis(length(med_breaks) - 1), legend.width=LEGENDWIDTH, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=LEGENDAXISCEX), legend.args=list(text=legend_text, side=3, font=2, line=0.5, cex=LEGENDMAINCEX))
        title(paste("MCD12Q2", metric_name, "MAD"), cex=2)

        #----------------------------------------
        # plot TS slope
        slope_r <- raster(s, 5)
        pvalue_r <- raster(s, 6)
        slope_r[pvalue_r >= 0.05] <- NA
        slope_qs <- quantile(slope_r, c(0, 0.02, 0.98, 1))
        neg_breaks <- c(-max(abs(slope_qs[c(1, 4)])), seq(-max(abs(slope_qs[2:3])), 0, len=128))
        pos_breaks <- c(seq(0, max(abs(slope_qs[2:3])), len=128), max(abs(slope_qs[c(1, 4)])))
        slope_breaks <- c(neg_breaks[1:(length(neg_breaks) - 1)], pos_breaks[2:length(pos_breaks)])
        slope_pal <- colorRampPalette(rev(brewer.pal(9, "RdBu")))
        par(mar=rep(1, 4), oma=c(0, 0, 2, 0), bg=WATERCOLOR, col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white", cex.main=1.5)
        plot(continent_shp, border=NA, col=LANDCOLOR, box=F)
        plot(slope_r, breaks=slope_breaks, col=slope_pal(length(slope_breaks) - 1), maxpixels=ncell(s), legend=F, xaxt="n", yaxt="n", bty="n", box=FALSE, add=T)
        legend_at <- signif(seq(slope_breaks[2], slope_breaks[length(slope_breaks) - 1], len=7), digits=2)
        legend_labels <- c(paste("<", legend_at[1], sep=""), as.character(legend_at[2:(length(legend_at) - 1)]), paste(">", legend_at[length(legend_at)], sep=""))
        plot(raster(matrix(c(legend_at[1], legend_at[length(legend_at)]))), legend.only=T, col=slope_pal(length(slope_breaks) - 1), legend.width=LEGENDWIDTH, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=LEGENDAXISCEX), legend.args=list(text=paste(legend_text, "/yr", sep=""), side=3, font=2, line=0.5, cex=LEGENDMAINCEX))
        title(paste("MCD12Q2", metric_name, "Sig. TS-slopes"), cex=2)

        dev.off()
    }

}
# dev.off() # close the continent file
