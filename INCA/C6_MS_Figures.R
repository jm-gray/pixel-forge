library(raster)
library(mblm)
library(rgdal)
library(viridis)
library(RColorBrewer)

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#### Functions
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

#--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
# makes a nice color ramp legend for rasters when using contrast enhanced breaks
PlotLegend <- function(breaks, viridis_colpal="magma", round_dig=0, num_labels=7, LEGENDAXISCEX=1, LEGENDMAINCEX=1, LEGENDWIDTH=1, LEGENDLABELCOL="black", ...){
    legend_at <- round(seq(breaks[2], breaks[length(breaks) - 1], len=num_labels), round_dig)
    legend_at_date <- legend_at
    legend_labels <- c(paste("<", legend_at_date[1]), as.character(legend_at_date[2:(length(legend_at_date) - 1)]), paste(">", legend_at_date[length(legend_at_date)]))
    plot(
        # raster(matrix(legend_at[1]:legend_at[length(legend_at)])),
        raster(matrix(seq(legend_at[1],legend_at[length(legend_at)], len=255))),
        legend.only=T,
        col=get(viridis_colpal)(length(breaks)-1),
        legend.width=LEGENDWIDTH,
        axis.args=list(
            at=legend_at,
            labels=legend_labels,
            cex.axis=LEGENDAXISCEX,
            col=LEGENDLABELCOL,
            col.axis=LEGENDLABELCOL
            ),
        legend.args=list(text="", side=3, font=2, line=0.5, cex=LEGENDMAINCEX),
        ...
    )
}

#--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
get_tile_raster <- function(tile, sds, dir=inca_dir) raster(dir(dir, pattern=paste(".*", tile, ".", sds, ".tif$", sep=""), full=T))

#--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
# Builds a VRT mosaic of files in file_list, possibly w/ band specification
BuildVRT <- function(file_list, out_file, band=NULL, vrtnodata=0, separate=F, overwrite=F){  
  sep_flag <- ifelse(separate, "-separate", "")
  overwrite_flag <- ifelse(overwrite, "-overwrite", "")
  if(!is.null(band)){
    sys_cmd <- paste("gdalbuildvrt", sep_flag, paste("-vrtnodata", vrtnodata), paste("-b", band), out_file, overwrite_flag, paste(file_list, collapse=" "))
    
  }else{
    sys_cmd <- paste("gdalbuildvrt", sep_flag, paste("-vrtnodata", vrtnodata), out_file, overwrite_flag, paste(file_list, collapse=" "))
  }
  system(sys_cmd)
  return(out_file)
}

#--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
# Returns the path to the matching INCA tile file
GetTile <- function(tile, metric, data_dir) dir(data_dir, pattern=paste(tile, "\\.", metric, "\\.tif$", sep=""), full=T)

#### Constants/definitions/names
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
inca_dir <- "/rsstu/users/j/jmgray2/SEAL/INCA/INCAglobaloutput"
fig_out_dir <- "/rsstu/users/j/jmgray2/SEAL/INCA/"

all_tiles <- c("h08v04", "h09v04", "h10v04", "h11v04", "h12v04", "h13v04", "h08v05", "h09v05", "h10v05", "h11v05", "h12v05", "h08v06", "h09v06", "h10v06", "h00v08", "h00v09", "h00v10", "h01v08", "h01v09", "h01v10", "h01v11", "h02v06", "h02v08", "h02v09", "h02v10", "h02v11", "h03v06", "h03v07", "h03v09", "h03v10", "h03v11", "h04v09", "h04v10", "h04v11", "h05v10", "h05v11", "h05v13", "h06v03", "h06v11", "h07v03", "h07v05", "h07v06", "h07v07", "h08v03", "h08v07", "h08v08", "h08v09", "h09v02", "h09v03", "h09v07", "h09v08", "h09v09", "h10v02", "h10v03", "h10v07", "h10v08", "h10v09", "h10v10", "h10v11", "h11v02", "h11v03", "h11v06", "h11v07", "h11v08", "h11v09", "h11v10", "h11v11", "h11v12", "h12v01", "h12v02", "h12v03", "h12v07", "h12v08", "h12v09", "h12v10", "h12v11", "h12v12", "h12v13", "h13v01", "h13v02", "h13v03", "h13v08", "h13v09", "h13v10", "h13v11", "h13v12", "h13v13", "h13v14", "h14v01", "h14v02", "h14v03", "h14v04", "h14v09", "h14v10", "h14v11", "h14v14", "h14v16", "h14v17", "h15v01", "h15v02", "h15v03", "h15v05", "h15v07", "h15v11", "h15v14", "h15v15", "h15v16", "h15v17", "h16v00", "h16v01", "h16v02", "h16v05", "h16v06", "h16v07", "h16v08", "h16v09", "h16v12", "h16v14", "h16v16", "h16v17", "h17v00", "h17v01", "h17v02", "h17v03", "h17v04", "h17v05", "h17v06", "h17v07", "h17v08", "h17v10", "h17v12", "h17v13", "h17v15", "h17v16", "h17v17", "h18v00", "h18v01", "h18v02", "h18v03", "h18v04", "h18v05", "h18v06", "h18v07", "h18v08", "h18v09", "h18v14", "h18v15", "h18v16", "h18v17", "h19v00", "h19v01", "h19v02", "h19v03", "h19v04", "h19v05", "h19v06", "h19v07", "h19v08", "h19v09", "h19v10", "h19v11", "h19v12", "h19v15", "h19v16", "h19v17", "h20v01", "h20v02", "h20v03", "h20v04", "h20v05", "h20v06", "h20v07", "h20v08", "h20v09", "h20v10", "h20v11", "h20v12", "h20v13", "h20v15", "h20v16", "h20v17", "h21v01", "h21v02", "h21v03", "h21v04", "h21v05", "h21v06", "h21v07", "h21v08", "h21v09", "h21v10", "h21v11", "h21v13", "h21v15", "h21v16", "h21v17", "h22v01", "h22v02", "h22v03", "h22v04", "h22v05", "h22v06", "h22v07", "h22v08", "h22v09", "h22v10", "h22v11", "h22v13", "h22v14", "h22v15", "h22v16", "h23v01", "h23v02", "h23v03", "h23v04", "h23v05", "h23v06", "h23v07", "h23v08", "h23v09", "h23v10", "h23v11", "h23v15", "h23v16", "h24v02", "h24v03", "h24v04", "h24v05", "h24v06", "h24v07", "h24v12", "h24v15", "h25v02", "h25v03", "h25v04", "h25v05", "h25v06", "h25v07", "h25v08", "h25v09", "h26v02", "h26v03", "h26v04", "h26v05", "h26v06", "h26v07", "h26v08", "h27v03", "h27v04", "h27v05", "h27v06", "h27v07", "h27v08", "h27v09", "h27v10", "h27v11", "h27v12", "h27v14", "h28v03", "h28v04", "h28v05", "h28v06", "h28v07", "h28v08", "h28v09", "h28v10", "h28v11", "h28v12", "h28v13", "h28v14", "h29v03", "h29v05", "h29v06", "h29v07", "h29v08", "h29v09", "h29v10", "h29v11", "h29v12", "h29v13", "h30v05", "h30v06", "h30v07", "h30v08", "h30v09", "h30v10", "h30v11", "h30v12", "h30v13", "h31v06", "h31v07", "h31v08", "h31v09", "h31v10", "h31v11", "h31v12", "h31v13", "h32v07", "h32v08", "h32v09", "h32v10", "h32v11", "h32v12", "h33v07", "h33v08", "h33v09", "h33v10", "h33v11", "h34v07", "h34v08", "h34v09", "h34v10", "h35v08", "h35v09", "h35v10")
ganguly_tiles <- c("h23v02", "h12v03", "h12v04", "h20v07", "h20v08", "h27v05")

doy_sds <- c("Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy")
other_sds <- c("EVI_Minimum", "EVI_Amplitude", "EVI_Area", "GSL")

inca_names <- c("median", "MAD", "mean", "sd", "ts_slope", "ts_pvalue", "num_vals")

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#### Ganguly plots
# Reproduce Ganguly 2010 plot
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

# First the DOY metrics
do_sds_stretch <- FALSE # T=each SDS gets 2/98 linear stretch, F=all SDS plotted 1:365
default_breaks <- c(-1e6, seq(1, 365, len=254), 1e6) # 255 colors
ncell_to_plot <- (2400^2)
output_file <- file.path(fig_out_dir, "MCD12Q2_C6_Ganguly_DOY_Metrics_same_ramp.pdf")
# output_file <- file.path(fig_out_dir, "MCD12Q2_C6_Ganguly_DOY_Metrics.pdf")
viridis_color_pal <- "magma"
NAcol <- rgb(0.3, 0.3, 0.3)

# some calculations for laying out the plot
raster_size_in <- 2
plot_mai <- rep(0, 4)
plot_omi <- c(0.5, 0.5, 0.5, 0.5)
pdf_height <- (length(doy_sds) * raster_size_in) + plot_omi[1] + plot_omi[3] + (sum(plot_mai[c(1, 3)]) * length(doy_sds))
pdf_width <- (length(ganguly_tiles) * raster_size_in) + plot_omi[2] + plot_omi[4] + (sum(plot_mai[c(2, 4)]) * length(ganguly_tiles))

# these indices control which plots get labels/legends
y_lab_plots <- seq(1, ((length(ganguly_tiles) + 1) * length(doy_sds)), by=(length(ganguly_tiles) + 1))
x_lab_plots <- 1:length(ganguly_tiles)

pdf(file=output_file, width=pdf_width, height=pdf_height)
# nf <- layout(matrix(1:(length(ganguly_tiles) * length(doy_sds)), nrow=length(doy_sds), byrow=T))
nf <- layout(matrix(1:((length(ganguly_tiles) + 1) * length(doy_sds)), nrow=length(doy_sds), byrow=T))
par(mar=rep(0.1, 4), oma=c(1, 2, 2, 1))
label_counter <- 1
breaks_list <- list()
for(sds in doy_sds){
    if(do_sds_stretch){
        # first, we generate 2/98% breaks for this particular SDS
        all_tiles_raster <- lapply(ganguly_tiles, get_tile_raster, sds=sds)
        tmp_vals <- lapply(all_tiles_raster, values)
        all_vals <- do.call(rbind, tmp_vals)
        all_vals_qs <- quantile(all_vals, c(0, 0.02, 0.98, 1), na.rm=T)
        this_sds_breaks <- c(all_vals_qs[1], seq(all_vals_qs[2], all_vals_qs[3], len=254), all_vals_qs[4])
        breaks_list[[which(doy_sds == sds)]] <- this_sds_breaks
        # print(paste("Got the breaks for", sds, ":", all_vals_qs))
        rm(all_vals, tmp_vals, all_tiles_raster)
    }else{
        # use defaults
        this_sds_breaks <- default_breaks
    }

    # loop through the tiles and make the plots
    for(tile in ganguly_tiles){
        # get the INCA output tile file, and grab the median layer
        in_file <- dir(inca_dir, pattern=paste(".*", tile, ".", sds, ".tif$", sep=""), full=T)
        inca_s <- stack(in_file)
        names(inca_s) <- inca_names
        inca_median <- inca_s$median
        # handle NA values
        NAvalue(inca_median) <- 32767
        inca_median[inca_median == -32768] <- NA
        plot(inca_median, breaks=this_sds_breaks, col=get(viridis_color_pal)(length(default_breaks) - 1), maxpixels=ncell_to_plot, legend= F, xaxt="n", yaxt="n", colNA=NAcol, axes=F, box=F)

        # add top tile label
        if(label_counter %in% x_lab_plots) mtext(tile, side=3, line=0, outer=F)
        if(label_counter %in% y_lab_plots) mtext(sds, side=2, line=0.5, outer=F)
        label_counter <- label_counter + 1
    }

    # plot the legend
    plot(1:10, type="n", xaxt="n", yaxt="n", axes=F, box=F)
    PlotLegend(breaks=this_sds_breaks, viridis_colpal = viridis_color_pal, smallplot=c(0, 0.1, 0.2, 0.8))
    label_counter <- label_counter + 1
}
dev.off()


#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# Now the derivative and other metrics
do_sds_stretch <- TRUE # T=each SDS gets 2/98 linear stretch, F=all SDS plotted 1:365
default_breaks <- c(-1e6, seq(1, 365, len=254), 1e6) # 255 colors
ncell_to_plot <- (2400^2)
output_file <- file.path(fig_out_dir, "MCD12Q2_C6_Ganguly_Other_Metrics.pdf")
viridis_color_pal <- "viridis"
NAcol <- rgb(0.3, 0.3, 0.3)

# some calculations for laying out the plot
raster_size_in <- 2
plot_mai <- rep(0, 4)
plot_omi <- c(0.5, 0.5, 0.5, 0.5)
pdf_height <- (length(other_sds) * raster_size_in) + plot_omi[1] + plot_omi[3] + (sum(plot_mai[c(1, 3)]) * length(other_sds))
pdf_width <- (length(ganguly_tiles) * raster_size_in) + plot_omi[2] + plot_omi[4] + (sum(plot_mai[c(2, 4)]) * length(ganguly_tiles))

# these indices control which plots get labels/legends
y_lab_plots <- seq(1, ((length(ganguly_tiles) + 1) * length(other_sds)), by=(length(ganguly_tiles) + 1))
x_lab_plots <- 1:length(ganguly_tiles)

pdf(file=output_file, width=pdf_width, height=pdf_height)
# nf <- layout(matrix(1:(length(ganguly_tiles) * length(doy_sds)), nrow=length(doy_sds), byrow=T))
nf <- layout(matrix(1:((length(ganguly_tiles) + 1) * length(other_sds)), nrow=length(other_sds), byrow=T))
par(mar=rep(0.1, 4), oma=c(1, 2, 2, 1))
label_counter <- 1
breaks_list <- list()
for(sds in other_sds){
    if(do_sds_stretch){
        # first, we generate 2/98% breaks for this particular SDS
        all_tiles_raster <- lapply(ganguly_tiles, get_tile_raster, sds=sds)
        tmp_vals <- lapply(all_tiles_raster, values)
        all_vals <- do.call(rbind, tmp_vals)
        all_vals_qs <- quantile(all_vals, c(0, 0.02, 0.98, 1), na.rm=T)
        this_sds_breaks <- c(all_vals_qs[1], seq(all_vals_qs[2], all_vals_qs[3], len=254), all_vals_qs[4])
        breaks_list[[which(other_sds == sds)]] <- this_sds_breaks
        # print(paste("Got the breaks for", sds, ":", all_vals_qs))
        rm(all_vals, tmp_vals, all_tiles_raster)
    }else{
        # use defaults
        this_sds_breaks <- default_breaks
    }

    # loop through the tiles and make the plots
    for(tile in ganguly_tiles){
        # get the INCA output tile file, and grab the median layer
        in_file <- dir(inca_dir, pattern=paste(".*", tile, ".", sds, ".tif$", sep=""), full=T)
        inca_s <- stack(in_file)
        names(inca_s) <- inca_names
        inca_median <- inca_s$median
        # handle NA values
        NAvalue(inca_median) <- 32767
        inca_median[inca_median == -32768] <- NA
        plot(inca_median, breaks=this_sds_breaks, col=get(viridis_color_pal)(length(default_breaks) - 1), maxpixels=ncell_to_plot, legend= F, xaxt="n", yaxt="n", colNA=NAcol, axes=F, box=F)

        # add top tile label
        if(label_counter %in% x_lab_plots) mtext(tile, side=3, line=0, outer=F)
        if(label_counter %in% y_lab_plots) mtext(sds, side=2, line=0.5, outer=F)
        label_counter <- label_counter + 1
    }

    # plot the legend
    plot(1:10, type="n", xaxt="n", yaxt="n", axes=F, box=F)
    if(sds == "EVI_Minimum") round_dig <- 2
    if(sds == "EVI_Amplitude") round_dig <- 2
    if(sds == "EVI_Area") round_dig <- 1
    if(sds == "GSL") round_dig <- 0
    PlotLegend(breaks=this_sds_breaks, viridis_colpal = viridis_color_pal, smallplot=c(0, 0.1, 0.2, 0.8), round_dig = round_dig)
    label_counter <- label_counter + 1
}
dev.off()

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#### MODIS tile grid plot
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

# reproject the oceans shapefile
# ogr2ogr -t_srs "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m" ~/Desktop/MCD12Q2_C6_Figures/oceans_sin ~/Desktop/ne_10m_ocean/ne_10m_ocean.shp
# # reproject the continents shapefile
# ogr2ogr -t_srs "+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m" ~/Desktop/MCD12Q2_C6_Figures/continents_sin ~/Desktop/continents/continent.shp

oceans_shp <- shapefile("~/Desktop/MCD12Q2_C6_Figures/oceans_sin/ne_10m_ocean.shp")
continents_shp <- shapefile("~/Desktop/MCD12Q2_C6_Figures/continents_sin/continent.shp")
modgrid_shp <- shapefile("~/Desktop/MCD12Q2_C6_Figures/modis_grid/modis_sinusoidal_grid_world.shp")
modgrid_shp$tile <- paste(paste("h", formatC(as.integer(modgrid_shp$h), flag="0", width=2), sep=""), paste("v", formatC(as.integer(modgrid_shp$v), flag="0", width=2), sep=""), sep="")

ganguly_tiles <- c("h23v02", "h12v03", "h12v04", "h20v07", "h20v08", "h27v05")
ganguly_tiles_shp <- modgrid_shp[modgrid_shp$tile %in% ganguly_tiles,]

# for labeling
x_centers <- coordinates(modgrid_shp[modgrid_shp$v == "0",])[,1]
y_centers <- coordinates(modgrid_shp[modgrid_shp$h == "0",])[,2]
h_labels <- sort(unique(as.integer(modgrid_shp$h)))
v_labels <- sort(unique(as.integer(modgrid_shp$v)))

land_color <- rgb(0.3, 0.3, 0.3)
ocean_color <- rgb(0.7, 0.7, 0.7)
label_x_offset <- diff(x_centers)[1]
label_y_offset <- diff(y_centers)[1]
pdf(file="~/Desktop/MCD12Q2_C6_Figures/modgrid.pdf", height=10, width=18)
par(mar=rep(0,4), oma=c(0, 2, 2, 0))
plot(oceans_shp, col=ocean_color, border=NA)
plot(continents_shp, col=land_color, border=NA, add=T)
plot(modgrid_shp, col=NA, border="black", lwd=1, add=T)
plot(ganguly_tiles_shp, col=NA, border="red", lwd=2.5, add=T, lty=1)
text(x=x_centers, y=rep(y_centers[length(y_centers)], length(x_centers)) + label_y_offset, labels=h_labels)
text(x=rep(x_centers[1], length(y_centers)) - label_x_offset, y=y_centers, labels=rev(v_labels))
dev.off()

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#### Latitudinal band plot
# Make latitudinal band plots
# NOTE: INCAsample.R does the data grab/munge to produce all_data 
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

load("~/Desktop/INCA_sample_ten_percent_all_data.Rdata")
# above file is on HPC: /rsstu/users/j/jmgray2/SEAL/INCA/INCAsample_ten_percent
lat_min <- -90
lat_max <- 90
lat_breaks <- seq(lat_min, lat_max, by=1)
center_lats <- lat_breaks[1:(length(lat_breaks) - 1)] + (diff(lat_breaks) / 2)

GetCenterLat <- function(x, lat_breaks, center_lats) center_lats[findInterval(x, lat_breaks, all.inside=T)]

all_data[, center_lat_1deg:=GetCenterLat(lat, lat_breaks = lat_breaks, center_lats = center_lats)]
all_data_med_1deg <- all_data[, lapply(.SD, median, na.rm=T), by=center_lat_1deg, .SDcols=8:17]
all_data_mad_1deg <- all_data[, lapply(.SD, mad, na.rm=T), by=center_lat_1deg, .SDcols=8:17]

MakeTransp <- function(thecolor, alpha=0.5) rgb(matrix(c(col2rgb(thecolor)), ncol=3), max=255, alpha=alpha*255)
PlotLatPheno <- function(med_data, mad_data, phenometric, plot_col=rgb(0.5, 0.5, 0.5), poly_transp=0.75, label_y_axt=FALSE, ylim=NULL, ...){
    med_vals <- med_data[, .(get(phenometric), center_lat_1deg)][order(center_lat_1deg)][, V1]
    mad_vals <- mad_data[, .(get(phenometric), center_lat_1deg)][order(center_lat_1deg)][, V1]
    plot_lats <- sort(med_data[, center_lat_1deg])
    poly_x <- c(med_vals - mad_vals, rev(med_vals + mad_vals))
    poly_y <- c(plot_lats, rev(plot_lats))
    if(is.null(ylim)) ylim <- range(plot_lats)
    plot(med_vals, plot_lats, type="n", xlim=range(poly_x), xlab="", ylab="", xaxt="n", yaxt="n", ylim=ylim)
    polygon(poly_x, poly_y, col=MakeTransp(plot_col, poly_transp), border=NA)
    points(med_vals, plot_lats, xlim=range(poly_x), type="l", col=plot_col, ...)
    title(phenometric)
    axis(side=1)
    if(label_y_axt){
        axis(side=2)
    }else{
        axis(side=2, labels=F)
    }
    axis(side=4, labels=F)
}

library(viridis)
# plot_cols <- c(magma(9)[1:7], viridis(4)[1:3])
plot_cols <- rep(rgb(0.4, 0.4, 0.4), 10)
metrics_to_plot <- names(all_data)[8:17]
ylim <- c(-56, 80)
# pdf(file="~/Desktop/MCD12Q2_C6_Figures/latitudinal_medians.pdf", height=8, width=12)
pdf(file="~/Desktop/MCD12Q2_C6_Figures/latitudinal_medians_bw.pdf", height=8, width=12)
layout(matrix(1:10, nrow=2, byrow=2))
par(mar=c(1, 0.25, 1.5, 0.25), oma=c(1, 3, 1, 1), tcl=0.35, mgp=c(3, 0, 0))
i <- 1
for(this_metric in metrics_to_plot){
    if(i %in% c(1, 6)) label_y_axt <- TRUE
    PlotLatPheno(all_data_med_1deg, all_data_mad_1deg, this_metric, lwd=2, plot_col=plot_cols[i], 
    label_y_axt=label_y_axt, ylim=ylim)
    if(i %in% c(1, 6)) mtext("Latitude", side=2, line=1)
    label_y_axt <- FALSE
    i <- i + 1
}
dev.off()

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#### Global phenometric plots
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
data_dir <- "/rsstu/users/j/jmgray2/SEAL/INCA/INCAglobaloutput"
out_dir <- "/rsstu/users/j/jmgray2/SEAL/INCA/INCAglobaloutput/continent_mosaics"

# output projection, bounds, resolution, naming, etc.
out_prefix <- "Global"
out_res <- 2500
out_suffix <- "Sin_2500m"
vrt_no_data_out <- "32767"
out_proj4 <- "'+proj=sinu +R=6371007.181 +nadgrids=@null +wktext +unit=m'"
doy_sds <- c("Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy")
other_sds <- c("EVI_Minimum", "EVI_Amplitude", "EVI_Area", "GSL")
all_sds <- c(doy_sds, other_sds)

# create global mosaics, then merge them into a single multiband raster, then
# reproject, resample, etc. with gdalwarp
tiles_to_mosaic <- gsub(".*(h[0-9]{2}v[0-9]{2}).*", "\\1", dir(data_dir, pattern=".*Peak.tif"))
for(this_metric in all_sds){
    print(paste("Doing:", this_metric))
    out_vrt_file <- file.path(out_dir, paste(paste(out_prefix, this_metric, sep="_"), ".vrt", sep=""))
    mosaic_files <- unlist(lapply(tiles_to_mosaic, GetTile, metric=this_metric, data_dir=data_dir))
    out_vrt_file <- BuildVRT(mosaic_files, out_file=out_vrt_file, vrtnodata=32767, overwrite=T)
    out_warp_file <- file.path(out_dir, paste(paste(out_prefix, this_metric, out_suffix, sep="_"), ".tif", sep=""))
    gdalwarp_cmd <- paste("gdalwarp -overwrite -r med -srcnodata 32767 -dstnodata 32767 -t_srs", out_proj4, "-tr", out_res, out_res, out_vrt_file, out_warp_file)
    system(gdalwarp_cmd)
}

#--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
# Let's do the actual plotting...
oceans_shp <- shapefile("~/Desktop/MCD12Q2_C6_Figures/oceans_sin/ne_10m_ocean.shp")
continents_shp <- shapefile("~/Desktop/MCD12Q2_C6_Figures/continents_sin/continent.shp")
modgrid_shp <- shapefile("~/Desktop/MCD12Q2_C6_Figures/modis_grid/modis_sinusoidal_grid_world.shp")
modgrid_shp$tile <- paste(paste("h", formatC(as.integer(modgrid_shp$h), flag="0", width=2), sep=""), paste("v", formatC(as.integer(modgrid_shp$v), flag="0", width=2), sep=""), sep="")
s <- stack("~/Desktop/global_mosaics/Global_GSL_Sin_2500m.tif")
# return(c(round(med), MAD, round(avg), stddev, ts.slope, ts.pvalue, num_vals))

# GSL median
qs <- quantile(raster(s, 1), c(0, 0.02, 0.98, 1))
breaks <- c(qs[1], seq(qs[2], qs[3], len=254), qs[4])
land_color <- rgb(0.3, 0.3, 0.3)
ocean_color <- rgb(0.7, 0.7, 0.7)
pdf(file="~/Desktop/MCD12Q2_C6_Figures/global_GSL_median.pdf", height=10, width=20)
par(mar=c(5, 0, 0, 0), oma=rep(0, 4))
plot(oceans_shp, col=ocean_color, border=NA)
plot(continents_shp, col=land_color, border=NA, add=T)
plot(raster(s, 1), colNA=NA, legend=F, add=T, breaks=breaks, col=magma(length(breaks) - 1), maxpixels=ncell(raster(s, 1)))
PlotLegend(breaks=breaks, viridis_colpal = "magma", round_dig=0, horiz=T)
dev.off()

# GSL MAD
qs <- quantile(raster(s, 2), c(0, 0.02, 0.98, 1))
breaks <- c(qs[1], seq(qs[2], qs[3], len=254), qs[4])
land_color <- rgb(0.3, 0.3, 0.3)
ocean_color <- rgb(0.7, 0.7, 0.7)
pdf(file="~/Desktop/MCD12Q2_C6_Figures/global_GSL_MAD.pdf", height=10, width=20)
par(mar=c(5, 0, 0, 0), oma=rep(0, 4))
plot(oceans_shp, col=ocean_color, border=NA)
plot(continents_shp, col=land_color, border=NA, add=T)
plot(raster(s, 2), colNA=NA, legend=F, add=T, breaks=breaks, col=viridis(length(breaks) - 1), maxpixels=ncell(raster(s, 2)))
PlotLegend(breaks=breaks, viridis_colpal = "viridis", round_dig=0, horiz=T)
dev.off()

# GSL trend
gsl_ts <- raster(s, 5)
gsl_p <- raster(s, 6)
gsl_ts[gsl_p > 0.05] <- NA

qs <- quantile(gsl_ts, c(0, 0.02, 0.98, 1))
ncols <- 255
neg_breaks <- seq(-max(abs(qs[2:3])), 0, len=round(ncols / 2))
neg_breaks <- c(qs[1], neg_breaks[-length(neg_breaks)])
pos_breaks <- seq(0, max(abs(qs[2:3])), len=round(ncols / 2))
pos_breaks <- c(pos_breaks[-1], qs[4])
breaks <- c(neg_breaks, pos_breaks)

# pal <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))
# pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
pal <- colorRampPalette(rev(brewer.pal(11, "PiYG")))

land_color <- rgb(0.3, 0.3, 0.3)
ocean_color <- rgb(0.7, 0.7, 0.7)
pdf(file="~/Desktop/MCD12Q2_C6_Figures/global_GSL_TSslope_v3.pdf", height=10, width=20)
par(mar=c(5, 0, 0, 0), oma=rep(0, 4))
plot(oceans_shp, col=ocean_color, border=NA)
plot(continents_shp, col=land_color, border=NA, add=T)
plot(gsl_ts, colNA=NA, legend=F, add=T, breaks=breaks, col=pal(length(breaks) - 1), maxpixels=ncell(gsl_ts))
PlotLegend(breaks=breaks, viridis_colpal = "pal", round_dig=0, horiz=T)
dev.off()


#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#### MCD12Q2 - NPN analysis
#### Derived from AGU_2018_figures.R
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
library(data.table)
library(mblm)
library(RColorBrewer)
load("~/Google Drive File Stream/My Drive/Projects/INCA/AGU_analysis/npnDT_all.Rdata")

#------------------------------------------------------------
nrow(npnDT_all) # 388,470 phenophase observations
length(npnDT_all[,unique(site_id)]) # 6056 sites

#------------------------------------------------------------
# how many site years of data for each phenophase?
tmp <- npnDT_all[,.(mean(first_yes_doy, na.rm=T)), by=.(site_id, ObsYear, phenophase_name)]
nrow(tmp) # 97,599 site-year-phenophase observations

plot_col <- brewer.pal(5, "Greys")[3]
png(file="~/Desktop/MCD12Q2_C6_Figures/USNPN_N_by_phenophase.png", width=2400, height=1250)
# par(mar=c(42, 9, 0, 0), oma=rep(1, 4), col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white", bg=NA)
par(mar=c(42, 9, 0, 0), oma=rep(1, 4))
barplot(tmp[,.N,by=phenophase_name][order(-N)][N>500,N], names=tmp[,.N,by=phenophase_name][order(-N)][N>500,phenophase_name], las=2, ylim=c(0, 9500), col=plot_col, border=NA, cex.names=4.5, cex.axis=4.5, lwd.ticks=5)
# caption note: only for categories with at least 500 site-years of a phenometric
# grid(nx=NA, ny=NULL)
box(lwd=5)
dev.off()

#------------------------------------------------------------
# how many site years by phenophase/SFT?
names(npnDT_all)
tmp <- npnDT_all[,.(mean(first_yes_doy, na.rm=T)), by=.(site_id, ObsYear, phenophase_name, species_functional_type)]
trash <- tmp[,.N,by=.(phenophase_name, species_functional_type)][N>100]
# for checking the matrix
# trash[phenophase_name == "Full flowering" & species_functional_type == "Deciduous broadleaf"]
# trash[phenophase_name == "All leaves fallen" & species_functional_type == "Deciduous broadleaf"]
# trash[phenophase_name == "Leaves" & species_functional_type == "Forb"]
m <- as.matrix(dcast(trash, species_functional_type~phenophase_name, value.var="N"))
pheno_names <- names(m[1,2:dim(m)[2]])
sft_names <- m[,1]
set.seed(1)
bar_cols <- sample(brewer.pal(length(sft_names), "Paired"))
# bar_cols <- sample(magma(length(sft_names)))
bp_m <- matrix(as.numeric(m[,2:dim(m)[2]]), ncol=length(2:dim(m)[2]), byrow=F)
# bp_m[which(sft_names == "Deciduous broadleaf"), which(pheno_names == "Full flowering")]
# bp_m[which(sft_names == "Deciduous broadleaf"), which(pheno_names == "All leaves fallen")]
# bp_m[which(sft_names == "Forb"), which(pheno_names == "Leaves")]
bp_m[is.na(bp_m)] <- 0
the_mask <- colSums(bp_m) > 500
bp_m <- bp_m[,the_mask]
pheno_names <- pheno_names[the_mask]

png(file="~/Desktop/MCD12Q2_C6_Figures/USNPN_N_by_phenophase_sft.png", width=2400, height=1250)
# par(mar=c(42, 9, 0, 0), oma=rep(1, 4), col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white", bg=NA)
par(mar=c(42, 9, 0, 0), oma=rep(1, 4))
barplot(bp_m[,order(-colSums(bp_m))], beside=F, las=2, names=pheno_names[order(-colSums(bp_m))], col=bar_cols, border=NA, cex.names=4.5, cex.axis=3.5, lwd.ticks=5, ylim=c(0, 13000))
legend("topright", legend=sft_names, fill=bar_cols, cex=4, border=NA, ncol=2, bty="n")
box(lwd=5)
dev.off()

#------------------------------------------------------------
# a faster single median Theil-Sen method than what is available in MBLM
# Modified from here: https://r.789695.n4.nabble.com/help-speeding-up-simple-Theil-regression-function-td4646923.html
np.lm.alt <-function(X, Y, ...){
        # Ch 9.2: Slope est. (X) for Theil statistic
        dat <- data.frame("X"=X, "Y"=Y)
        combos <- combn(nrow(dat), 2)
        i.s <- combos[1,]
        j.s <- combos[2,]
       
        Y.num <- dat[j.s,"Y"] - dat[i.s,"Y"]
        X.dom <- dat[j.s,"X"] - dat[i.s,"X"]
        slopes <- Y.num / X.dom
        med_slope <- median(slopes, na.rm=T)
        
        # Ch 9.4: Intercept est. for Theil statistic
        intercepts <- dat[,"Y"] - med_slope*dat[,"X"]
        med_intercept <- median(intercepts, na.rm=T)
        # out <- data.frame(Intercept, X)
        # return(out)
        list("slopes"=slopes, "intercepts"=intercepts, "coefficients"=c(med_intercept, med_slope))
}

#------------------------------------------------------------
# Calculate correlations and regressions between USNPN and MCD12Q2
DoNPNAnalysis <- function(DT, filter_exps, mcd_var, mid_duration=FALSE, min_duration=0, out_frac=0, first_yes=T, site_sd_thresh=30, makeplot=TRUE, out_file=NULL, line_cols=NULL, one_to_one_col="white", color_points=FALSE){
    # apply the filtering expression(s)
    old_N <- nrow(DT)
    n_filter <- rep(NA, length(filter_exps)) # track the reduction in rows for each expression in filter_exps
    i <- 1
    DT_filt <- DT # initialize
    for(filter_exp in filter_exps){
        DT_filt <- DT_filt[eval(filter_exp),]
        n_filter[i] <- old_N - nrow(DT_filt)
        paste(n_filter[i], "rows filtered to", nrow(DT_filt), "for filter", i)

        old_N <- nrow(DT_filt)
        i <- i + 1
    }
    
    # Calculate USNPN phenophase values for all individual_id + ObsYear;
    # measure the duration of the recorded phenophase, then:
    # 1) screen under the assumption that phenophase records with unrealistically short durations are not reliable
    # 2) create phenophase duration midpoint data (e.g. maybe Peak = midpoint of first and last "Yes, Leaves")
    DT_filt[, npn_phenophase_duration := last_yes_doy - first_yes_doy]
    DT_filt <- DT_filt[npn_phenophase_duration >= min_duration] # filter by min_duration
    dim(DT_filt)
    
    # aggregate to site
    if(first_yes){
        DT_filt <- DT_filt[, .(lat=latitude[1], long=longitude[1], sft=species_functional_type[1], mcd12=mcd12q1_2006[1], nlcd_prop=nlcd_mcd_agree[1], shd=site_shannon[1], npn_site=as.numeric(median(first_yes_doy, na.rm=T)), npn_site_lag=as.numeric(median(last_no_doy, na.rm=T)), npn_site_mid=median(first_yes_doy + (npn_phenophase_duration / 2), na.rm=T), var_npn=var(first_yes_doy, na.rm=T), mcd_site=as.numeric(median(get(mcd_var), na.rm=T)), lag_var=as.numeric(median(first_yes_doy - last_no_doy, na.rm=T))), by=.(site_id, ObsYear)]
    }else{
        DT_filt <- DT_filt[, .(lat=latitude[1], long=longitude[1], sft=species_functional_type[1], mcd12=mcd12q1_2006[1], nlcd_prop=nlcd_mcd_agree[1], shd=site_shannon[1], npn_site=as.numeric(median(last_yes_doy, na.rm=T)), npn_site_lag=as.numeric(median(first_no_doy, na.rm=T)), npn_site_mid=median(first_yes_doy + (npn_phenophase_duration / 2), na.rm=T), var_npn=var(last_yes_doy, na.rm=T), mcd_site=as.numeric(median(get(mcd_var), na.rm=T)), lag_var=as.numeric(median(last_yes_doy - first_no_doy, na.rm=T))), by=.(site_id, ObsYear)]
    }

    # eliminate high variance site-years
    # DT_filt <- DT_filt[sqrt(var_npn) <= site_sd_thresh | is.na(var_npn)]
    
    # extract regression variables    
    # lag <- DT_filt[, med_lag]
    # npn_obs <- DT_filt[, npn_site]
    if(mid_duration){
        npn_mid <- DT_filt[, npn_site_mid]
        npn_var <- DT_filt[, var_npn]
        mcd_obs <- DT_filt[, mcd_site]
        lc_weights <- DT_filt[, nlcd_prop]
    }else{
        npn_mid <- DT_filt[, (npn_site + npn_site_lag) / 2]
        npn_var <- DT_filt[, var_npn]
        mcd_obs <- DT_filt[, mcd_site]
        lc_weights <- DT_filt[, nlcd_prop]
    }
    
    # do a linear regression
    lm1 <- lm(npn_mid ~ mcd_obs, weights=lc_weights)
    lm1_confint <- confint(lm1)
    lm1_rsq <- summary(lm1)$r.squared

    # do a Theil-Sen regression
    tmp1 <- npn_mid[!is.na(npn_mid) & !is.na(mcd_obs)]
    tmp2 <- mcd_obs[!is.na(npn_mid) & !is.na(mcd_obs)]
    ts1 <- mblm(tmp1~tmp2)

    # get the approximate 95% CI
    # confint.mblm() produces limits that do not contain the estimate! So we make our own based on:
    # https://www.ucl.ac.uk/child-health/short-courses-events/about-statistical-courses/research-methods-and-statistics/chapter-8-content-8#:~:text=It%20is%20possible%20to%20similarly,confidence%20interval%20for%20the%20median.
    lower_ts_rank <- round((length(ts1$slopes) / 2) - (1.96 * sqrt(length(ts1$slopes)) / 2))
    upper_ts_rank <- round(1 + (length(ts1$slopes) / 2) + (1.96 * sqrt(length(ts1$slopes)) / 2))
    ts1_slope_confint <- sort(ts1$slopes)[c(lower_ts_rank, upper_ts_rank)]
    ts1_int_confint <- sort(ts1$intercepts)[c(lower_ts_rank, upper_ts_rank)]
    ts1_confint <- matrix(c(ts1_int_confint, ts1_slope_confint), nrow=2, byrow=T)
    # ts1_confint <- confint(ts1)
    # ts1_confint <- confint.mblm(ts1)
    ts1_rsq <- cor(predict(ts1), npn_mid[!is.na(npn_mid) & !is.na(mcd_obs)], use="complete.obs")^2

    # retrieve the USNPN phenophase_name from the filtering expression
    # this is necessary for printing diagnostics and plotting
    # NOTE: assumes that the phenophase_name spec is in the first filtering exp
    filter_l <- unlist(strsplit(as.character(filter_exps[1]), split="&"))
    npn_name <- gsub(".*\"(.*)\".*", "\\1", filter_l[grep("phenophase_name", filter_l)])
    if(makeplot){
        xlab <- paste("MCD12Q2", mcd_var)
        
        if(first_yes){
            ylab <- paste(npn_name, "First Yes", sep=", ")
        }else{
            ylab <- paste(npn_name, "Last Yes", sep=", ")
        }
        
        if(mid_duration){
            ylab <- paste("Midpoint of", npn_name, sep=" ")
        }

        xlim <- ylim <- range(c(mcd_obs, npn_mid), na.rm=T)
        # par(mar=c(2.5, 2.5, 1.5, 1.5), tcl=0.35)
        makeTransp <- function(thecol, alpha=0.7) rgb(t(matrix(col2rgb(thecol))), max=255, alpha=alpha*255)
        igbp_colors <- c("#51B6F5", "#218A21", "#31CD31", "#9ACD31", "#97FA97", "#8FBB8F", "#BB8F8F", "#F5DEB3", "#DBEB9D", "#FFD600", "#EFB766", "#4682B2", "#FAED73", "#FF0000", "#999355", "#F5F5DC", "#BDBDBD", "#000000")
        igbp_names <- c("water", "enf", "ebf", "dnf", "dbf", "mixed", "closed shrubs", "open shrubs", "woody savannas", "savannas", "grasslands", "perm wetlands", "croplands", "urban", "crop/natural mosaic", "snow and ice", "barren/sparse veg", "unclassified")
        igbp_colors <- sapply(igbp_colors, makeTransp, alpha=1)
        
        # ols_col <- ts_col <- rgb(0.4, 0.4, 0.4)
        ols_col <- magma(10)[5]
        ts_col <- magma(10)[7]
        
        if(color_points){
            plot(mcd_obs, npn_mid, pch=16, col=igbp_colors[DT_filt[, mcd12] + 1], xlim=xlim, ylim=ylim, cex=1.5, xlab="", ylab="", xaxt="n", yaxt="n")
            legend_names <- c("1:1 line", "OLS", "Theil-Sen", igbp_names[DT_filt[, sort(unique(mcd12))] + 1])
            legend_cols <- c(1, ols_col, ts_col, igbp_colors[DT_filt[, sort(unique(mcd12))] + 1])
            legend_lty <- legend_lwd <- rep(NA, length(legend_names))
            legend_lty[1:3] <- 3:1
            legend_lwd[1:3] <- 3
            legend_pch <- rep(16, length(legend_cols))
            legend_pch[1:3] <- NA
            legend("bottomright", legend=legend_names, lty=legend_lty, col=legend_cols, lwd=legend_lwd, pch=legend_pch, bty="n", bg="white")
        }else{
            plot(mcd_obs, npn_mid, pch=16, col=rgb(0.3, 0.3, 0.3, 0.5), xlim=xlim, ylim=ylim, cex=1.5, xlab="", ylab="", xaxt="n", yaxt="n")
            legend_names <- c("1:1 line", "OLS", "Theil-Sen")
            legend_cols <- c(1, ols_col, ts_col)
            legend_lty <- legend_lwd <- rep(NA, length(legend_names))
            legend_lty <- 3:1
            legend_lwd <- rep(3, 3)
            legend("bottomright", legend=legend_names, lty=legend_lty, col=legend_cols, lwd=legend_lwd, bty="n", bg="white")
        }
        
        abline(a=0, b=1, col=1, lty=3, lwd=2)
        abline(lm1, col=ols_col, lty=2, lwd=3)
        abline(ts1, col=ts_col, lty=1, lwd=3)
        cex_axis <- 1.2
        cex_labels <- 1.2
        axis(side=1, padj=-1.4, cex.axis=cex_axis, lwd.tick=1.5)
        axis(side=2, padj=1.4, cex.axis=cex_axis, lwd.tick=1.5)
        axis(side=3, padj=-1.4, labels=FALSE, lwd.tick=1.5)
        axis(side=4, padj=-1.4, labels=FALSE, lwd.tick=1.5)
        mtext(xlab, side=1, line=1.3, cex=cex_labels)
        mtext(ylab, side=2, line=1.3, cex=cex_labels)
    }

    # print some diagnostics
    me <- mean(mcd_obs - npn_mid, na.rm=T) # calculate mean and mean absolute error
    mae <- mean(abs(mcd_obs - npn_mid), na.rm=T)
    
    if(first_yes){
        npn_str <- npn_name
    }else{
        npn_str <- paste(npn_name, "(last yes)")
    }
    if(mid_duration) npn_str <- paste("midpoint of", npn_name)

    ols_intercept <- paste(round(lm1$coefficients[1], 1), paste("(", paste(round(lm1_confint[1,], 1), collapse=", "), ")", sep=""))
    ols_slope <- paste(round(lm1$coefficients[2], 1), paste("(", paste(round(lm1_confint[2,], 1), collapse=", "), ")", sep=""))
    ts_intercept <- paste(round(ts1$coefficients[1], 1), paste("(", paste(round(ts1_confint[1,], 1), collapse=", "), ")", sep=""))
    ts_slope <- paste(round(ts1$coefficients[2], 1), paste("(", paste(round(ts1_confint[2,], 1), collapse=", "), ")", sep=""))

    print(cat(paste(
        "\\multirow{2}{*}{", mcd_var, "} & ",
        "\\multirow{2}{*}{", npn_str, "} & ",
        "\\multirow{2}{*}{", nrow(DT_filt), "} & ",
        "\\multirow{2}{*}{", round(mae, 1), "} & ",
        "\\multirow{2}{*}{", round(me, 1), "} & ",
        "OLS & ",
        paste(ols_intercept, " & ", sep=""),
        paste(ols_slope, " & ", sep=""),
        paste(round(lm1_rsq, 2), " \\\\", sep=""),
    sep="")))
    print(cat(paste(
        "{} & ",
        "{} & ",
        "{} & ",
        "{} & ",
        "{} & ",
        "TS & ",
        paste(ts_intercept, " & ", sep=""),
        paste(ts_slope, " & ", sep=""),
        paste(round(ts1_rsq, 2), " \\\\", sep=""),
    sep="")))
        
    return(list(num_site_year=nrow(DT_filt), mcd_npn_me=me, mcd_npn_mae=mae, ols_mod=lm1, ols_conf=lm1_confint, ols_rsq=lm1_rsq, ts_mod=ts1, ts_conf=ts1_confint, ts_rsq=ts1_rsq, n_filter=n_filter))
}



# use DoNPNAnalysis() to quantify various MCD12Q2~USNPN relationships
all_forest_types <- c("Evergreen conifer", "Pine", "Evergreen broadleaf", "Deciduous conifer", "Deciduous broadleaf", "Drought deciduous broadleaf", "Semi-evergreen broadleaf")
enf_types <- c("Evergreen conifer", "Pine")
ebf_types <- c("Evergreen broadleaf", "Semi-evergreen broadleaf")
dbf_types <- c("Deciduous broadleaf", "Drought deciduous broadleaf", "Semi-evergreen broadleaf")
dnf_types <- c("Deciduous conifer")
grass_types <- c("Graminoid", "Forb")
# line_cols <- brewer.pal(7, "RdPu")[3:5]
sd_thresh <- 1e6

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
### All Forest Types
#------------------------------
## Greenup
# First Yes, Leaves ~ Greenup, DBF sfts, decid forest MCD12; first_yes_doy < 200 & Greenup > 50
filter_exp1 <- expression(phenophase_name == "Leaves" & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression(first_yes_doy < 200 & Greenup > 50)
filter_exps <- list(filter_exp1, filter_exp2)
leaves_fy_gup_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="Greenup", first_yes=T, makeplot=T, site_sd_thresh=sd_thresh)

# First Yes, Breaking leaf buds ~ Greenup, DBF sfts, decid forest MCD12, first_yes_doy < 200 & Greenup > 50
filter_exp1 <- expression(phenophase_name == "Breaking leaf buds" & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression(first_yes_doy < 200 & Greenup > 50)
filter_exps <- list(filter_exp1, filter_exp2)
# filter_exps <- filter_exp1 # check the necessity of filter expression 2
breakingleafbuds_fy_gup_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="Greenup", first_yes=T, makeplot=T, site_sd_thresh=sd_thresh)

# First Yes, Emerging leaves buds ~ Greenup, DBF sfts, decid forest MCD12, first_yes_doy < 200 & Greenup > 50
filter_exp1 <- expression(phenophase_name == "Emerging leaves" & species_functional_type %in% all_forest_types & mcd12q1_2006 %in% c(1:5, 6, 8))
filter_exp2 <- expression(first_yes_doy < 200 & Greenup > 50)
filter_exps <- list(filter_exp1, filter_exp2)
emergingleaves_fy_gup_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="Greenup", first_yes=T, makeplot=T, site_sd_thresh=sd_thresh)

#------------------------------
## MidGreenup
# Leaves, first_yes ~ MidGreenup, DBF sfts, decid forest MCD12, first_yes_doy < 200 & Greenup > 50
filter_exp1 <- expression(phenophase_name == "Leaves" & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression(first_yes_doy < 200 & MidGreenup > 50)
filter_exps <- list(filter_exp1, filter_exp2)
leaves_fy_midgup_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="MidGreenup", first_yes=T, makeplot=T, site_sd_thresh=sd_thresh)

# Midpoint Increasing leaf size ~ MidGreenup, mid_duration >= 30
filter_exp1 <- expression(phenophase_name == "Increasing leaf size" & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression(first_yes_doy < 200 & MidGreenup > 50)
filter_exps <- list(filter_exp1, filter_exp2)
leaves_fy_midgup_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="MidGreenup", mid_duration=TRUE, min_duration=30, makeplot=T, site_sd_thresh=sd_thresh)

#------------------------------
# Maturity
# Last yes, Increasing leaf size ~ Maturity, duration >= 30, DBF sfts, decid forest; MCD12 last_yes_doy < 200 & Maturity > 100
filter_exp1 <- expression(phenophase_name == "Increasing leaf size" & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression(last_yes_doy < 200 & Maturity > 100)
filter_exps <- list(filter_exp1, filter_exp2)
increasingleafsize_ly_maturity_dbf_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="Maturity", min_duration=30, first_yes=F, makeplot=T, site_sd_thresh=sd_thresh)

#------------------------------
# Peak
# Midpoint Leaves ~ Peak, duration >= 90, DBF sfts, decid forest
filter_exp1 <- expression(phenophase_name == "Leaves" & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression(last_yes_doy > 0)
filter_exps <- list(filter_exp1, filter_exp2)
increasingleafsize_ly_maturity_dbf_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="Peak", mid_duration=TRUE, min_duration=90, first_yes=F, makeplot=T, site_sd_thresh=sd_thresh)

#------------------------------
# Senescence
# First yes, Colored leaves ~ Senescence, DBF sfts, decid forest MCD12, first_yes_doy > 50
filter_exp1 <- expression(phenophase_name == "Colored leaves" & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression(first_yes_doy > 50)
filter_exps <- list(filter_exp1, filter_exp2)
coloredleaves_fy_sen_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="Senescence", first_yes=T, makeplot=T, site_sd_thresh=sd_thresh)

# First yes, Falling leaves ~ Senescence, DBF sfts, decid forest MCD12, first_yes_doy > 50
filter_exp1 <- expression(phenophase_name == "Falling leaves" & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression(first_yes_doy > 50)
filter_exps <- list(filter_exp1, filter_exp2)
fallingleaves_fy_sen_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="Senescence", first_yes=T, makeplot=T, site_sd_thresh=sd_thresh)

#------------------------------
# MidGreendown
# midpoint Colored leaves ~ MidGreendown, duration >=30, DBF sfts, decid forest MCD12
filter_exp1 <- expression(phenophase_name == "Colored leaves"  & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression(first_yes_doy > 0)
filter_exps <- list(filter_exp1, filter_exp2)
coloredleaves_fy_midgreendown_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="MidGreendown", mid_duration=TRUE, min_duration=30, first_yes=T, makeplot=T, site_sd_thresh=sd_thresh)

# midpoint Falling leaves ~ MidGreendown, duration >=30, DBF sfts, decid forest MCD12
filter_exp1 <- expression(phenophase_name == "Falling leaves"  & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression(first_yes_doy > 0)
filter_exps <- list(filter_exp1, filter_exp2)
coloredleaves_fy_midgreendown_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="MidGreendown", mid_duration=TRUE, min_duration=30, first_yes=T, makeplot=T, site_sd_thresh=sd_thresh)

# # Falling leaves, first_yes, MidGreendown, all forest; last_yes > 200; Dorm > 200
# # filter_exp <- expression(phenophase_name == "Falling leaves" & species_functional_type %in% all_forest_types & mcd12q1_2006 %in% c(1:5, 6, 8))
# filter_exp1 <- expression(phenophase_name == "Falling leaves" & species_functional_type %in% all_forest_types & mcd12q1_2006 %in% c(1:5, 6, 8))
# filter_exps <- filter_exp1
# fallingleaves_fy_midgreendown_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="MidGreendown", first_yes=T, makeplot=T, site_sd_thresh=sd_thresh)

# # All leaves colored, first_yes, MidGreendown, all forest;
# # filter_exp <- expression(phenophase_name == "All leaves colored" & species_functional_type %in% all_forest_types & mcd12q1_2006 %in% c(1:5, 6, 8) & first_yes_doy > 200)
# filter_exp1 <- expression(phenophase_name == "All leaves colored" & species_functional_type %in% all_forest_types & mcd12q1_2006 %in% c(1:5, 6, 8))
# filter_exp2 <- expression(first_yes_doy > 200)
# filter_exps <- list(filter_exp1, filter_exp2)
# allleavescolored_fy_midgreendown_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="MidGreendown", first_yes=T, makeplot=T, site_sd_thresh=sd_thresh)

#------------------------------
# Dormancy
# Last Yes, Leaves ~ Dormancy, duration>=90, DBF sfts, decid forest MCD12, last_yes_doy > 150 & Dormancy > 200
filter_exp1 <- expression(phenophase_name == "Leaves" & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression(last_yes_doy > 150 & Dormancy > 200)
filter_exps <- list(filter_exp1, filter_exp2)
leaves_ly_dor_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="Dormancy", min_duration=90, first_yes=F, makeplot=T, site_sd_thresh=sd_thresh)

# Last Yes, Falling leaves ~ Dormancy, duration>=30, DBF sfts, decid forest MCD12, last_yes_doy > 150 & Dormancy > 200
filter_exp1 <- expression(phenophase_name == "Falling leaves" & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression(last_yes_doy > 150 & Dormancy > 200)
filter_exps <- list(filter_exp1, filter_exp2)
leaves_ly_dor_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="Dormancy", min_duration=30, first_yes=F, makeplot=T, site_sd_thresh=sd_thresh)

# First Yes, All leaves colored ~ Dormancy, DBF sfts, decid forest MCD12, first_yes_doy > 150
filter_exp1 <- expression(phenophase_name == "All leaves colored" & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression(first_yes_doy > 150)
filter_exps <- list(filter_exp1, filter_exp2)
allleavescolored_fy_dor_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="Dormancy", first_yes=T, makeplot=T, site_sd_thresh=sd_thresh)

# First yes, All leaves fallen ~ Dormancy, DBF sfts, decid forest MCD12, first_yes_doy > 150
filter_exp1 <- expression(phenophase_name == "All leaves fallen" & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression(first_yes_doy > 150)
filter_exps <- list(filter_exp1, filter_exp2)
allleavesfallen_fy_dor_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="Dormancy", first_yes=T, makeplot=T, site_sd_thresh=sd_thresh)


#-------------------------------------------------
# output figure
png(file="~/Desktop/MCD12Q2_C6_Figures/USNPN_MCD12Q2_scatterplot.png", width=10, height=10, units="in", res=300)
layout(matrix(1:4, nrow=2, byrow=T))
par(mar=c(2.7, 2.7, 1.5, 1.5), tcl=0.25)
# First Yes, Leaves ~ Greenup, DBF sfts, decid forest MCD12, first_yes_doy < 200 & Greenup > 50
filter_exp1 <- expression(phenophase_name == "Leaves" & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression (first_yes_doy < 200 & Greenup > 50)
filter_exps <- list(filter_exp1, filter_exp2)
leaves_fy_gup_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="Greenup", first_yes=T, makeplot=T, site_sd_thresh=sd_thresh)
text(par()$usr[1], par()$usr[4], labels="a", adj=c(-0.5, 1.5), cex=1.85)

# First Yes, Leaves ~ MidGreenup, DBF sfts, decid forest MCD12, first_yes_doy < 200 & Greenup > 50
filter_exp1 <- expression(phenophase_name == "Leaves" & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression(first_yes_doy < 200 & MidGreenup > 50)
filter_exps <- list(filter_exp1, filter_exp2)
leaves_fy_midgup_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="MidGreenup", first_yes=T, makeplot=T, site_sd_thresh=sd_thresh)
text(par()$usr[1], par()$usr[4], labels="b", adj=c(-0.5, 1.5), cex=1.85)

# Senescence
# First yes, Colored leaves ~ Senescence, DBF sfts, decid forest MCD12, first_yes_doy > 50
filter_exp1 <- expression(phenophase_name == "Colored leaves" & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression(first_yes_doy > 50)
filter_exps <- list(filter_exp1, filter_exp2)
coloredleaves_fy_sen_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="Senescence", first_yes=T, makeplot=T, site_sd_thresh=sd_thresh)
text(par()$usr[1], par()$usr[4], labels="c", adj=c(-0.5, 1.5), cex=1.85)

# First Yes, All leaves fallen ~ Dormancy, all forest, lasst_yes_doy > 200
filter_exp1 <- expression(phenophase_name == "All leaves fallen" & species_functional_type %in% dbf_types & mcd12q1_2006 %in% c(3:5, 6, 8))
filter_exp2 <- expression(last_yes_doy > 200)
filter_exps <- list(filter_exp1, filter_exp2)
allleavesfallen_fy_dor_all_forest <- DoNPNAnalysis(npnDT_all, filter_exps, mcd_var="Dormancy", first_yes=T, makeplot=T, site_sd_thresh=sd_thresh)
text(par()$usr[1], par()$usr[4], labels="d", adj=c(-0.5, 1.5), cex=1.85)
dev.off()

