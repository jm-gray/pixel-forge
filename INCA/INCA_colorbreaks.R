library(viridis)
library(scico)
library(rgdal)
library(RColorBrewer)
data_dir <- "/Volumes/users/j/jmgray2/SEAL/INCA/INCAproducts"

# BG col: #333333
PlotLegend <- function(breaks, colpal, LEGENDAXISCEX=1, LEGENDMAINCEX=1, LEGENDWIDTH=1, LEGENDLABELCOL="white"){
    legend_at <- round(seq(breaks[2], breaks[length(breaks) - 1], len=7))
    legend_at_date <- legend_at
    legend_labels <- c(paste("<", legend_at_date[1]), as.character(legend_at_date[2:(length(legend_at_date) - 1)]), paste(">", legend_at_date[length(legend_at_date)]))
    plot(
        raster(matrix(legend_at[1]:legend_at[length(legend_at)])),
        legend.only=T,
        col=colpal(length(breaks)-1),
        legend.width=LEGENDWIDTH,
        axis.args=list(at=legend_at, labels=legend_labels, cex.axis=LEGENDAXISCEX, col=LEGENDLABELCOL, col.axis=LEGENDLABELCOL),
        legend.args=list(text="", side=3, font=2, line=0.5, cex=LEGENDMAINCEX)
    )
}

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# MidGup
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Median MidGup
r <- raster(file.path(data_dir, "CONUS_MidGup_Median_NAD83_02deg.tif"))
v <- values(r)
qs <- quantile(v, c(0, 0.02, 0.98, 1), na.rm=T)
num_cols <- 35
midgup_med_breaks <- c(qs[1], seq(qs[2], qs[3], len=num_cols - 1), qs[4])
midgup_med_breaks <- unique(round(midgup_med_breaks))
midgup_med_cols <- rev(plasma(length(midgup_med_breaks) - 1))
colpal <- colorRampPalette(rev(plasma(length(midgup_med_breaks) - 1)))
par(mar=rep(0, 4), oma=rep(0.5, 4), bg=rgb(0.2, 0.2, 0.2))
plot(r, breaks=midgup_med_breaks, col=rev(plasma(length(midgup_med_breaks)-1)), maxpixels=ncell(r), colNA=rgb(0.2,0.2,0.2), legend=F, axes=F, box=F)
PlotLegend(midgup_med_breaks, colpal, LEGENDWIDTH=2)

# MAD MidGup
r <- raster(file.path(data_dir, "CONUS_MidGup_MAD_NAD83_02deg.tif"))
v <- values(r)
qs <- quantile(v, c(0, 0.05, 0.95, 1), na.rm=T)
num_cols <- 35
midgup_mad_breaks <- c(qs[1], seq(qs[2], qs[3], len=num_cols - 1), qs[4])
midgup_mad_cols <- rev(viridis(length(midgup_mad_breaks) - 1))
colpal <- colorRampPalette(rev(viridis(length(midgup_mad_breaks) - 1)))
par(mar=rep(0, 4), oma=rep(0.5, 4), bg=rgb(0.2, 0.2, 0.2))
plot(r, breaks=midgup_mad_breaks, col=rev(viridis(length(midgup_mad_breaks)-1)), maxpixels=ncell(r), colNA=rgb(0.2,0.2,0.2), legend=F, axes=F, box=F)
PlotLegend(midgup_mad_breaks, colpal, LEGENDWIDTH=2)

# TSslope MidGup
r <- raster(file.path(data_dir, "CONUS_MidGup_TSslope_NAD83_02deg.tif"))
v <- values(r)
qs <- quantile(v, c(0, 0.02, 0.98, 1), na.rm=T)
num_cols <- 35
midgup_tsslope_breaks <- c(-90, seq(-2, 2, len=num_cols - 1), 90)
pal <- colorRampPalette(rev(scico(100, pal="tokyo")))
alt_pal <- colorRampPalette(rev(brewer.pal(9, "RdBu")))
midgup_tsslope_cols <- pal(length(midgup_tsslope_breaks) - 1)
midgup_tsslope_cols_alt <- alt_pal(length(midgup_tsslope_breaks) - 1)
par(mar=rep(0, 4), oma=rep(0.5, 4), bg=rgb(0.2, 0.2, 0.2))
plot(r, breaks=midgup_tsslope_breaks, col=pal(length(midgup_tsslope_breaks)-1), maxpixels=ncell(r), colNA=rgb(0.2,0.2,0.2), legend=F, axes=F, box=F)
PlotLegend(midgup_tsslope_breaks, pal, LEGENDWIDTH=2)

plot(r, breaks=midgup_tsslope_breaks, col=alt_pal(length(midgup_tsslope_breaks)-1), maxpixels=ncell(r), colNA=rgb(0.2,0.2,0.2), legend=F, axes=F, box=F)
PlotLegend(midgup_tsslope_breaks, alt_pal, LEGENDWIDTH=2)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# MidGdown
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Median MidGdown
r <- raster(file.path(data_dir, "CONUS_MidGdown_Median_NAD83_02deg.tif"))
v <- values(r)
qs <- quantile(v, c(0, 0.02, 0.98, 1), na.rm=T)
num_cols <- 35
midgdown_med_breaks <- c(qs[1], seq(qs[2], qs[3], len=num_cols - 1), qs[4])
midgdown_med_breaks <- unique(round(midgdown_med_breaks))
midgdown_med_cols <- rev(plasma(length(midgdown_med_breaks) - 1))
colpal <- colorRampPalette(rev(plasma(length(midgdown_med_breaks) - 1)))
par(mar=rep(0, 4), oma=rep(0.5, 4), bg=rgb(0.2, 0.2, 0.2))
plot(r, breaks=midgdown_med_breaks, col=rev(plasma(length(midgdown_med_breaks)-1)), maxpixels=ncell(r), colNA=rgb(0.2,0.2,0.2), legend=F, axes=F, box=F)
PlotLegend(midgdown_med_breaks, colpal, LEGENDWIDTH=2)
 
# MAD MidGdown
r <- raster(file.path(data_dir, "CONUS_MidGdown_MAD_NAD83_02deg.tif"))
v <- values(r)
qs <- quantile(v, c(0, 0.05, 0.95, 1), na.rm=T)
num_cols <- 35
midgdown_mad_breaks <- c(qs[1], seq(qs[2], qs[3], len=num_cols - 1), qs[4])
midgdown_mad_cols <- rev(viridis(length(midgdown_mad_breaks) - 1))
colpal <- colorRampPalette(rev(viridis(length(midgdown_mad_breaks) - 1)))
par(mar=rep(0, 4), oma=rep(0.5, 4), bg=rgb(0.2, 0.2, 0.2))
plot(r, breaks=midgdown_mad_breaks, col=rev(viridis(length(midgdown_mad_breaks)-1)), maxpixels=ncell(r), colNA=rgb(0.2,0.2,0.2), legend=F, axes=F, box=F)
PlotLegend(midgdown_mad_breaks, colpal, LEGENDWIDTH=2)

# TSslope MidGdown
r <- raster(file.path(data_dir, "CONUS_MidGdown_TSslope_NAD83_02deg.tif"))
v <- values(r)
qs <- quantile(v, c(0, 0.02, 0.98, 1), na.rm=T)
num_cols <- 35
midgdown_tsslope_breaks <- c(-90, seq(-2, 2, len=num_cols - 1), 90)
pal <- colorRampPalette(rev(scico(100, pal="tokyo")))
alt_pal <- colorRampPalette(rev(brewer.pal(9, "RdBu")))
midgdown_tsslope_cols <- pal(length(midgdown_tsslope_breaks) - 1)
midgdown_tsslope_cols_alt <- alt_pal(length(midgdown_tsslope_breaks) - 1)
par(mar=rep(0, 4), oma=rep(0.5, 4), bg=rgb(0.2, 0.2, 0.2))
plot(r, breaks=midgdown_tsslope_breaks, col=pal(length(midgdown_tsslope_breaks)-1), maxpixels=ncell(r), colNA=rgb(0.2,0.2,0.2), legend=F, axes=F, box=F)
PlotLegend(midgdown_tsslope_breaks, pal, LEGENDWIDTH=2)

plot(r, breaks=midgdown_tsslope_breaks, col=alt_pal(length(midgdown_tsslope_breaks)-1), maxpixels=ncell(r), colNA=rgb(0.2,0.2,0.2), legend=F, axes=F, box=F)
PlotLegend(midgdown_tsslope_breaks, alt_pal, LEGENDWIDTH=2)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# EVIarea
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Median EVIArea
r <- raster(file.path(data_dir, "CONUS_EVIarea_Median_NAD83_02deg.tif"))
v <- values(r)
qs <- quantile(v, c(0, 0.02, 0.98, 1), na.rm=T)
num_cols <- 35
eviarea_med_breaks <- c(qs[1], seq(qs[2], qs[3], len=num_cols - 1), qs[4])
eviarea_med_cols <- rev(plasma(length(eviarea_med_breaks) - 1))
colpal <- colorRampPalette(rev(plasma(length(eviarea_med_breaks) - 1)))
par(mar=rep(0, 4), oma=rep(0.5, 4), bg=rgb(0.2, 0.2, 0.2))
plot(r, breaks=eviarea_med_breaks, col=rev(plasma(length(eviarea_med_breaks)-1)), maxpixels=ncell(r), colNA=rgb(0.2,0.2,0.2), legend=F, axes=F, box=F)
PlotLegend(eviarea_med_breaks, colpal, LEGENDWIDTH=2)
 
# MAD MidGdown
r <- raster(file.path(data_dir, "CONUS_EVIarea_MAD_NAD83_02deg.tif"))
v <- values(r)
qs <- quantile(v, c(0, 0.05, 0.95, 1), na.rm=T)
num_cols <- 35
eviarea_mad_breaks <- c(qs[1], seq(qs[2], qs[3], len=num_cols - 1), qs[4])
eviarea_mad_cols <- rev(viridis(length(eviarea_mad_breaks) - 1))
colpal <- colorRampPalette(rev(viridis(length(eviarea_mad_breaks) - 1)))
par(mar=rep(0, 4), oma=rep(0.5, 4), bg=rgb(0.2, 0.2, 0.2))
plot(r, breaks=eviarea_mad_breaks, col=rev(viridis(length(eviarea_mad_breaks)-1)), maxpixels=ncell(r), colNA=rgb(0.2,0.2,0.2), legend=F, axes=F, box=F)
PlotLegend(eviarea_mad_breaks, colpal, LEGENDWIDTH=2)

# TSslope MidGdown
r <- raster(file.path(data_dir, "CONUS_EVIarea_TSslope_NAD83_02deg.tif"))
v <- values(r)
qs <- quantile(v, c(0, 0.02, 0.98, 1), na.rm=T)
num_cols <- 35
eviarea_tsslope_breaks <- c(-6, seq(-1.5, 1.5, len=num_cols - 1), 6)
pal <- colorRampPalette(rev(scico(100, pal="tokyo")))
alt_pal <- colorRampPalette(rev(brewer.pal(9, "RdBu")))
eviarea_tsslope_cols <- pal(length(eviarea_tsslope_breaks) - 1)
eviarea_tsslope_cols_alt <- alt_pal(length(eviarea_tsslope_breaks) - 1)
par(mar=rep(0, 4), oma=rep(0.5, 4), bg=rgb(0.2, 0.2, 0.2))
plot(r, breaks=eviarea_tsslope_breaks, col=pal(length(eviarea_tsslope_breaks)-1), maxpixels=ncell(r), colNA=rgb(0.2,0.2,0.2), legend=F, axes=F, box=F)
PlotLegend(eviarea_tsslope_breaks, pal, LEGENDWIDTH=2)

plot(r, breaks=eviarea_tsslope_breaks, col=alt_pal(length(eviarea_tsslope_breaks)-1), maxpixels=ncell(r), colNA=rgb(0.2,0.2,0.2), legend=F, axes=F, box=F)
PlotLegend(eviarea_tsslope_breaks, alt_pal, LEGENDWIDTH=2)



setwd("/Volumes/users/j/jmgray2/SEAL/INCA/INCAproducts")
save(midgup_med_breaks, midgup_med_cols, midgup_mad_breaks, midgup_mad_cols, midgup_tsslope_breaks, midgup_tsslope_cols, midgup_tsslope_cols_alt, midgdown_med_breaks, midgdown_med_cols, midgdown_mad_breaks, midgdown_mad_cols, midgdown_tsslope_breaks, midgdown_tsslope_cols, midgdown_tsslope_cols_alt, eviarea_med_breaks, eviarea_med_cols, eviarea_mad_breaks, eviarea_mad_cols, eviarea_tsslope_breaks, eviarea_tsslope_cols, eviarea_tsslope_cols_alt, file="sld_data_midgup_midgdown_eviarea.Rdata")
