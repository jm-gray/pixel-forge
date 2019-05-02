library(viridis)
library(scico)
library(rgdal)
library(RColorBrewer)

# BG col: #333333

#######################
# Median MidGup
r <- raster("CONUS_MidGup_Median_NAD83_004deg.tif")
r_25 <- raster("CONUS_MidGup_Median_NAD83_02deg.tif")
v <- values(r)
qs <- quantile(v, c(0, 0.02, 0.98, 1), na.rm=T)
v_25 <- values(r_25)
qs_25 <- quantile(v_25, c(0, 0.02, 0.98, 1), na.rm=T)
num_cols <- 35
med_breaks <- c(qs[1], seq(qs[2], qs[3], len=num_cols - 1), qs[4])
med_breaks <- unique(round(med_breaks))
med_cols <- rev(plasma(length(breaks) - 1))
plot(r, breaks=breaks, col=rev(plasma(length(breaks)-1)), maxpixels=5e6, colNA=rgb(0.2,0.2,0.2))

#######################
# MAD MidGup
r <- raster("CONUS_MidGup_MAD_NAD83_004deg.tif")
v <- values(r)
qs <- quantile(v, c(0, 0.05, 0.95, 1), na.rm=T)
r_25 <- raster("CONUS_MidGup_MAD_NAD83_02deg.tif")
v_25 <- values(r_25)
qs_25 <- quantile(v_25, c(0, 0.05, 0.95, 1), na.rm=T)
num_cols <- 35
mad_breaks <- c(min(c(qs_25, qs)), seq(qs[2], qs[3], len=num_cols - 1), max(c(qs_25, qs)))
plot(r, breaks=mad_breaks, col=rev(viridis(length(mad_breaks)-1)), maxpixels=5e6, colNA=rgb(0.2,0.2,0.2))
mad_cols <- rev(viridis(length(mad_breaks) - 1))

#######################
# MAD TSslope
r <- raster("CONUS_MidGup_TSslope_NAD83_004deg.tif")
v <- values(r)
qs <- quantile(v, c(0, 0.02, 0.98, 1), na.rm=T)
r_25 <- raster("CONUS_MidGup_TSslope_NAD83_02deg.tif")
v_25 <- values(r_25)
qs_25 <- quantile(v_25, c(0, 0.02, 0.98, 1), na.rm=T)
num_cols <- 35
tsslope_breaks <- c(-90, seq(-2, 2, len=num_cols - 1), 90)
pal <- colorRampPalette(rev(scico(100, pal="tokyo")))
alt_pal <- colorRampPalette(rev(brewer.pal(9, "RdBu")))
tsslope_cols <- pal(length(tsslope_breaks) - 1)
tsslope_cols_alt <- alt_pal(length(tsslope_breaks) - 1)
plot(r, breaks=tsslope_breaks, col=pal(length(tsslope_breaks) - 1), maxpixels=5e6, colNA=rgb(0.2,0.2,0.2))

# pal <- colorRampPalette(rev(brewer.pal(9, "RdBu")))
# plot(r, breaks=tsslope_breaks, col=pal(length(tsslope_breaks) - 1), maxpixels=5e6, colNA=rgb(0.2,0.2,0.2))


midgup_med_breaks <- med_breaks
midgup_med_colorss <- med_cols
midgup_mad_breaks <- mad_breaks
midgup_mad_colors <- mad_cols
midgup_tsslope_breaks <- tsslope_breaks
midgup_tsslope_colors <- tsslope_cols
alt_midgup_tsslope_colors <- tsslope_cols_alt

save(midgup_med_breaks, midgup_med_colorss, midgup_mad_breaks, midgup_mad_colors, midgup_tsslope_breaks, midgup_tsslope_colors, alt_midgup_tsslope_colors, file="sld_data_midgup.Rdata")
getwd()
