library(raster)
library(viridis)

mcd12q2c6_infile <- "~/Desktop/MCD12Q2.A2013001.h12v03.006.2018241135442.hdf"
mcd12q2c6_infile <- normalizePath(mcd12q2c6_infile)
product_year_date <- as.Date(gsub(".*A([0-9]{7}).*", "\\1", basename(mcd12q2c6_infile)), format="%Y%j")
doy_offset <- as.integer(product_year_date - as.Date("1970-1-1"))

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
GetSDS <- function(file_path, sds){
  # valid SDS: Greenup, MidGreenup, Peak, Senescence, MidGreendown, Dormancy, EVI_Minimum, EVI_Amplitude, NumCycles, QA_Detailed, QA_Overall
  return(paste("HDF4_EOS:EOS_GRID:\"", file_path, "\":MCD12Q2:", sds, sep = ""))
}

pheno_s <- stack(GetSDS(mcd12q2c6_infile, "Greenup"), GetSDS(mcd12q2c6_infile, "MidGreenup"), GetSDS(mcd12q2c6_infile, "Peak"), GetSDS(mcd12q2c6_infile, "Senescence"), GetSDS(mcd12q2c6_infile, "MidGreendown"), GetSDS(mcd12q2c6_infile, "Dormancy"), GetSDS(mcd12q2c6_infile, "EVI_Minimum"))

mcd12q2c6_na <- 32767
gup_r <- raster(GetSDS(mcd12q2c6_infile, "Greenup"))
NAvalue(gup_r) <- mcd12q2c6_na
gup_r <- gup_r - doy_offset
midgup_r <- raster(GetSDS(mcd12q2c6_infile, "MidGreenup"))
NAvalue(midgup_r) <- mcd12q2c6_na
midgup_r <- midgup_r - doy_offset

layout(matrix(1:2, nrow=1))
par(mar=rep(0, 4))
gup_v <- values(gup_r)
qs <- quantile(gup_v, c(0, 0.02, 0.98, 1), na.rm=T)
breaks <- unique(round(c(qs[1], seq(qs[2], qs[3], len=254), qs[4])))
plot(gup_r, breaks=breaks, col=inferno(length(breaks) - 1))

gup_diff_r <- midgup_r - gup_r
gup_diff_v <- values(gup_diff_r)
qs <- quantile(gup_diff_v, c(0, 0.02, 0.98, 1), na.rm=T)
breaks <- unique(round(c(qs[1], seq(qs[2], qs[3], len=254), qs[4])))
plot(gup_diff_r, breaks=breaks, col=inferno(length(breaks) - 1))

plot(gup_diff_r, breaks=c(qs[1], 40, qs[4]), col=c("grey", "red"), maxpixels=2.5e6)
myc <- click(gup_diff_r, xy=T, cell=T)
cells_to_query <- myc$cell[myc$value > 40]

r <- raster(matrix(1:9, nrow=3, byrow=T))

s <- stack(raster(matrix(1:9, nrow=3, byrow=T)), raster(matrix(10:18, nrow=3, byrow=T)), raster(matrix(19:27, nrow=3, byrow=T)))
s[c(2,4)]
plot(s)

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
in_dir <- "/rsstu/users/j/jmgray2/SEAL/C6SNOW"
in_files <- dir(in_dir, pattern="C6Snow.*.tif", full=T)
in_dates <- as.Date(gsub("C6Snow.h[0-9]{2}v[0-9]{2}.([0-9]{4}.[0-9]{3}).tif", "\\1", basename(in_files)), format="%Y.%j")
start_date <- as.Date("2012-1-1")
end_date <- as.Date("2014-12-31")
in_files <- in_files[in_dates >= start_date & in_dates <= end_date]
in_dates <- in_dates[in_dates >= start_date & in_dates <= end_date]
in_files <- in_files[order(in_dates)]
cells_to_query <- c(5634048, 5326931, 869289, 857296, 612809, 234598, 333006, 71791, 11796)

snow_stack <- stack(in_files)
system.time(tmp <- snow_stack[cells_to_query])
tmp_df <- data.frame(tmp)
names(tmp_df) <- paste(c("EVI2", "NDSI", "NDMI", "A2QA", "A2SNOW", "LST", "LSTQA"), rep(strftime(in_dates, format="%Y%j"), each=7), sep="_")
tmp_df$cell <- cells_to_query

save(tmp_df, file=file.path(in_dir, "SnowSamples_h12v03.Rdata"))

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
load("~/Desktop/SnowSamplesFullData.Rdata")
mcd12q2c6_infile <- "~/Desktop/MCD12Q2.A2013001.h12v03.006.2018241135442.hdf"
mcd12q2c6_infile <- normalizePath(mcd12q2c6_infile)
pheno_s <- stack(GetSDS(mcd12q2c6_infile, "Greenup"), GetSDS(mcd12q2c6_infile, "MidGreenup"), GetSDS(mcd12q2c6_infile, "Peak"), GetSDS(mcd12q2c6_infile, "Senescence"), GetSDS(mcd12q2c6_infile, "MidGreendown"), GetSDS(mcd12q2c6_infile, "Dormancy"), GetSDS(mcd12q2c6_infile, "EVI_Minimum"))
pheno_s <- subset(pheno_s, seq(1, 14, by=2))
NAvalue(pheno_s) <- 32767
# doy_offset <- as.integer(product_year_date - as.Date("1970-1-1"))

evi2_inds <- seq(1, dim(tmp_df)[2] - 1, by=7)
ndsi_inds <- seq(2, dim(tmp_df)[2] - 1, by=7)
ndmi_inds <- seq(3, dim(tmp_df)[2] - 1, by=7)
a2qa_inds <- seq(4, dim(tmp_df)[2] - 1, by=7)
a2snow_inds <- seq(5, dim(tmp_df)[2] - 1, by=7)
lst_inds <- seq(6, dim(tmp_df)[2] - 1, by=7)
lstqa_inds <- seq(7, dim(tmp_df)[2] - 1, by=7)

# tmp_df[which_pix, a2snow_inds]

# Plotting
PlotSnowSeries <- function(which_pix){
    # which_pix <- 1
    which_cell <- cells_to_query[which_pix]
    evi2_scale <- 1e4
    lst_scale <- 10
    ndsi_thresh <- 0.2
    bg_quant <- 0.05
    qual_cex <- c(1.5, 1.25, 1, 0.5)
    GetTranspCol <- function(the_col, alpha=0.5) rgb(t(col2rgb(the_col)), max=255, alpha=alpha*255)
    qual_cols <- rev(sapply(magma(10)[c(3, 5, 7, 9)], GetTranspCol))
    snow_col <- rgb(0, 0.3, 1)
    grey_col <- rgb(0.3, 0.3, 0.3, 0.5)

    evi2 <- tmp_df[which_pix, evi2_inds] / evi2_scale
    ndsi <- tmp_df[which_pix, ndsi_inds] / evi2_scale
    lst <- tmp_df[which_pix, lst_inds] / lst_scale
    a2qa <- as.integer(tmp_df[which_pix, a2qa_inds])
    snow <- as.logical(tmp_df[which_pix, a2snow_inds])
    snow_ndsi <- ndsi >= ndsi_thresh
    q2_vals <- pheno_s[which_cell]
    q2_pheno_dates <- as.Date(q2_vals[1:6], origin=as.Date("1970-1-1"))
    q2_bg <- q2_vals[7] / evi2_scale

    cexs <- rep(1.25, length(evi2))
    cexs[snow | snow_ndsi] <- 0.75
    pchs <- rep(16, length(evi2))
    pchs[snow & snow_ndsi] <- 8
    pchs[(snow & !snow_ndsi) | (snow & is.na(snow_ndsi))] <- 3
    pchs[(!snow & snow_ndsi) | (is.na(snow) & snow_ndsi)] <- 4
    cols <- qual_cols[a2qa + 1]
    cols[snow | snow_ndsi] <- snow_col
    pchs[snow | snow_ndsi]

    bg_val <- quantile(evi2[!snow & !snow_ndsi & !is.na(snow) & !is.na(snow_ndsi)], bg_quant)

    layout(matrix(1:2, nrow=2))
    par(mar=c(1, 4, 1, 1))
    plot(in_dates, evi2, xlab="", ylab="EVI2", type="p", pch=pchs, cex=cexs, col=cols)
    abline(h=bg_val, lwd=1.5, lty=2, col=grey_col)
    abline(h=q2_bg, lwd=1.5, lty=2, col="black")
    abline(v=q2_pheno_dates)
    legend("topleft", legend=c("A2QA=0", "A2QA=1", "A2QA=2", "A2QA=3", "NDSI Snow", "A2 Snow", "Both Snow", "5% BG", "Q2 Min"), pch=c(rep(16, 4), 4, 3, 8, NA, NA), col=c(qual_cols, rep(snow_col, 3), grey_col, "black"), pt.cex=1.6, bty="n", lty=c(rep(NA, 7), 2, 2))
    plot(in_dates, lst, xlab="", ylab="LST", type="p", pch=16, col=snow_col)
    abline(h=273.15, lwd=1.5, lty=2, col=grey_col)
    abline(v=q2_pheno_dates)
}

PlotSnowSeries(2)

# site_sin_coords <- data.frame(xyFromCell(raster(pheno_s, 1), cells_to_query))
# coordinates(site_sin_coords) <- ~x+y
# projection(site_sin_coords) <- projection(raster(pheno_s, 1))
# site_latlong_coords <- spTransform(site_sin_coords, CRS("+proj=longlat +ellps=GRS80"))
# site_data <- data.frame(site_id=1:length(cells_to_query), mod_cell=cells_to_query)
# site_latlong <- SpatialPointsDataFrame(site_latlong_coords, data=site_data)
# library(rgdal)
# writeOGR(site_latlong, ".", "SnowTestSites", driver="ESRI Shapefile")
# setwd("~/Desktop")


gup_diff_r <- midgup_r - gup_r
gup_diff_v <- values(gup_diff_r)
qs <- quantile(gup_diff_v, c(0, 0.02, 0.98, 1), na.rm=T)
breaks <- unique(round(c(qs[1], seq(qs[2], qs[3], len=254), qs[4])))
plot(gup_diff_r, breaks=breaks, col=inferno(length(breaks) - 1), legend=F, maxpixels=2.5e6)
PlotLegend(breaks, colorRampPalette(inferno(length(breaks) - 1)), legend_at = NA, labels = NA, num_breaks = 5)
points(xyFromCell(gup_diff_r, cells_to_query), col="red", pch=1, cex=1)

plot(xyFromCell(gup_diff_r, cells_to_query), col="red", pch=1, cex=1)
text(xyFromCell(gup_diff_r, cells_to_query), text=1:9)

