library(raster)
library(data.table)
library(viridis)
library(RColorBrewer)

#--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
GetSDS <- function(file_path, sds=NULL){
  # gets the SDS names for MCD12Q2 EOS HDFs
  # valid sds options: Greenup, MidGreenup, Peak, Senescence, MidGreendown, Dormancy, EVI_Minimum, EVI_Amplitude, NumCycles, QA_Detailed, QA_Overall
  # if sds is not provided, filenames for all SDS are returned
  all_sds <- c("NumCycles", "Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy", "EVI_Minimum", "EVI_Amplitude", "EVI_Area", "QA_Overall", "QA_Detailed")
  if(is.null(sds)){
    return(paste("HDF4_EOS:EOS_GRID:\"", file_path, "\":MCD12Q2:", all_sds, sep = ""))
  }else{
    return(paste("HDF4_EOS:EOS_GRID:\"", file_path, "\":MCD12Q2:", sds, sep = ""))
  }
}

#--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
GetSDSQ1 <- function(file_path, sds="LC_Type1"){
    # returns the proper EOS-HDF name for a particular MCD12Q1 C6 SDS
    return(paste("HDF4_EOS:EOS_GRID:\"", file_path, "\":MCD12Q1:", sds, sep = ""))
}

#--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
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

#--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
inca_dir <- "/rsstu/users/j/jmgray2/SEAL/INCA/INCAglobaloutput"
output_dir <- "/rsstu/users/j/jmgray2/SEAL/INCA/C6_Figures"
mcd12q2c6_dir <- "/rsstu/users/j/jmgray2/SEAL/INCA/MCD12Q2C6/MCD12Q2"
lc_dir <- "/rsstu/users/j/jmgray2/SEAL/MCD12Q1/2016/001"
europe_tiles <- apply(expand.grid(h=17:19, v=2:5), 1, function(x) paste(paste("h", formatC(x[1], width=2, flag="0"), sep=""), paste("v", formatC(x[2], width=2, flag="0"), sep=""), sep=""))

#--------
# Get LC data
europe_lc_files <- sapply(europe_tiles, function(tile) dir(lc_dir, pattern=paste("MCD12.*", tile, ".*hdf$", sep=""), full=T))
out_vrt_lc <- file.path(output_dir, "Europe_LC.vrt")
BuildVRT(GetSDSQ1(europe_lc_files), out_file=out_vrt_lc, band=1, vrtnodata = -9999)
europe_lc_r <- raster(out_vrt_lc)

#--------
# get the median and MAD INCA data for Europe
metric <- "Greenup"
# metric <- "EVI_Area"
europe_inca_files <- sapply(europe_tiles, function(x, metric) dir(inca_dir, pattern=paste(x, metric, sep=".*"), full=T), metric=metric)
# return(c(round(med), MAD, round(avg), stddev, ts.slope, ts.pvalue, ann.median.anoms, ann.median.anoms.mad))
out_vrt_median <- file.path(output_dir, paste("Europe_INCAmedian_", metric, ".vrt", sep=""))
BuildVRT(europe_inca_files, out_file=out_vrt_median, band=1, vrtnodata = -9999)
out_vrt_mad <- file.path(output_dir, paste("Europe_INCAmad_", metric, ".vrt", sep=""))
BuildVRT(europe_inca_files, out_file=out_vrt_mad, band=2, vrtnodata = -9999)

europe_median_r <- raster(out_vrt_median)
europe_mad_r <- raster(out_vrt_mad)
europe_mad_r[europe_mad_r == 0] <- NA

#--------
# get the MCD12Q2 C6 2018 data for Europe
europe_2018_files <- sapply(europe_tiles, function(tile) dir(file.path(mcd12q2c6_dir, "2018"), pattern=paste("MCD12.*", tile, ".*hdf$", sep=""), full=T))

# out_europe_midgreenup_2018_vrt <- file.path(output_dir, paste("Europe_"))
out_vrt_2018 <- file.path(output_dir, paste("Europe_2018_", metric, ".vrt", sep=""))
BuildVRT(GetSDS(europe_2018_files, sds=metric), band=1, vrtnodata=-9999, out_file=out_vrt_2018)

europe_2018_r <- raster(out_vrt_2018)
NAvalue(europe_2018_r) <- 32767

if(metric %in% c("EVI_Minimum", "EVI_Amplitude", "EVI_Area")){
    europe_2018_anomaly_mad_r <- (europe_2018_r - europe_median_r) / europe_mad_r
}else{
    doy_offset <- as.integer(as.Date("2018-1-1") - as.Date("1970-1-1"))
    europe_2018_doy_r <- europe_2018_r - doy_offset
    europe_2018_anomaly_mad_r <- (europe_2018_doy_r - europe_median_r) / europe_mad_r
}

#--------
# plot the anomaly as multiples of MAD across Europe
europe_shp <- shapefile("/rsstu/users/j/jmgray2/SEAL/INCA/continents_sin/continent.shp")
europe_shp <- europe_shp[europe_shp$CONTINENT %in% c("Europe", "Africa"),]

land_color <- rgb(0.4, 0.4, 0.4)
water_color <- rgb(0.2, 0.2, 0.2)

out_pdf <- file.path(output_dir, paste("Europe_2018AnomalyMAD_", metric, ".pdf", sep=""))

qs <- quantile(europe_2018_anomaly_mad_r, c(0, 0.02, 0.98, 1), na.rm=T)
extreme_value <- max(c(abs(qs[c(1:4)]), 3))
# breaks <- unique(round(c(qs[1], seq(qs[2], qs[3], len=254), qs[4])))
mad_breaks <- c(-extreme_value, seq(-3, 3, len=254), extreme_value)
if(metric %in% c("EVI_Minimum", "EVI_Amplitude", "EVI_Area")){
    anomaly_pal <- colorRampPalette(brewer.pal(11, "BrBG"))
}else{
    anomaly_pal <- colorRampPalette(c(rev(brewer.pal(7, "Blues")), brewer.pal(7, "Reds")))
}


fig_height <- 19
fig_width <- 18
dpi <- 300
maxpixels <- fig_height * dpi * fig_width * dpi
pdf(file=out_pdf, width=fig_width, height=fig_height)
par(mar=rep(1, 4), oma=c(0, 0, 0, 6))
# plot(europe_shp, col=land_color, border=NA, bg=water_color, xaxt="n", yaxt="n", xlab="", ylab="", xlim=extent(europe_2018_anomaly_mad_r)[1:2], ylim=extent(europe_2018_anomaly_mad_r)[3:4])
# plot(europe_2018_anomaly_mad_r, col=anomaly_pal(length(mad_breaks)), breaks=mad_breaks, legend=F, xaxt="n", yaxt="n", xlab="", ylab="", maxpixels=maxpixels, colNA=rgb(0.2, 0.2, 0.2))
ylim <- c(3778000, 7735000)
plot(europe_2018_anomaly_mad_r, maxpixels=1, legend=F, xaxt="n", yaxt="n", xlab="", ylab="", colNA=water_color, col=water_color, bg=water_color, ylim=ylim)
plot(europe_shp, col=land_color, border=NA, bg=water_color, add=T)
plot(europe_2018_anomaly_mad_r, col=anomaly_pal(length(mad_breaks)), breaks=mad_breaks, legend=F, maxpixels=maxpixels, colNA=NA, add=TRUE)
legend_at <- seq(-3, 3, by=1)
legend_labels <- c(paste("<=", legend_at[1]), legend_at[2:(length(legend_at) - 1)], paste(">=", legend_at[length(legend_at)]))
plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=anomaly_pal(length(mad_breaks) - 1), legend.width=1.5, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=2), legend.args=list(text="", side=3, font=2, line=0.5, cex=3))
# plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=plasma(length(med_breaks) - 1), legend.width=LEGENDWIDTH, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=LEGENDAXISCEX), legend.args=list(text="DOY", side=3, font=2, line=0.5, cex=LEGENDMAINCEX))
dev.off()

#--------
# Check LC-specific patterns
pheno_sample_dir <- "/rsstu/users/j/jmgray2/SEAL/INCA/GlobalPhenoSample"
europe_rdata_files <- sapply(europe_tiles, function(x) dir(pheno_sample_dir, pattern=x, full=T))
load(europe_rdata_files[1])
DTeurope <- DTfinal
for(rdata_file in europe_rdata_files[2:length(europe_rdata_files)]){
    load(rdata_file)
    DTeurope <- rbind(DTeurope, DTfinal)
}

# as.integer(as.Date(paste(year, "-1-1", sep="")) - as.Date("1970-1-1"))
# DTeurope[, value - as.integer(as.Date(paste(year, "-1-1", sep="")) - as.Date("1970-1-1"))]
# findInterval(c(30, 31, 70, 69), seq(30, 70, by=2), all.inside=T)
doy_phenometrics <- paste(rep(c("Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy"), each=2), 1:2, sep="")
DTeurope[, value_conv := value]
DTeurope[phenometric %in% doy_phenometrics, value_conv := value - as.integer(as.Date(paste(year, "-1-1", sep="")) - as.Date("1970-1-1"))]


DTmed <- DTeurope[year %in% 2001:2016, .(x=x[1], y=y[1], lat=lat[1], lon=lon[1], igbp=igbp[1], med_val=median(value_conv, na.rm=T), mad_val=mad(value_conv, na.rm=T)), by=.(tile, cell_num, phenometric)]
DTmed$val_2018 <- DTeurope[year == 2018, value_conv]
DTmed[, anomaly2018 := val_2018 - med_val]
DTmed[, anomaly2018mad := anomaly2018 / mad_val]
save(DTmed, file="Europe_DT_median_2018.Rdata")

the_metric <- "Greenup1"
lat_intervals <- seq(30, 70, by=0.5)
# tmp <- DTmed[phenometric == the_metric & igbp %in% c(12, 14), median(anomaly2018mad, na.rm=T), by=findInterval(lat, lat_intervals, all.inside=T)]
tmp <- DTmed[phenometric == the_metric, median(anomaly2018mad, na.rm=T), by=findInterval(lat, lat_intervals, all.inside=T)]
tmp[, middle_lat := lat_intervals[findInterval] + diff(lat_intervals)[1]/2]
tmp <- tmp[order(middle_lat)]
plot(tmp[,.(V1,middle_lat)], type="l")
abline(v=0, lty=3)


# the_metric <- "Greenup1"
which_igbps <- c(1:5)
igbp_list <- list(1:5, c(12, 14), 10)
igbp_names <- c("Forests", "Crops/mosaic", "Grasslands")
col_list <- c(brewer.pal(5, "Greens")[4], brewer.pal(5, "Oranges")[4], brewer.pal(5, "Purples")[4])
lat_intervals <- seq(30, 70, by=2)

plot(NA, NA, type="n", xlab=paste(the_metric, "Anomaly (days)"), ylab="Latitude", ylim=c(35, 70), xlim=c(-21, 21))
abline(v=0, lty=2)

i <- 1
for(which_igbps in igbp_list){
    DT_lat_all <- DTeurope[
        igbp %in% which_igbps & phenometric == the_metric & year %in% 2001:2016,
        median(value - as.integer(as.Date(paste(year, "-1-1", sep="")) - as.Date("1970-1-1")), na.rm=T),
        by=findInterval(lat, lat_intervals, all.inside=T)]
    DT_lat_all[, middle_lat := seq(30, 70, by=2)[findInterval] + 1]
    DT_lat_all <- DT_lat_all[order(middle_lat)]

    DT_lat_2018 <- DTeurope[
        igbp %in% which_igbps & phenometric == the_metric & year == 2018,
        median(value - as.integer(as.Date(paste(year, "-1-1", sep="")) - as.Date("1970-1-1")), na.rm=T),
        by=findInterval(lat, lat_intervals, all.inside=T)]
    DT_lat_2018[, middle_lat := seq(30, 70, by=2)[findInterval] + 1]
    DT_lat_2018 <- DT_lat_2018[order(middle_lat)]

    # plot the anomalies as DOYS
    points(DT_lat_2018[, V1] - DT_lat_all[, V1], DT_lat_all[, middle_lat], type="l", lwd=2, col=col_list[i])
    i <- i + 1
}
legend("topright", legend=igbp_names, col=col_list, lwd=2)

the_metric <- "EVI_Amplitude1"
which_igbps <- c(1:5)
igbp_list <- list(1:5, c(12, 14), 10)
igbp_names <- c("Forests", "Crops/mosaic", "Grasslands")
col_list <- c(brewer.pal(5, "Greens")[4], brewer.pal(5, "Oranges")[4], brewer.pal(5, "Purples")[4])
lat_intervals <- seq(30, 70, by=2)

plot(NA, NA, type="n", xlab=paste(the_metric, "Anomaly"), ylab="Latitude", ylim=c(35, 70), xlim=c(-21, 21))
abline(v=0, lty=2)

i <- 1
for(which_igbps in igbp_list){
    DT_lat_all <- DTeurope[
        igbp %in% which_igbps & phenometric == the_metric & year %in% 2001:2016,
        median(value, na.rm=T),
        by=findInterval(lat, lat_intervals, all.inside=T)]
    DT_lat_all[, middle_lat := seq(30, 70, by=2)[findInterval] + 1]
    DT_lat_all <- DT_lat_all[order(middle_lat)]

    DT_lat_2018 <- DTeurope[
        igbp %in% which_igbps & phenometric == the_metric & year == 2018,
        median(value, na.rm=T),
        by=findInterval(lat, lat_intervals, all.inside=T)]
    DT_lat_2018[, middle_lat := seq(30, 70, by=2)[findInterval] + 1]
    DT_lat_2018 <- DT_lat_2018[order(middle_lat)]

    # plot the anomalies as DOYS
    points(DT_lat_2018[, V1] - DT_lat_all[, V1], DT_lat_all[, middle_lat], type="l", lwd=2, col=col_list[i])
    i <- i + 1
}
legend("topright", legend=igbp_names, col=col_list, lwd=2)


MakeTransp <- function(thecolor, alpha=0.5) rgb(matrix(c(col2rgb(thecolor)), ncol=3), max=255, alpha=alpha*255)

GetLatAnomAndPlot <- function(metric, igbps=NULL, anomaly_year=2018, other_years=NULL, lat_step=2, conv_doy=TRUE, plot=TRUE, plot_col="#3182BD"){
    lat_intervals <- seq(30, 70, by=lat_step)
    if(is.null(igbps)) igbps <- 1:14
    if(is.null(other_years)) other_years <- 2001:2018

    if(conv_doy){
        DTsub <- DTeurope[
            igbp %in% igbps & phenometric == metric, 
            .(year, lat, igbp, value - as.integer(as.Date(paste(year, "-1-1", sep="")) - as.Date("1970-1-1")))
        ]
    }else{
        DTsub <- DTeurope[
            igbp %in% igbps & phenometric == metric, 
            .(year, lat, igbp, value)
        ]
    }
    setnames(DTsub, c("year", "lat", "igbp", "value"))

    # get latitudinal band averages and s.d. across all years
    DT_lat_all <- DTsub[
        year %in% other_years,
        .(mean(value, na.rm=T), sd(value, na.rm=T)),
        by=findInterval(lat, lat_intervals, all.inside=T)
    ]
    setnames(DT_lat_all, c("int", "pheno_mean", "pheno_sd"))
    DT_lat_all[, middle_lat := seq(30, 70, by=2)[int] + 1]
    DT_lat_all <- DT_lat_all[order(middle_lat)]

    # get latitudinal band averages and s.d. in anomaly year
    DT_lat_anom <- DTsub[
        year %in% anomaly_year,
        .(mean(value, na.rm=T), sd(value, na.rm=T)),
        by=findInterval(lat, lat_intervals, all.inside=T)
    ]
    setnames(DT_lat_anom, c("int", "pheno_mean", "pheno_sd"))
    DT_lat_anom[, middle_lat := seq(30, 70, by=2)[int] + 1]
    DT_lat_anom <- DT_lat_anom[order(middle_lat)]
    
    # plot the anomalies as DOYS

    # points(DT_lat_2018[, V1] - DT_lat_all[, V1], DT_lat_all[, middle_lat], type="l", lwd=2, col=col_list[i])

    


}

DT <- data.table(x=rnorm(10), y=rnorm(10), z=sample(LETTERS[1:3], 10, rep=T))
DT[, mean_x:=mean(x), by=z]
brewer.pal(5, "Blues")[4]