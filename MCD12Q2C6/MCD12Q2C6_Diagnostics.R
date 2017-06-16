#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
library(raster)
library(RColorBrewer)
library(rgdal)
library(parallel)
library(argparse)
library(vioplot)

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
PhenoNormals <- function(x, years=NULL){
  # function takes a vector of DOY's (x) and input years (years) and returns:
  # x's median, MAD then the rounded median anomalies of x for all years, then that anomaly divided
  # by the MAD (how many median absolute deviations from the median is the year?)

  if(is.null(years)) years <- 1:length(x) # if no years are given, assume they are in the right order and have no gaps
  med <- median(x, na.rm=T)
  MAD <- mad(x, na.rm=T)
  
  # calculate annual median anomalies: raw and as multiples of MAD
  ann.median.anoms <- round(x - med)
  ann.median.anoms.mad <- (x - med) / MAD

  return(c(round(med), MAD, ann.median.anoms, ann.median.anoms.mad))
}

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
PlotTile <- function(r, lwmask, cutoffs=c(0, 365), round_digs=0, pal=colorRampPalette(rev(brewer.pal(11, "Spectral"))), title=NA, pdf_out=NULL, plot_height=12, plotNEW=F, garish=F, continents=NA,...){
  MAXPIXELS <- 5.7e6
  WATERCOLOR <- rgb(0.7, 0.7, 0.7)
  LANDCOLOR <- rgb(0.2, 0.2, 0.2)
  LEGENDAXISCEX <- 1
  LEGENDMAINCEX <- 1
  LEGENDWIDTH <- 3.5

  breaks <- round(c(minValue(r), seq(cutoffs[1], cutoffs[2]), maxValue(r)), round_digs)
  breaks <- unique(breaks)

  plot(lwmask, breaks=c(-1, 0.5, 1.5, 10), col=c(WATERCOLOR, LANDCOLOR, WATERCOLOR), xaxt="n", yaxt="n", legend=F, bty="n", box=FALSE, maxpixels=MAXPIXELS)
  plot(r, breaks=breaks, col=pal(length(breaks) - 1), maxpixels=MAXPIXELS, legend=F, xaxt="n", yaxt="n", bty="n", box=FALSE, add=T, ...)
  legend_at <- round(seq(breaks[2], breaks[length(breaks) - 1], len=7), round_digs)
  legend_at_date <- legend_at
  legend_labels <- c(paste("<", legend_at_date[1]), as.character(legend_at_date[2:(length(legend_at_date) - 1)]), paste(">", legend_at_date[length(legend_at_date)]))
  plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=pal(length(breaks)-1), legend.width=LEGENDWIDTH, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=LEGENDAXISCEX), legend.args=list(text="", side=3, font=2, line=0.5, cex=LEGENDMAINCEX))
  title(title)
}



#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
cl <- makeCluster(8)
tile <- "h12v04"
data_dir <- "/projectnb/modislc/data/mcd12_out/phen_out/c6"
output_dir <- "/projectnb/modislc/users/joshgray/C6_Diagnostics"
lwmask_dir <- "/projectnb/modislc/data/mcd12_in/c6/ancillary_layers/C6_LW_Mask/lw_mask_500m"
doy_metrics <- c("Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy")
other_metrics <- c("NumCycles", "QA_Overall", "EVI_Amplitude", "EVI_Area", "EVI_Minimum")
years <- 2001:2014

# 2=land, 1=water
tile_lwmask <- raster(dir(lwmask_dir, pattern=paste("map.*", tile, "$", sep=""), full=T))

for(metric in doy_metrics){
  in_files <- Sys.glob(file.path(data_dir, "*", "001", metric, paste(metric, "*", tile, "*[0-9]", sep="")))
  pheno_s <- stack(in_files)
  pheno_s <- subset(pheno_s, seq(1, nlayers(pheno_s), by=2))
  pheno_s[(pheno_s == 32767) | (pheno_s == -32768)] <- NA
  tmp_lw <- raster(pheno_s, 1)
  values(tmp_lw) <- values(tile_lwmask)
  pheno_s[tmp_lw == 1] <- NA
  doy_offset <- unlist(lapply(years, function(x) as.numeric(as.Date(paste(x, "-1-1", sep=""), format="%Y-%j") - as.Date("1970-1-1"))))
  pheno_s <- pheno_s - doy_offset
  pheno_stack_v <- values(pheno_s)
  system.time(pheno_output <- parApply(cl, pheno_stack_v, 1, PhenoNormals))

  #------------------------------------------
  # plotting
  #------------------------------------------
  #-------
  # violin plot: distribution of metric across all years
  layout(matrix(1:2, nrow=2))
  par(mar=c(2, 4, 1, 1))

  spectral_pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  vio_colors <- spectral_pal(length(years))
  plot(NA, ylim=c(-30, 365), xlim=c(0, length(years) + 1), xlab="", ylab="DOY", xaxt="n")
  for(i in 1:length(years)) vioplot(pheno_stack_v[,i][!is.na(pheno_stack_v[,i])], add=T, col=vio_colors[i], at=i)
  axis(1, at=1:length(years), labels=years)
  title=paste("Annual tile-level distributions:", metric)

  # calculate missing fraction among land values
  lw_v <- values(tmp_lw)
  land_missing_frac <- apply(pheno_stack_v, 2, function(x, lw_vals=lw_v) {x <- x[lw_vals == 2]; sum(is.na(x))}) / sum(lw_v == 2)
  plot(1:length(years), land_missing_frac, xlim=c(0, length(years) + 1), xlab="", ylab="Land NA Fraction", xaxt="n", type="h", lwd=6, col=spectral_pal(5)[1])
  axis(1, at=1:length(years), labels=years)
  title=paste("Annual tile-level missing (land-only):", metric)

  #-------
  # plot Median and MAD for all years
  median_pal <- colorRampPalette(brewer.pal(9, "PuBuGn"))
  tmp_r <- raster(pheno_s)
  median_r <- tmp_r
  values(median_r) <- t(pheno_output[1,])
  mad_r <- tmp_r
  values(mad_r) <- t(pheno_output[2,])

  # layout(matrix(1:2, nrow=1))
  par(mfrow=c(1, 2), mar=c(1, 1, 1, 3))

  # plot the median raster
  qs <- quantile(median_r, c(0.02, 0.98))
  PlotTile(median_r, tmp_lw, cutoffs=qs, pal=spectral_pal, title=paste("Median", metric))

  qs <- quantile(mad_r, c(0.02, 0.98))
  PlotTile(mad_r, tmp_lw, cutoffs=qs, pal=median_pal, round_digs=1, title=paste("MAD", metric))

  #-------
  # plot annual anomalies and QA/QC values
  par(mfrow=c(3, 5), mar=c(1, 1, 1, 3))
  year <- years[1]
  for(i in 3:(3 + length(years))){
    values(tmp_r) <- t(pheno_output[i,])
    PlotTile(mad_r, tmp_lw, cutoffs=qs, pal=median_pal, round_digs=1, title=year)
    year <- year + 1
  }


}
