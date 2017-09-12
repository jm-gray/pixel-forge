#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
library(raster)
library(RColorBrewer)
library(rgdal)
library(parallel)
library(argparse)
# library(vioplot)

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
}

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
# /projectnb/landsat/users/dsm/eval_modis_lc_061917/MCD12I6
PlotTileAllMetrics <- function(tile, cl, years=2001:2008, doy_metrics=c("Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy"), data_dir="/projectnb/landsat/users/dsm/eval_modis_lc_061917/MCD12I6", output_dir="/projectnb/modislc/users/joshgray/C6_Diagnostics", lwmask_dir="/projectnb/modislc/data/mcd12_in/c6/ancillary_layers/C6_LW_Mask/lw_mask_500m"){
  # 2=land, 1=water
  print(paste("Doing", tile))
  tile_lwmask <- raster(dir(lwmask_dir, pattern=paste("map.*", tile, "$", sep=""), full=T))
  # overall_qa_files <- Sys.glob(file.path(data_dir, "*", "001", "QA_Overall", paste("QA_Overall", "*", tile, "*[0-9]", sep="")))
  # example_geom_r <- raster(Sys.glob(file.path(data_dir, "2001", "001", "Greenup", paste("Greenup", "*", tile, "*[0-9]", sep=""))))
  overall_qa_files <- Sys.glob(file.path(data_dir, "QA_Overall", paste("QA_Overall", "*", tile, "*[0-9]", sep="")))
  example_geom_r <- raster(Sys.glob(file.path(data_dir, "Greenup", paste("Greenup", "*", tile, "*[0-9]", sep="")))[1])
  projection(tile_lwmask) <- projection(example_geom_r)
  extent(tile_lwmask) <- extent(example_geom_r)

  # make the output PDF
  pdf(file=file.path(output_dir, paste("C6_Diagnostic_", tile, ".pdf", sep="")), width=12, height=7) # NOTE: change this!

  # make plot of overall QA
  layout(matrix(1:15, nrow=3, byrow=T))
  par(mar=rep(1,4))
  i <- 1
  for(overall_qa_file in overall_qa_files){
    r <- raster(overall_qa_file)
    r[tile_lwmask == 1] <- NA
    r[(r == 32767) | (r == -32768)] <- NA
    qa_pal <- colorRampPalette(brewer.pal(5, "PuRd")[2:5])
    PlotTile(r, tile_lwmask, breaks=seq(-0.5, 3.5, by=1), pal=qa_pal, title=paste("Overall QA", years[i]), plot_legend=F)
    i <- i + 1
  }
  plot(1:10, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  legend("top", legend=c("0","1","2","3"), fill=qa_pal(4), bg="white", horiz=T, title="Overall QA Value (first cycle)")

  # Make a plots of NumCycles
  layout(matrix(1:15, nrow=3, byrow=T))
  par(mar=rep(1,4))
  # num_cycles_files <- Sys.glob(file.path(data_dir, "*", "001", "NumCycles", paste("NumCycles", "*", tile, "*[0-9]", sep="")))
  num_cycles_files <- Sys.glob(file.path(data_dir, "NumCycles", paste("NumCycles", "*", tile, "*[0-9]", sep="")))
  i <- 1
  for(num_cycles_file in num_cycles_files){
    r <- raster(num_cycles_file)
    r[tile_lwmask == 1] <- NA
    r[(r == 32767) | (r == -32768)] <- NA
    num_cycles_pal <- colorRampPalette(brewer.pal(5, "YlGnBu"))
    PlotTile(r, tile_lwmask, breaks=c(0.5, 1.5, 2.5, 3.5, 100), pal=num_cycles_pal, title=paste("Num Cycles", years[i]), plot_legend=F)
    i <- i + 1
  }
  plot(1:10, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  legend("top", legend=c("1","2","3",">3"), fill=num_cycles_pal(4), bg="white", horiz=T, title="Num Veg Cycles")

  # Make a plots of EVI_Amplitude
  layout(matrix(1:15, nrow=3, byrow=T))
  par(mar=rep(1,4))
  # evi_amp_files <- Sys.glob(file.path(data_dir, "*", "001", "EVI_Amplitude", paste("EVI_Amplitude", "*", tile, "*[0-9]", sep="")))
  evi_amp_files <- Sys.glob(file.path(data_dir, "EVI_Amplitude", paste("EVI_Amplitude", "*", tile, "*[0-9]", sep="")))
  evi_amp_s <- stack(evi_amp_files)
  evi_amp_s <- subset(evi_amp_s, seq(1, nlayers(evi_amp_s), by=2))
  evi_amp_s[tile_lwmask == 1] <- NA
  evi_amp_s[(evi_amp_s == 32767) | (evi_amp_s == -32768)] <- NA
  qs <- quantile(evi_amp_s, c(0.02, 0.98), na.rm=T)
  cutoffs <- c(min(qs[,1]), max(qs[,2]))
  for(i in 1:nlayers(evi_amp_s)){
    r <- raster(evi_amp_s, i)
    evi_amp_pal <- colorRampPalette(brewer.pal(5, "Greens"))
    PlotTile(r, tile_lwmask, cutoffs=cutoffs, pal=evi_amp_pal, title=paste("EVI_Amplitude", years[i]), plot_legend=F)
  }
  plot(1:10, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  PlotTile(r, tile_lwmask, cutoffs=cutoffs, pal=evi_amp_pal, title="", legend_only=T, round=3, scale=1e4, legend_title="EVI Amplitude")

  # Make a plots of EVI_Area
  layout(matrix(1:15, nrow=3, byrow=T))
  par(mar=rep(1,4))
  # evi_area_files <- Sys.glob(file.path(data_dir, "*", "001", "EVI_Area", paste("EVI_Area", "*", tile, "*[0-9]", sep="")))
  evi_area_files <- Sys.glob(file.path(data_dir, "EVI_Area", paste("EVI_Area", "*", tile, "*[0-9]", sep="")))
  evi_area_s <- stack(evi_area_files)
  evi_area_s <- subset(evi_area_s, seq(1, nlayers(evi_area_s), by=2))
  evi_area_s[tile_lwmask == 1] <- NA
  evi_area_s[(evi_area_s == 32767) | (evi_area_s == -32768)] <- NA
  qs <- quantile(evi_area_s, c(0.02, 0.98), na.rm=T)
  cutoffs <- c(min(qs[,1]), max(qs[,2]))
  for(i in 1:nlayers(evi_area_s)){
    r <- raster(evi_area_s, i)
    evi_area_pal <- colorRampPalette(brewer.pal(5, "Blues"))
    PlotTile(r, tile_lwmask, cutoffs=cutoffs, pal=evi_area_pal, title=paste("EVI_Area", years[i]), plot_legend=F)
  }
  plot(1:10, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  PlotTile(r, tile_lwmask, cutoffs=cutoffs, pal=evi_area_pal, title="", legend_only=T, round=3, scale=10, legend_title="EVI Area")

  # Make a plots of EVI_Minimum
  layout(matrix(1:15, nrow=3, byrow=T))
  par(mar=rep(1,4))
  # evi_min_files <- Sys.glob(file.path(data_dir, "*", "001", "EVI_Minimum", paste("EVI_Minimum", "*", tile, "*[0-9]", sep="")))
  evi_min_files <- Sys.glob(file.path(data_dir, "EVI_Minimum", paste("EVI_Minimum", "*", tile, "*[0-9]", sep="")))
  evi_min_s <- stack(evi_min_files)
  evi_min_s <- subset(evi_min_s, seq(1, nlayers(evi_min_s), by=2))
  evi_min_s[tile_lwmask == 1] <- NA
  evi_min_s[(evi_min_s == 32767) | (evi_min_s == -32768)] <- NA
  qs <- quantile(evi_min_s, c(0.02, 0.98), na.rm=T)
  cutoffs <- c(min(qs[,1]), max(qs[,2]))
  for(i in 1:nlayers(evi_min_s)){
    r <- raster(evi_min_s, i)
    evi_min_pal <- colorRampPalette(brewer.pal(5, "Reds"))
    PlotTile(r, tile_lwmask, cutoffs=cutoffs, pal=evi_min_pal, title=paste("EVI_Minimum", years[i]), plot_legend=F)
  }
  plot(1:10, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  PlotTile(r, tile_lwmask, cutoffs=cutoffs, pal=evi_min_pal, title="", legend_only=T, round=3, scale=1e4, legend_title="EVI Minimum")

  # loop through each DOY metrics and plot: violin plot of each year's values, missing fraction for each year, median raster, MAD raster,
  for(metric in doy_metrics){
    # in_files <- Sys.glob(file.path(data_dir, "*", "001", metric, paste(metric, "*", tile, "*[0-9]", sep="")))
    in_files <- Sys.glob(file.path(data_dir, metric, paste(metric, "*", tile, "*[0-9]", sep="")))
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
    # boxplots: distribution of metric across all years
    layout(matrix(1:2, nrow=2))
    par(mar=c(2, 4, 3, 1))

    spectral_pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
    vio_colors <- spectral_pal(length(years))
    # get the proper ylim
    for(i in 1:length(years)){
      qs <- quantile(pheno_stack_v[, i], c(0.25, 0.75), na.rm=T)
      iqr <- qs[2] - qs[1]
      if(i == 1){
        ylim <- c(qs[1] - (1.5 * iqr), qs[2] + (1.5 * iqr))
      }else{
        ylim[1] <- min(ylim[1], qs[1] - (1.55 * iqr))
        ylim[2] <- max(ylim[2], qs[2] + (1.55 * iqr))
      }
    }
    plot(NA, ylim=ylim, xlim=c(0, length(years) + 1), xlab="", ylab="DOY", xaxt="n", bty="n")
    # for(i in 1:length(years)) vioplot(pheno_stack_v[,i][!is.na(pheno_stack_v[,i])], add=T, col=vio_colors[i], at=i)
    for(i in 1:length(years)) boxplot(pheno_stack_v[,i][!is.na(pheno_stack_v[,i])], add=T, col=vio_colors[i], at=i, outline=F, yaxt="n")
    axis(1, at=1:length(years), labels=years)
    title(paste("Annual tile-level distributions:", metric))

    # calculate missing fraction among land values
    lw_v <- values(tmp_lw)
    land_missing_frac <- apply(pheno_stack_v, 2, function(x, lw_vals=lw_v) {x <- x[lw_vals == 2]; sum(is.na(x))}) / sum(lw_v == 2)
    plot(1:length(years), land_missing_frac, xlim=c(0, length(years) + 1), xlab="", ylab="Land NA Fraction", xaxt="n", type="h", lwd=6, col=spectral_pal(5)[1])
    axis(1, at=1:length(years), labels=years)
    title(paste("Annual tile-level missing (land-only):", metric))

    #-------
    # plot Median and MAD for all years
    median_pal <- colorRampPalette(brewer.pal(9, "PuBuGn"))
    tmp_r <- raster(pheno_s)
    median_r <- tmp_r
    values(median_r) <- t(pheno_output[1,])
    mad_r <- tmp_r
    values(mad_r) <- t(pheno_output[2,])

    layout(matrix(1:2, nrow=1))
    # par(mfrow=c(1, 2), mar=c(1, 1, 1, 3))
    # layout(matrix(1:3, nrow=1))
    par(mar=c(1, 1, 1, 3), oma=rep(1, 4))

    # plot the median raster
    qs <- quantile(median_r, c(0.02, 0.98))
    PlotTile(median_r, tmp_lw, cutoffs=qs, pal=spectral_pal, title=paste("Median", metric))

    qs <- quantile(mad_r, c(0.02, 0.98))
    PlotTile(mad_r, tmp_lw, cutoffs=qs, pal=median_pal, round_digs=1, title=paste("MAD", metric))

    # #-------
    median_anom_breaks <- c(-1 * 365, seq(-21, 21, len=254), 365)

    layout(matrix(1:15, nrow=3, byrow=T))
    par(mar=rep(1,4))
    for(i in 1:nlayers(pheno_s)){
      values(tmp_r) <- t(pheno_output[2 + i,])
      median_anom_pal <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))
      PlotTile(tmp_r, tile_lwmask, breaks=median_anom_breaks, pal=median_anom_pal, title=paste(metric, "DOY Anomaly", years[i]), plot_legend=F)
    }
    plot(1:10, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
    PlotTile(r, tile_lwmask, breaks=median_anom_breaks, pal=median_anom_pal, title="", legend_only=T, legend_title="DOY Anomaly")

  }
  dev.off()
}

#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
arg_parser <- ArgumentParser()
arg_parser$add_argument("-tile", type="character") # tile to process
arg_parser$add_argument("-cores", type="integer", default=8) # tile to process
args <- arg_parser$parse_args()

cl <- makeCluster(args$cores) # NOTE: change to command line arg
doy_metrics <- c("Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy")
PlotTileAllMetrics(tile=args$tile, cl=cl, doy_metrics=doy_metrics)


# # submission shell script:
# #!/bin/bash
#
# echo Submitting tile $1
# R --vanilla < /projectnb/modislc/users/joshgray/C6_Diagnostics/MCD12Q2C6_Diagnostics.R --args -tile $1

# to submit all tiles:
tiles <- scan("/projectnb/modislc/users/joshgray/MCD12Q2C6/gltiles.txt", what=character(), quiet=T)
for(tile in tiles){
  sys_cmd <- paste("qsub -V -l h_rt=02:00:00 -pe omp 8 /projectnb/modislc/users/joshgray/C6_Diagnostics/run_diagnostics.sh", tile)
  system(sys_cmd)
}
