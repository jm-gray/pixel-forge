.libPaths("/home/jmgray2/R/x86_64-unknown-linux-gnu-library/3.1")
library(raster)
library(RColorBrewer)
library(argparse)

PlotINCASummary <- function(input_file, metric_name="INCA", LWMASKDIR="/share/jmgray2/MODIS/LWMASK500"){
  # plotting for INCA output
  MAXPIXELS <- 2.5e6
  PVALUETHRESH <- 0.05
  MEDIANCUTOFFS <- NULL
  MEDIANCUTOFFSQUANTS <- c(0.02, 0.98)
  SLOPECUTOFFS <- NULL
  SLOPECUTOFFSQUANTS <- c(0.05, 0.95)
  MADCUTOFFS <- c(0, 21)
  MADCUTOFFSQUANTS <- c(0, 0.98)
  WATERCOLOR <- rgb(0.7, 0.7, 0.7)
  LANDCOLOR <- rgb(0.2, 0.2, 0.2)
  LEGENDAXISCEX <- 1
  LEGENDMAINCEX <- 1
  LEGENDWIDTH <- 3.5

  # the color palettes
  div_pal <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))
  median_pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  mad_pal <- colorRampPalette(brewer.pal(9, "PuBuGn"))

  ###################################
  # input_file <- "~/Desktop/INCA_summary_h11v04_halfspring.tif"
  tile <- gsub(".*(h[0-9][0-9]v[0-9][0-9]).*", "\\1", basename(input_file))
  lwmask_file <- dir(LWMASKDIR, pattern=paste(".*", tile, ".*", "bin$", sep=""), full=T)

  #DEBUG
  print(input_file)
  print(lwmask_file)

  s <- stack(input_file)
  lwmask <- raster(lwmask_file)
  r_med <- raster(s, 1)
  r_med[lwmask != 1] <- NA
  r_mad <- raster(s, 2)
  r_mad[lwmask != 1] <- NA
  r_slope <- raster(s, 5)
  r_slope[lwmask != 1] <- NA
  r_pvalue <- raster(s, 6)
  r_pvalue[lwmask != 1] <- NA
  r_slope[r_pvalue > PVALUETHRESH] <- NA
  r_land <- lwmask
  r_land[lwmask != 1] <- NA

  ###################################

  ###################################
  # set up breaks for median plotting
  if(is.null(MEDIANCUTOFFS)){
    qs_med <- quantile(r_med, c(0, MEDIANCUTOFFSQUANTS[1], MEDIANCUTOFFSQUANTS[2], 1))
  }else{
    qs_med <- c(minValue(r_med), MEDIANCUTOFFS[1], MEDIANCUTOFFS[2], maxValue(r_med))
  }
  median_breaks <- c(qs_med[1], seq(qs_med[2], qs_med[3], len=254), qs_med[4])
  median_breaks <- unique(round(median_breaks))

  # set up breaks for MAD plotting
  if(is.null(MADCUTOFFS)){
    qs_mad <- quantile(r_mad, c(0, MADCUTOFFSQUANTS[1], MADCUTOFFSQUANTS[2], 1))
  }else{
    qs_mad <- c(minValue(r_mad), MADCUTOFFS[1], MADCUTOFFS[2], maxValue(r_mad))
  }
  mad_breaks <- c(qs_mad[1], seq(qs_mad[2], qs_mad[3], len=254), qs_mad[4])
  mad_breaks <- unique(round(mad_breaks))

  # set up breaks for slope plotting
  if(is.null(SLOPECUTOFFS)){
    qs_slope <- quantile(r_slope, c(0, SLOPECUTOFFSQUANTS[1], SLOPECUTOFFSQUANTS[2], 1))
  }else{
    qs_slope <- c(minValue(r_slope), SLOPECUTOFFS[1], SLOPECUTOFFS[2], maxValue(r_slope))
  }
  # the following is necessary to ensure that the number of breaks are even on
  # the positive/negative side, thus ensuring 0 is the middle color
  slope_extreme <- max(abs(qs_slope[1]), qs_slope[4])
  slope_start <- max(abs(qs_slope[2]), qs_slope[3])
  slope_breaks <- c(-1 * slope_extreme, seq(-1 * slope_start, slope_start, len=254), slope_extreme)

  ###################################
  # setup the plot
  # quartz(h=7, w=21.5)
  # if(!is.null(out_pdf)){
  #   pdf(file=out_pdf, h=7, w=24)
  # }else{
  #   quartz(h=7, w=24)
  # }

  layout(matrix(1:3, nrow=1))
  par(mar=c(1, 1, 1, 3), oma=rep(1, 4))

  ###################################
  # plot the median DOY map
  plot(lwmask, breaks=c(-1, 0.5, 1.5, 10), col=c(WATERCOLOR, LANDCOLOR, WATERCOLOR), xaxt="n", yaxt="n", legend=F, bty="n", box=FALSE)
  plot(r_med, breaks=median_breaks, col=median_pal(length(median_breaks) - 1), maxpixels=MAXPIXELS, legend=F, xaxt="n", yaxt="n", bty="n", box=FALSE, add=T)
  legend_at <- round(seq(median_breaks[2], median_breaks[length(median_breaks) - 1], len=7))
  legend_at_date <- legend_at
  legend_labels <- c(paste("<", legend_at_date[1]), as.character(legend_at_date[2:(length(legend_at_date) - 1)]), paste(">", legend_at_date[length(legend_at_date)]))
  plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=median_pal(length(median_breaks)-1), legend.width=LEGENDWIDTH, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=LEGENDAXISCEX), legend.args=list(text="", side=3, font=2, line=0.5, cex=LEGENDMAINCEX))
  # title("Medsian")
  title(paste(tile, metric_name, "Median"))

  # plot the MAD map
  plot(lwmask, breaks=c(-1, 0.5, 1.5, 10), col=c(WATERCOLOR, LANDCOLOR, WATERCOLOR), xaxt="n", yaxt="n", legend=F, bty="n", box=FALSE)
  plot(r_mad, breaks=mad_breaks, col=mad_pal(length(mad_breaks) - 1), maxpixels=MAXPIXELS, legend=F, xaxt="n", yaxt="n", bty="n", box=FALSE, add=T)
  legend_at <- round(seq(mad_breaks[2], mad_breaks[length(mad_breaks) - 1], len=7))
  legend_at_date <- legend_at
  legend_labels <- c(paste("<", legend_at_date[1]), as.character(legend_at_date[2:(length(legend_at_date) - 1)]), paste(">", legend_at_date[length(legend_at_date)]))
  plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=mad_pal(length(mad_breaks)-1), legend.width=LEGENDWIDTH, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=LEGENDAXISCEX), legend.args=list(text="", side=3, font=2, line=0.5, cex=LEGENDMAINCEX))
  # title("Median Abs. Dev.")
  title(paste(tile, metric_name, "Med. Abs. Dev."))

  # plot the slope map
  plot(lwmask, breaks=c(-1, 0.5, 1.5, 10), col=c(WATERCOLOR, LANDCOLOR, WATERCOLOR), xaxt="n", yaxt="n", legend=F, bty="n", box=FALSE)
  # plot(r_land, col=SLOPEBGCOLOR, xaxt="n", yaxt="n", legend=F, bty="n", box=FALSE)
  plot(r_slope, breaks=slope_breaks, col=div_pal(length(slope_breaks) - 1), maxpixels=MAXPIXELS, legend=F, add=T, xaxt="n", yaxt="n", bty="n", box=FALSE)
  legend_at <- round(seq(slope_breaks[2], slope_breaks[length(slope_breaks) - 1], len=7))
  legend_at_date <- legend_at
  legend_labels <- c(paste("<", legend_at_date[1]), as.character(legend_at_date[2:(length(legend_at_date) - 1)]), paste(">", legend_at_date[length(legend_at_date)]))
  plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=div_pal(length(slope_breaks)-1), legend.width=LEGENDWIDTH, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=LEGENDAXISCEX), legend.args=list(text="", side=3, font=2, line=0.5, cex=LEGENDMAINCEX))
  title(paste(tile, metric_name, "Theil-Sen Slope (where p <=", PVALUETHRESH, ")"))

  # plot the annual median anomalies as multiples of the MAD
  MADANOMCUTOFFS <- c(-3, 3)
  mad_anoms <- subset(s, 21:34)
  madanom_extreme <- max(c(max(maxValue(mad_anoms)), abs(min(minValue(mad_anoms)))))
  madanom_breaks <- c(-1*madanom_extreme, seq(MADANOMCUTOFFS[1], MADANOMCUTOFFS[2], len=254), madanom_extreme)
  years <- 2001:2014

  # quartz(h=12, w=12)
  # layout(matrix(1:16, nrow=4, byrow=TRUE))
  layout(matrix(1:14, nrow=2, byrow=TRUE))
  par(mar=c(0.25, 0.25, 0.75, 0.25), oma=c(1,2,1,3.5))
  for(i in 1:nlayers(mad_anoms)){
    plot(lwmask, breaks=c(-1, 0.5, 1.5, 10), col=c(WATERCOLOR, LANDCOLOR, WATERCOLOR), xaxt="n", yaxt="n", legend=F, bty="n", box=FALSE)
    plot(raster(mad_anoms, i), breaks=madanom_breaks, col=div_pal(length(madanom_breaks) - 1), legend=FALSE, xaxt="n", yaxt="n", bty="n", box=FALSE, add=T)
    title(years[i])
    if(i == 1 | i == 8){
      mtext(paste(tile, metric_name), side=2, line=0.5)
    }
    if(i == 14){
      legend_at <- round(seq(madanom_breaks[2], madanom_breaks[length(madanom_breaks) - 1], len=7))
      legend_at_date <- legend_at
      legend_labels <- c(paste("<", legend_at_date[1]), as.character(legend_at_date[2:(length(legend_at_date) - 1)]), paste(">", legend_at_date[length(legend_at_date)]))
      plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=div_pal(length(madanom_breaks)-1), legend.width=1, axis.args=list(at=legend_at, labels=legend_labels, cex.axis=LEGENDAXISCEX), legend.args=list(text="", side=3, font=2, line=0.5, cex=LEGENDMAINCEX))
    }
  }
  # if(!is.null(out_pdf)) dev.off()
}

# Example:
# PlotINCASummary("~/Desktop/INCA_summary_h10v04_halfspring.tif", out_pdf="~/Desktop/test.pdf", LWMASKDIR="~/Desktop", metric_name="HalfSpring")

#===========================================================================
# parse the command line arguments
arg_parser <- ArgumentParser()
arg_parser$add_argument("-tile", type="character") # tile to process
arg_parser$add_argument("-output_dir", type="character", default="/share/jmgray2/INCA") # output directory
arg_parser$add_argument("-data_dir", type="character", default="/share/jmgray2/INCA/output") # input binary splined evi data directory
args <- arg_parser$parse_args()

metrics <- c("ogi", "halfspring", "halffall", "dormancy", "gsl")
# metrics <- c("halfspring")
out_pdf <- file.path(args$out_dir, paste("INCA_plot_", args$tile, ".pdf", sep=""))
pdf(file=out_pdf, h=7, w=24)
for(metric in metrics){
  in_file <- dir(args$data_dir, pattern=paste("INCA_summary_", args$tile, "_", metric, ".tif", sep=""), full=T)
  print(in_file)
  PlotINCASummary(in_file, metric_name=metric)
}
dev.off()
