library(raster)
library(RColorBrewer)

setwd("/projectnb/modislc/users/joshgray/MCD12Q2C6/CLM_output/CLM_extract")

#================================================================
# Annual mosaic for particular metric
#----------------------------------------------------------------
# make a list of RasterLayers for a particular year
year <- 2012
# phenophase <- "twentygup"
phenophase <- "fiftygup"
rasters <- lapply(Sys.glob(paste("*", phenophase, year, "*", sep="*")), raster)

# merge them into a single raster layer
m <- do.call(mosaic, c(rasters, fun="mean"))
m_5k <- aggregate(m, fact=10, fun=median)

# get some colors for plotting
pal <- colorRampPalette(brewer.pal(11, "Spectral"))

# make some breaks
jan1 <- as.numeric(as.Date(paste(year, "-1-1", sep="")) - as.Date("1970-1-1"))
mar1 <- as.numeric(as.Date(paste(year, "-3-1", sep="")) - as.Date("1970-1-1"))
apr1 <- as.numeric(as.Date(paste(year, "-4-1", sep="")) - as.Date("1970-1-1"))
may1 <- as.numeric(as.Date(paste(year, "-5-1", sep="")) - as.Date("1970-1-1"))
jun1 <- as.numeric(as.Date(paste(year, "-6-1", sep="")) - as.Date("1970-1-1"))
aug1 <- as.numeric(as.Date(paste(year, "-8-1", sep="")) - as.Date("1970-1-1"))
dec31 <- as.numeric(as.Date(paste(year, "-12-31", sep="")) - as.Date("1970-1-1"))
# breaks <- c(jan1, seq(mar1, jun1), dec31)
breaks <- c(jan1, seq(apr1, jun1), dec31) # good for N America 50% GUP
# breaks <- c(jan1, seq(may1, aug1), dec31)

# plot
num_labs <- 7 # number of labels on legend
par(mar=rep(1, 4))
# plot(m, col=pal(length(breaks) - 1), breaks=breaks, xaxt="n", yaxt="n", legend=F)
plot(m, col=pal(length(breaks) - 1), breaks=breaks, xaxt="n", yaxt="n", legend=F, maxpixels=5e6)

# add a reasonable legend
legend_at <- round(seq(breaks[2], breaks[length(breaks) - 1], len=7))
legend_at_date <- as.Date(legend_at, origin="1970-1-1")
legend_labels <- c(
					paste("<", legend_at_date[1]),
					as.character(legend_at_date[2:(length(legend_at_date) - 1)]),
					paste(">", legend_at_date[length(legend_at_date)])
				)

# pdf(h=14, w=18, file="/projectnb/modislc/users/joshgray/NA_mcd12q2_c5_2013.pdf")
# par(mar=rep(0.75, 4), oma=c(1, 1, 1, 5))
# plot(na_r, breaks=breaks, col=pal(length(breaks)-1), maxpixels=2.3e7, legend=F)
plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=pal(length(breaks)-1), axis.args=list(at=legend_at, labels=legend_labels))


# add a title
title(paste(phenophase, year, sep=":"))


#================================================================
# C6 NBAR input 2013/2014
#----------------------------------------------------------------
setwd("/projectnb/modislc/users/joshgray/MCD12Q2C6")

gup <- raster("extracted_nbarC6/twentygup_h12v04_2013.tif")
midgup <- raster("extracted_nbarC6/midgup_h12v04_2013.tif")
peak <- raster("extracted_nbarC6/peak_h12v04_2013.tif")
gup_rsq <- raster("extracted_nbarC6/gup_rsquared_h12v04_2013.tif")
gup_missing <- raster("extracted_nbarC6/gup_missing_h12v04_2013.tif")

pal <- colorRampPalette(brewer.pal(11, "Spectral"))
year <- 2013

#----------------------
# plot GUP
breaks <- DateBreaks(year, bottom_col=paste(year, "-4-1", sep=""), top_col=paste(year, "-6-1", sep=""))
PlotDateLegend(gup, breaks, pal, maxpixels=2e6)

# gup_coarse <- aggregate(gup, factor=8, fun=median)
# PlotDateLegend(gup_coarse, breaks, pal, maxpixels=2e6)

#----------------------
# plot MIDGUP
breaks <- DateBreaks(year, bottom_col=paste(year, "-5-1", sep=""), top_col=paste(year, "-6-15", sep=""))
PlotDateLegend(midgup, breaks, pal, maxpixels=2e6)

#----------------------
# plot peak
breaks <- DateBreaks(year, bottom_col=paste(year, "-6-7", sep=""), top_col=paste(year, "-8-24", sep=""))
PlotDateLegend(peak, breaks, pal, maxpixels=2e6)


#---------------------------------------------------------------------
PlotDateLegend <- function(r, breaks, pal, ...){
	plot(r, col=pal(length(breaks) - 1), breaks=breaks, xaxt="n", yaxt="n", legend=F, ...)

	# add a reasonable legend
	legend_at <- round(seq(breaks[2], breaks[length(breaks) - 1], len=7))
	legend_at_date <- as.Date(legend_at, origin="1970-1-1")
	legend_labels <- c(paste("<", legend_at_date[1]), as.character(legend_at_date[2:(length(legend_at_date) - 1)]), paste(">", legend_at_date[length(legend_at_date)]))
	plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=pal(length(breaks)-1), axis.args=list(at=legend_at, labels=legend_labels))
}

#
PlotDoyReadableLegend <- function(r, breaks, pal, ...){
	# par(bg="white")
	par(col.lab="white", col.axis="white", col.main="white", col.sub="white", fg="white", bg="black", mar=c(3.5, 3.5, 1, 1), oma=rep(1, 4))
	plot(r, col=pal(length(breaks) - 1), breaks=breaks, xaxt="n", yaxt="n", legend=F, ...)

	# add a reasonable legend
	legend_at <- round(seq(breaks[2], breaks[length(breaks) - 1], len=7))
	legend_at_date <- strftime(as.Date(paste("2001-",legend_at, sep=""),format="%Y-%j"),"%b-%d")
	# legend_at_date <- as.Date(legend_at, origin="1970-1-1")

	legend_labels <- c(paste("<", legend_at_date[1]), as.character(legend_at_date[2:(length(legend_at_date) - 1)]), paste(">", legend_at_date[length(legend_at_date)]))
	# plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=pal(length(breaks)-1), axis.args=list(at=legend_at, labels=legend_labels, col="white", lab.col="white"))
	plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=pal(length(breaks)-1), axis.args=list(at=legend_at, labels=legend_labels))
}


#---------------------------------------------------------------------
PlotDOYLegend <- function(r, breaks, pal, ...){
	plot(r, col=pal(length(breaks) - 1), breaks=breaks, xaxt="n", yaxt="n", legend=F, ...)

	# add a reasonable legend
	legend_at <- round(seq(breaks[2], breaks[length(breaks) - 1], len=7))
	# legend_at_date <- as.Date(legend_at, origin="1970-1-1")
	legend_labels <- c(paste("<", legend_at[1]), as.character(legend_at[2:(length(legend_at) - 1)]), paste(">", legend_at[length(legend_at)]))
	plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=pal(length(breaks)-1), axis.args=list(at=legend_at, labels=legend_labels))
}

#---------------------------------------------------------------------
PlotPheno <- function(r, breaks=NULL, lin_stretch=c(0.02, 0.98), max_cols=255, pal=colorRampPalette(brewer.pal(11, "Spectral")), ...){
	if(is.null(breaks)){
		quants <- quantile(values(r), c(0, lin_stretch, 1), na.rm=T)
		breaks <- c(quants[1], seq(quants[2], quants[3], len=max_cols - 1), quants[4])
		breaks <- unique(breaks)
	}

	plot(r, col=pal(length(breaks) - 1), breaks=breaks, xaxt="n", yaxt="n", legend=F, ...)

	# add a reasonable legend
	legend_at <- round(seq(breaks[2], breaks[length(breaks) - 1], len=7))
	# legend_at_date <- as.Date(legend_at, origin="1970-1-1")
	legend_labels <- c(paste("<", legend_at[1]), as.character(legend_at[2:(length(legend_at) - 1)]), paste(">", legend_at[length(legend_at)]))
	plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=pal(length(breaks)-1), axis.args=list(at=legend_at, labels=legend_labels))
}

#---------------------------------------------------------------------
PlotDiff <- function(r, pal=colorRampPalette(brewer.pal(9, "RdBu")), ...){
	v <- values(r)
	breaks <- c(min(v, na.rm=T), c(-1*rev(seq(3.5, 21.5, by=3)), seq(3.5, 21.5, by=3)), max(v, na.rm=T))
	breaks <- unique(round(breaks))

	plot(r, col=pal(length(breaks) - 1), breaks=breaks, xaxt="n", yaxt="n", legend=F, ...)

	# add a reasonable legend
	legend_at <- round(seq(breaks[2], breaks[length(breaks) - 1], len=7))
	# legend_at_date <- as.Date(legend_at, origin="1970-1-1")
	legend_labels <- c(paste("<", legend_at[1]), as.character(legend_at[2:(length(legend_at) - 1)]), paste(">", legend_at[length(legend_at)]))
	plot(raster(matrix(legend_at[1]:legend_at[length(legend_at)])), legend.only=T, col=pal(length(breaks)-1), axis.args=list(at=legend_at, labels=legend_labels))
}



#---------------------------------------------------------------------
DateBreaks <- function(year, bottom=paste(year, "-1-1", sep=""), top=paste(year, "-12-31", sep=""), bottom_col=paste(year, "-4-1", sep=""), top_col=paste(year, "-6-1", sep=""), origin="1970-1-1"){
	bottom <- as.numeric(as.Date(bottom) - as.Date(origin))
	top <- as.numeric(as.Date(top) - as.Date(origin))
	bottom_col <- as.numeric(as.Date(bottom_col) - as.Date(origin))
	top_col <- as.numeric(as.Date(top_col) - as.Date(origin))
	breaks <- c(bottom, seq(bottom_col, top_col), top)
	return(breaks)
}
