library(raster)
library(RColorBrewer)
library(parallel)


GetEVI2 <- function(mcd43a4_file_path){
	SDS_red <- paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":MOD_Grid_BRDF:Nadir_Reflectance_Band1", sep = "")
	r_red <- raster(SDS_red) # read in the red SDS as a raster
	SDS_nir <- paste("HDF4_EOS:EOS_GRID:\"", mcd43a4_file_path, "\":MOD_Grid_BRDF:Nadir_Reflectance_Band2", sep = "")
	r_nir <- raster(SDS_nir) # read in the red SDS as a raster
	EVI2 <- 2.5 * ((r_nir - r_red) / (r_nir + (2.4 * r_red) + 1)) # calculate EVI2
	return(EVI2)
}

DateFunc <- function(x){
	as.Date(paste(unlist(strsplit(x, split="/"))[8:9], collapse="-"), format="%Y-%j")
}

read_int_file <- function(in_file, dsize, nbands, s_flag=T, bsq_flag=F){
  f_size <- file.info(in_file)$size # get the total dimensions of file = nrow*ncol*dsize*nbands
  f <- file(in_file,"rb") # open the file
  tot_size <- (f_size / dsize) # dsize is typically 2 for integers
  temp <- readBin(f, integer(), n=tot_size, endian="little", size=dsize, signed=s_flag) # read all of the data into an array
  close(f) # close the file
  # re-order the temp array into a matrix; if bsq_flag then read by row, otherwise by column
  ifelse(bsq_flag, byr_flag<-FALSE, byr_flag<-TRUE)
  temp <- matrix(temp, ncol=nbands, byrow=byr_flag)
  return(temp) # return the data
}


cl <- makeCluster(16)
clusterExport(cl, c("GetEVI2"))
clusterEvalQ(cl, {library(raster)})

data_dir <- "/projectnb/modislc/data/mcd12_in/c5/mcd43a4"
tile <- "h11v04"
tile_files <- Sys.glob(file.path(data_dir, "*", "*", paste("*", tile, "*", sep="")))
dates <- do.call("c", lapply(tile_files, DateFunc)) # get the dates for each file

start_date <- as.Date("2010-1-1")
end_date <- as.Date("2016-5-1")

# read in all rasters that are between start_date and end_date, make into a stack, then extract values
tile_files_sub <- tile_files[(dates >= start_date) & (dates <= end_date)]
my_dates <- dates[(dates >= start_date) & (dates <= end_date)]
evi2_list <- parLapply(cl, tile_files_sub, GetEVI2)
evi2_stack <- stack(unlist(evi2_list))
evi2_stack_values <- values(evi2_stack)

# get the LW mask
lw_mask_file <- Sys.glob(file.path("/projectnb/modislc/data/mcd12_in/c5/ancillary_layers/C5_LW_Mask/lw_mask_500m", paste("*", tile, "*", sep="")))
lw_mask_values <- read_int_file(lw_mask_file, 1, 1,s_flag=F)
lw_mask <- raster(evi2_stack)
values(lw_mask) <- lw_mask_values

# The Anderson's ethanol factory location
# is_h h fwd tp] lat 42.257777  long -84.790458  =>  vert tile 4  horiz tile 11  line 1857.63  samp 1737.15
anderson_ethanol <- c(xFromCol(evi2_stack, 1737.15), yFromRow(evi2_stack, 1857.63))
sub_extent <- c(100, 100) * xres(tmp_r)
anderson_extent <- extent(anderson_ethanol[1] - round(sub_extent[1]/2), anderson_ethanol[1] + round(sub_extent[1]/2), anderson_ethanol[2] - round(sub_extent[2]/2), anderson_ethanol[2] + round(sub_extent[2]/2))

my_pal <- colorRampPalette(brewer.pal(11, "Spectral"))
evi_breaks <- c(-1e6, seq(0.1, 0.7, len=254), 1e6)
out_dir <- "/projectnb/modislc/users/joshgray/C3/AndersonMovie"
for(i in 1:nlayers(evi2_stack)){
	print(paste("Working on", i, "of", nlayers(evi2_stack)))
	out_file <- file.path(out_dir, paste("AndersonMovie", formatC(i, width=3, flag="0"), ".jpg", sep=""))
	jpeg(file=out_file, width=(2400*2), height=(2400), quality=75)
	layout(matrix(1:2, nrow=1))
	par(mar=rep(1,4), bg="black", oma=rep(1, 4))
	tmp_date <- my_dates[i]
	tmp_r <- raster(evi2_stack, i)
	tmp_r[lw_mask != 1] <- NA
	plot(tmp_r, col=my_pal(255), legend=F, breaks=evi_breaks, xaxt="n", yaxt="n", maxpixels=2.5e6)
	abline(v=anderson_ethanol[1], lty=1, col="black", lwd=4)
	abline(h=anderson_ethanol[2], lty=1, col="black", lwd=4)
	plot(anderson_extent, add=T, col="black", lwd=4)
	text(par()$usr[1], par()$usr[4], labels=dates[i], col="black", cex=5, adj=c(-0.1, 2.1))

	tmp_sub <- crop(tmp_r, anderson_extent)
	plot(tmp_sub, col=my_pal(255), legend=F, breaks=evi_breaks, xaxt="n", yaxt="n", maxpixels=2.5e6)
	# abline(v=anderson_ethanol[1], lty=1, col="white", lwd=1.5)
	# abline(h=anderson_ethanol[2], lty=1, col="white", lwd=1.5)
	points(anderson_ethanol[1], anderson_ethanol[2], col="black", pch=4)
	dev.off()
}

# make into a movie with ffmpeg:
# module load ffmpeg/git-2013-12-18-4a2570f
# ffmpeg -framerate 4 -i AndersonMovie%03d.jpg -pix_fmt yuv420p -vf scale="720:trunc(ow/a/2)*2" AndersonEthanolEVI.mp4
