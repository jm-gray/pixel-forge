#=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
# Apply the multispectral Kalman Filter over four years and plot each result
# bsub -q cnr -W 16:00 -n 24 -R "span[ptile=24]" -o /share/jmgray2/KF_multispec/kf_multifig.out.%J -e /share/jmgray2/KF_multispec/kf_multifig.err.%J R CMD BATCH --vanilla /share/jmgray2/KF_multispec/KF_multispec_figures.R

library(tools)
library(raster)
library(rgdal)
library(parallel)
library(RColorBrewer)
library(dlm)
library(reshape2)

fig_out_dir <- "/share/jmgray2/KF_multispec/output_figures"

source("/share/jmgray2/KF_multispec/KF_functions.R")
cl <- makeCluster(24)
clusterEvalQ(cl, {source("/share/jmgray2/KF_multispec/KF_functions.R"); library(dlm); library(reshape2)})

load("/share/jmgray2/KF_multispec/nebraska_multispec_workspace.Rdata")

landsat_years <- as.integer(strftime(landsat_dates, format="%Y"))
modis_years <- as.integer(strftime(modis_dates, format="%Y"))
num_years <- length(unique(c(landsat_years, modis_years)))

# blank_r <- tmp_r; values(blank_r) <- NA
# out_s_blue <- do.call(stack, replicate(length(unique(c(landsat_years, modis_years)))*365, blank_r))
# out_s_green <- do.call(stack, replicate(length(unique(c(landsat_years, modis_years)))*365, blank_r))
# out_s_red <- do.call(stack, replicate(length(unique(c(landsat_years, modis_years)))*365, blank_r))
# out_s_nir <- do.call(stack, replicate(length(unique(c(landsat_years, modis_years)))*365, blank_r))

blue_inds <- 1:1460
green_inds <- 1461:2920
red_inds <- 2921:4380
nir_inds <- 4381:5840

output_data <- matrix(NA, nrow=ncell(tmp_r), ncol=5840)

rows_to_do <- 1e4
row_seq <- seq(1, nrow(Y), by=rows_to_do)

for(i in row_seq){
	print(i)
	i_end <- min(i + rows_to_do - 1, nrow(Y))
	tmp <- parApply(cl, Y[i:i_end,], 1, FuseLandsatModisMultispec, landsat_dates=landsat_dates, modis_dates=modis_dates, landsat_sensor=landsat_sensor, multispectral_cdl_process_cov=multispectral_cdl_process_cov, cdl_types=cdl_types)
	# tmp <- apply(Y[i:i_end,], 1, FuseLandsatModisMultispec, landsat_dates=landsat_dates, modis_dates=modis_dates, landsat_sensor=landsat_sensor, multispectral_cdl_process_cov=multispectral_cdl_process_cov, cdl_types=cdl_types)
	# print(dim(tmp))
	# out_s_blue[i:i_end] <- c(t(tmp[blue_inds,]))
	# out_s_green[i:i_end] <- c(t(tmp[green_inds,]))
	# out_s_red[i:i_end] <- c(t(tmp[red_inds,]))
	# out_s_nir[i:i_end] <- c(t(tmp[nir_inds,]))
	output_data[i:i_end,] <- t(tmp)
}

# plot_dates <- sort(as.Date(paste(rep(sort(unique(c(landsat_years, modis_years))), each=365), rep(1:365, num_years), sep="-"), format="%Y-%j"))
# for(i in 1:nlayers(out_s_blue)){
# 	out_file <- file.path(fig_out_dir, paste(as.character(strftime(plot_dates[i], format="%Y%j")), "_kf_multispec.jpg", sep=""))
# 	jpeg(file=out_file, width=1214, height=772, quality=75)
# 	plot_s <- stack(raster(out_s_nir, i), raster(out_s_red, i), raster(out_s_green, i))
# 	plotRGB(plot_s, stretch="lin")
# 	dev.off()
# }
#

plot_dates <- sort(as.Date(paste(rep(sort(unique(c(landsat_years, modis_years))), each=365), rep(1:365, num_years), sep="-"), format="%Y-%j"))
for(i in 1:length(plot_dates)){
	out_file <- file.path(fig_out_dir, paste(as.character(strftime(plot_dates[i], format="%Y%j")), "_kf_multispec.jpg", sep=""))
	jpeg(file=out_file, width=1214, height=772, quality=75)
	tmp_red <- tmp_r
	tmp_blue <- tmp_r
	tmp_green <- tmp_r
	values(tmp_red) <- output_data[, min(nir_inds) + (i - 1)]
	values(tmp_blue) <- output_data[, min(red_inds) + (i - 1)]
	values(tmp_green) <- output_data[, min(green_inds) + (i - 1)]
	plot_s <- stack(tmp_red, tmp_blue, tmp_green)
	plotRGB(plot_s, stretch="lin")
	dev.off()
}
