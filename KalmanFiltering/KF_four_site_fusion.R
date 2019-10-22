#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# DLM fusion for KF four site
# Josh Gray, August 2019
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# prelims
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
library(raster)
library(rgdal)
library(dlm)
library(viridis)
library(parallel)
library(abind)
source(normalizePath("~/Documents/pixel-forge/KalmanFiltering/KF_functions.R"))
data_dir <- normalizePath("/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat")
output_dir <- normalizePath("/Volumes/users/j/jmgray2/SEAL/KF_four_site_landsat/output")

site_to_do <- "central_valley"

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Gather all image files, dates, and order by date
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# get ARD data
ard_files <- dir(file.path(data_dir, site_to_do, "cv_ard"), pattern=".*vrt$", full=T)
ard_dates <- as.Date(gsub(".*_([0-9]{8})_[0-9]{8}.*vrt$", "\\1", basename(ard_files)), format="%Y%m%d")
ard_files <- ard_files[order(ard_dates)]
ard_dates <- sort(ard_dates)

# get MODIS files and dates
modis_files <- dir(file.path(data_dir, site_to_do, "cv_mcd43"), pattern=".*vrt$", full=T)
modis_dates <- as.Date(gsub(".*_([0-9]{4}_[0-9]{3})_.*", "\\1", basename(modis_files)), format="%Y_%j")
modis_files <- modis_files[order(modis_dates)]
modis_dates <- sort(modis_dates)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Make a plot of the image dates
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
xlim <- as.Date(paste("2018", strftime(range(c(ard_dates, modis_dates)), format="%m"), c("1", "31"), sep="-"))
ylim <- c(0, 1)
pch <- 16
cexs <- c(1, 2)
cols <- magma(6)[c(2,4)]
plot(NA, type="n", xlim=xlim, ylim=ylim, xlab="", ylab="", xaxt="n", yaxt="n")
axis(side=1, at=seq(min(xlim), max(xlim), by="month"), labels=as.Date(seq(min(xlim), max(xlim), by="month")))
points(ard_dates, rep(0.5, length(ard_dates)), pch=pch, cex=cexs[2], col=cols[2])
points(modis_dates, rep(0.5, length(modis_dates)), pch=pch, cex=cexs[1], col=cols[1])
legend("top", legend=c("ARD", "MODIS"), pch=pch, col=rev(cols), bty="n", cex=1.5, horiz=T)


#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# create RMSE maps
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("GetRMSE"))
clusterEvalQ(cl, {library(raster); library(rgdal); source("/Users/jmgray2/Documents/pixel-forge/KalmanFiltering/KF_functions.R")})
num_rows_planet <- 1000
example_planet_r <- raster(planet_files[1])
row_seq <- seq(1, nrow(example_planet_r), by=num_rows_planet)

# Make Planet-OLI RMSE maps
oli_2_planet_match_bands <- 2:5
the_cell_map <- dir(file.path(output_dir, "CellMaps"), pattern="OLI", full=T)
system.time(tmp <- parLapply(cl, row_seq, GetRMSE, rows_to_do=num_rows_planet, match_files=oli_files, ref_files=planet_files, match_dates=oli_dates, ref_dates=planet_dates, ref_to_match_cell_map_file=the_cell_map, ref_bands=1:4, match_bands=oli_2_planet_match_bands, match_qa_band=NA, match_qa_good_vals=NA))
system.time(tmp_c <- do.call(rbind, tmp))
tmp_rmse_s <- stack(planet_files[1])
values(tmp_rmse_s) <- c(tmp_c)
out_file <- file.path(output_dir, "RMSEMaps", "Planet2OLI_RMSE.tif")
writeRaster(tmp_rmse_s, file=out_file)
rm(tmp, tmp_c, tmp_rmse_s)

# Make Planet-MSI RMSE maps
msi_2_planet_match_bands <- c(2:4, 8)
the_cell_map <- dir(file.path(output_dir, "CellMaps"), pattern="MSI", full=T)
system.time(tmp <- parLapply(cl, row_seq, GetRMSE, rows_to_do=num_rows_planet, match_files=msi_files, ref_files=planet_files, match_dates=msi_dates, ref_dates=planet_dates, ref_to_match_cell_map_file=the_cell_map, ref_bands=1:4, match_bands=msi_2_planet_match_bands, match_qa_band=NA, match_qa_good_vals=NA))
system.time(tmp_c <- do.call(rbind, tmp))
tmp_rmse_s <- stack(planet_files[1])
values(tmp_rmse_s) <- c(tmp_c)
out_file <- file.path(output_dir, "RMSEMaps", "Planet2MSI_RMSE.tif")
writeRaster(tmp_rmse_s, file=out_file)
rm(tmp, tmp_c, tmp_rmse_s)

# Make Planet-MODIS RMSE maps
modis_2_planet_match_bands <- c(3, 4, 1, 2)
the_cell_map <- dir(file.path(output_dir, "CellMaps"), pattern="MODIS", full=T)
system.time(tmp <- parLapply(cl, row_seq, GetRMSE, rows_to_do=num_rows_planet, match_files=modis_files, ref_files=planet_files, match_dates=modis_dates, ref_dates=planet_dates, ref_to_match_cell_map_file=the_cell_map, ref_bands=1:4, match_bands=modis_2_planet_match_bands, match_qa_band=NA, match_qa_good_vals=NA))
system.time(tmp_c <- do.call(rbind, tmp))
tmp_rmse_s <- stack(planet_files[1])
values(tmp_rmse_s) <- c(tmp_c)
out_file <- file.path(output_dir, "RMSEMaps", "Planet2MODIS_RMSE.tif")
writeRaster(tmp_rmse_s, file=out_file)
rm(tmp, tmp_c, tmp_rmse_s)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# Get state covariance/correlation
state_bands <- 2:7
# qa_band <- 7
# qa_good_vals <- c(0)
oli_state_cov <- array(NA, dim=c(length(state_bands), length(state_bands), length(oli_dates)))
oli_state_cor <- array(NA, dim=c(length(state_bands), length(state_bands), length(oli_dates)))
for(i in 1:length(oli_dates)){
    print(paste("Doing OLI date", oli_dates[i]))
    this_date <- oli_dates[i]
    oli_file <- oli_files[oli_dates == this_date]
    oli_s <- stack(oli_file)
    oli_v <- values(oli_s)
    # oli_qa <- oli_v[, rep(qa_band, length(state_bands))]
    oli_v <- oli_v[, state_bands]
    # oli_v[!(oli_qa %in% qa_good_vals)] <- NA
    oli_state_cov[,, i] <- cov(oli_v, use="pairwise.complete.obs")
    oli_state_cor[,, i] <- cor(oli_v, use="pairwise.complete.obs")
}

state_bands <- c(2:4, 8, 12:13)
# qa_band <- 7
# qa_good_vals <- c(0)
msi_state_cov <- array(NA, dim=c(length(state_bands), length(state_bands), length(msi_dates)))
msi_state_cor <- array(NA, dim=c(length(state_bands), length(state_bands), length(msi_dates)))
for(i in 1:length(msi_dates)){
    print(paste("Doing MSI date", oli_dates[i]))
    this_date <- msi_dates[i]
    msi_file <- msi_files[msi_dates == this_date]
    msi_s <- stack(msi_file)
    msi_v <- values(msi_s)
    # oli_qa <- oli_v[, rep(qa_band, length(state_bands))]
    msi_v <- msi_v[, state_bands]
    # oli_v[!(oli_qa %in% qa_good_vals)] <- NA
    msi_state_cov[,, i] <- cov(msi_v, use="pairwise.complete.obs")
    msi_state_cor[,, i] <- cor(msi_v, use="pairwise.complete.obs")
}

# fit a smoothing spline to the covariance values to get a time-varying
# version of band-to-band covariance that covers all of the output dates
# also visualize to make sure it's what we expect
# out_dates <- seq.Date(min(oli_dates), max(oli_dates), by="day")
combined_dates <- c(oli_dates, msi_dates)
# all_dates <- c(planet_dates, oli_dates, msi_dates, modis_dates)
# out_dates <- seq.Date(min(oli_dates), max(oli_dates), by="day")
pred_dates <- seq.Date(min(combined_dates), max(combined_dates), by="day")
smooth_combined_cov <- array(NA, dim=c(length(state_bands), length(state_bands), length(pred_dates)))
band_names <- c("Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2")
pt_cols <- c(viridis(5)[3], plasma(5)[3])
layout(matrix(1:36, nrow=6, byrow=T))
par(mar=c(2, 1, 4, 1), oma=c(4, 4, 1, 1))
for(i in 1:6){
    for(j in 1:6){
        oli_cov_vals <- oli_state_cov[i, j, ]
        msi_cov_vals <- msi_state_cov[i, j, ]
        the_cov_vals <- c(oli_state_cov[i, j, ], msi_state_cov[i, j, ])
        spl <- smooth.spline(as.integer(combined_dates)[!is.na(the_cov_vals)], the_cov_vals[!is.na(the_cov_vals)])
        smooth_cov <- predict(spl, as.integer(pred_dates))$y
        smooth_combined_cov[i, j, ] <- smooth_cov
        ylim <- range(c(oli_state_cov[i, j, ], msi_state_cov[i, j, ]), na.rm=T)
        plot(oli_dates, oli_state_cov[i, j, ], pch=16, xlim=range(pred_dates), ylim=ylim, col=pt_cols[1], xlab="", ylab="")
        points(msi_dates, msi_state_cov[i, j, ], pch=16, col=pt_cols[2][])
        points(pred_dates, smooth_cov, type="l", col=rgb(0.3, 0.3, 0.3))
        abline(h=0, lty=2)
        if(i == 1) title(band_names[j])
        if(j == 1) mtext(band_names[i], side=2, line=2)
        if(i == 6 & j == 6) legend("bottomleft", legend=c("OLI", "MSI"), pch=16, col=pt_cols, title="Covariance", ncol=2)
    }
}

# create a time varying state covariance X index matrix
state_cov_tv_mat <- matrix(0, nrow=length(state_bands), ncol=length(state_bands), byrow=T)
state_cov_tv_mat[upper.tri(state_cov_tv_mat, diag=T)] <- 1:((6*(6+1))/2)
tmp <- t(state_cov_tv_mat)
diag(tmp) <- 0
state_cov_tv_mat <- state_cov_tv_mat + tmp

# create X
X <- matrix(NA, nrow=length(pred_dates), ncol=max(state_cov_tv_mat))
for(i in 1:max(state_cov_tv_mat)){
    mat_inds <- which(state_cov_tv_mat == i, arr.ind=T)
    X[, i] <- smooth_combined_cov[mat_inds[1, 1], mat_inds[1, 2],]
}
# identical(X[,4], smooth_oli_cov[1,3,])
save(state_cov_tv_mat, X, file=file.path(output_dir, "tv_state_cov.Rdata"))

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# set up prototype DLM
data_files <- list(planet_files, oli_files, msi_files, modis_files)
data_dates <- list(planet_dates, oli_dates, msi_dates, modis_dates)
cell_maps <- c(NA, file.path(output_dir, "CellMaps", "Planet2OLI_cellmap.tif"), file.path(output_dir, "CellMaps", "Planet2MSI_cellmap.tif"), file.path(output_dir, "CellMaps", "Planet2MODIS_cellmap.tif"))
qa_bands <- c(NA, NA, NA, NA) # the qa band, if any, for each data set
oli_qa_weights <- matrix(c(0, 1, 1, 0, 2, 0, 3, 1, 4, 1, 5, 0), ncol=2, byrow=T)
msi_qa_weights <- matrix(c(0, 1, 1, 0, 2, 0, 3, 1, 4, 1, 5, 0), ncol=2, byrow=T)
modis_qa_weights <- matrix(c(0, 1, 1, 0.75, 2, 0.5, 3, 0.25, 4, 0), ncol=2, byrow=T)
qa_weights <- list(NA, oli_qa_weights, msi_qa_weights, modis_qa_weights) # 2-col matrices of QA band values and associated weights for each data set
# measure_states <- list(1:4, 1:6, 1:6, 1:6) # the ordering of measured system states for each data set
measure_states <- list(1:4, 2:7, c(2:4, 8, 12:13), c(3, 4, 1:2, 6:7)) # the ordering of measured system states for each data set
rmse_maps <- c(NA, file.path(output_dir, "RMSEMaps", "Planet2OLI_RMSE.tif"), file.path(output_dir, "RMSEMaps", "Planet2MSI_RMSE.tif"), file.path(output_dir, "RMSEMaps", "Planet2MODIS_RMSE.tif"))
rmse_states <- list(NA, c(1:4, NA, NA), c(1:4, NA, NA), c(1:4, NA, NA)) # the order of available RMSE information for each system state
rmse_na_rep <- c(NA, 4, 4, 4) # if a state has no RMSE information, this is the band we'll copy it from
meas_errors <- list(rep(100, 4), rep(100, 6), rep(100, 6), rep(100, 6)) # the default measurement errors for each state and dataset

start_row <- 1
rows_to_get <- 500
# system.time(the_data <- GetKFData(start_row=start_row, num_row=rows_to_get, data_files=data_files, data_dates=data_dates, cell_maps=cell_maps, rmse_maps=rmse_maps, rmse_states=rmse_states, rmse_na_rep=rmse_na_rep, meas_errors=meas_errors, qa_bands=qa_bands, qa_weights=qa_weights, measure_states=measure_states))
system.time(the_data <- GetKFData(start_row=start_row, num_row=rows_to_get, data_files=data_files, data_dates=data_dates, cell_maps=cell_maps, rmse_maps=rmse_maps, rmse_states=rmse_states, rmse_na_rep=rmse_na_rep, meas_errors=meas_errors, qa_bands=qa_bands, qa_weights=qa_weights, measure_states=measure_states, out_dates=pred_dates))
system.time(data_chunk <- do.call(abind, the_data))
dim(data_chunk)

# Sanity check: plot some pixels
tmp_s <- stack(planet_files[3])
plot(raster(tmp_s, 4), col=plasma(255))
# myc <- click(tmp_s, cell=T)
plotRGB(tmp_s, 3, 2, 1, stretch="lin")
points(xyFromCell(tmp_s,myc$cell), col="white", pch=16)
text(xyFromCell(tmp_s,myc$cell), col="white", pch=16, labels=1:nrow(myc), pos=3)

the_pixel <- myc$cell[3]
band_names <- c("Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2")
planet_array_start <- 1
oli_array_start <- 5
msi_array_start <- 11
modis_array_start <- 17
plot_dates <- sort(unique(c(planet_dates, oli_dates, msi_dates, modis_dates)))
pchs <- 1:4
plot_cols <- plasma(5)[1:4]
# plot_cols <- c(brewer.pal(5, "Greys")[4], brewer.pal(5, "Blues")[4], brewer.pal(5, "Reds")[4], brewer.pal(5, "Greens")[4])
layout(matrix(1:6, nrow=6))
par(mar=c(2, 4, 1, 1))
for(the_band in 1:6){
    if(the_band <= 4){
        ylim <- range(data_chunk[the_pixel, , c(the_band - 1 + planet_array_start, the_band - 1 + oli_array_start, the_band - 1 + msi_array_start, the_band - 1 + modis_array_start)], na.rm=T)
    }else{
        ylim <- range(data_chunk[the_pixel, , c(the_band - 1 + oli_array_start, the_band - 1 + msi_array_start, the_band - 1 + modis_array_start)], na.rm=T)
    }
    # plot(plot_dates, rep(NA, length(out_dates)), type="n",  xlab="", ylab=band_names[the_band], ylim=ylim)
    plot(plot_dates, rep(NA, length(plot_dates)), type="n",  xlab="", ylab=band_names[the_band], ylim=ylim)
    if(the_band <= 4) points(plot_dates, data_chunk[the_pixel, 2:(length(plot_dates) + 1), the_band - 1 + planet_array_start], pch=pchs[1], col=plot_cols[1])
    points(plot_dates, data_chunk[the_pixel, 2:(length(plot_dates) + 1), the_band - 1 + oli_array_start], pch=pchs[2], col=plot_cols[2])
    points(plot_dates, data_chunk[the_pixel, 2:(length(plot_dates) + 1), the_band - 1 + msi_array_start], pch=pchs[3], col=plot_cols[3])
    points(plot_dates, data_chunk[the_pixel, 2:(length(plot_dates) + 1), the_band - 1 + modis_array_start], pch=pchs[4], col=plot_cols[4])
    if(the_band == 1) legend("topright", legend=c("Planet", "OLI", "MSI", "MODIS"), pch=pchs, col=plot_cols)
}

# out_dates <- as.Date(sort(unique(unlist(plot_dates))), origin="1970-1-1")
pixel_to_do <- myc$cell[1]
num_states <- 6
sensors <- list(1:4, 1:6, 1:6, 1:6)
state_cov <- 1e4
default_meas_cov <- 100
bihar_dlm <- MakeMultiDLM(num_states=num_states, sensors=sensors)
system.time(kf_result <- DoBiharKF(data_chunk[pixel_to_do,,], bihar_dlm, sensor=sensors, state_cov=state_cov, default_meas_cov=default_meas_cov))

# plot the results
sensor_vec <- unlist(sensors)
library(RColorBrewer)
true_signal_color <- brewer.pal(5, "Greys")[3]
main_cols <- c(brewer.pal(5, "Blues")[3], brewer.pal(5, "Greens")[3], brewer.pal(5, "Reds")[3], brewer.pal(5, "Purples")[3], brewer.pal(5, "Oranges")[3], brewer.pal(5, "RdPu")[3])
err_cols <- c(brewer.pal(5, "Blues")[2], brewer.pal(5, "Greens")[2], brewer.pal(5, "Reds")[2], brewer.pal(5, "Purples")[2], brewer.pal(5, "Oranges")[2], brewer.pal(5, "RdPu")[2])

layout(matrix(1:num_states, nrow=num_states))
par(mar=c(2, 3, 1, 1))
for(i in 1:num_states) PlotForecast(kf_result[, i], kf_result[, i + num_states], t=pred_dates, signal=split(data_chunk[pixel_to_do, 2:dim(data_chunk)[2], sensor_vec == i], col(as.matrix(data_chunk[pixel_to_do, 2:dim(data_chunk)[2], sensor_vec == i]))), colmain=main_cols[i], colerr=err_cols[i])

# make KF with time-varying state covariance
load(file.path(output_dir, "tv_state_cov.Rdata")) # load the X and state cov index matrices
bihar_dlm_tv <- MakeMultiDLM(num_states=num_states, sensors=sensors)
# X(bihar_dlm_tv) <- X * 0.005
X(bihar_dlm_tv) <- abs(X)^0.5 * sign(X)
JW(bihar_dlm_tv) <- state_cov_tv_mat
system.time(kf_result <- DoBiharKF(data_chunk[pixel_to_do,,], bihar_dlm_tv, sensor=sensors, state_cov=state_cov, default_meas_cov=default_meas_cov))
sensor_vec <- unlist(sensors)

true_signal_color <- brewer.pal(5, "Greys")[3]
main_cols <- c(brewer.pal(5, "Blues")[3], brewer.pal(5, "Greens")[3], brewer.pal(5, "Reds")[3], brewer.pal(5, "Purples")[3], brewer.pal(5, "Oranges")[3], brewer.pal(5, "RdPu")[3])
err_cols <- c(brewer.pal(5, "Blues")[2], brewer.pal(5, "Greens")[2], brewer.pal(5, "Reds")[2], brewer.pal(5, "Purples")[2], brewer.pal(5, "Oranges")[2], brewer.pal(5, "RdPu")[2])
layout(matrix(1:num_states, nrow=num_states))
par(mar=c(2, 5, 1, 1))
# for(i in 1:num_states) PlotForecast(kf_result[, i], kf_result[, i + num_states], t=pred_dates, signal=split(data_chunk[pixel_to_do, 2:dim(data_chunk)[2], sensor_vec == i], col(as.matrix(data_chunk[pixel_to_do, 2:dim(data_chunk)[2], sensor_vec == i]))), colmain=main_cols[i], colerr=err_cols[i])
# for(i in 1:num_states) PlotForecast(kf_result[, i], kf_result[, i + num_states], t=pred_dates, signal=split(data_chunk[pixel_to_do, 2:dim(data_chunk)[2], sensor_vec == i], col(as.matrix(data_chunk[pixel_to_do, 2:dim(data_chunk)[2], sensor_vec == i]))), colmain=main_cols[i], colerr=err_cols[i], ylab=band_names[i])
for(i in 1:num_states) PlotForecast(kf_result[, i], kf_result[, i + num_states], t=pred_dates, signal=split(data_chunk[pixel_to_do, 2:dim(data_chunk)[2], GetSensorInd(sensors, i)], col(as.matrix(data_chunk[pixel_to_do, 2:dim(data_chunk)[2], GetSensorInd(sensors, i)]))), colmain=main_cols[i], colerr=err_cols[i], ylab=band_names[i])
legend("topright", legend=c("Planet", "OLI", "MSI", "MODIS"), pch=1:4, col="#636363", horiz=T, bty="o")


#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# do the fusion
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("DoBiharKF"))
clusterEvalQ(cl, {library(dlm)})
# system.time(chunk_kf_results <- parApply(cl, data_chunk, 1, DoBiharKF, the_dlm=bihar_dlm, sensors=sensors, state_cov=state_cov, default_meas_cov=default_meas_cov))
system.time(chunk_kf_results <- parApply(cl, data_chunk, 1, DoBiharKF, the_dlm=bihar_dlm_tv, sensors=sensors, state_cov=state_cov, default_meas_cov=default_meas_cov))
# convert the chunk output to an array that has dimensions: [dates, state est+se, pixels]
chunk_kf_result_mat <- array(chunk_kf_results, dim=c(length(out_dates), num_states * 2, dim(chunk_kf_results)[2]))

# write the results for this chunk
example_planet_r <- raster(planet_files[[1]])
out_dir <- normalizePath("~/Desktop/HeumannOutput")
out_root <- "Heumann_v0_"
for(i in 1:length(out_dates)){
    this_date <- out_dates[i]
    data_to_write <- t(chunk_kf_result_mat[i,,])
    out_file <- file.path(out_dir, paste(out_root, strftime(this_date, format="%Y_%j"), sep=""))
    WriteENVIMultiband(t(chunk_kf_result_mat[i,,]), out_file=out_file, example_r=example_planet_r)
}

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# visualization

in_files <- dir(out_dir, pattern="Heumann_v0_2018_[0-9]{3}$", full=T)
plot_out_dir <- normalizePath("~/Desktop/HeumannOutput/plots432")
i <- 1
for(in_file in in_files){
    print(paste("Doing", i, "of", length(in_files)))
    tmp_s <- stack(in_file)
    out_file <- file.path(plot_out_dir, paste(basename(in_file), ".jpeg", sep=""))
    jpeg(file=out_file, width=494, height=494, quality=100)
    plotRGB(tmp_s, 4, 3, 2, stretch="lin")
    dev.off()
    i <- i + 1
}

# ffmpeg -framerate 10 -pattern_type glob -i '*.jpeg' -pix_fmt yuv420p heumann_example_432.mp4

in_files <- dir(out_dir, pattern="Heumann_v0_2018_[0-9]{3}$", full=T)
plot_out_dir <- normalizePath("~/Desktop/HeumannOutput/plots543")
i <- 1
for(in_file in in_files){
    print(paste("Doing", i, "of", length(in_files)))
    tmp_s <- stack(in_file)
    out_file <- file.path(plot_out_dir, paste(basename(in_file), ".jpeg", sep=""))
    jpeg(file=out_file, width=494, height=494, quality=100)
    plotRGB(tmp_s, 5, 4, 3, stretch="lin")
    dev.off()
    i <- i + 1
}

# ffmpeg -framerate 10 -pattern_type glob -i '*.jpeg' -pix_fmt yuv420p heumann_example_543.mp4


