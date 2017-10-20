load("~/Google Drive/Projects/KalmanFiltering/KF_fusion_data_new/nebraska_multispec_workspace.Rdata")
source("~/Documents/pixel-forge/KalmanFiltering/KF_functions.R")
library(dlm)
library(RColorBrewer)
library(parallel)

#-------------------------------------------------------------------------
expand_x <- function(x, landsat_dates, modis_dates, pred_dates, nbands=4, landsat_band_start, landsat_band_end, modis_band_start, modis_band_end){
  X <- matrix(NA, nrow=nbands * 2, ncol=length(pred_dates)) # setup mtarix
  modis_inds <- match(pred_dates, modis_dates) # band-specific modis indices for each date
  landsat_inds <- match(pred_dates, landsat_dates) # same for landsat
  # assign the values to the matrix for each band, matched to the proper date
  for(i in 1:nbands){
    X[i,] <- x[landsat_band_start[i]:landsat_band_end[i]][landsat_inds]
    X[i + nbands,] <- x[modis_band_start[i]:modis_band_end[i]][modis_inds]
  }
  return(X)
}

#-------------------------------------------------------------------------
plotdlmresults <- function(Y, dlm_result, pred_dates, plot_band=1, nbands=4){
  plot_cols <- c("#DE2D26", "#3182BD", "#636363", "#CCCCCC") # x1, x2, filt, filt error
  ylim <- range(c(Y[plot_band,], Y[plot_band + nbands,], dlm_result[,plot_band]), na.rm=T)
  plot(pred_dates, Y[plot_band,], type="n", ylim=ylim)
  points(pred_dates, Y[plot_band,], col=plot_cols[1], pch=1)
  points(pred_dates, Y[plot_band + nbands,], col=plot_cols[2], pch=2)
  points(pred_dates, dlm_result[,plot_band], type="l", lty=2, lwd=2, col=plot_cols[3])
}

#-------------------------------------------------------------------------
# determine indices for each landsat and modis band in a row of the data
x_skip <- 5 # first value is the landcover code (CDL), next 4 are band-specific MODIS-Landsat RMSE values for that pixel
nbands <- 4 # we have blue, green, red, and nir surface reflectance from Landsat and MODIS
landsat_band_start <- x_skip + seq(1, length(landsat_dates) * nbands, by=length(landsat_dates))
landsat_band_end <- landsat_band_start + length(landsat_dates) - 1
modis_band_start <- seq(landsat_band_end[nbands], ncol(Y) - 1, by=length(modis_dates)) + 1
modis_band_end <- modis_band_start + length(modis_dates) - 1

#-------------------------------------------------------------------------
# create the data matrix for all unique landsat and modis dates (no interpolation)
pred_dates <- sort(unique(c(landsat_dates, modis_dates)))
pixel_number <- 1
Y_tmp <- expand_x(Y[pixel_number,], landsat_dates, modis_dates, pred_dates, nbands=4, landsat_band_start, landsat_band_end, modis_band_start, modis_band_end)

#-------------------------------------------------------------------------
# simple KF with arbitrary errors
mydlm <- MakeMultiDLM(4, 2)
mydlm$m0 <- Y_tmp[1:4, min(which(!is.na(Y_tmp[1, ])))] # set initial condition to first non-missing Landsat observation (fragile: perhaps not all Landsat band obs are non-missing!)
landsat_obs_error <- 100 # assume arbitrary constant landsat observation error
diag(mydlm$V) <- c(rep(landsat_obs_error, 4), Y[pixel_number, 2:5]) # modis-landsat long-term RMSE is in columns 2:5 of Y
# and the process error matrix
process_error <- 10 # assume arbitrary constant process error
diag(mydlm$W) <- rep(process_error, 4)
m_filt <- dlmFilter(t(Y_tmp), mydlm) # filtering
m_smooth <- dlmSmooth(t(Y_tmp), mydlm) # smoothing

layout(matrix(1:4, nrow=2, byrow=T))
par(mar=rep(1, 4))
for(i in 1:4) plotdlmresults(Y_tmp, dropFirst(m_smooth$s), pred_dates=pred_dates, plot_band=i) # for smoothing distributions

#-------------------------------------------------------------------------
# time-varying errors
mydlm <- MakeMultiDLM(4, 2)
mydlm$m0 <- Y_tmp[1:4, min(which(!is.na(Y_tmp[1, ])))] # set initial condition to first non-missing Landsat observation (fragile: perhaps not all Landsat band obs are non-missing!)
diag(mydlm$V) <- c(rep(landsat_obs_error, 4), x[2:5]) # modis-landsat long-term RMSE is in columns 2:5 of Y
mydlm$JV <- diag(c(1:nbands, rep(0, nbands))) # allow the Landsat observation error to be time-varying, use long-term MODIS-landsat RMSE for MODIS obs error
# prototype the X matrix which holds time-varying values. First nbands will hold Landsat time-varying error
X <- matrix(NA, nrow=ncol(Y_tmp), ncol=nbands)
# create a landsat error multiplier vector
etm_mult_error <- 0.05 # multiplicative error for ETM+ (Landsat 7) obs
tm_mult_error <- 0.07 # multiplicative error for TM (Landsat 5) obs
landsat_error <- rep(NA, ncol(Y_tmp))
landsat_error[pred_dates %in% landsat_dates][landsat_sensor == "LE7"] <- etm_mult_error
landsat_error[pred_dates %in% landsat_dates][landsat_sensor == "LT5"] <- tm_mult_error
# populated the time-varying error values for each band
for(i in 1:nbands) X[, i] <- Y_tmp[i, ] * landsat_error
mydlm$X <- X

# and the process error matrix
process_error <- 50 # assume arbitrary constant process error
diag(mydlm$W) <- rep(process_error, 4)

m_smooth <- dlmSmooth(t(Y_tmp), mydlm) # smoothing
layout(matrix(1:4, nrow=2, byrow=T))
par(mar=rep(1, 4))
for(i in 1:4) plotdlmresults(Y_tmp, dropFirst(m_smooth$s), pred_dates=pred_dates, plot_band=i) # for smoothing distributions



#-------------------------------------------------------------------------
# recreate Landsat data for a particular date
blue_r <- raster(tmp_r)
green_r <- raster(tmp_r)
red_r <- raster(tmp_r)
nir_r <- raster(tmp_r)
match_date <- all_dates[96]
landsat_ind <- which(landsat_dates == match_date)
# which(landsat_date == all_dates[310])
values(blue_r) <- Y[, c(landsat_band_start[1]:landsat_band_end[1])[landsat_ind]]
values(green_r) <- Y[, c(landsat_band_start[2]:landsat_band_end[2])[landsat_ind]]
values(red_r) <- Y[, c(landsat_band_start[3]:landsat_band_end[3])[landsat_ind]]
values(nir_r) <- Y[, c(landsat_band_start[4]:landsat_band_end[4])[landsat_ind]]
s <- stack(blue_r, green_r, red_r, nir_r)
plotRGB(s, stretch="lin")
points(xyFromCell(tmp_r, pixel_number))
