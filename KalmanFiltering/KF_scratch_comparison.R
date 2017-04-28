#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Jack knife evaluation of KF error

set.seed(42)
t <- 1:100
# trueSignal <- predict(smooth.spline(1:100, cumsum(rnorm(100))))$y # random walk
trueSignal <- 10 * sin(((2 * pi) / 200) * t) + 1
zA <- trueSignal + rnorm(length(trueSignal), mean=0, sd=1)
zB <- trueSignal + rnorm(length(trueSignal), mean=0, sd=2)
zC <- trueSignal + rnorm(length(trueSignal), mean=0, sd=3)
z <- matrix(c(zA, zB, zC), nrow=3, byrow=T)

A <- matrix(1, nrow=1) # process transition
B <- matrix(1, nrow=1) # control
Q <- matrix(0.1, nrow=1) # process noise
H <- matrix(c(1, 1, 1), nrow=3, byrow=T) # observation
R <- matrix(c(1, 0, 0, 0, 4, 0, 0, 0, 9), nrow=3, byrow=T) # obs noise
kf <- list(A=A, B=B, Q=Q, H=H, R=R)

# external control
uA <- c(0, diff(trueSignal)) # append 0 to head b/c 1st step is X0 -> X1
u <- matrix(uA, nrow=1, byrow=T)

# run the KF
X <- list(x=1, P=1) # initialize
kf_result <- Kfilter(X, z, u, kf)
filterX <- unlist(lapply(kf_result, "[[", 1)) # note, this only works for univariate x!
filterP <- unlist(lapply(kf_result, "[[", 2))


# obs_index <- 1 # row in z where the "true" observations are held; the ones from which the resid is calc


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# real data
landsat_data_dir <- "~/Desktop/KF_fusion_data_new/landsat"
landsat_evi2_files <- dir(landsat_data_dir, pattern="evi2.tif", full=T)
landsat_evi2_files <- landsat_evi2_files[order(GetLandsatDate(landsat_evi2_files))]
landsat_dates <- GetLandsatDate(landsat_evi2_files)
landsat_doys <- as.integer(strftime(landsat_dates, format="%j"))

modis_data_dir <- "~/Desktop/KF_fusion_data_new/modis"
modis_evi2_files <- dir(modis_data_dir, pattern="evi2.tif", full=T)
modis_dates <- as.Date(unlist(lapply(modis_evi2_files, GetModisDate)), origin="1970-1-1")
modis_evi2_files <- modis_evi2_files[order(modis_dates)]
modis_doys <- as.integer(strftime(modis_dates, format="%j"))

rmse_sub <- raster("~/Desktop/KF_fusion_data_new/rmse_sub.tif")

year_to_do <- 2008
rows_to_read <- 1e3
tmp_r <- raster(landsat_evi2_files[1])
is <- 1:(ceiling(nrow(tmp_r) / rows_to_read))

i <- 1
start_row <- ((i - 1) * rows_to_read) + 1
nrows <- min(rows_to_read, (nrow(tmp_r) - ((i - 1) * rows_to_read))) # last block may have less rows
landsat_v <- GetValuesGDAL(landsat_evi2_files, start_row, nrows)
modis_v <- GetValuesGDAL(modis_evi2_files, start_row, nrows)
modis_snow_v <- GetValuesGDAL(modis_snow_files, start_row, nrows)
modis_v[modis_snow_v == 1] <- NA # screen out all snow
rmse_v <- getValues(rmse_sub, start_row, nrows)
V <- cbind(rmse_v, landsat_v, modis_v)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
KF_fusion <- function(z, x_dates, y_dates, kf, u, X0, rmse_col=1, smooth=T, use_rmse=T){
  rmse_tmp <- z[rmse_col]
  z <- z[-rmse_col]
  z_annual <- matrix(NA, nrow=2, ncol=365)
  x_doys <- as.integer(strftime(x_dates, format="%j"))
  y_doys <- as.integer(strftime(y_dates, format="%j"))
  z_annual[1, x_doys] <- z[1:length(x_dates)]
  z_annual[2, y_doys] <- z[(length(x_dates) + 1):length(z)]

  # adjust rmse for modis obs
  if(use_rmse) kf$R[2, 2] <- rmse_tmp^2

  if(smooth){
    kf_result <- Ksmooth(X0, z_annual, u, kf)
  }else{
    kf_result <- Kfilter(X0, z_annual, u, kf)
  }

  filterX <- unlist(lapply(kf_result, "[[", 1)) # note, this only works for univariate x!
  filterP <- unlist(lapply(kf_result, "[[", 2))
  return(list(filterX=filterX, filterP=filterP))
}

# EVI2, landsat + modis obs
A <- matrix(1, nrow=1) # process transition
B <- matrix(1, nrow=1) # control
Q <- 0.00005 * matrix(1, nrow=1) # process noise
H <- matrix(c(1, 1), nrow=2, byrow=T) # observation
R <- 0.001 * matrix(c(1, 0, 0, 1), nrow=2) # obs noise
kf <- list(A=A, B=B, Q=Q, H=H, R=R)

uA <- rep(0, 365) # no external forcing
u <- matrix(uA, nrow=1, byrow=T)

# initial conditions
X0 <- list(x=0, P=1) # initialize

x <- V[1,] / 1e4
tmp_result <- KF_fusion(x, landsat_dates, modis_dates, kf, u, X0)
x <- x[-1] # remove rmse column for plotting simplicity
filterX <- tmp_result$filterX
filterP <- tmp_result$filterP

# plot some results
upper <- filterX + (2 * sqrt(filterP))
lower <- filterX - (2 * sqrt(filterP))
mycols <- c(brewer.pal(5, "Greys")[4], brewer.pal(5, "Reds")[4], brewer.pal(5, "Blues")[4], brewer.pal(5, "Greens")[4], brewer.pal(5, "Purples")[4])
plot(1:365, filterX, ylim=c(min(c(x, filterX, lower), na.rm=T), max(c(x, filterX, upper), na.rm=T)), type="n", xlab="", ylab="EVI2")
polygon(x=c(1:365, 365:1), y=c(upper, rev(lower)), border=NA, col=rgb(0.8, 0.8, 0.8))
points(as.integer(strftime(landsat_dates, format="%j")), x[1:length(landsat_dates)], col=mycols[2], pch=1)
points(as.integer(strftime(modis_dates, format="%j")), x[(length(landsat_dates) + 1): length(x)], col=mycols[3], pch=2)
# points(1:365, z_annual[2,], col=mycols[3], pch=2)
points(1:365, filterX, col=mycols[1], type="l", lty=2)
legend("topleft", legend=c("landsat", "modis", "KF"), col=c(mycols[2], mycols[3], mycols[1]), lty=c(NA, NA, 2), lwd=c(NA, NA, 1), pch=c(1, 2, NA))


#----------------------------
# Compare by_hand with dlm implementation
A <- matrix(1, nrow=1) # process transition
B <- matrix(1, nrow=1) # control
Q <- matrix(1, nrow=1) # process noise
H <- matrix(c(1, 1), nrow=2, byrow=T) # observation
R <- matrix(c(1, 0, 0, 1), nrow=2) # obs noise
kf <- list(A=A, B=B, Q=Q, H=H, R=R)
u <- matrix(rep(0, 365), nrow=1)
X0 <- list(x=0, P=1) # initialize

# build the dlm KF model
tmp_dlm <- MakeMultiDLM(num_states=1, sensors=2)
tmp_dlm$GG <- A
tmp_dlm$W <- Q
tmp_dlm$FF <- H
tmp_dlm$V <- R
tmp_dlm$m0 <- X0$x
tmp_dlm$C0 <- X0$P

# get the data
x <- V[1, -rmse_col] / 1e4
v_annual <- matrix(NA, nrow=2, ncol=365)
landsat_doys <- as.integer(strftime(landsat_dates, format="%j"))
modis_doys <- as.integer(strftime(modis_dates, format="%j"))
v_annual[1, landsat_doys] <- x[1:length(landsat_dates)]
v_annual[2, modis_doys] <- x[(length(landsat_dates) + 1):length(x)]

# do the KF "by hand"
kf_bh_result <- Kfilter(X0, v_annual, u, kf)
kf_bh_X <- unlist(lapply(kf_bh_result, "[", c("x")), rec=T)
kf_bh_P <- unlist(lapply(kf_bh_result, "[", c("P")), rec=T)

# do the KF w/ dlm
kf_dlm_result <- dlmFilter(t(v_annual), tmp_dlm)
kf_dlm_X <- dropFirst(kf_dlm_result$m)
kf_dlm_P <- sqrt(unlist(dropFirst(dlmSvd2var(filt$U.C, filt$D.C))))

plot(1:365, rep(NA, 365), type="n", ylim=c(min(kf_dlm_X, kf_bh_X, v_annual, na.rm=T), max(kf_dlm_X, kf_bh_X, v_annual, na.rm=T)))
points(1:365, kf_bh_X, type="l", col=mycols[2], lwd=2)
points(1:365, kf_dlm_X, type="l", col=mycols[3])
points(1:365, v_annual[1, ], pch=1)
points(1:365, v_annual[2, ], pch=2)
legend("topleft", legend=c("JG", "dlm", "landsat", "modis"), col=c(mycols[2], mycols[3], 1, 1), pch=c(NA, NA, 1, 2), lty=c(1, 1, NA, NA))



#######################################
# Results are identical! So my code works, at least for filtering of 1 state and 1 sensor
# 1 state 1 sensor
A <- matrix(1, nrow=1) # process transition
B <- matrix(0, nrow=1) # control
Q <- matrix(0.01, nrow=1) # process noise
# H <- matrix(1, nrow=1) # observation
H <- matrix(c(1, 1), nrow=2, byrow=T) # observation
# R <- matrix(1, nrow=1) # obs noise
# R <- matrix(c(1, 0, 0, 1), nrow=2, byrow=T) # obs noise
R <- matrix(c(0.06, 0, 0, 0.16), nrow=2, byrow=T) # obs noise
kf <- list(A=A, B=B, Q=Q, H=H, R=R)
u <- matrix(rep(0, 365), nrow=1)
X0 <- list(x=0, P=1) # initialize

# build the dlm KF model
# tmp_dlm <- MakeMultiDLM(num_states=1, sensors=1)
tmp_dlm <- MakeMultiDLM(num_states=1, sensors=2)
tmp_dlm$GG <- A
tmp_dlm$W <- Q
tmp_dlm$FF <- H
tmp_dlm$V <- R
tmp_dlm$m0 <- X0$x
tmp_dlm$C0 <- X0$P

# get the data
x <- V[1, -rmse_col] / 1e4
# v_annual <- matrix(NA, nrow=1, ncol=365)
v_annual <- matrix(NA, nrow=2, ncol=365)
landsat_doys <- as.integer(strftime(landsat_dates, format="%j"))
modis_doys <- as.integer(strftime(modis_dates, format="%j"))
v_annual[1, landsat_doys] <- x[1:length(landsat_dates)]
v_annual[2, modis_doys] <- x[(length(landsat_dates) + 1):length(x)]

# do the KF "by hand"
kf_bh_result <- Kfilter(X0, v_annual, u, kf)
kf_bh_X <- unlist(lapply(kf_bh_result, "[", c("x")), rec=T)
kf_bh_P <- unlist(lapply(kf_bh_result, "[", c("P")), rec=T)

ks_bh_result <- Ksmooth(X0, v_annual, u, kf)
ks_bh_X <- unlist(lapply(ks_bh_result, "[", c("x")), rec=T)
ks_bh_P <- unlist(lapply(ks_bh_result, "[", c("P")), rec=T)

# do the KF w/ dlm
kf_dlm_result <- dlmFilter(t(v_annual), tmp_dlm)
kf_dlm_X <- dropFirst(kf_dlm_result$m)
kf_dlm_P <- unlist(dropFirst(dlmSvd2var(kf_dlm_result$U.C, kf_dlm_result$D.C)))

ks_dlm_result <- dlmSmooth(t(v_annual), tmp_dlm)
ks_dlm_X <- dropFirst(ks_dlm_result$s)
ks_dlm_P <- unlist(dropFirst(dlmSvd2var(ks_dlm_result$U.S, ks_dlm_result$D.S)))


layout(matrix(1:2, nrow=2))
par(mar=rep(1, 4))
plot(1:365, rep(NA, 365), type="n", ylim=c(min(kf_dlm_X, kf_bh_X, v_annual, na.rm=T), max(kf_dlm_X, kf_bh_X, v_annual, na.rm=T)))
points(1:365, kf_bh_X, type="l", col=mycols[2], lwd=2)
points(1:365, kf_dlm_X, type="l", col=mycols[3], lty=2)
points(1:365, v_annual[1, ], pch=1)
legend("topleft", legend=c("JG", "dlm", "landsat", "modis"), col=c(mycols[2], mycols[3], 1, 1), pch=c(NA, NA, 1, 2), lty=c(1, 1, NA, NA))

plot(1:365, rep(NA, 365), type="n", ylim=c(min(ks_dlm_X, ks_bh_X, v_annual, na.rm=T), max(ks_dlm_X, ks_bh_X, v_annual, na.rm=T)))
points(1:365, ks_bh_X, type="l", col=mycols[2], lwd=2)
points(1:365, ks_dlm_X, type="l", col=mycols[3], lty=2)
points(1:365, v_annual[1, ], pch=1)
points(1:365, v_annual[2, ], pch=2)
# legend("topleft", legend=c("JG", "dlm", "landsat", "modis"), col=c(mycols[2], mycols[3], 1, 1), pch=c(NA, NA, 1, 2), lty=c(1, 1, NA, NA))
