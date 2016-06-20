#==========================================================================
# Prelims
#==========================================================================
library(dlm)
library(animation)

#-----------------------------------------------------
ShiftSample <- function(y, num_sample, shift=0, noise_frac=0, miss_fraction=0, seed=42, retInds=F){
	# Function to shift, sample, and add missing data to time series
	# takes even samples, assumes uniform multiplicative noise
	# do the shifting
	if(shift < 0){
		y <- c(y[abs(shift):length(y)], rep(y[length(y)], abs(shift) - 1))
	}else if(shift > 0){
		y <- c(rep(y[1], shift + 1), y[2:(length(y) - shift)])
	}

	# do the sampling
	set.seed(seed)
	sample_inds <- round(seq(1, length(y), length=num_sample))
	y <- y[sample_inds]

	# add noise
	y <- y * runif(num_sample, 1 - noise_frac, 1 + noise_frac)

	# add missing data
	y[sample(1:length(y), round(num_sample * miss_fraction))] <- NA

	ifelse(retInds, return(list(y, sample_inds)), return(y))
}

#-----------------------------------------------------
CrossingDates <- function(y, upfracs=c(0.2, 0.5, 1), downfracs=c(0.5, 0.2)){
	# retrieve relative threshold crossing dates
	peak_ind <- which(y == max(y, na.rm=T))
	crossing_inds <- rep(NA, length(upfracs) + length(downfracs))
	i <- 1
	up_y <- y[1:peak_ind]
	down_y <- y[(peak_ind + 1):length(y)]
	if(length(upfracs) > 0){
		min_y <- min(up_y, na.rm=T)
		max_y <- max(up_y, na.rm=T)
		amp <- max_y - min_y
		for(upfrac in upfracs){
			crossing_inds[i] <- min(which(up_y >= min_y + (amp * upfrac)))
			i <- i + 1
		}
	}
	if(length(downfracs) > 0){
		min_y <- min(down_y, na.rm=T)
		max_y <- max(down_y, na.rm=T)
		amp <- max_y - min_y
		for(downfrac in downfracs){
			crossing_inds[i] <- min(which(down_y <= min_y + (amp * downfrac))) + peak_ind
			i <- i + 1
		}
	}
	return(crossing_inds)
}

#-----------------------------------------------------
GetKFPheno <- function(x, m, x_inds, lin_interp=F, ret_plot=F, Ksmooth=T, upfracs=c(0.2), downfracs=c()){
	# Kalman Filter the data and estimate phenology
	kf_filt <- dlmFilter(x, m) # get the Kalman Filter object
	if(Ksmooth){
		kf_smooth <- dlmSmooth(kf_filt) # get a Kalman Smoother object
		smoothed <- dropFirst(kf_smooth$s) # get the Kalman Smoothed estimate of x
		v <- unlist(dlmSvd2var(kf_smooth$U.S, kf_smooth$D.S))
		pl <- smoothed + qnorm(0.05, sd=sqrt(v[-1]))
		pu <- smoothed - qnorm(0.05, sd=sqrt(v[-1]))
	}else{
		v <- unlist(dlmSvd2var(kf_filt$U.C, kf_filt$D.C))
		smoothed <- dropFirst(kf_filt$m)
		pl <- smoothed + qnorm(0.05, sd=sqrt(v[-1]))
		pu <- smoothed - qnorm(0.05, sd=sqrt(v[-1]))
	}

	# interpolate gaps w/ either a spline or linearly
	ifelse(lin_interp, x_daily <- approx(x_inds, smoothed, xout=1:365)$y, x_daily <- predict(smooth.spline(x_inds, smoothed), 1:365)[[2]])
	# get the SOS estimate from the interpolated data
	# sos_estimate <- CrossingDates(x_daily, upfracs=0.2, downfracs=c())
	sos_estimate <- CrossingDates(x_daily, upfracs=upfracs, downfracs=downfracs)
	# return the data, with or without an object to plot
	if(ret_plot){
		# # get the 90% CI
		# pl <- dropFirst(kf_smooth$s) + qnorm(0.05, sd=sqrt(v[-1]))
		# pu <- dropFirst(kf_smooth$s) - qnorm(0.05, sd=sqrt(v[-1]))
		plot_obj <- list(x=x_inds, y=smoothed, pl=pl, pu=pu)
		return(list(sos_estimate, plot_obj)) # return a list w/ the SOS est and the plotting object
	}else{
		return(sos_estimate) # return just the SOS estimate
	}
}

#==========================================================================
# Processing
#==========================================================================

# Create some fake data using the GenLog model: symmetric about DOY 183
mygenlog <- function(t, A=0, K=0.5, B=0.04, M=100, S=0.8) A + ((K - A) / (1 + (S * exp(-B * (t - M)))) ^ (1 / S)) # GenLog model
doy <- 1:365
sim_evi <- mygenlog(1:183, A=0.02, K=0.45, B=0.07, M=110, S=0.3) # calculate gen logistic function for doy 1:183
sim_evi <- c(sim_evi, rev(sim_evi[1:182])) # mirror for senescence

#-----------------------------------------------------
# create sample data.frame of shifted, sampled, noisy signal
NUM_TO_DO <- 1e4
SAMPLES_PER_YEAR <- 52
sim_evi_df <- data.frame(matrix(NA, nrow=NUM_TO_DO, ncol=SAMPLES_PER_YEAR)) # initialize the data.frame
shifts <- round(rnorm(NUM_TO_DO, sd=14)) # get a sample of different seasonal shifts
shift_inds <- ShiftSample(sim_evi, SAMPLES_PER_YEAR, retInds=T)[[2]] # get the sample indices
for(i in 1:NUM_TO_DO){
	sim_evi_df[i,] <- ShiftSample(sim_evi, SAMPLES_PER_YEAR, shift=shifts[i], noise_frac=0.1, seed=i)
}

# calculate the mean and variance over all of the simulated t.s.
mean_evi <- colMeans(sim_evi_df, na.rm=T)
sd_evi <- apply(sim_evi_df, 2, sd, na.rm=T)
# var_evi <- apply(sim_evi_df, 2, "var", na.rm=T)

#-----------------------------------------------------
# Create a dlm KF to estimate EVI evolution
# m0 <- mean_evi[1] # prior estimate
m0 <- sim_evi[shift_inds[1]] # prior estimate
C0 <- max(sd_evi[1], 0.01) # prior error; need some min value or dlm can fail
FF <- 1 # observation model
V <- 0.005 # observation error; NOTE: why does this have to be so low to create realistic bounds in the KF estimates?
GG <- 1 # process model; ignored b/c time-variant parameters are given
W <- 1 # process model; ignored b/c time-variant parameters are given
JFF <- 0 # observation model is time-invariant
JV <- 0 # observation error is time-invariant
JGG <- 1 # process model is time-variant, those values are in col 1 of X
JW <- 2 # process error is time-variant, those values are in col 2 of X
X <- array(1, dim=c(SAMPLES_PER_YEAR, 2))
# calculate the process model from
rates <- sim_evi[shift_inds][2:length(shift_inds)] / sim_evi[shift_inds][1:(length(shift_inds) - 1)] # calculate rate of change from the mean EVI
# mean_smooth_evi <- predict(smooth.spline(shift_inds[!is.na(mean_evi)], mean_evi[!is.na(mean_evi)], shift_inds))[[2]]
# rates <- mean_smooth_evi[2:length(mean_smooth_evi)] / mean_smooth_evi[1:(length(mean_smooth_evi) - 1)] # calculate rate of change from the mean EVI
X[, 1] <- c(rates, rates[length(rates)]) # use rates to define time-varying process model in dlm
# NOTE: is this right? is it the residual process error, or the process parameter uncertainty?
X[, 2] <- sd_evi
sim_dlm <- dlm(m0=m0, C0=C0, FF=FF, V=V, GG=GG, W=W, JFF=JFF, JV=JV, JGG=JGG, JW=JW, X=X) # create the dlm

#-----------------------------------------------------
# Kalman Filtering examples:
x <- sim_evi_df[sample(1:NUM_TO_DO, 1),]
dlmFilt <- dlmFilter(x, sim_dlm) # filter the series with the model
v <- unlist(dlmSvd2var(dlmFilt$U.C, dlmFilt$D.C))
pl <- dropFirst(dlmFilt$m) + qnorm(0.05, sd=sqrt(v[-1]))
pu <- dropFirst(dlmFilt$m) - qnorm(0.05, sd=sqrt(v[-1]))
plot(1:length(x), x, col="seagreen")
lines(dropFirst(dlmFilt$m), col="brown")
lines(pl, col="brown", lty=2)
lines(pu, col="brown", lty=2)

# Kalman Smoothing
# my_dlm_smooth <- dlmSmooth(z, my_dlm) # smooth the series with the model
dlmSmooth <- dlmSmooth(dlmFilt) # smooth the series with the model
v <- unlist(dlmSvd2var(dlmSmooth$U.S, dlmSmooth$D.S))
pl <- dropFirst(dlmSmooth$s) + qnorm(0.05, sd=sqrt(v[-1]))
pu <- dropFirst(dlmSmooth$s) - qnorm(0.05, sd=sqrt(v[-1]))
plot(1:length(x), x, col="black", ylim=c(0, 0.8))
lines(dropFirst(dlmSmooth$s), col="darkgreen", lwd=2)
lines(pl, col="darkgreen", lty=2)
lines(pu, col="darkgreen", lty=2)

#-----------------------------------------------------
# determine the "achievable error" due to sampling and noise (no shifts)
mean_true_sos <- CrossingDates(sim_evi, upfracs=0.2, downfracs=c())
noshift_error <- rep(NA, NUM_TO_DO)
for(i in 1:NUM_TO_DO){
	tmp_ts <- ShiftSample(sim_evi, SAMPLES_PER_YEAR, shift=0, noise_frac=0.1, seed=i)
	noshift_error[i] <- GetKFPheno(tmp_ts, sim_dlm, shift_inds) - mean_true_sos
}
mean(abs(noshift_error)) # 0.5738, but should round up to +/- one day; 97% of the data are +/- 1 day

#-----------------------------------------------------
# forecast phenology for an individual time series through time and plot it
N <- sample(1:1e2, 1)
upfracs <- 1; downfracs <- c()
ylim <- quantile(sim_evi_df[N,], c(0, 1), na.rm=T) * c(0.8, 1.2)
phenometric_ests <- rep(NA, length(shift_inds))
mean_true_phenometric <- CrossingDates(sim_evi, upfracs=upfracs, downfracs=downfracs) # the true mean SOS timing from simulated signal
true_phenometric <- mean_true_phenometric + shifts[N] # the sample's true SOS, taking into account the shift

saveGIF(
	for(i in 1:length(shift_inds)){
		x <- sim_evi_df[N,]
		x[i:length(x)] <- NA
		# tmp <- GetKFSOS(x, sim_dlm, shift_inds, ret_plot=T)
		tmp <- GetKFPheno(x, sim_dlm, shift_inds, upfracs=c(1), downfracs=c(), ret_plot=T)
		phenometric_ests[i] <- tmp[[1]]

		# plot
		# par(mar=c(3,3,1,1), oma=rep(1,4), ask=T)
		par(mar=c(3,3,1,1), oma=rep(1,4))
		layout(matrix(1:2, nrow=1))
		plot(tmp[[2]]$x, tmp[[2]]$y, type="n", ylim=ylim, xlab="", ylab="")
		mtext("EVI", side=2, line=2, cex=1.25)
		mtext("DOY", side=1, line=2.2, cex=1.25)
		polygon(x=c(shift_inds, rev(shift_inds)), y=c(tmp[[2]]$pu, rev(tmp[[2]]$pl)), col="lightgrey", border=NA)
		points(shift_inds, x, pch=16, col=1)
		points(tmp[[2]]$x, tmp[[2]]$y, pch=4, col=4, type="b")
		abline(v=tmp[[1]], lwd=2, col=2, lty=2)
		text(tmp[[1]], ylim[2], labels=tmp[[1]], cex=1.7, col=2, adj=c(-0.1, 0.7))
		legend("topright", legend=c("Data", "Forecast", "PhenoMetric Est", "95% CI"), pch=c(16, 4, NA, 15), col=c(1, 4, 2, "lightgrey"), lty=c(NA, NA, 2, NA), bty="o", bg="white", pt.cex=c(1,1,1,2))
		title("KF PhenoMetric Forecast",cex=1.75)

		# plot error
		error_ylim <- c(-30, 30)
		# error_ylim <- c(min(phenometric_ests) - true_phenometric, max(phenometric_ests - true_phenometric)) + c(-2, 2)
		# error_ylim <- c(min(c(phenometric_ests[1] - true_phenometric, 0)), max(c(phenometric_ests[1] - true_phenometric, 0))) + c(-1, 2)
		# plot(shift_inds, rep(1, length(shift_inds)), type="n", ylim=error_ylim)
		# points(shift_inds[1:i], sos_ests[1:i] - true_sos, type="b")
		plot(shift_inds - true_phenometric, rep(1, length(shift_inds)), type="n", ylim=error_ylim)
		polygon(x=c(rep(shift_inds[1] - true_phenometric, 2), rep(shift_inds[length(shift_inds)] - true_phenometric, 2)), y=c(1, -1, -1, 1), border=NA, col="pink", density=NA)
		points(shift_inds[1:i] - true_phenometric, phenometric_ests[1:i] - true_phenometric, type="b", pch=16, col=1, lwd=2)
		abline(v=0, lty=3, lwd=1.5)
		abline(h=0, lty=3, lwd=1.5)
		# abline(h=c(-1, 1), lty=2, col=2)
		legend("topright", legend=c("Error", "Sample/Noise Error"), col=c(1, "pink"), pch=c(NA, 15), lty=c(1, NA), bty="o", bg="white", pt.cex=c(1, 2))
		mtext("Estimate - True PhenoMetric (DOY)", side=2, line=2, cex=1.25)
		mtext("Days from True PhenoMetric", side=1, line=2.2, cex=1.25)
		title("PhenoMetric Error Evolution",cex=1.75)
	},
movie.name="KF_Pheno_forecast.gif", interval=0.33, ani.height=600, ani.width=1200)

#-----------------------------------------------------
# forecast phenology for an individual time series through time and plot it
# create an oversize array to store the error-by-day
error_doy_add <- 400
error_mat <- matrix(NA, nrow=NUM_TO_DO, ncol=(error_doy_add * 2) + 1)

mean_true_sos <- CrossingDates(sim_evi, upfracs=0.2, downfracs=c()) # the true mean SOS timing from simulated signal
system.time( # takes about 20 minutes to do 1e4
	# for(N in 1:1000){
	for(N in 1:NUM_TO_DO){
		true_sos <- mean_true_sos + shifts[N] # the sample's true SOS, taking into account the shift
		tmp_sos_errors <- rep(NA, length(shift_inds))
		for(i in 1:length(shift_inds)){
			x <- sim_evi_df[N,]
			x[i:length(x)] <- NA
			tmp_sos_errors[i] <- GetKFPheno(x, sim_dlm, shift_inds, upfracs=0.2, downfracs=c()) - true_sos
		}
		error_mat[N, shift_inds - true_sos + error_doy_add] <- tmp_sos_errors
		# doy_diff <- shift_inds - true_sos
	}
)
sd_error <- apply(abs(error_mat), 2, sd, na.rm=T)
mean_absolute_error <- colMeans(abs(error_mat), na.rm=T)

plot(1:dim(error_mat)[2] - error_doy_add, colMeans(abs(error_mat), na.rm=T), xlim=c(-175, 210), xlab="", ylab="", type="n")
polygon(x=c(-175, 210, 210, -175), y=c(1, 1, -1, -1), border=NA, col="pink", density=NA)
# polygon(
# 	x=c(1:dim(error_mat)[2] - error_doy_add, rev(1:dim(error_mat)[2] - error_doy_add)),
# 	y=c(mean_absolute_error + (2 * sd_error), rev(mean_absolute_error - (2 * sd_error))),
# 	col="lightgrey", border=0
# )
points(1:dim(error_mat)[2] - error_doy_add, mean_absolute_error, xlim=c(-175, 210), pch=16)

abline(v=0, lty=2)
abline(h=0, lty=2)
legend("topright", legend=c("Error", "Sample/Noise Error"), col=c(1, "pink"), pch=c(16, 15), bty="o", bg="white", pt.cex=c(1, 2))
mtext("abs(Estimate - True SOS) (DOY)", side=2, line=2, cex=1.25)
mtext("Days from SOS", side=1, line=2.2, cex=1.25)
title("Mean Absolute Forecast Error",cex=1.75)


#-----------------------------------------------------
# forecast phenology for an individual time series through time and plot it
# create an oversize array to store the error-by-day
error_doy_add <- 400
sos_error_mat <- matrix(NA, nrow=NUM_TO_DO, ncol=(error_doy_add * 2) + 1)
halfspring_error_mat <- matrix(NA, nrow=NUM_TO_DO, ncol=(error_doy_add * 2) + 1)
peak_error_mat <- matrix(NA, nrow=NUM_TO_DO, ncol=(error_doy_add * 2) + 1)
halfautumn_error_mat <- matrix(NA, nrow=NUM_TO_DO, ncol=(error_doy_add * 2) + 1)
eos_error_mat <- matrix(NA, nrow=NUM_TO_DO, ncol=(error_doy_add * 2) + 1)

mean_true_sos <- CrossingDates(sim_evi, upfracs=0.2, downfracs=c()) # the true mean SOS timing from simulated signal
mean_true_halfspring <- CrossingDates(sim_evi, upfracs=0.5, downfracs=c())
mean_true_peak <- CrossingDates(sim_evi, upfracs=1, downfracs=c())
mean_true_halfautumn <- CrossingDates(sim_evi, upfracs=c(), downfracs=c(0.5))
mean_true_eos <- CrossingDates(sim_evi, upfracs=c(), downfracs=c(0.2))
system.time(
	for(N in 1:1e2){
	# for(N in 1:NUM_TO_DO){
		true_sos <- mean_true_sos + shifts[N] # the sample's true SOS, taking into account the shift
		true_halfspring <- mean_true_halfspring + shifts[N]
		true_peak <- mean_true_peak + shifts[N]
		true_halfautumn <- mean_true_halfautumn + shifts[N]
		true_eos <- mean_true_eos + shifts[N]

		tmp_sos_errors <- rep(NA, length(shift_inds))
		tmp_halfspring_errors <- tmp_sos_errors
		tmp_peak_errors <- tmp_sos_errors
		tmp_halfautumn_errors <- tmp_sos_errors
		tmp_eos_errors <- tmp_sos_errors

		for(i in 1:length(shift_inds)){
			x <- sim_evi_df[N,]
			x[i:length(x)] <- NA

			tmp_sos_errors[i] <- GetKFPheno(x, sim_dlm, shift_inds, upfracs=c(0.2), downfrac=c()) - true_sos
			tmp_halfspring_errors[i] <- GetKFPheno(x, sim_dlm, shift_inds, upfracs=c(0.5), downfrac=c()) - true_halfspring
			tmp_peak_errors[i] <- GetKFPheno(x, sim_dlm, shift_inds, upfracs=c(1), downfrac=c()) - true_peak
			tmp_halfautumn_errors[i] <- GetKFPheno(x, sim_dlm, shift_inds, upfracs=c(), downfrac=c(0.5)) - true_halfautumn
			tmp_eos_errors[i] <- GetKFPheno(x, sim_dlm, shift_inds, upfracs=c(), downfrac=c(0.2)) - true_eos
		}
		sos_error_mat[N, shift_inds - true_sos + error_doy_add] <- tmp_sos_errors
		halfspring_error_mat[N, shift_inds - true_halfspring + error_doy_add] <- tmp_halfspring_errors
		peak_error_mat[N, shift_inds - true_peak + error_doy_add] <- tmp_peak_errors
		halfautumn_error_mat[N, shift_inds - true_halfautumn + error_doy_add] <- tmp_halfautumn_errors
		eos_error_mat[N, shift_inds - true_eos + error_doy_add] <- tmp_eos_errors
		# doy_diff <- shift_inds - true_sos
	}
)
mean_absolute_sos_error <- colMeans(abs(sos_error_mat), na.rm=T)
mean_absolute_halfspring_error <- colMeans(abs(halfspring_error_mat), na.rm=T)
tmp <- peak_error_mat
tmp[is.infinite(tmp)] <- NA
# mean_absolute_peak_error <- colMeans(abs(peak_error_mat), na.rm=T)
mean_absolute_peak_error <- colMeans(abs(tmp), na.rm=T)
mean_absolute_halfautumn_error <- colMeans(abs(halfautumn_error_mat), na.rm=T)
mean_absolute_eos_error <- colMeans(abs(eos_error_mat), na.rm=T)

# plot(1:dim(error_mat)[2] - error_doy_add, colMeans(abs(error_mat), na.rm=T), xlim=c(-175, 210), xlab="", ylab="", type="n")
plot(1:dim(sos_error_mat)[2] - error_doy_add, colMeans(abs(sos_error_mat), na.rm=T), xlim=c(-250, 210), xlab="", ylab="", type="n", ylim=c(0, 45))
polygon(x=c(-250, 210, 210, -250), y=c(1, 1, -1, -1), border=NA, col="pink", density=NA)
# polygon(
# 	x=c(1:dim(error_mat)[2] - error_doy_add, rev(1:dim(error_mat)[2] - error_doy_add)),
# 	y=c(mean_absolute_error + (2 * sd_error), rev(mean_absolute_error - (2 * sd_error))),
# 	col="lightgrey", border=0
# )
# points(1:dim(error_mat)[2] - error_doy_add, mean_absolute_error, xlim=c(-175, 210), pch=16)
points(1:dim(sos_error_mat)[2] - error_doy_add, mean_absolute_sos_error, pch=16, col=1)
points(1:dim(halfspring_error_mat)[2] - error_doy_add, mean_absolute_halfspring_error, pch=16, col=2)
points(1:dim(peak_error_mat)[2] - error_doy_add, mean_absolute_peak_error, pch=16, col=3)
points(1:dim(halfautumn_error_mat)[2] - error_doy_add, mean_absolute_halfautumn_error, pch=16, col=4)
points(1:dim(eos_error_mat)[2] - error_doy_add, mean_absolute_eos_error, pch=16, col=5)

abline(v=0, lty=2)
abline(h=0, lty=2)
# legend("topright", legend=c("Error", "Sample/Noise Error"), col=c(1, "pink"), pch=c(16, 15), bty="o", bg="white", pt.cex=c(1, 2))
legend("topright", legend=c("SOS MAE(t)", "1/2 Spring MAE(t)", "Peak MAE(t)", "1/2 Autumn MAE(t)", "EOS MAE(t)", "Sample/Noise Error"), col=c(1:5, "pink"), pch=c(rep(16, 5), 15), bty="o", bg="white", pt.cex=c(rep(1, 5), 2))
mtext("abs(Estimate - True SOS) (DOY)", side=2, line=2, cex=1.25)
mtext("t (Days from event)", side=1, line=2.2, cex=1.25)
title("Mean Absolute Forecast Error",cex=1.75)


as.numeric(cut(0:error_doy_add, seq(0, error_doy_add + 7, by=7), include.lowest=T))
as.numeric(cut(error_doy_add:0, seq(0, error_doy_add + 7, by=7), include.lowest=T))
