# load the dlm library and the Kalman Filtering functions
library(dlm)
source("https://raw.github.com/jm-gray/pixel-forge/master/KalmanFiltering/KF_functions.R")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# example: single state, single sensor
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
N <- 500
y <- cumsum(rnorm(N)) # create random walk data

dlm1 <- MakeMultiDLM()
W(dlm1) <- 10; V(dlm1) <- 100
# W(dlm1) <- 0.1; V(dlm1) <- 100 # favor the model
# W(dlm1) <- 1000; V(dlm1) <- 10 # favor the observations

filt <- dlmFilter(y, dlm1)
filt_m <- dropFirst(filt$m)
filt_se <- sqrt(unlist(dropFirst(dlmSvd2var(filt$U.C, filt$D.C))))
PlotForecast(filt_m, filt_se, signal=y)

# NOTE: I do not entirely understand this fitting procedure, and it doesn't yield
# results that make sense to me yet. I'll leave it here for the sake of a demo
buildfun <- function(x, dlm){
  W(dlm) <- exp(x[1])
  V(dlm) <- exp(x[2])
  return(dlm)
}
fit <- dlmMLE(y, parm=c(100,100), build=buildfun, dlm=dlm1)
dlm1_fit <- buildfun(fit$par, dlm1)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# example: 3 states, 1 sensor
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
N <- 200
num_states <- 3
y <- do.call(cbind, lapply(1:num_states, function(x, N) cumsum(rnorm(N)), N=N))

dlm2 <- MakeMultiDLM(num_states=3)
diag(W(dlm2)) <- 10; diag(V(dlm2)) <- 100

filt <- dlmFilter(y, dlm2)
layout(matrix(1:num_states, ncol=1)); par(mar=rep(1, 4))
for(i in 1:num_states){
  filt_m <- dropFirst(filt$m)[, i]
  filt_se <- sqrt(unlist(lapply(dropFirst(dlmSvd2var(filt$U.C, filt$D.C)), function(x, i) return(x[i, i]), i=i)))
  PlotForecast(filt_m, filt_se, signal=y[,i])
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# example: 3 states each observed by 2 different sensors
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
N <- 250
num_states <- 3
sensors <- 2
y <- do.call(cbind, lapply(1:num_states, function(x, N) cumsum(rnorm(N)), N=N)) # the true signal
y1 <- y + array(rnorm(prod(dim(y)), sd=0.5), dim=dim(y)) # y + noise
y2 <- y + array(rnorm(prod(dim(y))), dim=dim(y)) # y + more noise
y <- cbind(y1, y2)
y[sample(1:prod(dim(y)), 400)] <- NA # add some missing data

dlm3 <- MakeMultiDLM(num_states=num_states, sensors=sensors)
diag(W(dlm3)) <- 10; diag(V(dlm3)) <- 100

filt <- dlmFilter(y, dlm3)
layout(matrix(1:num_states, ncol=1)); par(mar=rep(1, 4))
if(is.list(sensors)){
  sensor_list <- unlist(sensors)
}else{
  sensor_list <- rep(1:num_states, sensors)
}
for(i in 1:num_states){
  filt_m <- dropFirst(filt$m)[, i]
  filt_se <- sqrt(unlist(lapply(dropFirst(dlmSvd2var(filt$U.C, filt$D.C)), function(x, i) return(x[i, i]), i=i)))
  PlotForecast(filt_m, filt_se, signal=split(y[, sensor_list == i], col(y[, sensor_list == i])))
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# example: 3 states, 3 sensors, but not every state observed by every sensor
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
N <- 250
num_states <- 3
sensors <- list(1:3, 1:2, 2:3)
y <- do.call(cbind, lapply(1:num_states, function(x, N) cumsum(rnorm(N)), N=N)) # the true signal
y1 <- y + array(rnorm(prod(dim(y)), sd=0.5), dim=dim(y)) # y + noise
y2 <- y[,1:2] + array(rnorm(prod(dim(y[,1:2]))), dim=dim(y[,1:2])) # y + more noise
y3 <- y[,2:3] + array(rnorm(prod(dim(y[,2:3])), sd=1.5), dim=dim(y[,2:3])) # y + even more noise
y <- cbind(y1, y2, y3)

dlm4 <- MakeMultiDLM(num_states=num_states, sensors=sensors)
diag(W(dlm4)) <- 10; diag(V(dlm4)) <- 100

filt <- dlmFilter(y, dlm4)
layout(matrix(1:num_states, ncol=1)); par(mar=rep(1, 4))
if(is.list(sensors)){
  sensor_list <- unlist(sensors)
}else{
  sensor_list <- rep(1:num_states, sensors)
}
for(i in 1:num_states){
  filt_m <- dropFirst(filt$m)[, i]
  filt_se <- sqrt(unlist(lapply(dropFirst(dlmSvd2var(filt$U.C, filt$D.C)), function(x, i) return(x[i, i]), i=i)))
  PlotForecast(filt_m, filt_se, signal=split(y[, sensor_list == i], col(y[, sensor_list == i])))
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# example: 1 state, 1 sensor, time-varying
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# first, create some fake EVI data at weekly resolution for several years
params <- list(midpoint=c(180, 4), gup_duration=c(120, 30), gdown_duration=c(120, 15), gup_rate=c(7, sd=3.5), gdown_rate=c(20, 3.5))
start_year <- 2001
numyears <- 4
sensor1_sd <- 0.04
sensor1_samples_per_year <- 52
fakeEVI <- SimSVI(numyears, params)
dates <- as.Date(unlist(lapply(start_year:(start_year + numyears - 1), function(x) as.Date(paste(x, 1:365, sep="-"), format="%Y-%j"))), origin="1970-1-1")
sample_dates <- as.Date(unlist(lapply(start_year:(start_year + numyears - 1), function(x) as.Date(paste(x, round(seq(1,365,len=sensor1_samples_per_year)), sep="-"), format="%Y-%j"))), origin="1970-1-1")
y1 <- fakeEVI[match(sample_dates, dates)] + rnorm(length(sample_dates), sd=sensor1_sd)
# doys <- as.integer(strftime(sample_dates, format="%j"))

# define the process model prior as the mean rate of change over a thousand random time series
tmp <- matrix(SimSVI(1e3, params), nrow=1e3, byrow=T) # generate a lot of EVI curves
prior_mean <- unlist(lapply(unique(as.integer(strftime(sample_dates,format="%j"))), function(x) mean(tmp[,x], na.rm=T)))
prior_var <- unlist(lapply(unique(as.integer(strftime(sample_dates,format="%j"))), function(x) var(tmp[,x], na.rm=T)))
prior_rate <- prior_mean[2:length(prior_mean)] / prior_mean[1:(length(prior_mean) - 1)] # calculate rate of mean curve
prior_rate <- c(prior_rate, prior_rate[length(prior_rate)]) # append the last value so we preserve the proper size

dlm5 <- MakeMultiDLM()
JGG(dlm5) <- 1 # store the time-varying values of GG in col 1 of X
JW(dlm5) <- 2 # store the time-varying values of W in col 2 of X
X(dlm5) <- matrix(c(rep(prior_rate, numyears), rep(prior_var, numyears)), ncol=2)
V(dlm5) <- sensor1_sd^2
m0(dlm5) <- prior_mean[1]
C0(dlm5) <- as.numeric(prior_var[1])

filt <- dlmFilter(y1, dlm5)
filt_m <- dropFirst(filt$m)
filt_se <- sqrt(unlist(dropFirst(dlmSvd2var(filt$U.C, filt$D.C))))
PlotForecast(filt_m, filt_se, signal=y1)

# and a forecast:
y1_sub <- c(y1[1:(52*2)], y1[52:64], rep(NA, 39))
filt <- dlmFilter(y1_sub, dlm5)
filt_m <- dropFirst(filt$m)
filt_se <- sqrt(unlist(dropFirst(dlmSvd2var(filt$U.C, filt$D.C))))
PlotForecast(filt_m, filt_se, signal=y1_sub)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# example: 1 state, 2 sensors, time-varying
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# create two time series (borrowing from the above example)
y1 <- fakeEVI[match(sample_dates, dates)] + rnorm(length(sample_dates), sd=sensor1_sd)
y1[sample(1:length(y1), round(length(y1)/5))] <- NA
sensor2_sd <- 0.06
y2 <- fakeEVI[match(sample_dates, dates)] + rnorm(length(sample_dates), sd=sensor2_sd)
y2[sample(1:length(y2), round(length(y2)/5))] <- NA
y <- cbind(y1, y2)
# doys <- as.integer(strftime(sample_dates, format="%j"))

dlm6 <- MakeMultiDLM(1, 2)
JGG(dlm6) <- 1 # store the time-varying values of GG in col 1 of X
JW(dlm6) <- 2 # store the time-varying values of W in col 2 of X
X(dlm6) <- matrix(c(rep(prior_rate, numyears), rep(prior_var, numyears)), ncol=2)
diag(V(dlm6)) <- c(sensor1_sd^2, sensor2_sd^2)
m0(dlm6) <- prior_mean[1]
C0(dlm6) <- as.numeric(prior_var[1])

filt <- dlmFilter(y, dlm6)
filt_m <- dropFirst(filt$m)
filt_se <- sqrt(unlist(dropFirst(dlmSvd2var(filt$U.C, filt$D.C))))
PlotForecast(filt_m, filt_se, signal=split(y, col(y)), ylab="EVI2", xlab="")
