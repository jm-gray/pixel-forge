library(dlm)

#------------------------------------------------
MakeMultiDLM <- function(num_states=1, sensors=1, time_varying=FALSE){
  # NOTE: time_varying option doesn't work right now, but you would just need to create a dlm
  # of the proper state/sensor configuration, and then copy GG to JGG, FF to JFF, and so on
  # then fill in X column indices in the proper places. That is, the architecture of the KF
  # that results from this function is correct

  m0 <- rep(0, num_states) # initial state vector
  C0 <- diag(1e6, num_states)

  GG <- diag(num_states) # state evolution matrix
  W <- diag(1e6, num_states) # evolution error covariance

  # create observation matrix FF
  if(is.list(sensors)){
    sensors_list <- unlist(sensors)
  }else{
    sensors_list <- rep(1:num_states, sensors)
  }
  Frow <- function(state, num_states) c(rep(0, (state - 1)), 1, rep(0, num_states - state))
  FF <- do.call(rbind, lapply(sensors_list, Frow, num_states=num_states))

  V <- diag(1e6, dim(FF)[1]) # observation error covariance

  if(time_varying){
    JFF <- 0 # column index in X of time varying FF values
    JGG <- 0 # column index in X of time varying GG values
    JW <- 0 # column index in X of time varying W values
    JV <- 0 # column index in X of time varying V values
    X <- 0 # matrix holding time varying parameters in columns
    my_dlm <- dlm(m0=m0, C0=C0, GG=GG, W=W, FF=FF, V=V, JGG=JGG, JFF=JFF, JW=JW, JV=JV, X=X)
    return(my_dlm)
  }else{
    my_dlm <- dlm(m0=m0, C0=C0, GG=GG, W=W, FF=FF, V=V)
    return(my_dlm)
  }
}

#------------------------------------------------
PlotForecast <- function(filt_m, filt_se, signal, t=NULL, conf_level=0.05, ylim=NULL, ...){
  if(is.null(t)) t <- 1:length(filt_m)
  colmain <- "#3182BD"; colerr <- "#BDD7E7"; colsignal <- "#636363"
  lower <- filt_m + qnorm(conf_level, sd=filt_se)
  upper <- filt_m + qnorm(1 - conf_level, sd=filt_se)

  if(is.null(ylim)){
    ylim <- range(c(upper, lower), na.rm=T) * c(0.9, 1.1)
    if(diff(range(c(upper, lower), na.rm=T)) > (2 * diff(range(filt_m, na.rm=T)))){
      # ylim <- range(filt_m, na.rm=T) * c(0.55, 1.5)
      # ylim <- c(quantile(lower,c(0.1)), quantile(upper,c(0.9)))
      ylim <- range(filt_m, na.rm=T) + sd(filt_m)*c(-1, 1)
    }
  }

  # ylim <- c(min(lower, na.rm=T), max(upper, na.rm=T)) * c(0.9, 1.1)
  plot(t, filt_m, type="n", lwd=2, col=2, lty=2, ylim=ylim, ...)
  polygon(x=c(t, rev(t)), y=c(upper, rev(lower)), border=NA, col=colerr)
  points(t, filt_m, type="l", lwd=2, col=colmain, lty=1)
  if(is.list(signal)){
    for(i in 1:length(signal)){
      tmp_y <- signal[[i]]
      points(t, tmp_y, type="p", pch=i, cex=0.75, col=colsignal)
    }
  }else{
    points(t, signal, type="p", pch=1, cex=0.75, col=colsignal)
  }
}

#------------------------------------------------
ApplyAnnualKF <- function(y, dlm, period=52){
  tmp_dlm <- dlm
  filt_m <- c()
  filt_se <- c()
  for(i in seq(1, length(y), by=period)){
    tmp_y <- y[i:(i + period - 1)]
    filt <- dlmFilter(tmp_y, tmp_dlm)
    filt_m <- c(filt_m, dropFirst(filt$m))
    filt_se <- c(filt_se, sqrt(unlist(dropFirst(dlmSvd2var(filt$U.C, filt$D.C)))))
    m0(tmp_dlm) <- filt_m[length(filt_m)]
    C0(tmp_dlm) <- filt_se[length(filt_m)]^2
  }
  return(list(filt_m, filt_se))
}

#------------------------------------------------
doublesigmoid <- function(t, params) return(0.5 * (tanh(((t - params$M) - params$c1) / params$w1) - tanh(((t - params$M) - params$c2) / params$w2)))

#------------------------------------------------
SimSVI <- function(numyears, params, seed=42, quiet=T){
  set.seed(seed)
  svi <- c()
  t <- 1:365
  for(i in 1:numyears){
    M <- rnorm(1, mean=params$midpoint[1], sd=params$midpoint[2])
    c1 <- -0.5 * rnorm(1, mean=params$gup_duration[1], sd=params$gup_duration[2])
    c2 <- 0.5 * rnorm(1, mean=params$gdown_duration[1], sd=params$gdown_duration[2])
    w1 <- rnorm(1, mean=params$gup_rate[1], sd=params$gup_rate[2])
    w2 <- rnorm(1, mean=params$gdown_rate[1], sd=params$gdown_rate[2])
    svi <- c(svi, doublesigmoid(t, list(M=M, c1=c1, c2=c2, w1=abs(w1), w2=abs(w2))))
    if(!quiet) print(c(M=M, c1=c1, c2=c2, w1=w1, w2=w2))
  }
  return(svi)
}

#-------------------------------------------------------------------
#-------------------------------------------------------------------
# example: single state, single sensor
#-------------------------------------------------------------------
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
filt <- dlmFilter(y, dlm1_fit)

#-------------------------------------------------------------------
#-------------------------------------------------------------------
# example: 3 states, 1 sensor
#-------------------------------------------------------------------
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

#-------------------------------------------------------------------
#-------------------------------------------------------------------
# example: 3 states each observed by 2 different sensors
#-------------------------------------------------------------------
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

#-------------------------------------------------------------------
#-------------------------------------------------------------------
# example: 3 states, 3 sensors, but not every state observed by every sensor
#-------------------------------------------------------------------
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

#-------------------------------------------------------------------
#-------------------------------------------------------------------
# example: 1 state, 1 sensor, time-varying
#-------------------------------------------------------------------
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

#-------------------------------------------------------------------
#-------------------------------------------------------------------
# example: 1 state, 2 sensors, time-varying
#-------------------------------------------------------------------
# create two time series (borrowing from the above example)
y1 <- fakeEVI[match(sample_dates, dates)] + rnorm(length(sample_dates), sd=sensor1_sd)
y1[sample(1:length(y1), round(length(y1)/5))] <- NA
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
