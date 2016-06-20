# example of using dlmSeas()
# from: http://robjhyndman.com/talks/ABS3.R
f <- function(x, period=12) sin((2*pi*x)/period)

buildfun <- function(x, period=12){
  mod <- dlmModPoly(1) + dlmModSeas(period)
  V(mod) <- exp(x[1])
  diag(W(mod))[1:(period - 1)] <- exp(x[2:period])
  return(mod)
}

# plot.forecast <- function(fc, colmain="#2171B5", colerr="#BDD7E7", colx="#636363", ...){
plot.forecast <- function(fc, color="Grey", ...){
  mycols <- c("#636363", "#DE2D26", "#3182BD", "#31A354", "#756BB1", "#CCCCCC", "#FCAE91", "#BDD7E7", "#BAE4B3", "#CBC9E2")
  # mycols <- c(brewer.pal(5, "Greys")[[4]], brewer.pal(5, "Reds")[[4]], brewer.pal(5, "Blues")[[4]], brewer.pal(5, "Greens")[[4]], brewer.pal(5, "Purples")[[4]], brewer.pal(5, "Greys")[[2]], brewer.pal(5, "Reds")[[2]], brewer.pal(5, "Blues")[[2]], brewer.pal(5, "Greens")[[2]], brewer.pal(5, "Purples")[[2]])
  if(color=="Blue"){
    colmain <- mycols[3]
    colerr <- mycols[8]
    colx <- mycols[1]
  }else if(color=="Red"){
    colmain <- mycols[2]
    colerr <- mycols[7]
    colx <- mycols[1]
  }else if(color=="Green"){
    colmain <- mycols[4]
    colerr <- mycols[9]
    colx <- mycols[1]
  }else if(color=="Purple"){
    colmain <- mycols[5]
    colerr <- mycols[10]
    colx <- mycols[1]
  }else if(color=="Grey"){
    colmain <- mycols[1]
    colerr <- mycols[6]
    colx <- mycols[1]
  }

  ylim <- c(min(fc$lower, na.rm=T), max(fc$upper, na.rm=T)) * c(0.9, 1.1)
  plot(fc$t, fc$mean, type="n", lwd=2, col=2, lty=2, ylim=ylim, xlab="t", ylab="y", ...)
  # polygon(x=c(1:length(fc$upper), rev(1:length(fc$upper))), y=c(fc$upper, rev(fc$lower)), border=NA, col=colerr)
  polygon(x=c(fc$t, rev(fc$t)), y=c(fc$upper, rev(fc$lower)), border=NA, col=colerr)
  points(fc$t, fc$mean, type="l", lwd=2, col=colmain, lty=1)
  points(fc$t, fc$x, type="p", pch=16, cex=0.75, col=colx)
  legend("top", legend=c("data", "filtered", "95% CI"), col=c(colx, colmain, colerr), pch=c(16, NA, 15), lty=c(NA, 1, NA), pt.cex=c(1, 1, 1.5), bty="n", horiz=T)
}

#-------------------------------------------------------------
mygenlog <- function(t, A=0, K=0.5, B=0.04, M=100, S=0.8) A + ((K - A) / (1 + (S * exp(-B * (t - M)))) ^ (1 / S)) # GenLog model
doy <- 1:365
sim_evi <- mygenlog(1:183, A=0.02, K=0.45, B=0.07, M=110, S=0.3) # calculate gen logistic function for doy 1:183
sim_evi <- c(sim_evi, rev(sim_evi[1:182])) # mirror for senescence

# c1: midpoint of greenup
# c2: midpoint of greendown
# M: maturity/transition date
doublesigmoid <- function(t, params) return(0.5 * (tanh(((t - params$M) - params$c1) / params$w1) - tanh(((t - params$M) - params$c2) / params$w2)))
# t <- seq(1, 365, len=1e3)
# params <- list(c1=-64, c2=46, w1=0.4, w2=15, M=182)
# plot(t, doublesigmoid(t, params), type="l")

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

params <- list(
  midpoint=c(180, 4),
  gup_duration=c(120, 30),
  gdown_duration=c(120, 15),
  gup_rate=c(7, sd=3.5),
  gdown_rate=c(20, 3.5))

numyears <- 10
fakeEVI <- SimSVI(numyears, params)
plot(fakeEVI, type="l")
abline(v=seq(1, (365 * numyears) + 1, by=365), lty=3, col="grey")

#-------------------------------------------------------------
MakeSeasonalDLM <- function(period, numstates, sensors=1, statecov=F, default_process_variance=1e6, default_process_covariance=1e6, default_observation_variance=1e6, initial_state=0, initial_variance=1e6){
  zeros <- function(x) matrix(0, nrow=x, ncol=x) # create square matrix of zeros with dimension x
  circ <- function(x) rbind(cbind(rep(0, x - 1), diag(x - 1)), c(1, rep(0, x - 1))) # create a circulant matrix of period x
  Grow <- function(state, period, numstates) do.call("cbind", c(rep(list(zeros(period)), state - 1), list(circ(period)), rep(list(zeros(period)), numstates - state))) # create a row of the KF state transition matrix (G)
  Frow <- function(state, period, numstates) c(rep(0, (state - 1) * period), 1, rep(0, period - 1), rep(0, (numstates - state) * period))

  # create the KF state transition matrix (GG)
  GG <- do.call("rbind", lapply(1:numstates, Grow, period=period, numstates=numstates))

  # create the process error matrix (W) w/ or w/o covariance among seasonal factors of different states
  if(statecov){
    # encode default covariance between seasonal factors of different states
    W <- do.call(cbind, rep(list(do.call("rbind", rep(list(diag(period)), numstates))), numstates)) * default_process_covariance
    diag(W) <- default_process_variance # set seasonal factor variances
  }else{
    W <- diag(period * numstates) * default_process_variance # set seasonal factor variances
  }

  # create the observation matrix (FF)
  trash <- ifelse(is.list(sensors), sensor_list <- unlist(sensors), sensor_list <- rep(1:numstates, sensors))
  FF <- do.call("rbind", lapply(sensor_list, Frow, period=period, numstates=numstates))

  # create observation error matrix (V)
  V <- diag(dim(FF)[1]) * default_observation_variance

  # create the initial state vector (m0)
  m0 <- array(rep(initial_state, numstates * period))

  # create the initial state covariance matrix (C0)
  C0 <- diag(numstates * period) * initial_variance

  # create and return the dlm
  return(dlm(m0=m0, C0=C0, GG=GG, FF=FF, V=V, W=W))
}

#-------------------------------------------------------------
SeasonalDLMBuildFun <- function(x, period, numstates, sensors=1, statecov=F, default_process_variance=1e6, default_process_covariance=1e6, default_observation_variance=1e6, initial_state=0, initial_variance=1e6){
  tmp_dlm <- MakeSeasonalDLM(period=period, numstates=numstates, sensors=sensors, statecov=statecov, default_process_variance=default_process_variance, default_process_covariance=default_process_covariance, default_observation_variance=default_observation_variance, initial_state=initial_state, initial_variance=initial_variance)

  # update the process error matrix
  diag(W(tmp_dlm)) <- exp(x[1:(numstates * period)])

  # update the observation error matrix
  trash <- ifelse(is.list(sensors), sensornum <- length(unlist(sensors)), sensornum <- sensors)
  diag(V(tmp_dlm)) <- exp(x[(numstates * period + 1):(numstates * period + sensornum)])

  return(tmp_dlm)
}

#-------------------------------------------------------------
# period <- 12 # about 5 seconds to fit
period <- 26 # about 55 seconds to fit
# period <- 52
numstates <- 1
sensors <- 1

# create some fake data
# t <- 1:(period*9)
# y <- f(t, period=period) + rnorm(length(t), sd=0.1)
# y[y < 0] <- 0
sample_inds <- round(seq(1, (365 * numyears), by=365 / period))
y <- fakeEVI[sample_inds]
y_noise <- y + rnorm(length(y), sd=0.03)
plot(sample_inds, fakeEVI[sample_inds], col=1, pch=16, cex=0.5)
points(1:length(fakeEVI), fakeEVI, type="l", col="pink", lty=3)
abline(v=seq(1, (365 * numyears) + 1, by=365), lty=1, col="grey")


# fit the DLM
trash <- ifelse(is.list(sensors), initial_params <- c(rep(0, period), rep(0, length(unlist(sensors)))), initial_params <- c(rep(0, period), rep(0, sensors * numstates)))
# system.time(fit <- dlmMLE(y, parm=initial_params, build=SeasonalDLMBuildFun, period=period, numstates=numstates, sensors=sensors))
# mySeasonalDLM <- SeasonalDLMBuildFun(fit$par, period=period, numstates=numstates, sensors=sensors)
system.time(fit <- dlmMLE(y_noise, parm=initial_params, build=SeasonalDLMBuildFun, period=period, numstates=numstates, sensors=sensors))
mySeasonalDLM <- SeasonalDLMBuildFun(fit$par, period=period, numstates=numstates, sensors=sensors)


# Filt <- dlmFilter(y, mod=mySeasonalDLM)
# Smooth <- dlmSmooth(y, mod=mySeasonalDLM)
# Fore <- dlmForecast(Filt, nAhead=12*3)
# fsd <- sqrt(unlist(Fore$Q))
# pl <- Fore$f + qnorm(0.05, sd = fsd)
# pu <- Fore$f + qnorm(0.95, sd = fsd)
# fc <- list(mean=Fore$f, lower=pl, upper=pu, x=y, level=90)

oldW <- W(mySeasonalDLM)
# W(mySeasonalDLM) <- W(mySeasonalDLM) *1e2
# W(mySeasonalDLM) <- oldW

# oldV <- V(mySeasonalDLM)
# V(mySeasonalDLM) <- V(mySeasonalDLM) * 1e-1

oldm0 <- m0(mySeasonalDLM)
# m0(mySeasonalDLM) <- colMeans(matrix(y, ncol=26, byrow=T)) # start out with mean seasonal factors

W <- W(mySeasonalDLM); V <- V(mySeasonalDLM)
FF <- FF(mySeasonalDLM); GG <- GG(mySeasonalDLM)
m0 <- m0(mySeasonalDLM); C0 <- c(mySeasonalDLM)
save(y, y_noise, W, V, GG, FF, m0, C0, Filt, Smooth, file="~/Desktop/ForLuis.Rdata")

#----------------------------------------------------------------
# compute and plot the mean and 95% CI of the filtering distribution
mycols <- c(brewer.pal(5, "Greys")[[4]], brewer.pal(5, "Reds")[[4]], brewer.pal(5, "Blues")[[4]], brewer.pal(5, "Greens")[[4]], brewer.pal(5, "Purples")[[4]], brewer.pal(5, "Greys")[[2]], brewer.pal(5, "Reds")[[2]], brewer.pal(5, "Blues")[[2]], brewer.pal(5, "Greens")[[2]], brewer.pal(5, "Purples")[[2]])
# Filt <- dlmFilter(y, mod=mySeasonalDLM)
# y_noise <- y + rnorm(length(y), sd=0.03)
Filt <- dlmFilter(y_noise, mod=mySeasonalDLM)
# Fore <- dlmForecast(mySeasonalDLM, nAhead=26)
filt.m <- dropFirst(Filt$m[,1])
filt.se <- sqrt(unlist(lapply(dropFirst(dlmSvd2var(Filt$U.C, Filt$D.C)), function(x) x[1,1])))
filt.pl <- filt.m + qnorm(0.05, sd=filt.se)
filt.pu <- filt.m + qnorm(0.95, sd=filt.se)
# plotting
par(mfrow=c(2,1), mar=c(2,4,2,0), oma=rep(1,4))
plot.forecast(list(mean=filt.m, upper=filt.pu, lower=filt.pl, x=y_noise, t=sample_inds))
title("dlmFilter() results with dlmMLE() fitted model")
abline(v=seq(1, (365 * numyears) + 1, by=365), lty=3, col=mycols[6])
par(mar=c(3,4,1,0))
plot(sample_inds, filt.se, type="l", col=mycols[4], lwd=2, xlab="t", ylab="filtered s.e.")
abline(v=seq(1, (365 * numyears) + 1, by=365), lty=3, col=mycols[6])

#----------------------------------------------------------------
# compute and plot the mean and 95% CI of the smoothing distribution
# Smooth <- dlmSmooth(y, mod=mySeasonalDLM)
Smooth <- dlmSmooth(y_noise, mod=mySeasonalDLM)
smooth.m <- dropFirst(Smooth$s[,1])
smooth.se <- sqrt(unlist(lapply(dropFirst(dlmSvd2var(Smooth$U.S, Smooth$D.S)), function(x) x[1,1])))
smooth.pl <- smooth.m + qnorm(0.05, sd=smooth.se)
smooth.pu <- smooth.m + qnorm(0.95, sd=smooth.se)
# plotting
par(mfrow=c(2,1), mar=c(2,4,2,0), oma=rep(1,4))
plot.forecast(list(mean=smooth.m, upper=smooth.pu, lower=smooth.pl, x=y_noise, t=sample_inds))
title("dlmSmooth() results with dlmMLE() fitted model")
abline(v=seq(1, (365 * numyears) + 1, by=365), lty=3, col=mycols[6])
par(mar=c(3,4,1,0))
plot(sample_inds, smooth.se, type="l", col=mycols[4], lwd=2, xlab="t", ylab="smoothed s.e.")
abline(v=seq(1, (365 * numyears) + 1, by=365), lty=3, col=mycols[6])

#
# system.time(fit <- dlmMLE(y, parm=rep(0, period), build=buildfun, period=period))
# seas_dlm <- buildfun(fit$par, period=period)
#
# ySmooth <- dlmSmooth(y, mod = seas_dlm)
# x <- cbind(y, dropFirst(ySmooth$s[,c(1,2)]))
# plot(t, y, type="l", lwd=2, col=rgb(0.35, 0.35, 0.35)); points(t, x[,2] + x[,3], type="l", col=2, lty=2)
#
# Filt <- dlmFilter(y, mod=seas_dlm)
# Fore <- dlmForecast(Filt, nAhead=20)
# fsd <- sqrt(unlist(Fore$Q))
# pl <- Fore$f + qnorm(0.05, sd = fsd)
# pu <- Fore$f + qnorm(0.95, sd = fsd)
# fc <- list(mean=Fore$f, lower=pl, upper=pu, x=y, level=90)
#
# plot.forecast(fc, main="This is a Test")

num_extra_years <- 4
y_tmp <- c(y_noise, rep(NA, num_extra_years * period))
# sample_inds_tmp <- round(seq(1, (365 * (numyears + 1 + num_extra_years)), by=365 / period))
sample_inds_tmp <- round(seq(1, (365 * (numyears + num_extra_years)), by=365 / period))
Filt <- dlmFilter(y_tmp, mod=mySeasonalDLM)
# Fore <- dlmForecast(mySeasonalDLM, nAhead=26)
filt.m <- dropFirst(Filt$m[,1])
filt.se <- sqrt(unlist(lapply(dropFirst(dlmSvd2var(Filt$U.C, Filt$D.C)), function(x) x[1,1])))
filt.pl <- filt.m + qnorm(0.05, sd=filt.se)
filt.pu <- filt.m + qnorm(0.95, sd=filt.se)
# plotting
par(mfrow=c(2,1), mar=c(2,4,2,0), oma=rep(1,4))
plot.forecast(list(mean=filt.m, upper=filt.pu, lower=filt.pl, x=y_tmp, t=sample_inds_tmp))
title("dlmFilter() results with dlmMLE() fitted model")
abline(v=seq(1, (365 * numyears) + 1, by=365), lty=3, col=mycols[6])
par(mar=c(3,4,1,0))
plot(sample_inds, filt.se, type="l", col=mycols[4], lwd=2, xlab="t", ylab="filtered s.e.")
abline(v=seq(1, (365 * numyears) + 1, by=365), lty=3, col=mycols[6])
