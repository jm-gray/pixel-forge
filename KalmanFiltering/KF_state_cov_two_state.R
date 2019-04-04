library(dlm)
library(RColorBrewer)
source("/Users/jmgray2/Documents/pixel-forge/KalmanFiltering/KF_functions.R")

# create the data
set.seed(42)
N <- 100
true_A <- predict(smooth.spline(cumsum(rnorm(N)), spar=0.65))$y
true_B <- -1 * true_A + 2
sensor_a_A <- true_A + rnorm(N, mean=0, sd=1)
sensor_a_B <- true_B + rnorm(N, mean=0, sd=1)
sensor_b_A <- true_A + rnorm(N, mean=0, sd=1)
y <- cbind(sensor_a_A, sensor_a_B, sensor_b_A)
true_y <- cbind(true_A, true_B)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Experiment 1: all data
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
num_states <- 2
sensors <- list(1:2, 1)
dlm_ex1 <- MakeMultiDLM(num_states=num_states, sensors=sensors)
diag(W(dlm_ex1)) <- 0.07
diag(V(dlm_ex1)) <- 1
filt_ex1 <- dlmFilter(y, dlm_ex1)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Experiment 2: missing half of obs for state B
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
tmp_y <- y
tmp_y[51:100, 2] <- NA
dlm_ex2 <- dlm_ex1
filt_ex2 <- dlmFilter(tmp_y, dlm_ex2)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Experiment 3: missing half of obs for state B; w/ state cov
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
dlm_ex3 <- dlm_ex1
state_cov <- cov(cbind(true_A, true_B))
# W(dlm_ex3) <- state_cov # why doesn't this work?
W(dlm_ex3) <- state_cov %*% t(state_cov) # and this does...?
# W(dlm_ex3) <- W(dlm_ex3) * 0.001 # make smoother...arbitrarily
filt_ex3 <- dlmFilter(tmp_y, dlm_ex3)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Experiment 4: missing half of obs for state B; w/ obs cov
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
dlm_ex4 <- dlm_ex1 # copy the previous model
obs_cov <- cov(y)
V(dlm_ex4) <- obs_cov
filt_ex4 <- dlmFilter(tmp_y, dlm_ex4)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# plotting
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sensor_vec <- unlist(sensors)
true_signal_color <- brewer.pal(5, "Greys")[3]
main_cols <- c(brewer.pal(5, "Blues")[3], brewer.pal(5, "Reds")[3])
err_cols <- c(brewer.pal(5, "Blues")[2], brewer.pal(5, "Reds")[2])
ylims <- list(c(-3, 7), c(-5, 8))

layout(matrix(1:8, nrow=2))
par(mar=c(2, 3, 1, 1))
for(i in 1:num_states){
    filt_m_ex1 <- dropFirst(filt_ex1$m)[, i]
    filt_se_ex1 <- sqrt(unlist(lapply(dropFirst(dlmSvd2var(filt_ex1$U.C, filt_ex1$D.C)), function(x, i) return(x[i, i]), i=i)))
    PlotForecast(filt_m_ex1, filt_se_ex1, signal=split(y[, sensor_vec == i], col(as.matrix(y[, sensor_vec == i]))), colmain=main_cols[i], colerr=err_cols[i], ylim=ylims[[i]])
    points(true_y[, i], type="l", lty=2, lwd=2, col=true_signal_color)
    if(i == 1) title("Experiment 1")
    legend("top", legend=c("Sensor A", "Sensor B", "KF Filt", "KF 95%", "Truth"), pch=c(1, 2, NA, 15, NA), lty=c(NA, NA, 1, NA, 2), lwd=c(NA, NA, 2, NA, 2), col=c("#636363", "#636363", main_cols[i], err_cols[i], true_signal_color), bg="white", horiz=F, ncol=2)
}

for(i in 1:num_states){
    filt_m_ex2 <- dropFirst(filt_ex2$m)[, i]
    filt_se_ex2 <- sqrt(unlist(lapply(dropFirst(dlmSvd2var(filt_ex2$U.C, filt_ex2$D.C)), function(x, i) return(x[i, i]), i=i)))
    PlotForecast(filt_m_ex2, filt_se_ex2, signal=split(tmp_y[, sensor_vec == i], col(as.matrix(tmp_y[, sensor_vec == i]))), colmain=main_cols[i], colerr=err_cols[i], ylim=ylims[[i]])
    points(true_y[, i], type="l", lty=2, lwd=2, col=true_signal_color)
    if(i == 1) title("Experiment 2")
    legend("top", legend=c("Sensor A", "Sensor B", "KF Filt", "KF 95%", "Truth"), pch=c(1, 2, NA, 15, NA), lty=c(NA, NA, 1, NA, 2), lwd=c(NA, NA, 2, NA, 2), col=c("#636363", "#636363", main_cols[i], err_cols[i], true_signal_color), bg="white", horiz=F, ncol=2)
}

for(i in 1:num_states){
    filt_m_ex3 <- dropFirst(filt_ex3$m)[, i]
    filt_se_ex3 <- sqrt(unlist(lapply(dropFirst(dlmSvd2var(filt_ex3$U.C, filt_ex3$D.C)), function(x, i) return(x[i, i]), i=i)))
    PlotForecast(filt_m_ex3, filt_se_ex3, signal=split(tmp_y[, sensor_vec == i], col(as.matrix(tmp_y[, sensor_vec == i]))), colmain=main_cols[i], colerr=err_cols[i], ylim=ylims[[i]])
    points(true_y[, i], type="l", lty=2, lwd=2, col=true_signal_color)
    if(i == 1) title("Experiment 3")
    legend("top", legend=c("Sensor A", "Sensor B", "KF Filt", "KF 95%", "Truth"), pch=c(1, 2, NA, 15, NA), lty=c(NA, NA, 1, NA, 2), lwd=c(NA, NA, 2, NA, 2), col=c("#636363", "#636363", main_cols[i], err_cols[i], true_signal_color), bg="white", horiz=F, ncol=2)
}

for(i in 1:num_states){
    filt_m_ex4 <- dropFirst(filt_ex4$m)[, i]
    filt_se_ex4 <- sqrt(unlist(lapply(dropFirst(dlmSvd2var(filt_ex4$U.C, filt_ex4$D.C)), function(x, i) return(x[i, i]), i=i)))
    PlotForecast(filt_m_ex4, filt_se_ex4, signal=split(tmp_y[, sensor_vec == i], col(as.matrix(tmp_y[, sensor_vec == i]))), colmain=main_cols[i], colerr=err_cols[i], ylim=ylims[[i]])
    points(true_y[, i], type="l", lty=2, lwd=2, col=true_signal_color)
    if(i == 1) title("Experiment 4")
    legend("top", legend=c("Sensor A", "Sensor B", "KF Filt", "KF 95%", "Truth"), pch=c(1, 2, NA, 15, NA), lty=c(NA, NA, 1, NA, 2), lwd=c(NA, NA, 2, NA, 2), col=c("#636363", "#636363", main_cols[i], err_cols[i], true_signal_color), bg="white", horiz=F, ncol=2)
}
