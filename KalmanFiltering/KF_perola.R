KF_perola <- function(vi_in, qall, ut){
  # Kalman filter of vegeation data with the model:
  # X(t) = AX(t-1) + BU(t) + e(t)
  # Y(t) = CX(t) + w(q)
  # where:
  # X(t) is the vegetation index
  # U(t) is differentiated "normal" seasonal trajectory (external signal)
  # e(t) is model noise
  # Y(t) is vegetation index from satellite data (observable)
  # w(q) is measurement noise
  # A, B, C are in this case coeffcients
  # Settings are tailored for NDVI, birch forest Abisko
  # Parameters:
  # vi_in = matrix with vegetationindex, pixels x vi values.
  # qall = quality data in classes 1-3
  # ut = external signal
  # params = [A B C Re Rw_1 Rxx_1 xtt_1]
  # where: Re is model noise, Rw_1 measurement noise, Rxx_1, xxt_1
  # are initial values. Default is given below.

  # Quality data (for ndvi)
  # Weights for adjusting Rw depending in QA; highest QA to the right
  qw <- c(3, 0.5, 0)

  # define parameters
  A <- 1
  Bq <- c(0.8, 1, 1)
  Cq <- c(1.1, 1, 1)
  Re <- 0.5
  Rw_1 <- 8 # Measure variance

  nbr_pixels <- dim(vi_in)[1]
  N <- dim(vi_in)[2]

  vi_kalman <- matrix(0, nbr_pixels, N) # Kalman filtered
  er_kalman <- matrix(0, nbr_pixels, N) # Estimate errors, just for testing

  for(p in 1:nbr_pixels){
    # get signal and QA
    y <- vi_in[p,]
    qa <- qall[p,]

    # set initial values
    Rxx_1 <- 5 # error covariance
    xtt_1 <- 0 # Since starting in winter

    xsave <- rep(0, N) # for temporarily holding output

    # Kalman filter
    for(k in 2:N){
      # To get appropriate values according to QA for present date
      Rw <- qw[qa[k]] * Rw_1
      C <- Cq[qa[k]]
      B <- Bq[qa[k]]

      # Update
      Ryy <- C * Rxx_1 * t(C) + Rw # Just to split the Kalman gain below
      Kt <- Rxx_1 * t(C) * solve(Ryy) # Kalman gain
      xtt <- xtt_1 + (Kt * (y[k] - C * xtt_1))
      Rxx <- Rxx_1 - (Kt * C * Rxx_1)

      # Save the state
      xsave[k] <- xtt
      er_kalman[p, k] <- xtt_1 - xtt

      # Predict
      Rxx_1 = A * Rxx * t(A) + Re
      xtt_1 <- A * xtt + B * ut[p, k]
    }
    # Saving all seasons for Kalman filtered data and one pixel
    vi_kalman[p,] <- xsave
  }
  return(list(vi_kalman, er_kalman))
}
