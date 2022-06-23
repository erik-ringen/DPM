delta_theta <- function(X, A, b) {
  
  J <- ncol(X)
  med_X <- apply(X, 2, median)
  mad_X <- apply(X, 2, mad)
  
  theta_med <- rep(NA, J)
  theta_mad <- rep(NA, J)
  diff_theta <- rep(NA, J)
  
  for (j in 1:J) {
    theta_med[j] = -( sum(med_X[-j]*A[j,-j]) + b[j])/A[j,j]
    theta_mad[j] = -( sum((med_X[-j] + mad_X[-j])*A[j,-j]) + b[j])/A[j,j]
    diff_theta[j] = (theta_mad[j] - theta_med[j]) / mad_X[j]
  }
  
  return(diff_theta)
}
