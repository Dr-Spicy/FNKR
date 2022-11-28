#' The objective funion during coordiante descent
#' Pretty much a MSE on NW-est on train with lbd_gamma as the kernel bandwidth
#' @param x is the gamma in the optimization 
#' @param lbd_jp1 the the vector with j-th location as 1 and 0 elsewhere
#' @param tau_cand all of the candidate tau vector 
#' @param s the index of currently in use tau, used together with tau_cand
#' @param lbd_mjp_c the lambda minuse j prime vector in use currently, c is for current
#' @param y the response vector used in optimization, alway gonna be yTrain here
#' @param Dmatrix the functional predictor distance matrix, always gonna be Dtrain_all here
#' 
#' @export result The MSE metric


opt1 <- function(x,lbd_jp1, tau_cand, s,lbd_mjp_c,y, Dmatrix){
  
  # lbd_gamma is the lambda(j,gamma) and since diter(j) is the same for all lbd here, so no longer need to specify
  lbd_gamma = (x*lbd_jp1+(tau_cand[s]-x)*lbd_mjp_c)
  # Define a squared prediction error containing vector 
  result = rep(0, length(y))
  # Calculate the inside of exponential part with X only
  exp.partx = tensor::tensor(Dmatrix, lbd_gamma[1:P]^2, 3, 1)
  # Combine two exponential parts, and do exp
  exp_part = exp(-1/2*exp.partx)
  # Get the numerator of the NW estimator
  numerator = t(exp_part) %*% y
  # Get the denominator of the NW estimator
  denominator = t(exp_part) %*% rep(1, length(y))
  # Update the squared prediction error containing vector
  result = (y - numerator/denominator)^2
  # Return the MSE
  return(mean(result))
}
