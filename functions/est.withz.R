#' The NW-like estimation with mixed predictors as a vector having the same dim as y
#' @param y the response vector used in optimization, alway gonna be yTrain here
#' @param Dmatrix the functional predictor distance matrix, always gonna be Dtrain_all here
#' @param lbd the vector of kernel bandwidth 
#' 
#' @export result the NW-like estimation


est.withoutz <- function(y = yTrain, Dmatrix, lbd){
   
   # Calculate the inside of exponential part with X only
   exp.partx = tensor::tensor(Dmatrix, lbd[1:P]^2, 3, 1)
   # Combine two exponential parts, and do exp
   exp_part = exp(-1/2*(exp.partx ))
   # Get the numerator of the NW estimator
   numerator = t(exp_part) %*% y
   # Get the denominator of the NW estimator
   denominator = t(exp_part) %*% rep(1, length(y))
   # Get the export
   return(numerator/denominator)
}

#' The NW-like estimation with mixed predictors as a vector having the same dim as y
#' @param y the response vector used in optimization, alway gonna be yTrain here
#' @param Dmatrix the functional predictor distance matrix, always gonna be Dtrain_all here
#' @param DmatrixZ the vector predictor distance matrix, always gonna be DzTrain here
#' @param lbd the vector of kernel bandwidth 
#' 
#' @export result the NW-like estimation


est.withz <- function(y = yTrain, Dmatrix, DmatrixZ, lbd){
   
   # Calculate the inside of exponential part with X only
   exp.partx = tensor::tensor(Dmatrix, lbd[1:P]^2, 3, 1)
   # Calculate the inside of exponential part with Z only
   exp.partz = tensor::tensor(DmatrixZ,lbd[1:Q+P]^2,3,1)
   # Combine two exponential parts, and do exp
   exp_part = exp(-1/2*(exp.partx + exp.partz))
   # Get the numerator of the NW estimator
   numerator = t(exp_part) %*% y
   # Get the denominator of the NW estimator
   denominator = t(exp_part) %*% rep(1, length(y))
   # Get the export
   return(numerator/denominator)
}

#' The NW-like estimation with mixed predictors as a vector having the same dim as y (With Discrete)
#' @param y the response vector used in optimization, alway gonna be yTrain here
#' @param Dmatrix the functional predictor distance matrix, always gonna be Dtrain_all here
#' @param DmatrixZ the continuous vector predictor distance matrix, always gonna be DzTrain here
#' @param DmatrixZd the discrete vector predictor distance matrix, always gonna be DzdTrain here
#' @param lbd the vector of kernel bandwidth 
#' 
#' @export result the NW-like estimation


est.withz.d <- function(y = yTrain, Dmatrix, DmatrixZ, DmatrixZd, lbd){
   
   # Calculate the inside of exponential part with X only
   exp.partx = tensor::tensor(Dmatrix, lbd[1:P]^2, 3, 1)
   # Calculate the inside of exponential part with Z only
   exp.partz = tensor::tensor(DmatrixZ,lbd[1:Q+P]^2,3,1)
   # Calculate the inside of exponential part with discrete Z only
   exp.partzd = tensor::tensor(DmatrixZd,lbd[1:R+(Q+P)]^2,3,1)
   # Combine two exponential parts, and do exp
   exp_part = exp(-1/2*(exp.partx + exp.partz+ exp.partzd))
   # Get the numerator of the NW estimator
   numerator = t(exp_part) %*% y
   # Get the denominator of the NW estimator
   denominator = t(exp_part) %*% rep(1, length(y))
   # Get the export
   return(numerator/denominator)
}
