#' MSE of NW-est on tune with lbd as the kernel bandwidth vector
#' @param yout the response vector used in the outer summation of the MSE function, alway gonna be yTune here
#' @param Dmatrix the functional predictor distance matrix, always gonna be Dtune_all here
#' @param lbd the vector of the kernel bandwidth, always gonna be lbd_c here
#' 
#' @export result The MSE metric



mse.tune.withoutz <- function(yout = yTune){
   # Define a squared prediction error containing vector 
   result = rep(0, length(yout))
   # Update the squared prediction error containing vector
   result = (yout - est.withoutz(y = yTrain, Dmatrix = Dtune_all, lbd = lbd_c))^2
   # Return the MSE
   return(mean(result))
}

#' MSE of NW-est on tune with lbd as the kernel bandwidth vector
#' @param yout the response vector used in the outer summation of the MSE function, alway gonna be yTune here
#' @param Dmatrix the functional predictor distance matrix, always gonna be Dtune_all here
#' @param DmatrixZ the vector predictor distance matrix, always gonna be DzTune here
#' @param lbd the vector of the kernel bandwidth, always gonna be lbd_c here
#' 
#' @export result The MSE metric



mse.tune <- function(yout = yTune){
   # Define a squared prediction error containing vector 
   result = rep(0, length(yout))
   # Update the squared prediction error containing vector
   result = (yout - est.withz(y = yTrain, Dmatrix = Dtune_all, DmatrixZ = DzTune, lbd = lbd_c))^2
   # Return the MSE
   return(mean(result))
}

#' MSE of NW-est on tune with lbd as the kernel bandwidth vector (With discrete)
#' @param yout the response vector used in the outer summation of the MSE function, alway gonna be yTune here
#' @param Dmatrix the functional predictor distance matrix, always gonna be Dtune_all here
#' @param DmatrixZ the continuous vector predictor distance matrix, always gonna be DzTune here
#' @param DmatrixZd the discrete vector predictor distance matrix, always gonna be DzdTune here
#' @param lbd the vector of the kernel bandwidth, always gonna be lbd_c here
#' 
#' @export result The MSE metric



mse.tune.d <- function(yout = yTune){
   # Define a squared prediction error containing vector 
   result = rep(0, length(yout))
   # Update the squared prediction error containing vector
   result = (yout - est.withz.d(y = yTrain, Dmatrix = Dtune_all, DmatrixZ = DzTune, DmatrixZd = DzdTune, lbd = lbd_c))^2
   # Return the MSE
   return(mean(result))
}

#' MSE of NW-est on test with lbd as the kernel bandwidth vector
#' @param yout the response vector used in the outer summation of the MSE function, alway gonna be yTune here
#' @param Dmatrix the functional predictor distance matrix, always gonna be Dtest_all here
#' @param lbd the vector of the kernel bandwidth selected on tune set, always gonna be opted_lbd[repetationnumber, kiter, ] here
#' 
#' @export result The MSE metric


mse.test.withoutz <- function(yout = yTest){
   # Define a squared prediction error containing vector 
   result = rep(0, length(yout))
   # Update the squared prediction error containing vector
   result = (yout - est.withoutz(y = yTrain, Dmatrix = Dtest_all, lbd = conv_sln[which.min(MSE),]))^2
   # Return the MSE
   return(mean(result))
}

#' RMSE of NW-est on test with lbd as the kernel bandwidth vector
#' @param yout the response vector used in the outer summation of the MSE function, alway gonna be yTune here
#' @param Dmatrix the functional predictor distance matrix, always gonna be Dtest_all here
#' @param lbd the vector of the kernel bandwidth selected on tune set, always gonna be opted_lbd[repetationnumber, kiter, ] here
#' 
#' @export result The MSE metric


rmse.test.withoutz <- function(yout = yTest){
   # Define a squared prediction error containing vector 
   result = rep(0, length(yout))
   # Update the squared prediction error containing vector
   result = (yout - est.withoutz(y = yTrain, Dmatrix = Dtest_all, lbd = conv_sln[which.min(MSE),]))^2
   # Return the MSE
   return(sqrt(mean(result)))
}


#' MSE of NW-est on test with lbd as the kernel bandwidth vector
#' @param yout the response vector used in the outer summation of the MSE function, alway gonna be yTune here
#' @param Dmatrix the functional predictor distance matrix, always gonna be Dtest_all here
#' @param DmatrixZ the vector predictor distance matrix, always gonna be DzTest here
#' @param lbd the vector of the kernel bandwidth selected on tune set, always gonna be opted_lbd[repetationnumber, kiter, ] here
#' 
#' @export result The MSE metric



mse.test <- function(yout = yTest){
   # Define a squared prediction error containing vector 
   result = rep(0, length(yout))
   # Update the squared prediction error containing vector
   result = (yout - est.withz(y = yTrain, Dmatrix = Dtest_all, DmatrixZ = DzTest, lbd = conv_sln[which.min(MSE),]))^2
   # Return the MSE
   return(mean(result))
}

#' RMSE of NW-est on test with lbd as the kernel bandwidth vector
#' @param yout the response vector used in the outer summation of the MSE function, alway gonna be yTune here
#' @param Dmatrix the functional predictor distance matrix, always gonna be Dtest_all here
#' @param DmatrixZ the vector predictor distance matrix, always gonna be DzTest here
#' @param lbd the vector of the kernel bandwidth selected on tune set, always gonna be opted_lbd[repetationnumber, kiter, ] here
#' 
#' @export result The MSE metric



rmse.test <- function(yout = yTest){
   # Define a squared prediction error containing vector 
   result = rep(0, length(yout))
   # Update the squared prediction error containing vector
   result = (yout - est.withz(y = yTrain, Dmatrix = Dtest_all, DmatrixZ = DzTest, lbd = conv_sln[which.min(MSE),]))^2
   # Return the MSE
   return(sqrt(mean(result)))
}

#' MSE of NW-est on test with lbd as the kernel bandwidth vector (With discrete)
#' @param yout the response vector used in the outer summation of the MSE function, alway gonna be yTune here
#' @param Dmatrix the functional predictor distance matrix, always gonna be Dtest_all here
#' @param DmatrixZ the continuous vector predictor distance matrix, always gonna be DzTest here
#' @param DmatrixZd the discrete vector predictor distance matrix, always gonna be DzdTest here
#' @param lbd the vector of the kernel bandwidth selected on tune set, always gonna be opted_lbd[repetationnumber, kiter, ] here
#' 
#' @export result The MSE metric



mse.test.d <- function(yout = yTest){
   # Define a squared prediction error containing vector 
   result = rep(0, length(yout))
   # Update the squared prediction error containing vector
   result = (yout - est.withz.d(y = yTrain, Dmatrix = Dtest_all, DmatrixZ = DzTest, DmatrixZd = DzdTest, lbd = conv_sln[which.min(MSE),]))^2
   # Return the MSE
   return(mean(result))
}

#' RMSE of NW-est on test with lbd as the kernel bandwidth vector (With discrete)
#' @param yout the response vector used in the outer summation of the MSE function, alway gonna be yTune here
#' @param Dmatrix the functional predictor distance matrix, always gonna be Dtest_all here
#' @param DmatrixZ the continuous vector predictor distance matrix, always gonna be DzTest here
#' @param DmatrixZd the discrete vector predictor distance matrix, always gonna be DzdTest here
#' @param lbd the vector of the kernel bandwidth selected on tune set, always gonna be opted_lbd[repetationnumber, kiter, ] here
#' 
#' @export result The MSE metric



rmse.test.d <- function(yout = yTest){
   # Define a squared prediction error containing vector 
   result = rep(0, length(yout))
   # Update the squared prediction error containing vector
   result = (yout - est.withz.d(y = yTrain, Dmatrix = Dtest_all, DmatrixZ = DzTest, DmatrixZd = DzdTest, lbd = conv_sln[which.min(MSE),]))^2
   # Return the MSE
   return(sqrt(mean(result)))
}