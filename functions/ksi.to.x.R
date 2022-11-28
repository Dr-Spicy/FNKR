#### Define the functions to generate functional predictor from Ksi ####

eigFunct_sin <- function(s, diter){
   temp = sin(diter*pi*s/5) / sqrt(5)
   return(temp)
}

eigFunct_cos <- function(s, diter){
   temp = -cos(diter*pi*s/5) / sqrt(5)
   return(temp)
}

meanFunct <- function(s,diter, status){
   if (status == '0'){
      temp = rep(0,length(s))
   } else if (status == 'l+s'){
      temp = s + sin(s)
   } else if (status == 'l'){
      temp = s
   }
   return(temp)
}

gen_Xi <- function(N, K, diter, eigenvalue_sq){
   # Create FPC scores
   Ksi <- matrix(rnorm(N*K[diter],mean=0,sd = 1), N, K[diter])
   Ksi <- apply(Ksi, 2, scale)
   Ksi <- Ksi %*% diag(eigenvalue_sq[,diter])
   return(Ksi)
}


gen_X <- function(N, K, M, diter, Ksi_all, s, status){
   xTrue <- matrix(0, nrow = N, ncol = M)
   for (i in 1:N){
      xTrue[i,] = Ksi_all[[diter]][i,] %*% t(matrix(c(eigFunct_sin(s,diter),
                                                      eigFunct_cos(s,diter)),
                                                    ncol=K[diter])) + meanFunct(s,diter, status)
   }
   return(xTrue)
}

# gen_X_4 <- function(N, K, M, diter, Ksi_all, s, status){
#   xTrue <- matrix(0, nrow = N, ncol = M)
#   for (i in 1:N){
#     xTrue[i,] = Ksi_all[[diter]][i,] %*% t(matrix(c(eigFunct_sin(s,diter*2),eigFunct_sin(s,(diter*2-1)),eigFunct_cos(s,diter*2),eigFunct_cos(s,(diter*2-1))),ncol=K[diter])) + meanFunct(s,diter, status)
#   }
#   return(xTrue)
# }

gen_X_4 <- function(N, K, M, diter, Ksi_all, s, status){
   xTrue <- matrix(0, nrow = N, ncol = M)
   xTrue = Ksi_all[[diter]] %*% t(matrix(c(eigFunct_sin(s,2),eigFunct_sin(s,1),eigFunct_cos(s,2),eigFunct_cos(s,1)),ncol=K[diter])) + t(matrix(rep(meanFunct(s,diter, status),N),  ncol=N))
   
   return(xTrue)
}



## This function works with fda package in R
# K is a vector of the number of all predictors' generating basis 
gen_X_K <- function(N, M, diter, Ksi_all, s, status){
   xTrue <- matrix(0, nrow = N, ncol = M)
   
   xTrue = Ksi_all[[diter]] %*% t(eval.basis(times,generate.basis[[diter]])) + t(matrix(rep(meanFunct(times, diter, status), N), ncol=N))
   
   return(xTrue)
}
