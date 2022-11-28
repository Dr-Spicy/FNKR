
set.seed(12)


# env
Run.env = 'Statlab'

#### Define the simulation parameters ####

# Tell me the model index to run
model.index <- 'FAM'
# Tell me the number of repetition
repeatlimit <- 1000
# Set the number of subjects (N) and the size of train, tune, test
N <- 800; n1 <- .25*N; n2 <- .25*N; n3 <- .5*N;
# number of functional predictors (P) and number of non-fucntional predictors (Q)
P <- 5; Q <- 5; R <- 0;
# number of the grid points per f predictor to sample from (M) 
M <- 101;
#### Define the functional predictor's parameters ####

# number of obs per subject curve
sparsity_lvl <- seq(M,M,1);
# set the variance of eigen components of all the functional predictors
FPCv <- c(16,4,1,0.25); me <- sqrt(sum(FPCv)/me.ratio);
# mean function of the functional predictor (l+s:linear+sine)
meanfstat = "l+s";
# Set the continuum of the functional curves are defined on
sc <- seq(0,10,length.out = M)
# number of eigen components per func pred to generate, K is a P-dim vector (potentially different # for each func pred)
numberofK <- 4; K <- rep(numberofK,P)
# define the weights
weights_aa=c(1,1,1,1)

#### load the functions and modules ####
if (Run.env == 'PC') {
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/true.model.R")
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/ksi.to.x.R")
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/ksi.to.yTrue.R")
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/FPCAmaxtime.R")
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/pred.FPCs.R")
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/get.C.mat.R")
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/est.withz.R")
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/opt3.R")
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/mse.prdc.R")
} else if (Run.env == 'Statlab') {
  source("~/Functional Kernel Regression/functions/true.model.R")
  source("~/Functional Kernel Regression/functions/ksi.to.x.R")
  source("~/Functional Kernel Regression/functions/ksi.to.yTrue.R")
  source("~/Functional Kernel Regression/functions/FPCAmaxtime.R")
  source("~/Functional Kernel Regression/functions/pred.FPCs.R")
  source("~/Functional Kernel Regression/functions/get.C.mat.R")
  source("~/Functional Kernel Regression/functions/est.withz.R")
  source("~/Functional Kernel Regression/functions/opt3.R")
  source("~/Functional Kernel Regression/functions/mse.prdc.R")
}

#### Define the results containers ####
mu.x1 = rep(0, repeatlimit); mu.x2 = rep(0, repeatlimit); 
mu.z1 = rep(0, repeatlimit); mu.z2 = rep(0, repeatlimit);
repetationnumber=1

for (repetationnumber in 1:repeatlimit) {
  #### Define all the needed matrices to store intermediate results####
  xTrue_all = list(); xNoisy_all = list(); xSparse_all = list(); Ksi_all = list();
  Zmatrix = matrix(0, nrow = N, ncol = Q); eigenvalue_sq = matrix(0,nrow = numberofK, ncol = P);
  ## Define the beta matrix
  betamatrix = array(0,dim = c(P,M))
  
  #### Generate the FPC scores, Xtrues, then sparsify, eventually and/or ME ####
  ## Work on the all the predictors one by one
  
  for (diter in 1:P){
    # Generate eigvenvalues for the current predictor
    eigenvalue_sq[,diter] <- sqrt(FPCv)
    # Create FPC scores
    Ksi_all[[diter]] <- gen_Xi(N,K,diter, eigenvalue_sq)
    # Generate the functional X
    xTrue_all[[diter]] <- gen_X_4(N, K, M, diter, Ksi_all, sc, status = meanfstat)
  }
  
  #### generate the beta's ####
  betamatrix[1,] = beta1(sc)
  betamatrix[2,] = beta2(sc)
  betamatrix[3,] = beta3(sc)
  betamatrix[4,] = beta4(sc)
  betamatrix[5,] = beta5(sc)
  
  
  #### Generate the non-functional predictors as a standardized Q-variate normal dist ####
  Zmatrix <- MASS::mvrnorm(n=N, mu = rep(0,Q), Sigma = diag(rep(1,Q)))
  Zmatrix <- apply(Zmatrix, 2, scale)
  Zmatrix <- Zmatrix %*% diag(rep(1,Q))
  Z1 = Zmatrix[,1]
  Z2 = Zmatrix[,2]
  Z3 = Zmatrix[,3]
  Z4 = Zmatrix[,4]
  Z5 = Zmatrix[,5]
  
  
  
  #### Generate the true scalar response yTrue according to each true model ####
  yTrue <- Generate.yTrue(model.index, Ksi_all)
  
  yTrue.no.x1 = rep(0,N)
  ## 1st functional predictor
  diter = 1
  # apply the defined linear fucntion on ksi_11 and ksi_13
  temp = matrix(rep(apply(xTrue_all[[diter]], 1, mean), M), nrow = N, ncol = M) %*% betamatrix[diter,] * (max(sc)- min(sc))/(M-1)
  # update the yTrue
  yTrue.no.x1 = yTrue.no.x1 + temp
  
  # ## 2th predictor
  diter = 2
  # apply the defined quadratic fucntion on ksi_11 and ksi_13
  temp = xTrue_all[[diter]] %*% betamatrix[diter,] * (max(sc)- min(sc))/(M-1)
  # update the yTrue
  yTrue.no.x1 = yTrue.no.x1 + temp
  
  ## 1st scalar predictor
  temp = weights_aa[3] * sqrt(2) * Z1
  # update the yTrue
  yTrue.no.x1 = yTrue.no.x1 + temp 
  
  ## 2nd scalar predictor
  temp = weights_aa[4] * 1 *Z2^2
  # update the yTrue
  yTrue.no.x1 = yTrue.no.x1 + temp 
  
  mu.x1[repetationnumber] = mean((yTrue - yTrue.no.x1)**2)
  
  yTrue.no.x2 = rep(0,N)
  ## 1st functional predictor
  diter = 1
  # apply the defined linear fucntion on ksi_11 and ksi_13
  temp = xTrue_all[[diter]] %*% betamatrix[diter,] * (max(sc)- min(sc))/(M-1)
  # update the yTrue
  yTrue.no.x2 = yTrue.no.x2 + temp
  
  # ## 2th predictor
  diter = 2
  # apply the defined quadratic fucntion on ksi_11 and ksi_13
  temp = matrix(rep(apply(xTrue_all[[diter]], 1, mean), M), nrow = N, ncol = M) %*% betamatrix[diter,] * (max(sc)- min(sc))/(M-1)
  # update the yTrue
  yTrue.no.x2 = yTrue.no.x2 + temp
  
  ## 1st scalar predictor
  temp = weights_aa[3] * sqrt(2) * Z1
  # update the yTrue
  yTrue.no.x2 = yTrue.no.x2 + temp 
  
  ## 2nd scalar predictor
  temp = weights_aa[4] * 1 *Z2^2
  # update the yTrue
  yTrue.no.x2 = yTrue.no.x2 + temp 
  
  mu.x2[repetationnumber] = mean((yTrue - yTrue.no.x2)**2)
  
  yTrue.no.z1 = rep(0,N)
  ## 1st functional predictor
  diter = 1
  # apply the defined linear fucntion on ksi_11 and ksi_13
  temp = xTrue_all[[diter]] %*% betamatrix[diter,] * (max(sc)- min(sc))/(M-1)
  # update the yTrue
  yTrue.no.z1 = yTrue.no.z1 + temp
  
  # ## 2th predictor
  diter = 2
  # apply the defined quadratic fucntion on ksi_11 and ksi_13
  temp = xTrue_all[[diter]] %*% betamatrix[diter,] * (max(sc)- min(sc))/(M-1)
  # update the yTrue
  yTrue.no.z1 = yTrue.no.z1 + temp
  
  ## 1st scalar predictor
  temp = 0
  # update the yTrue
  yTrue.no.z1 = yTrue.no.z1 + temp 
  
  ## 2nd scalar predictor
  temp = weights_aa[4] * 1 *Z2^2
  # update the yTrue
  yTrue.no.z1 = yTrue.no.z1 + temp 
  
  mu.z1[repetationnumber] = mean((yTrue - yTrue.no.z1)**2)
  
  yTrue.no.z2 = rep(0,N)
  ## 1st functional predictor
  diter = 1
  # apply the defined linear fucntion on ksi_11 and ksi_13
  temp = xTrue_all[[diter]] %*% betamatrix[diter,] * (max(sc)- min(sc))/(M-1)
  # update the yTrue
  yTrue.no.z2 = yTrue.no.z2 + temp
  
  # ## 2th predictor
  diter = 2
  # apply the defined quadratic fucntion on ksi_11 and ksi_13
  temp = xTrue_all[[diter]] %*% betamatrix[diter,] * (max(sc)- min(sc))/(M-1)
  # update the yTrue
  yTrue.no.z2 = yTrue.no.z2 + temp
  
  ## 1st scalar predictor
  temp = weights_aa[3] * sqrt(2) * Z1
  # update the yTrue
  yTrue.no.z2 = yTrue.no.z2 + temp 
  
  ## 2nd scalar predictor
  temp = 0
  # update the yTrue
  yTrue.no.z2 = yTrue.no.z2 + temp 
  
  mu.z2[repetationnumber] = mean((yTrue - yTrue.no.z2)**2)
}
mux1 = mean(mu.x1)
mux2 = mean(mu.x2)
muz1 = mean(mu.z1)
muz2 = mean(mu.z2)

max = max(mux1, mux2, muz1, muz2)
Pesx1 = sqrt(mux1/max)
Pesx2 = sqrt(mux2/max)

Pesz1 = sqrt(muz1/max)
Pesz2 = sqrt(muz2/max)
