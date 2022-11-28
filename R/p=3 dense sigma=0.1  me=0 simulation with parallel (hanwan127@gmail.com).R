library(data.table)
library(fdapace)
library(doParallel)
library(foreach)
library(ggplot2)
library(nloptr)


## to unregister the cl
registerDoSEQ()
cl <- makeCluster(8)
registerDoParallel(cl)
# stopCluster(cl)

#### generating a toy dense functional dataset from scratch####

# Set the number of subjects (N) and the
# number of predictors (P)
# number of measurements per predictor(M)
N <- 200;  # when N=500, PACE time = 657s
P <- 3;
M <- 101;
sigma <- 0.1;
set.seed(100)

# Define the continuum
s <- seq(0,10,length.out = M)


# Define the number of eigencomponents
K_tot <- 10

# Define the mean and all possible eigencomponents
meanFunct1 <- function(s) s
meanFunct2 <- function(s) s + sin(s)
meanFunct3 <- function(s) s + cos(s)

eigFunct1 <- function(s) cos(s*pi/5) / sqrt(5)
eigFunct2 <- function(s) sin(s*pi/5) / sqrt(5)
eigFunct3 <- function(s) cos(2*s*pi/5) / sqrt(5)
eigFunct4 <- function(s) sin(2*s*pi/5) / sqrt(5)
eigFunct5 <- function(s) cos(4*s*pi/5) / sqrt(5)
eigFunct6 <- function(s) sin(4*s*pi/5) / sqrt(5)
eigFunct7 <- function(s) cos(6*s*pi/5) / sqrt(5)
eigFunct8 <- function(s) sin(6*s*pi/5) / sqrt(5)
eigFunct9 <- function(s) cos(8*s*pi/5) / sqrt(5)
eigFunct10 <- function(s) sin(8*s*pi/5) / sqrt(5)

meanFunct <- matrix(c(meanFunct1(s),meanFunct2(s),meanFunct3(s)), ncol=3)

eigFunct <- matrix(c(eigFunct1(s),eigFunct2(s),eigFunct3(s),eigFunct4(s),
                     eigFunct5(s),eigFunct6(s),eigFunct7(s),eigFunct8(s),
                     eigFunct9(s),eigFunct10(s)), ncol=10)

#### Define the needed hyperparameters ####
# number of eigen components per predictor
K <- c(3,3,3)
#

#### Define all the needed matrices ####
xTrue_all = list()
Ksi_all = list()
Dtrain_all = list()
Dtune_all = list()
Dtest_all = list()


#### Now, work on the all the predictors one by one ####

start_time = Sys.time()

for (diter in 1:P){
  # Create FPC scores
  lambda_sq = sort(sample(c(1,2,3), K[diter], replace = TRUE), decreasing = T)
  Ksi <- matrix(rnorm(N*K[diter],mean=0,sd=lambda_sq), N, K[diter])
  Ksi <- apply(Ksi, 2, scale)
  Ksi <- Ksi %*% diag(lambda_sq)
  Ksi_all[[diter]] <- Ksi
  
  xTrue <- Ksi %*% t(eigFunct[,seq(diter:(diter+2))]) + t(matrix(rep(meanFunct[,diter],N), nrow=M))
  xTrue_all[[diter]] <- xTrue
  
  ### separate the simulated data into training, tuning and testing
  training = MakeFPCAInputs(IDs = rep(1:(.6*N), each=M), tVec=rep(s,(.6*N)), yVec = t(xTrue[1:(.6*N),]))
  tuning = MakeFPCAInputs(IDs = rep((.6*N+1):(.8*N), each=M), tVec=rep(s,(.2*N)), yVec = t(xTrue[(.6*N+1):(.8*N),]))
  testing = MakeFPCAInputs(IDs = rep((.8*N+1):N, each=M), tVec=rep(s,(.2*N)), yVec = t(xTrue[(.8*N+1):N,]))
  
  ####Running FPCA on training set ####
  
  # Do FPCA on this sparse sample with 99% FVE cut
  # Notice that sparse FPCA will smooth the data internally (Yao et al., 2005)
  # Smoothing is the main computational cost behind sparse FPCA
  FPCAdense <- FPCA(training$Ly, training$Lt,
                    list(plot = F, dataType ='Dense', nRegGrid = 1001, useBinnedData = 'OFF', error = F))
  
  # Get the central S matrix
  temp <- t(FPCAdense$phi) %*% FPCAdense$phi
  S <- temp %*% diag(1/FPCAdense$lambda) %*% temp
  
  #### Predict the first K* FPC scores according to conditional expectation on the tuning ####
  pred_zeta <- predict(FPCAdense, tuning$Ly, tuning$Lt, K=FPCAdense$selectK)
  
  #### Predict the first K* FPC scores according to conditional expectation on the tuning ####
  pred_eta <- predict(FPCAdense, testing$Ly, testing$Lt, K=FPCAdense$selectK)
  
  #### Calculate all the pair-wise distances between all training samples ####
  Dtrain <- matrix(0, nrow = length(training$Ly), ncol = length(training$Ly))
  # start_time = Sys.time()
  Dtrain = foreach (i = 1:ncol(Dtrain), .combine=cbind) %:%
    foreach (j = 1:nrow(Dtrain), .combine=rbind) %dopar%{
      temp = FPCAdense$xiEst[j,] - FPCAdense$xiEst[i,]
      Dtrain[j,i] = t(temp) %*% S %*% temp
    }
  # End_time = Sys.time()
  # End_time - start_time
  
  Dtrain = Dtrain/((M-1)/10)^2
  
  
  #### Calculate all the pair-wise distances between all training samples and all tuning samples####
  Dtune <- matrix(nrow = length(training$Ly), ncol = length(tuning$Ly))
  # start_time = Sys.time()
  Dtune = foreach (i = 1:ncol(Dtune), .combine=cbind) %:%
    foreach (j = 1:nrow(Dtune), .combine=rbind) %dopar%{
      temp = FPCAdense$xiEst[j,] - pred_zeta[i,]
      Dtune[j,i] = t(temp) %*% S %*% temp
    }
  # End_time = Sys.time()
  # End_time - start_time
  
  Dtune = Dtune/((M-1)/10)^2
  
  #### Calculate all the pair-wise distances between all training samples and all testing samples####
  Dtest <- matrix(nrow = length(training$Ly), ncol = length(testing$Ly))
  #start_time = Sys.time()
  Dtest = foreach (i = 1:ncol(Dtest), .combine=cbind) %:%
    foreach (j = 1:nrow(Dtest), .combine=rbind) %dopar%{
      temp = FPCAdense$xiEst[j,] - pred_eta[i,]
      Dtest[j,i] = t(temp) %*% S %*% temp
    }
  #End_time = Sys.time()
  #End_time - start_time
  
  Dtest = Dtest/((M-1)/10)^2
  
  Dtrain_all[[diter]] = Dtrain
  Dtune_all[[diter]] = Dtune
  Dtest_all[[diter]] = Dtest
  
}

End_time = Sys.time()
End_time - start_time

#### save the objects ####
save(Dtrain_all,file = 'Dtrain_all_3_D_nome')
save(Dtune_all,file = 'Dtune_all_3_D_nome')
save(Dtest_all,file = 'Dtest_all_3_D_nome')


#### Compute the Y thru a non-linear relationship by defined nl functions ####
# m1 is exp then find the average
# m2 is qudratic then find the average
# m3 is cubic then find the average

yTrue <- apply(exp(xTrue_all[[1]]),1,sum)/M + apply(xTrue_all[[2]]^2,1,sum)/M 

y <- yTrue +  rnorm(N,0,sigma)

#### Find the optimal total budget over the tuning sets ####

pi_cand = c(1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1, 5, 1e1, 5e1, 1e2, 5e2, 1e3)

S = length(pi_cand)

# perform the modified coordinate descent algo on the tuning sets
foreach (s = 1:S) %:%
  foreach (jp = 1:P) %dopar% {
    # Initialize with equal smoothing
    lbd_jp = rep(pi_cand[s]/P, P)
    # update the initial vector as the current solution
    lbd_jp_c = lbd_jp
    # define the two tool vectors
    lbd_jp1 = rep(0,P)
    lbd_jp1[jp] = 1
    lbd_mjp_c = lbd_jp_c
    lbd_mjp_c[jp] = 0 
    lbd_mjp_c = lbd_mjp_c/sum(lbd_mjp_c)
    
    
  }


### Do the optimization over gamma, denoted as z here ###

# define the objective function
eval_f0_tune <- function (z,lbd_jp1, pi_cand,s,lbd_mjp_c,jp,
                          Dtune_all,){
  lambdav = (z*lbd_jp1+(pi_cand[s]-z)*lbd_mjp_c)
  
  
  return()
}
