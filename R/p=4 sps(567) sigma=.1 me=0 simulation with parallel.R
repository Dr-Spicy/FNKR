library(data.table)
library(fdapace)
library(doParallel)
library(foreach)
library(ggplot2)
library(nloptr)
library(Rmpfr)

source("def_func.R")


## to unregister the cl
registerDoSEQ()
cl <- makeCluster(8)
registerDoParallel(cl)
# stopCluster(cl)

#### generating a toy dense functional dataset from scratch####

# Set the number of subjects (N) and the
# number of predictors (P)
# number of measurements per predictor(M)
N <- 300;  # when N=500, PACE time = 657s
P <- 4;
M <- 101;
me <- 0;
sigma = 1;
rpt_tol = 1e-2
set.seed(100)

# Define the continuum
s <- seq(0,10,length.out = M)


# Define the number of eigencomponents
K_tot <- 12



#### Define the needed hyperparameters and fucntions ####
# number of eigen components per predictor
K <- c(3,3,3,3)

#### Define all the needed matrices ####
xTrue_all = list()
xNoisy_all = list()
Ksi_all = list()
Dtrain_all = list()
Dtune_all = list()
Dtest_all = list()
eigFunct = list()
yTrue = 0

# sample number of measurement per subjects
num_m = sample(x=matrix(rep(c(6,8,10),N),nrow= N), size = N)

# sample ti's for each sample based on the # of measurement on [0,10] as a list
sampled_t = list()
for (i in 1:N){
  sampled_t[[i]] = sort(runif(num_m[i],0,10)) # ascending order
}

### have to make sure the the range of training is larger than the range of testing and tuning
# find the min, max 
min_t = min(unlist(sampled_t))
max_t = max(unlist(sampled_t))
# find the index of the lists containes min and max
inx_mint = which(sapply(sampled_t,function(x) min_t %in% x))
inx_maxt = which(sapply(sampled_t,function(x) max_t %in% x))
# save the 2 extreme lists as temps
temp1 <- sampled_t[[inx_mint]]
temp2 <- sampled_t[[inx_maxt]]
# replace the 2 extreme lists by the first 2 lists
sampled_t[[inx_mint]] <- sampled_t[[1]]
sampled_t[[inx_maxt]] <- sampled_t[[2]]
# repalce the first 2 lists by the 2 temps
sampled_t[[1]] <- temp1
sampled_t[[2]] <- temp2

#### Now, work on the all the predictors one by one ####
## Compute the Y thru a non-linear relationship by defined nl functions ##
# m1 is exp then find the average
# m2 is qudratic then find the average
# m3 is cubic then find the average

## 1st predictor
diter = 1

# Create FPC scores
lambda_sq = sort(sample(c(1,2,3), K[diter], replace = TRUE), decreasing = T)
Ksi <- matrix(rnorm(N*K[diter],mean=0,sd=lambda_sq/10), N, K[diter])
Ksi <- apply(Ksi, 2, scale)
Ksi <- Ksi %*% diag(lambda_sq)
Ksi_all[[diter]] <- Ksi
xTrue <- list()
for (i in 1:N){
  xTrue[[i]] = Ksi[i,] %*% t(matrix(c(eigFunct1(sampled_t[[i]]),
                                      eigFunct2(sampled_t[[i]]),eigFunct3(sampled_t[[i]])), 
                             ncol=K[diter]))+ meanFunct1(sampled_t[[i]])
}
xTrue_all[[diter]] <- xTrue

# Add measurement errors in X on X_true
xNoisy <- xTrue
for (i in 1:N){
  xNoisy[[i]] = xNoisy[[i]] + rnorm(xTrue[[i]],sd = me)
}
xNoisy_all[[diter]] <- xNoisy

## m1 is sin then find the average
# apply the defined fucntion on this X
temp = sapply(xNoisy_all[[diter]],sin)
# convert list to matrix with fillers = 0
temp = list.as.matrix(temp, byrow = T, filler = 0)
# get the mean over row for each sample
temp1  = apply(temp,1,mean)
# update the yTrue
yTrue = yTrue + t(temp1)

## 2nd predictor
diter = 2

# Create FPC scores
lambda_sq = sort(sample(c(1,2,3), K[diter], replace = TRUE), decreasing = T)
Ksi <- matrix(rnorm(N*K[diter],mean=0,sd=lambda_sq/10), N, K[diter])
Ksi <- apply(Ksi, 2, scale)
Ksi <- Ksi %*% diag(lambda_sq)
Ksi_all[[diter]] <- Ksi
xTrue <- list()
for (i in 1:N){
  xTrue[[i]] = Ksi[i,] %*% t(matrix(c(eigFunct4(sampled_t[[i]]),
                                      eigFunct5(sampled_t[[i]])), 
                                    ncol=K[diter]))+ meanFunct2(sampled_t[[i]])
}
xTrue_all[[diter]] <- xTrue

# Add measurement errors in X on X_true
xNoisy <- xTrue
for (i in 1:N){
  xNoisy[[i]] = xNoisy[[i]] + rnorm(xTrue[[i]],sd = me)
}
xNoisy_all[[diter]] <- xNoisy

# m2 is qudratic then find the average
# apply the defined fucntion on this X
temp = sapply(xNoisy_all[[diter]],quadratic)
# convert list to matrix with fillers = 0
temp = list.as.matrix(temp, byrow = T, filler = 0)
# get the mean over row for each sample
temp1  = apply(temp,1,mean)
# update the yTrue
yTrue <- yTrue + t(temp1)

## 3rd predictor
diter = 3

# Create FPC scores
lambda_sq = sort(sample(c(1,2,3), K[diter], replace = TRUE), decreasing = T)
Ksi <- matrix(rnorm(N*K[diter],mean=0,sd=lambda_sq/10), N, K[diter])
Ksi <- apply(Ksi, 2, scale)
Ksi <- Ksi %*% diag(lambda_sq)
Ksi_all[[diter]] <- Ksi
xTrue <- list()
for (i in 1:N){
  xTrue[[i]] = Ksi[i,] %*% t(matrix(c(eigFunct7(sampled_t[[i]]),
                                      eigFunct8(sampled_t[[i]]),eigFunct9(sampled_t[[i]])), 
                                    ncol=K[diter]))+ meanFunct3(sampled_t[[i]])
}
xTrue_all[[diter]] <- xTrue
# Add measurement errors in X on X_true
xNoisy <- xTrue
for (i in 1:N){
  xNoisy[[i]] = xNoisy[[i]] + rnorm(xTrue[[i]],sd = me)
}
xNoisy_all[[diter]] <- xNoisy

# m3 is cos then find the average
# apply the defined fucntion on this X
temp = sapply(xNoisy_all[[diter]],cos)
# convert list to matrix with fillers = 0
temp = list.as.matrix(temp, byrow = T, filler = 0)
# get the mean over row for each sample
temp1  = apply(temp,1,mean)
# update the yTrue (nothing b/c x3 is not important)
# yTrue <- yTrue + t(temp1)

## 4th predictor
diter = 4

# Create FPC scores
lambda_sq = sort(sample(c(1,2,3), K[diter], replace = TRUE), decreasing = T)
Ksi <- matrix(rnorm(N*K[diter],mean=0,sd=lambda_sq/10), N, K[diter])
Ksi <- apply(Ksi, 2, scale)
Ksi <- Ksi %*% diag(lambda_sq)
Ksi_all[[diter]] <- Ksi
xTrue <- list()
for (i in 1:N){
  xTrue[[i]] = Ksi[i,] %*% t(matrix(c(eigFunct10(sampled_t[[i]]),
                                      eigFunct11(sampled_t[[i]]),eigFunct12(sampled_t[[i]])), 
                                    ncol=K[diter]))+ meanFunct4(sampled_t[[i]])
}
xTrue_all[[diter]] <- xTrue
# Add measurement errors in X on X_true
xNoisy <- xTrue
for (i in 1:N){
  xNoisy[[i]] = xNoisy[[i]] + rnorm(xTrue[[i]],sd = me)
}
xNoisy_all[[diter]] <- xNoisy

# m4 is abs then find the average
# apply the defined fucntion on this X
temp = sapply(xNoisy_all[[diter]],abs)
# convert list to matrix with fillers = 0
temp = list.as.matrix(temp, byrow = T, filler = 0)
# get the mean over row for each sample
temp1  = apply(temp,1,mean)
# update the yTrue (nothing b/c x4 is not important)
# yTrue <- yTrue + t(temp1)



#### add error term to yTrue as y ####
y <- yTrue +  rnorm(N,0,sigma)



#### Now, find all the distance matrices one by one ####

start_time = Sys.time()

for (diter in 1:P){
  ### separate the simulated data into training, tuning and testing
  training = list(x= xTrue_all[[diter]][c(1:(.6*N))], 
                  t= sampled_t[c(1:(.6*N))], y = y[1:(.6*N)])
  tuning = list(x= xTrue_all[[diter]][c((.6*N+1):(.8*N))], 
                t= sampled_t[c((.6*N+1):(.8*N))], y = y[(.6*N+1):(.8*N)])
  testing = list(x= xTrue_all[[diter]][c((.8*N+1):N)], 
                 t= sampled_t[c((.8*N+1):N)], y = y[(.8*N+1):N])

  ####Running FPCA on training set ####

  # Do FPCA on this sparse sample with 99% FVE cut
  # Notice that sparse FPCA will smooth the data internally (Yao et al., 2005)
  # Smoothing is the main computational cost behind sparse FPCA
  FPCAsparse <- FPCA(training$x, training$t,
                    list(plot = F, FVEthreshold = 0.90, nRegGrid = 301, useBinnedData = 'OFF', error = F))

  # Get the central S matrix
  temp <- t(FPCAsparse$phi) %*% FPCAsparse$phi
  S <- temp %*% diag(1/FPCAsparse$lambda) %*% temp

  #### Predict the first K* FPC scores according to conditional expectation on the tuning ####
  pred_zeta <- predict(FPCAsparse, tuning$x, tuning$t, K=FPCAsparse$selectK)

  #### Predict the first K* FPC scores according to conditional expectation on the tuning ####
  pred_eta <- predict(FPCAsparse, testing$x, testing$t, K=FPCAsparse$selectK)

  #### Calculate all the pair-wise distances between all training samples ####
  Dtrain <- matrix(0, nrow = length(training$x), ncol = length(training$x))
  # start_time = Sys.time()
  Dtrain = foreach (i = 1:ncol(Dtrain), .combine=cbind) %:%
    foreach (j = 1:nrow(Dtrain), .combine=rbind) %dopar%{
      temp = FPCAsparse$xiEst[j,] - FPCAsparse$xiEst[i,]
      Dtrain[j,i] = t(temp) %*% S %*% temp
    }
  # End_time = Sys.time()
  # End_time - start_time

  Dtrain = Dtrain/((M-1)/10)^2


  #### Calculate all the pair-wise distances between all training samples and all tuning samples####
  Dtune <- matrix(nrow = length(training$x), ncol = length(tuning$x))
  # start_time = Sys.time()
  Dtune = foreach (i = 1:ncol(Dtune), .combine=cbind) %:%
    foreach (j = 1:nrow(Dtune), .combine=rbind) %dopar%{
      temp = FPCAsparse$xiEst[j,] - pred_zeta[i,]
      Dtune[j,i] = t(temp) %*% S %*% temp
    }
  # End_time = Sys.time()
  # End_time - start_time

  Dtune = Dtune/((M-1)/10)^2

  #### Calculate all the pair-wise distances between all training samples and all testing samples####
  Dtest <- matrix(nrow = length(training$x), ncol = length(testing$x))
  #start_time = Sys.time()
  Dtest = foreach (i = 1:ncol(Dtest), .combine=cbind) %:%
    foreach (j = 1:nrow(Dtest), .combine=rbind) %dopar%{
      temp = FPCAsparse$xiEst[j,] - pred_eta[i,]
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
save(Dtrain_all,file = 'Dtrain_all_4_Sps_sig.1_nome.Rdata')
save(Dtune_all,file = 'Dtune_all_4_Sps_sig.1_nome.Rdata')
save(Dtest_all,file = 'Dtest_all_4_Sps_sig.1_nome.Rdata')





#### Find the optimal total budget over the tuning sets ####

pi_cand = seq(0, 10, length.out = 11)

S = length(pi_cand)

MSE = rep(0,S)
zopt = rep(0,S)
conv_sln = matrix(0,nrow = S, ncol = P)

s=3;jp=1



for (s in 2:S){
  # Initialize with equal smoothing
  lbd_0 = rep(pi_cand[s]/P, P)
  # update the initial vector as the current solution
  lbd_c = lbd_0
  # define a vector to record the MSE for the current s
  mem = c(1e15,1e12)
  
  repeat{
    for (jp in 1:P){
      
      # define the two tool vectors
      lbd_jp1 = rep(0,P)
      lbd_jp1[jp] = 1
      lbd_mjp_c = lbd_c
      lbd_mjp_c[jp] = 0 
      lbd_mjp_c = lbd_mjp_c/sum(lbd_mjp_c)
      
      ### Do the optimization over gamma, denoted as x here ###
      res <- nloptr(x0 = 0,
                    eval_f = opt,
                    lb = 0, ub = pi_cand[s],
                    opts = list("algorithm" = 'NLOPT_LN_COBYLA', 'print_level' = 1,
                                "xtol_rel" = 1e-7),
                    lbd_jp1 = lbd_jp1, pi_cand = pi_cand, s = s,lbd_mjp_c =lbd_mjp_c, y = y,
                    jp = jp, lbd = lbd_c, Dmatrix = Dtrain_all)
      
      zopt = res$solution
      
      # update the j'th component of the current solution with zopt
      # lbd_c = zopt * lbd_jp1 + (pi_cand[s] - zopt) * lbd_mjp_c
      lbd_c[jp] = zopt * lbd_jp1[jp] + (pi_cand[s] - zopt) * lbd_mjp_c[jp]
    }
    
    # compute the MSE at the current s
    temp1 = mse(y = y, Dmatrix = Dtune_all, lbd = lbd_c)
    # append the current MSE to mem
    mem = c(mem,temp1)
    # check the convergence condition 
    mem_conv = abs((mem[length(mem)]-mem[(length(mem)-1)])/mem[length(mem)])
    
    if (mem_conv < rpt_tol) break
  }
      
  conv_sln[s,] = lbd_c
}








# perform the modified coordinate descent algo on the tuning sets
foreach (s = 1:S) %:%
  foreach (jp = 1:P, .packages = 'nloptr') %dopar% {
    # Initialize with equal smoothing
    lbd_0 = rep(pi_cand[s]/P, P)
    # update the initial vector as the current solution
    lbd_c = lbd_0
    # define the two tool vectors
    lbd_jp1 = rep(0,P)
    lbd_jp1[jp] = 1
    lbd_mjp_c = lbd_jp_c
    lbd_mjp_c[jp] = 0 
    lbd_mjp_c = lbd_mjp_c/sum(lbd_mjp_c)
    
    ### Do the optimization over gamma, denoted as x here ###
    res <- nloptr(x0 = 0,
                  eval_f = opt,
                  lb = 0, ub = pi_cand[s],
                  opts = list("algorithm" = 'NLOPT_LN_COBYLA', 'print_level' = 1,
                              "xtol_rel" = 1e-8),
                  lbd_jp1 = lbd_jp1, pi_cand = pi_cand,s = s,lbd_mjp_c =lbd_mjp_c, y = y,
                  jp = jp, lbd = lbd_c, Dmatrix = Dtrain_all)
    
    zopt = res$solution

    # find the lbd_gamma first
  }    
    
    # opt until convergence
    repeat{
      res <- nloptr(x0 = 0.01*pi_cand[s],
                    eval_f = opt,
                    lb = 0, ub = pi_cand[s],
                    opts = list("algorithm" = 'NLOPT_LN_COBYLA', 'print_level' = 1,
                                "xtol_rel" = 1e-8),
                    lbd_jp1 = lbd_jp1, pi_cand = pi_cand,s = s,lbd_mjp_c =lbd_mjp_c, y = y,
                    jp = jp, lbd = lbd_c, Dmatrix = Dtrain_all)
      zopt = res$solution
      
      # update the j'th component of the current solution with zopt
      lbd_c[jp] = zopt * lbd_jp1[jp] + (pi_cand[s] - zopt) * lbd_mjp_c[jp]
      
      # update the corresponding lbd_mjp_c
      lbd_mjp_c = lbd_jp_c
      lbd_mjp_c[jp] = 0 
      lbd_mjp_c = lbd_mjp_c/sum(lbd_mjp_c)
      
      if (ip >5) break
    }

  
    


#### trail codes for optimization
ipp = 2
ip = 1



z = 0.5 * pi_cand[s]

lbd_gamma = (z*lbd_jp1+(pi_cand[s]-z)*lbd_mjp_c)



exp_part = exp(-2*temp)

obj = vector()

 

End_time = Sys.time()
End_time - start_time




lbd_gamma = findlbdgamma(z,lbd_jp1, pi_cand,s,lbd_mjp_c)


for (s in 1:S) {
  for (jp in 2:2)  {
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
    
    ### Do the optimization over gamma, denoted as z here ###
    res <- nloptr(x0 = 0,
                  eval_f = obj_f,
                  lb = 0, ub = pi_cand[s],
                  opts = list("algorithm" = 'NLOPT_LN_COBYLA', 'print_level' = 1,
                              "xtol_rel" = 1e-8),
                  lbd_jp1 = lbd_jp1, pi_cand = pi_cand,s = s,lbd_mjp_c =lbd_mjp_c,
                  jp = jp ,n1 = n1, n2 = n2, lbd_gamma = lbd_gamma,Dmatrix = Dtune_all)
    
    zopt[s] = res$solution
    OPT[s] = res$objective
    # find the lbd_gamma first
  }    
}



