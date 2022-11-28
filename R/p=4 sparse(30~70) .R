library(data.table)
library(fdapace)
library(doParallel)
library(foreach)
library(ggplot2)
library(nloptr)
library(tidyverse)

source("def_func_copy.R")


## to unregister the cl
registerDoSEQ()
cl <- makeCluster(3)
registerDoParallel(cl)
# stopCluster(cl)

#### generating a toy dense functional dataset from scratch####

# Set the number of subjects (N) and the
# number of predictors (P)
# number of measurements per predictor(M)
N <- 40;
n1 = .6*N; n2 = .2*N; n3 = .2*N;
P <- 4;
M <- 101;
me <- 0.1;

sparsity_lvl = seq(30,70,by=10)
sigma = 1;
rpt_tol = 1e-4



set.seed(123)

# Define the continuum
s <- seq(0,10,length.out = M)


#### Define the needed hyperparameters and fucntions ####
# number of eigen components per predictor
K <- rep(2,P)

#### Define all the needed matrices ####
xTrue_all = list()
xSparse_all = list()
Ksi_all = list()
Dtrain_all = array(0, dim = c(n1,n1,P))
Dtune_all = array(0, dim = c(n1,n2,P))
Dtest_all = array(0, dim = c(n1,n3,P))
Dtrain_all_true = array(0, dim = c(n1,n1,P))
Dtune_all_true = array(0, dim = c(n1,n2,P))
Dtest_all_true = array(0, dim = c(n1,n3,P))
Dtrain_all_residual_perc = array(0, dim = c(n1,n1,P))
Dtune_all_residual_perc = array(0, dim = c(n1,n2,P))
Dtest_all_residual_perc = array(0, dim = c(n1,n3,P))
eigFunct = list()
eigenvalue_sq = matrix(0,nrow = 2, ncol = P)
yTrue = 0


training = list()
tuning = list()
testing = list()
FPCAsparse = list()
S = list()
pred_zeta = list()
pred_eta = list()

#### Generate the FPC scores, Xtrues, then sparsify, eventually and/or ME ####
## Work on the all the predictors one by one

for (diter in 1:P){
  # Generate eigvenvalues for the current predictor
  eigenvalue_sq[,diter] = sort(sample(seq(1,6), K[diter], replace = FALSE),decreasing = T)
  # Create FPC scores
  Ksi_all[[diter]] <- gen_Xi(N,K,diter, eigenvalue_sq)
  # Generate the functional X
  xTrue_all[[diter]] <- gen_X(N, K, M, diter, Ksi_all, s)
  # Spasify the functional X
  xSparse_all[[diter]] <- Sparsify(samp = xTrue_all[[diter]], pts = s, sparsity = sparsity_lvl)
  
  # Add measurement errors in X on X_true
  xSparse_all[[diter]]$xNoisy = lapply(xSparse_all[[diter]]$Ly, function(x) x + rnorm(length(x), mean = 0, sd = me))
}

#### Calculate yTrue and yNoisy ####
# f1 is sin then find the average
# f2 is qudratic then find the average
# f3 is cubic then find the average

## 1st predictor
diter = 1

## f1 is sin then find the average
# apply the defined fucntion on this X
temp = t(sapply(xSparse_all[[diter]]$xNoisy,sin))
# # convert list to matrix with fillers = 0
temp = list.as.matrix(temp, byrow = T, filler = 0)
# get the mean over row for each sample
temp1  = apply(temp,1,mean)
# update the yTrue
yTrue = yTrue + t(temp1)

## 2nd predictor
diter = 2

# m2 is qudratic then find the average
# apply the defined fucntion on this X
temp = t(sapply(xSparse_all[[diter]]$xNoisy,quadratic))
# # convert list to matrix with fillers = 0
temp = list.as.matrix(temp, byrow = T, filler = 0)
# get the mean over row for each sample
temp1  = apply(temp,1,mean)
# update the yTrue
yTrue <- yTrue + t(temp1)

## 3rd predictor
diter = 3

# m3 is cos then find the average
# apply the defined fucntion on this X
temp = t(sapply(xSparse_all[[diter]]$xNoisy,cos))
# # convert list to matrix with fillers = 0
temp = list.as.matrix(temp, byrow = T, filler = 1)
# get the mean over row for each sample
temp1  = apply(temp,1,mean)
# update the yTrue (nothing b/c x3 is not important)
# yTrue <- yTrue + t(temp1)

# ## 4th predictor
# diter = 4
#
# # m4 is abs then find the average
# # apply the defined fucntion on this X
# temp = t(sapply(xSparse_all[[diter]],abs))
# # # convert list to matrix with fillers = 0
# # temp = list.as.matrix(temp, byrow = T, filler = 0)
# # get the mean over row for each sample
# temp1  = apply(temp,1,mean)
# # update the yTrue (nothing b/c x4 is not important)
# # yTrue <- yTrue + t(temp1)



#### add error term to yTrue as y ####
y <- yTrue +  rnorm(N,0,sigma)



#### Now, find all the distance matrices one by one ####
diter = 1
start_time = Sys.time()

for (diter in 1:P){
  ### separate the simulated data into training, tuning and testing
  training[[diter]] = list(x= xSparse_all[[diter]]$xNoisy[c(1:n1)],
                           t= xSparse_all[[diter]]$Lt[c(1:n1)], y = y[1:n1], xi = Ksi_all[[diter]][c(1:n1),])
  tuning[[diter]] = list(x= xSparse_all[[diter]]$xNoisy[c((n1+1):(n1+n2))],
                         t= xSparse_all[[diter]]$Lt[c((n1+1):(n1+n2))], y = y[(n1+1):(n1+n2)],xi = Ksi_all[[diter]][(n1+1):(n1+n2),])
  testing[[diter]] = list(x= xSparse_all[[diter]]$xNoisy[c((n1+n2+1):N)],
                          t= xSparse_all[[diter]]$Lt[c((n1+n2+1):N)], y = y[(n1+n2+1):N],xi = Ksi_all[[diter]][(n1+n2+1):N,])
  
  ####Running FPCA on training set ####
  
  # Do FPCA on this sparse sample with 99% FVE cut
  # Notice that sparse FPCA will smooth the data internally (Yao et al., 2005)
  # Smoothing is the main computational cost behind sparse FPCA
  FPCAsparse[[diter]] <- FPCA(training[[diter]]$x, training[[diter]]$t,
                              list(plot = F, dataType = 'Sparse',FVEthreshold = 0.98, kernel = 'epan',
                                   nRegGrid = M, useBinnedData = 'OFF', useBinnedCov = T, methodBwCov = 'GMeanAndGCV'))
  
  # Get the central S matrix
  temp <- t(FPCAsparse[[diter]]$phi) %*% FPCAsparse[[diter]]$phi
  S[[diter]] <- temp %*% diag(1/FPCAsparse[[diter]]$lambda) %*% temp
  
  #### Predict the first K* FPC scores according to conditional expectation on the tuning ####
  pred_zeta[[diter]] <- predict(FPCAsparse[[diter]], tuning[[diter]]$x, tuning[[diter]]$t, K=FPCAsparse[[diter]]$selectK)
  
  #### Predict the first K* FPC scores according to conditional expectation on the tuning ####
  pred_eta[[diter]] <- predict(FPCAsparse[[diter]], testing[[diter]]$x, testing[[diter]]$t, K=FPCAsparse[[diter]]$selectK)
  
  #### Calculate all the pair-wise distances between all training samples ####
  Dtrain <- matrix(0, nrow = length(training[[diter]]$x), ncol = length(training[[diter]]$x))
  # start_time = Sys.time()
  Dtrain = foreach (i = 1:ncol(Dtrain), .combine=cbind) %:%
    foreach (j = 1:nrow(Dtrain), .combine=rbind) %dopar%{
      temp = FPCAsparse[[diter]]$xiEst[j,] - FPCAsparse[[diter]]$xiEst[i,]
      Dtrain[j,i] = t(temp) %*% S[[diter]] %*% temp
    }
  # End_time = Sys.time()
  # End_time - start_time
  
  Dtrain = Dtrain/((M-1)/10)^2
  
  #### verify the correctness of distances by finding the true distances ####
  
  Dtrain_true <- matrix(0, nrow = length(training[[diter]]$x), ncol = length(training[[diter]]$x))
  start_time = Sys.time()
  Dtrain_true = foreach (i = 1:ncol(Dtrain_true), .combine=cbind) %:%
    foreach (j = 1:nrow(Dtrain_true), .combine=rbind) %dopar%{
      temp = Ksi_all[[diter]][j,] - Ksi_all[[diter]][i,]
      Dtrain_true[j,i] = t(temp) %*% diag(1/eigenvalue_sq[,diter]^2) %*% temp
    }
  End_time = Sys.time()
  End_time - start_time
  print(Sys.time())
  Dtrain_residual_perc <- abs(Dtrain - Dtrain_true)/Dtrain_true
  quantile(Dtrain_residual_perc, na.rm = T)
  
  
  #### Calculate all the pair-wise distances between all training samples and all tuning samples####
  Dtune <- matrix(nrow = length(training[[diter]]$x), ncol = length(tuning[[diter]]$x))
  # start_time = Sys.time()
  Dtune = foreach (i = 1:ncol(Dtune), .combine=cbind) %:%
    foreach (j = 1:nrow(Dtune), .combine=rbind) %dopar%{
      temp = FPCAsparse[[diter]]$xiEst[j,] - pred_zeta[[diter]][i,]
      Dtune[j,i] = t(temp) %*% S[[diter]] %*% temp
    }
  # End_time = Sys.time()
  # End_time - start_time
  
  Dtune = Dtune/((M-1)/10)^2
  
  #### verify the correctness of distances by finding the true Dtune ####
  
  Dtune_true <- matrix(nrow = length(training[[diter]]$x), ncol = length(tuning[[diter]]$x))
  start_time = Sys.time()
  Dtune_true = foreach (i = 1:ncol(Dtune_true), .combine=cbind) %:%
    foreach (j = 1:nrow(Dtune_true), .combine=rbind) %dopar%{
      temp = Ksi_all[[diter]][j,] - Ksi_all[[diter]][(n1+i),]
      Dtune_true[j,i] = t(temp) %*% diag(1/eigenvalue_sq[,diter]^2) %*% temp
    }
  End_time = Sys.time()
  End_time - start_time
  print(Sys.time())
  Dtune_residual_perc <- abs(Dtune - Dtune_true)/Dtune_true
  quantile(Dtune_residual_perc, na.rm = T)
  
  
  
  #### Calculate all the pair-wise distances between all training samples and all testing samples####
  Dtest <- matrix(nrow = length(training[[diter]]$x), ncol = length(testing[[diter]]$x))
  #start_time = Sys.time()
  Dtest = foreach (i = 1:ncol(Dtest), .combine=cbind) %:%
    foreach (j = 1:nrow(Dtest), .combine=rbind) %dopar%{
      temp = FPCAsparse[[diter]]$xiEst[j,] - pred_eta[[diter]][i,]
      Dtest[j,i] = t(temp) %*% S[[diter]] %*% temp
    }
  #End_time = Sys.time()
  #End_time - start_time
  
  Dtest = Dtest/((M-1)/10)^2
  
  #### verify the correctness of distances by finding the true Dtest ####
  
  Dtest_true <- matrix(nrow = length(training[[diter]]$x), ncol = length(testing[[diter]]$x))
  start_time = Sys.time()
  Dtest_true = foreach (i = 1:ncol(Dtest_true), .combine=cbind) %:%
    foreach (j = 1:nrow(Dtest_true), .combine=rbind) %dopar%{
      temp = Ksi_all[[diter]][j,] - Ksi_all[[diter]][(n1+n2+i),]
      Dtest_true[j,i] = t(temp) %*% diag(1/eigenvalue_sq[,diter]^2) %*% temp
    }
  End_time = Sys.time()
  End_time - start_time
  print(Sys.time())
  Dtest_residual_perc <- abs(Dtest - Dtest_true)/Dtest_true
  quantile(Dtest_residual_perc, na.rm = T)
  
  
  Dtrain_all[[diter]] = Dtrain
  Dtrain_all_true[[diter]] = Dtrain_true
  Dtrain_all_residual_perc[[diter]] = Dtrain_residual_perc
  Dtune_all[[diter]] = Dtune
  Dtune_all_true[[diter]] = Dtune_true
  Dtune_all_residual_perc[[diter]] = Dtune_residual_perc
  Dtest_all[[diter]] = Dtest
  Dtest_all_true[[diter]] = Dtest_true
  Dtest_all_residual_perc[[diter]] = Dtest_residual_perc
  
}

End_time = Sys.time()
End_time - start_time

#### save the objects ####
save(Dtrain_all, file = '~/Functional Kernel Regression/FKRS/Feb 23 p=3 sparse(30-70) results/Dtrain_all.Rdata')
save(Dtrain_all_true, file = '~/Functional Kernel Regression/FKRS/Feb 23 p=3 sparse(30-70) results/Dtrain_all_true.Rdata')
save(Dtrain_all_residual_perc, file = '~/Functional Kernel Regression/FKRS/Feb 23 p=3 sparse(30-70) results/Dtrain_all_residual_perc.Rdata')

save(Dtune_all, file = '~/Functional Kernel Regression/FKRS/Feb 23 p=3 sparse(30-70) results/Dtune_all.Rdata')
save(Dtune_all_true, file = '~/Functional Kernel Regression/FKRS/Feb 23 p=3 sparse(30-70) results/Dtune_all_true.Rdata')
save(Dtune_all_residual_perc, file = '~/Functional Kernel Regression/FKRS/Feb 23 p=3 sparse(30-70) results/Dtune_all_residual_perc.Rdata')

save(Dtest_all, file = '~/Functional Kernel Regression/FKRS/Feb 23 p=3 sparse(30-70) results/Dtest_all.Rdata')
save(Dtest_all_true, file = '~/Functional Kernel Regression/FKRS/Feb 23 p=3 sparse(30-70) results/Dtest_all_true.Rdata')
save(Dtest_all_residual_perc, file = '~/Functional Kernel Regression/FKRS/Feb 23 p=3 sparse(30-70) results/Dtest_all_residual_perc.Rdata')







#### Find the optimal total budget over the tuning sets ####

tau_cand = seq(0, 20, length.out = 101)

S = length(tau_cand)

MSE = rep(0,S)
zopt = rep(0,S)
conv_sln = matrix(0,nrow = S, ncol = P)

s=21;jp=1



for (s in 2:S){
  # Initialize with equal smoothing
  lbd_0 = rep(tau_cand[s]/P, P)
  # update the initial vector as the current solution
  lbd_c = lbd_0
  # define a vector to record the MSE for the current s
  mem = c(1e15,1e12)
  # define a vector to record the converging res$objective
  convg_obj = rep(0,50)
  
  #
  repeatnumber = 1
  
  repeat{
    repeatnumber = repeatnumber + 1
    message('Current repetation number is:', (repeatnumber-1))
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
                    lb = 0, ub = tau_cand[s],
                    opts = list("algorithm" = 'NLOPT_GN_DIRECT', 'print_level' = 1,
                                "xtol_rel" = 1e-4),
                    lbd_jp1 = lbd_jp1, tau_cand = tau_cand, s = s,lbd_mjp_c =lbd_mjp_c, y = y,
                    jp = jp, lbd = lbd_c, Dmatrix = Dtrain_all)
      
      zopt = res$solution
      
      # update the j'th component of the current solution with zopt
      lbd_c = zopt * lbd_jp1 + (tau_cand[s] - zopt) * lbd_mjp_c
    }
    
    convg_obj[repeatnumber] = res$objective
    
    
    if (abs(convg_obj[repeatnumber]- convg_obj[repeatnumber-1]) < rpt_tol) {
      print(paste('Convergence achieved at the', (repeatnumber-1), 'th iteration.'))
      break
    }
  }
  
  
  # compute the MSE at the current s over tune
  temp1 = mse_tune(y = y, lbd = lbd_c)
  # append the current MSE to mem
  mem = c(mem,temp1)
  # check the convergence condition
  mem_conv = abs((mem[length(mem)]-mem[(length(mem)-1)])/mem[length(mem)])
  
  conv_sln[s,] = lbd_c
}








# perform the modified coordinate descent algo on the tuning sets
foreach (s = 1:S) %:%
  foreach (jp = 1:P, .packages = 'nloptr') %dopar% {
    # Initialize with equal smoothing
    lbd_0 = rep(tau_cand[s]/P, P)
    # update the initial vector as the current solution
    lbd_c = lbd_0
    # define the two tool vectors
    lbd_jp1 = rep(0,P)
    lbd_jp1[jp] = 1
    lbd_mjp_c = lbd_c
    lbd_mjp_c[jp] = 0
    lbd_mjp_c = lbd_mjp_c/sum(lbd_mjp_c)
    
    ### Do the optimization over gamma, denoted as x here ###
    res <- nloptr(x0 = 0,
                  eval_f = opt,
                  lb = 0, ub = tau_cand[s],
                  opts = list("algorithm" = 'NLOPT_LN_COBYLA', 'print_level' = 1,
                              "xtol_rel" = 1e-6),
                  lbd_jp1 = lbd_jp1, tau_cand = tau_cand,s = s,lbd_mjp_c =lbd_mjp_c, y = y,
                  jp = jp, lbd = lbd_c, Dmatrix = Dtrain_all)
    
    zopt = res$solution
    
    # find the lbd_gamma first
  }

# opt until convergence
repeat{
  res <- nloptr(x0 = 0.01*tau_cand[s],
                eval_f = opt,
                lb = 0, ub = tau_cand[s],
                opts = list("algorithm" = 'NLOPT_LN_COBYLA', 'print_level' = 1,
                            "xtol_rel" = 1e-8),
                lbd_jp1 = lbd_jp1, tau_cand = tau_cand,s = s,lbd_mjp_c =lbd_mjp_c, y = y,
                jp = jp, lbd = lbd_c, Dmatrix = Dtrain_all)
  zopt = res$solution
  
  # update the j'th component of the current solution with zopt
  lbd_c[jp] = zopt * lbd_jp1[jp] + (tau_cand[s] - zopt) * lbd_mjp_c[jp]
  
  # update the corresponding lbd_mjp_c
  lbd_mjp_c = lbd_jp_c
  lbd_mjp_c[jp] = 0
  lbd_mjp_c = lbd_mjp_c/sum(lbd_mjp_c)
  
  if (ip >5) break
}





#### trail codes for optimization
ipp = 2
ip = 1



z = 0.5 * tau_cand[s]

lbd_gamma = (z*lbd_jp1+(tau_cand[s]-z)*lbd_mjp_c)



exp_part = exp(-2*temp)

obj = vector()



End_time = Sys.time()
End_time - start_time




lbd_gamma = findlbdgamma(z,lbd_jp1, tau_cand,s,lbd_mjp_c)


for (s in 1:S) {
  for (jp in 1:P)  {
    # Initialize with equal smoothing
    lbd_jp = rep(tau_cand[s]/P, P)
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
                  lb = 0, ub = tau_cand[s],
                  opts = list("algorithm" = 'NLOPT_LN_COBYLA', 'print_level' = 1,
                              "xtol_rel" = 1e-8),
                  lbd_jp1 = lbd_jp1, tau_cand = tau_cand,s = s,lbd_mjp_c =lbd_mjp_c,
                  jp = jp ,n1 = n1, n2 = n2, lbd_gamma = lbd_gamma,Dmatrix = Dtune_all)
    
    zopt[s] = res$solution
    OPT[s] = res$objective
    # find the lbd_gamma first
  }
}

#### verify the correctness of distances by finding the true distances ####

Dtrain_true <- matrix(nrow = length(training$x), ncol = length(training$x))
start_time = Sys.time()
Dtrain_true = foreach (i = 1:ncol(Dtrain_true), .combine=cbind) %:%
  foreach (j = 1:nrow(Dtrain_true), .combine=rbind) %dopar%{
    temp = Ksi_all[[1]][j,] - Ksi_all[[1]][i,]
    Dtrain_true[j,i] = t(temp) %*% diag(1/eigenvalue_sq[,diter]^2) %*% temp
  }
End_time = Sys.time()
End_time - start_time
print(Sys.time())
Dtrain_residual_perc <- abs(Dtrain - Dtrain_true)/Dtrain_true
quantile(Dtrain_residual_perc, na.rm = T)


rbenchmark::benchmark(opt = {opt(x,lbd_jp1, tau_cand, s,lbd_mjp_c,y, lbd, Dmatrix=Dtrain_all)},
                      opt1 = {opt1(x,lbd_jp1, tau_cand, s,lbd_mjp_c,y, lbd, Dmatrix=Dtrain_all)},
                      opt2 = {opt2(x,lbd_jp1, tau_cand, s,lbd_mjp_c,y, lbd, Dmatrix=Dtrain_all)}, replications = 5000)
