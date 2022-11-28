library(data.table)
library(fdapace)
library(doParallel)
library(foreach)
library(ggplot2)
library(nloptr)
library(doSNOW)

source("def_func_copy.R")


ncore = 12

## to unregister the cl
# stopCluster(cl)
registerDoSEQ()
cl <- makeCluster(ncore)
registerDoParallel(cl, cores = ncore)
# stopCluster(cl)

#### generating a toy dense functional dataset from scratch####

# Set the number of subjects (N) and the
# number of predictors (P)
# number of measurements per predictor(M)
repeatlimit <- 1
N <- 1000;
n1 = .6*N; n2 = .2*N; n3 = .2*N;
P <- 6;
M <- 101;
me <- 0.1;

sparsity_lvl = seq(15,30,by=5)
sigma = 0.1;
rpt_tol = 5e-3

set.seed(1120)

# Define the continuum
sc <- seq(0,10,length.out = M)

# Define the tau candidates
tau_cand = seq(0,2, length.out = 21)

S = length(tau_cand)

# Define the results containers
opted_tau = rep(NA, repeatlimit)
opted_lbd = array(NA, dim = c(repeatlimit,P))
tune_MSE_all = array(NA, dim = c(repeatlimit, S))
lbd_c_all = array(NA, dim=c(repeatlimit,S,P))

repetationnumber = 0

for (repetationnumber in 1:repeatlimit){
  rep_st_time = Sys.time()
  message('Current simulation index is:', repetationnumber)
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
  C = list()
  pred_zeta = list()
  pred_eta = list()
  
  #### Generate the FPC scores, Xtrues, then sparsify, eventually and/or ME ####
  ## Work on the all the predictors one by one
  
  for (diter in 1:P){
    # Generate eigvenvalues for the current predictor
    eigenvalue_sq[,diter] = sort(sample(seq(1,5), K[diter], replace = FALSE),decreasing = T)
    # Create FPC scores
    Ksi_all[[diter]] <- gen_Xi(N,K,diter, eigenvalue_sq)
    # Generate the functional X
    xTrue_all[[diter]] <- gen_X(N, K, M, diter, Ksi_all, sc)
    # Spasify the functional X
    xSparse_all[[diter]] <- Sparsify(samp = xTrue_all[[diter]], pts = sc, sparsity = sparsity_lvl)
    
    # Add measurement errors in X on X_true
    xSparse_all[[diter]]$xNoisy = lapply(xSparse_all[[diter]]$Ly, function(x) x + rnorm(length(x), mean = 0, sd = me))
    
    # Standardize the xNoisy
    mean = mean(unlist(xSparse_all[[diter]]$xNoisy))
    sd = sd(unlist(xSparse_all[[diter]]$xNoisy))
    xSparse_all[[diter]]$xNoisy = lapply(xSparse_all[[diter]]$xNoisy, function(x) (x - mean)/sd)
    
  }
  
  #### Calculate yTrue and yNoisy ####
  # f1 is sin then find the average
  # f2 is qudratic then find the average
  # f3 is cubic then find the average
  
  ## 1st predictor
  diter = 1
  
  ## f1 is abs then find the average
  # apply the defined fucntion on this X
  temp = t(sapply(xSparse_all[[diter]]$xNoisy,f3))
  # get the mean over row for each sample
  temp1  = unlist(lapply(temp, mean))
  # update the yTrue
  yTrue = yTrue + t(temp1)
  
  ## 2nd predictor
  diter = 2
  
  # m2 is qudratic then find the average
  # apply the defined fucntion on this X
  temp = t(sapply(xSparse_all[[diter]]$xNoisy,f1))
  # get the mean over row for each sample
  temp1  = unlist(lapply(temp, mean))
  # update the yTrue
  yTrue <- yTrue + t(temp1)
  
  ## 3rd predictor
  diter = 3
  
  # m3 is cubic then find the average
  # apply the defined fucntion on this X
  temp = t(sapply(xSparse_all[[diter]]$xNoisy,f1))
  # get the mean over row for each sample
  temp1  = unlist(lapply(temp, median))
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
  
  # start_time = Sys.time()
  
  for (diter in 1:P){
    ### separate the simulated data into training, tuning and testing
    training[[diter]] = list(x= xSparse_all[[diter]]$xNoisy[c(1:n1)],
                             t= xSparse_all[[diter]]$Lt[c(1:n1)], y = y[1:n1], xi = Ksi_all[[diter]][c(1:n1),])
    tuning[[diter]] = list(x= xSparse_all[[diter]]$xNoisy[c((n1+1):(n1+n2))],
                           t= xSparse_all[[diter]]$Lt[c((n1+1):(n1+n2))], y = y[(n1+1):(n1+n2)],xi = Ksi_all[[diter]][(n1+1):(n1+n2),])
    testing[[diter]] = list(x= xSparse_all[[diter]]$xNoisy[c((n1+n2+1):N)],
                            t= xSparse_all[[diter]]$Lt[c((n1+n2+1):N)], y = y[(n1+n2+1):N],xi = Ksi_all[[diter]][(n1+n2+1):N,])
    
    
  }
  
  # End_time = Sys.time()
  # End_time - start_time
  
  ####Running FPCA on training set ####
  
  # Do FPCA on this sparse sample with 99% FVE cut
  # Notice that sparse FPCA will smooth the data internally (Yao et al., 2005)
  # Smoothing is the main computational cost behind sparse FPCA
  clusterExport(cl, c('M','training','FPCA'))
  
  FPCAsparse<- parLapply(cl, 1:P, fun = function(diter, train = training){FPCA(training[[diter]]$x, training[[diter]]$t,
                                                                               list(plot = F, dataType = 'Sparse',FVEthreshold = 0.975, kernel = 'epan',
                                                                                    nRegGrid = M, useBinnedData = 'OFF', useBinnedCov = T, methodBwCov = 'GMeanAndGCV'))})
  
  # Get the central C matrices
  clusterExport(cl, c('FPCAsparse'))
  
  temp <- parLapply(cl,1:P, fun = function(diter, FPCAobj = FPCAsparse){t(FPCAobj[[diter]]$phi) %*% FPCAobj[[diter]]$phi})
  
  clusterExport(cl, c('FPCAsparse','temp'))
  
  C <- parLapply(cl,1:P, fun = function(diter, FPCAobj = FPCAsparse, te = temp){te[[diter]] %*% diag(1/FPCAobj[[diter]]$lambda) %*% te[[diter]]})
  
  #### Predict the first K* FPC scores according to conditional expectation on the tuning ####
  clusterExport(cl, c('FPCAsparse','tuning'))
  pred_zeta <- parLapply(cl,1:P, fun = function(diter, FPCAobj = FPCAsparse, tem = tuning){
    predict(FPCAobj[[diter]], tem[[diter]]$x, tem[[diter]]$t, K=FPCAobj[[diter]]$selectK)
  })
  
  #### Predict the first K* FPC scores according to conditional expectation on the tuning ####
  clusterExport(cl, c('FPCAsparse','testing'))
  pred_eta <- parLapply(cl,1:P, fun = function(diter, FPCAobj = FPCAsparse, tem = testing){
    predict(FPCAobj[[diter]], tem[[diter]]$x, tem[[diter]]$t, K=FPCAobj[[diter]]$selectK)
  })
  
  # stopCluster(cl)
  # 
  # cl <- makeCluster(ncore)
  # registerDoParallel(cl, cores = ncore)
  
  #### Calculate all the pair-wise squared distances between all training samples ####
  Dtrain_all = foreach (diter = 1:P, .combine='cbind') %dopar% {
    
    
    Dtrain <- matrix(0, nrow = n1, ncol = n1)
    
    for (i in 1:ncol(Dtrain)){
      for (j in 1:nrow(Dtrain)){
        temp = FPCAsparse[[diter]]$xiEst[j,] - FPCAsparse[[diter]]$xiEst[i,]
        Dtrain[j,i] = (t(temp) %*% C[[diter]] %*% temp) /((M-1)/10)^2
      }
    }
    Dtrain
  }
  
  dim(Dtrain_all) = c(n1,n1,P)
  
  #### Calculate all the pair-wise distances between all training samples and all tuning samples####
  Dtune_all = foreach(diter = 1:P, .combine='cbind') %dopar% {
    
    Dtune <- matrix(0,nrow = n1, ncol = n2)
    
    for (i in 1:ncol(Dtune)){
      for (j in 1:nrow(Dtune)){
        temp =  FPCAsparse[[diter]]$xiEst[j,] - pred_zeta[[diter]][i,]
        Dtune[j,i] = (t(temp) %*% C[[diter]] %*% temp) /((M-1)/10)^2
      }
    }
    Dtune
  }
  
  dim(Dtune_all) = c(n1,n2,P)
  
  #### Calculate all the pair-wise distances between all training samples and all testing samples####
  Dtest_all = foreach(diter = 1:P, .combine='cbind') %dopar% {
    
    Dtest <- matrix(0,nrow = n1, ncol = n3)
    
    for (i in 1:ncol(Dtest)){
      for (j in 1:nrow(Dtest)){
        temp =  FPCAsparse[[diter]]$xiEst[j,] - pred_eta[[diter]][i,]
        Dtest[j,i] = (t(temp) %*% C[[diter]] %*% temp) /((M-1)/10)^2
      }
    }
    Dtest
  }
  
  dim(Dtest_all) = c(n1,n3,P)
  
  #### save the objects ####
  # save(Dtrain_all, file = '~/Functional Kernel Regression/FKRS/Mar 6 p=6 sparse(10-30) results/Dtrain_all.Rdata')
  # save(Dtune_all, file = '~/Functional Kernel Regression/FKRS/Mar 6 p=6 sparse(10-30) results/Dtune_all.Rdata')
  # save(Dtest_all, file = '~/Functional Kernel Regression/FKRS/Mar 6 p=6 sparse(10-30) results/Dtest_all.Rdata')
  # save(Ksi_all, file = '~/Functional Kernel Regression/FKRS/Mar 6 p=6 sparse(10-30) results/Ksi_all.Rdata')
  # save(pred_zeta, file = '~/Functional Kernel Regression/FKRS/Mar 6 p=6 sparse(10-30) results/pred_zeta.Rdata')
  # save(pred_eta, file = '~/Functional Kernel Regression/FKRS/Mar 6 p=6 sparse(10-30) results/pred_eta.Rdata')
  # save(FPCAsparse, file = '~/Functional Kernel Regression/FKRS/Mar 6 p=6 sparse(10-30) results/FPCAsparse.Rdata')
  # save(xSparse_all, file = '~/Functional Kernel Regression/FKRS/Mar 6 p=6 sparse(10-30) results/xSparse_all.Rdata')
  # save(y, file = '~/Functional Kernel Regression/FKRS/Mar 6 p=6 sparse(10-30) results/y.Rdata')
  
  
  #### Find the optimal total budget over the tuning sets ####
  
  
  MSE = rep(0,S)
  zopt = rep(0,S)
  conv_sln = matrix(0,nrow = S, ncol = P)
  
  
  for (s in 1:S){
    # Initialize with equal smoothing
    lbd_0 = rep(tau_cand[s]/P, P)
    # update the initial vector as the current solution
    lbd_c = lbd_0
    # define a vector to record the MSE for the current s
    mem = rep(0,S)
    # define a vector to record the converging res$objective
    convg_obj = rep(0,S)
    
    message(tau_cand[s],' is currently set as the total budget')
    
    
    iterationnumber = 1
    
    repeat{
      iterationnumber = iterationnumber + 1
      message('Current iteration number is:', (iterationnumber-1), ' for budget = ', tau_cand[s])
      for (jp in 1:P){
        message(jp,'-th component of lambda is updating')
        # define the two tool vectors
        lbd_jp1 = rep(0,P)
        lbd_jp1[jp] = 1
        lbd_mjp_c = lbd_c
        lbd_mjp_c[jp] = 0
        if (sum(lbd_mjp_c) == 0){
          lbd_mjp_c = rep(0,P)
        } else {
          lbd_mjp_c = lbd_mjp_c/sum(lbd_mjp_c)
        }
        
        ### Do the optimization over gamma, denoted as x here ###
        res <- nloptr(x0 = runif(1)*tau_cand[s],
                      eval_f = opt1,
                      lb = 0, ub = tau_cand[s],
                      opts = list("algorithm" = 'NLOPT_LN_COBYLA', 'print_level' = 0,
                                  "xtol_rel" = 1e-3),
                      lbd_jp1 = lbd_jp1, tau_cand = tau_cand, s = s,lbd_mjp_c =lbd_mjp_c,
                      y = y, lbd = lbd_c, Dmatrix = Dtrain_all)
        
        zopt = res$solution
        
        # update thecurrent solution with zopt
        lbd_c = zopt * lbd_jp1 + (tau_cand[s] - zopt) * lbd_mjp_c
      }
      print(paste('Current slution is',lbd_c,'at', (iterationnumber-1), 'th iteration.'))
      
      convg_obj[iterationnumber] = res$objective
      
      
      if (abs(convg_obj[iterationnumber]- convg_obj[iterationnumber-1]) < rpt_tol) {
        print(paste('Convergence achieved at the', (iterationnumber-1), 'th iteration.'))
        break
      }
    }
    
    # store the current lbd vector
    conv_sln[s,] = lbd_c
    
    # compute the MSE at the current s over tune
    MSE[s] = mse_tune(y = y, lbd = lbd_c)
  }
  
  
  # save the tuning and optimized lambda
  
  # save(MSE, file = '~/Functional Kernel Regression/FKRS/Mar 6 p=6 sparse(10-30) results/MSE.Rdata')
  # save(conv_sln, file = '~/Functional Kernel Regression/FKRS/Mar 6 p=6 sparse(10-30) results/conv_sln.Rdata')
  
  # find out the index corresponding to the smallest MSE
  # which.min(MSE)
  
  # the selected tau is
  # print(paste('The selected tau by MSE over tuning is', tau_cand[which.min(MSE)]))
  # 1.35
  opted_tau[repetationnumber] = tau_cand[which.min(MSE)]
  
  tune_MSE_all[repetationnumber,] = MSE
  
  ### The optimized lambda
  # print(paste('The optimized lambda is', round(conv_sln[which.min(MSE),],4)))
  # round(conv_sln[which.min(MSE),],4)
  # 0.2016 0.0000 1.6984 0.0000 0.0000 0.0000
  opted_lbd[repetationnumber,] =  conv_sln[which.min(MSE),]
  
  lbd_c_all[repetationnumber,,] = conv_sln
  
  message('Current simulation of index:', repetationnumber, ' is done using time', (Sys.time()-rep_st_time) )
  
}


save(opted_tau, file = '~/Functional Kernel Regression/FKRS/0310 x1q x2l results/opted_tau.Rdata')
save(opted_lbd, file = '~/Functional Kernel Regression/FKRS/0310 x1q x2l results/opted_lbd.Rdata')
save(tune_MSE_all, file = '~/Functional Kernel Regression/FKRS/0310 x1q x2l results/tune_MSE_all.Rdata')
save(lbd_c_all, file = '~/Functional Kernel Regression/FKRS/0310 x1q x2l results/lbd_c_all.Rdata')

