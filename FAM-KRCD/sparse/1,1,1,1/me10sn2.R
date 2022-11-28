library(data.table)
library(fdapace)
library(doParallel)
library(foreach)
library(nloptr)
library(doSNOW)
library(fda.usc)

set.seed(271967)

# env
Run.env = 'Statlab'

#### Define the simulation parameters ####

# Tell me the model index to run
model.index <- 'FAM'
# Tell me the number of repetition
repeatlimit <- 5
# Set the number of subjects (N) and the size of train, tune, test
N <- 400; n1 <- .5*N; n2 <- .5*N; n3 <- 1*N;
# number of functional predictors (P) and number of non-fucntional predictors (Q)
P <- 5; Q <- 5; R <-0;
# number of the grid points per f predictor to sample from (M) 
M <- 51;
# set the variance of eigen components of all the functional predictors
FPCv <- c(16,4,1,0.25);
# Decide the S/N ratio of the measurement error (as to X) as me
me.ratio <- 10;    me <- sqrt(sum(FPCv)/me.ratio);
# Decide the S/N ratio of the regression model error (as to Y) as noise
noise.ratio <- 2; 
# define the weights
weights_aa=c(1,1,1,1)
# Set the method to choose the smoothing bandwidth when estimating the mean func and Cov func during FPCA
methodforbwCov = 'GCV'; methodforbwMu = 'GCV';
# method to calculate the distance, Dv use the inverse estimated var, D1 use just 1
method.for.Dist = 'Dv'
# min and max number of eigen components per function predictor to estimate during FPCA
min.K.used <- 3;  max.K.used <- 4; number.of.K.used <- max.K.used - min.K.used + 1;
# Set the convergence condition of the optimization during the coordinate descent algo
rpt_tol <- 1e-2
# Define the tau candidates as the searching grid during coordinate descent algo
if (method.for.Dist == 'D1'){
  tau_cand <- seq(0.6, 1.9, 0.05)
} else if (method.for.Dist == 'Dv') {
  tau_cand <- seq(5.5, 9.0, 0.1)} # tau_cand <- c(seq(7.6,8.4,0.2), seq(8.6,10, 0.1), seq(10.2, 11, .2))}


# Get the length of the searching grid of tau
S <- length(tau_cand)

#### Define the functional predictor's parameters ####

# number of obs per subject curve
sparsity_lvl <- seq(6,18,1);
# set the variance of eigen components of all the functional predictors
FPCv <- c(16,4,1,0.25);
# mean function of the functional predictor (l+s:linear+sine)
meanfstat = "l+s";
# Set the continuum of the functional curves are defined on
sc <- seq(0,10,length.out = M)
# number of eigen components per func pred to generate, K is a P-dim vector (potentially different # for each func pred)
numberofK <- 4; K <- rep(numberofK,P)

#### load the functions and modules ####
if (Run.env == 'PC') {
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/true.model.R")
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/ksi.to.x.R")
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/ksi.to.yTrue.R")
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/FPCAmaxtime.R")
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/pred.FPCs.R")
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/get.C.mat.R")
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/est.withz.R")
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/opt2.R")
  source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/mse.prdc.R")
} else if (Run.env == 'Statlab') {
  source("~/Functional Kernel Regression/functions/true.model.R")
  source("~/Functional Kernel Regression/functions/ksi.to.x.R")
  source("~/Functional Kernel Regression/functions/ksi.to.yTrue.R")
  source("~/Functional Kernel Regression/functions/FPCAmaxtime.R")
  source("~/Functional Kernel Regression/functions/pred.FPCs.R")
  source("~/Functional Kernel Regression/functions/get.C.mat.R")
  source("~/Functional Kernel Regression/functions/est.withz.R")
  source("~/Functional Kernel Regression/functions/opt2.R")
  source("~/Functional Kernel Regression/functions/mse.prdc.R")
}

#### Define the results containers ####
# (left contains results from all possible K choices during FPCA, right contains the best result out of the left one)

# Record all the converged lambda vector for each possible tau value
lbd_c_all = array(0, dim=c(repeatlimit,number.of.K.used,S,(P+Q)));     lbd_c_all.chosenone = array(0, dim=c(repeatlimit,S,(P+Q)));
# Record the MSE from tune set, used to pick the best K choice
tune_MSE_all = array(0, dim = c(repeatlimit,number.of.K.used, S)); tune_MSE_all.chosenone = array(0, dim = c(repeatlimit, S));
# Record the optimized tau value  
opted_tau = array(0, dim = c(repeatlimit,number.of.K.used));       opted_tau.chosenone = array(0, dim = c(repeatlimit));
# Record the optimized lambda vector corresponding to optimized tau
opted_lbd = array(0, dim = c(repeatlimit,number.of.K.used, (P+Q)));    opted_lbd.chosenone = array(0, dim = c(repeatlimit, (P+Q)));
# Record the MSE from test set, used to compare performance
final_eval = array(0, dim = c(repeatlimit,number.of.K.used));      final_eval.chosenone = array(0, dim = c(repeatlimit));

K_result = array(0, dim = c(repeatlimit,number.of.K.used));      K_result.chosenone = array(0, dim = c(repeatlimit));
csm = array(0, dim = c(repeatlimit,number.of.K.used));              csm.chosenone = array(0, dim = c(repeatlimit));    csm.all = array(0, dim = c(repeatlimit));

# Record the distances between all pairs of functional curves in train and all vector predictors in train
# Those distance matrices are 3-D tensors, 1st dim is length of the train, 2nd dim is length of train/tune/test, 3rd dim is the length of lambda vector
Dtrain_all = array(0, dim = c(n1,n1,P)); DzTrain = array(0, dim = c(n1,n1,Q))
Dtune_all = array(0, dim = c(n1,n2,P)); DzTune = array(0, dim = c(n1,n2,Q))
Dtest_all = array(0, dim = c(n1,n3,P)); DzTest = array(0, dim = c(n1,n3,Q))


#### Prep for parallel computing ####

# number of cores to call
ncore = P
# Create & Register cluster 
registerDoSEQ()
cl <- makeCluster(ncore, type = "SOCK")
registerDoParallel(cl, cores = ncore)


#### Start the Simu loop of repetation ####

for (repetationnumber in 1:repeatlimit){
  
  rep_st_time = Sys.time()
  
  message('Current simulation index is:', repetationnumber)
  
  #### Define all the needed matrices to store intermediate results####
  xTrue_all = list(); xNoisy_all = list(); xSparse_all = list(); Ksi_all = list();
  eigFunct = list(); eigenvalue_sq = matrix(0,nrow = numberofK, ncol = P);
  training = list(); tuning = list(); testing = list();
  FPCAsparse = list(); C = list(); pred_zeta = list(); pred_eta = list();
  Zmatrix = matrix(0, nrow = N, ncol = Q); yTrue = 0; 
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
    # Spasify the functional X
    xSparse_all[[diter]] <- Sparsify(samp = xTrue_all[[diter]], pts = sc, sparsity = sparsity_lvl)
    # Add measurement errors in X on X_true
    xSparse_all[[diter]]$xNoisy <- lapply(xSparse_all[[diter]]$Ly, function(x) x + rnorm(length(x), mean = 0, sd = me))
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
  
  #### add regression error to yTrue as y ####
  noise <- sqrt(8/noise.ratio)
  y <- yTrue +  rnorm(N,0,noise)
  
  #### Partition every functional predictor's ksi, t, and X to the train, tune and test set ####
  condition = TRUE
  while (condition == TRUE){
    
    ### sample out the index of the training set
    id_train = sample(seq(1,N),n1)
    
    # sample out the index of the tuning set
    id_tune = setdiff(seq(1,N),id_train)
    
    # get the index of testing set by set substraction
    id_test = seq(1,N)
    
    # generate the indicator of each predictor's status of range
    temp = rep(TRUE,P)
    
    for (diter in 1:P){
      ### separate the simulated data into training, tuning and testing
      training[[diter]] = list(x= xSparse_all[[diter]]$xNoisy[id_train], 
                               t= xSparse_all[[diter]]$Lt[id_train],  xi = Ksi_all[[diter]][id_train,])
      tuning[[diter]] = list(x= xSparse_all[[diter]]$xNoisy[id_tune], 
                             t= xSparse_all[[diter]]$Lt[id_tune], xi = Ksi_all[[diter]][id_tune,])
      testing[[diter]] = list(x= xSparse_all[[diter]]$xNoisy[id_test], 
                              t= xSparse_all[[diter]]$Lt[id_test], xi = Ksi_all[[diter]][id_test,])
      
      # get the lower range of training's time
      min_t_train = range(unlist(training[[diter]]$t))[1]
      
      # get the upper range of training's time
      max_t_train = range(unlist(training[[diter]]$t))[2]
      
      # verify the lower range of training's time equals the min of T
      cond1  = min_t_train > min(sc)
      
      # verify the upper range of training's time equals the min of T
      cond2  = max_t_train < max(sc)
      
      temp[diter] = cond1 | cond2 
    }
    if(sum(temp)>0) {
      condition = T
    } else {
      condition = F
    }
  }
  
  #### partition Zmatrix to zTrain, zTune, zTest in the same fashion as the functional predictors ####
  zTrain <- Zmatrix[id_train,];   zTune <- Zmatrix[id_tune,];   zTest <- Zmatrix[id_test,];
  
  #### partition y to yTrain, yTune, yTest in the same fashion as the functional predictors ####
  yTrain <- y[id_train];  yTune <- y[id_tune];   yTest <- y[id_test];
  
  
  #### Calculate the 3 distance matrices standardized by the median distance within Z, non-functional predictors ####
  for (diter in 1:Q){
    
    # find the difference b/t M1 and M2, then calculate their L1 distance by dividing the corresponding variance
    temp = zTrain[,diter] %*% t(rep(1,n1))
    DzTrain[,,diter] =  ((temp - t(temp))/sqrt(1))**2
    DzTrain[,,diter] = DzTrain[,,diter] / median(DzTrain[,,diter])
    
    temp1 = zTrain[,diter] %*% t(rep(1,n2))
    temp2 = rep(1,n1) %*% t(zTune[,diter])
    DzTune[,,diter] = (temp1 - temp2/sqrt(1))**2
    DzTune[,,diter] = DzTune[,,diter] / median(DzTune[,,diter])
    
    temp1 = zTrain[,diter] %*% t(rep(1,n3))
    temp2 = rep(1,n1) %*% t(zTest[,diter])
    DzTest[,,diter] = (temp1 - temp2/sqrt(1))**2
    DzTest[,,diter] = DzTest[,,diter] / median(DzTest[,,diter])
  }
  
  
  #### Running FPCA on training set ####
  # Smoothing is the main computational cost behind sparse FPCA
  
  FPCAsparse <- foreach(diter = 1:P, .packages = 'fdapace') %dopar% {
    zs = FPCA(training[[diter]]$x, training[[diter]]$t,
              list(plot = F, dataType = 'Sparse', methodSelectK = max.K.used, kernel = 'gauss',
                   nRegGrid = M, useBinnedData = 'OFF', useBinnedCov = T, maxK = K[diter], lean = T,
                   methodBwMu = methodforbwMu, methodBwCov = methodforbwCov))
    zs
  } 
  # Print the FPCA time
  FPCAmaxtime()
  
  ### Use FPCA on train result to predict FPC scores on tune and test ###
  
  # Predict the first K* FPC scores according to conditional expectation on the tuning 
  pred_zeta = pred.FPCs(FPCAsparse,tuning)
  
  # Predict the first K* FPC scores according to conditional expectation on the testing 
  pred_eta = pred.FPCs(FPCAsparse,testing)
  
  
  #### star the k-loop to go through all K candidate during FPCA ####
  
  for (kiter in min.K.used:max.K.used){
    # Specify the current working K
    K.in.use = kiter
    K_result[repetationnumber, kiter-min.K.used] = K.in.use 
    
    #### Calculate all the pair-wise squared distances between all training samples ####
    
    # Get the central C matrices
    C = get.C.mat(method.for.Dist, FPCAsparse, K.in.use)
    
    #### Calculate all the pair-wise squared distances between all training samples ####
    Dtrain1 = foreach (diter = 1:P, .combine='cbind') %dopar% {
      
      Dtrain <- matrix(0, nrow = n1, ncol = n1)
      
      for (i in 1:ncol(Dtrain)){
        for (j in 1:nrow(Dtrain)){
          temp = FPCAsparse[[diter]]$xiEst[j,1:K.in.use] - FPCAsparse[[diter]]$xiEst[i,1:K.in.use]
          Dtrain[j,i] = (t(temp) %*% C[[diter]] %*% temp) 
        }
      }
      Dtrain = Dtrain / median(Dtrain)
      Dtrain
    }
    
    dim(Dtrain1) = c(n1,n1,P)
    Dtrain_all = Dtrain1
    
    
    #### Calculate all the pair-wise distances between all training samples and all tuning samples####
    Dtune1 = foreach(diter = 1:P, .combine='cbind') %dopar% {
      
      Dtune <- matrix(0,nrow = n1, ncol = n2)
      
      for (i in 1:ncol(Dtune)){
        for (j in 1:nrow(Dtune)){
          temp =  FPCAsparse[[diter]]$xiEst[j,1:K.in.use] - pred_zeta[[diter]][i,1:K.in.use]
          Dtune[j,i] = (t(temp) %*% C[[diter]] %*% temp)
        }
      }
      Dtune = Dtune / median(Dtune)
      Dtune
    }
    
    dim(Dtune1) = c(n1,n2,P)
    Dtune_all = Dtune1
    
    #### Calculate all the pair-wise distances between all training samples and all testing samples####
    Dtest1 = foreach(diter = 1:P, .combine='cbind') %dopar% {
      
      Dtest <- matrix(0,nrow = n1, ncol = n3)
      
      for (i in 1:ncol(Dtest)){
        for (j in 1:nrow(Dtest)){
          temp =  FPCAsparse[[diter]]$xiEst[j,1:K.in.use] - pred_eta[[diter]][i,1:K.in.use]
          Dtest[j,i] = (t(temp) %*% C[[diter]] %*% temp)
        }
      }
      Dtest = Dtest / median(Dtest)
      Dtest
    }
    
    dim(Dtest1) = c(n1,n3,P)
    Dtest_all = Dtest1
    
    #### Find the optimal total budget over the tuning sets ####
    CD_st_time = Sys.time()
    
    MSE = rep(0,S);   zopt = rep(0,S);
    conv_sln = matrix(0,nrow = S, ncol = P+Q)
    
    opt_res = foreach (s = 1:S, .packages = 'nloptr', .combine = 'rbind')%dopar%{
      # Initialize with equal smoothing
      lbd_0 = rep(tau_cand[s]/(P+Q), P+Q)
      # update the initial vector as the current solution
      lbd_c = lbd_0
      # define a vector to record the MSE for the current s
      mem = rep(0,S)
      # define a vector to record the converging res$objective
      convg_lbd = rep(0,(P+Q))
      
      message(tau_cand[s],' is currently set as the total budget')
      
      
      iterationnumber = 1
      repeat{
        iterationnumber = iterationnumber + 1
        message('Current iteration number is:', (iterationnumber-1), ' for budget = ', tau_cand[s])
        for (jp in 1:(P+Q)){
          message(jp,'-th component of lambda is updating')
          # define the two tool vectors
          lbd_jp1 = rep(0,P+Q)
          lbd_jp1[jp] = 1
          lbd_mjp_c = lbd_c
          lbd_mjp_c[jp] = 0
          if (sum(lbd_mjp_c) == 0){
            lbd_mjp_c = rep(0,P+Q)
          } else {
            lbd_mjp_c = lbd_mjp_c/sum(lbd_mjp_c)
          }
          
          ### Do the optimization over gamma, denoted as x here ###
          res <- nloptr(x0 = 0,
                        eval_f = opt2,
                        lb = 0, ub = tau_cand[s],
                        opts = list("algorithm" = 'NLOPT_GN_DIRECT_NOSCAL', 'print_level' = 0,
                                    "xtol_rel" = 1e-4, maxeval = 2e2),
                        lbd_jp1 = lbd_jp1, tau_cand = tau_cand, s = s,lbd_mjp_c =lbd_mjp_c,
                        y = yTrain, Dmatrix = Dtrain_all, DmatrixZ = DzTrain)
          
          zopt = res$solution
          
          # update thecurrent solution with zopt
          lbd_c = zopt * lbd_jp1 + (tau_cand[s] - zopt) * lbd_mjp_c
        }
        print(paste('Current slution is',lbd_c,'at', (iterationnumber-1), 'th iteration.'))
        
        convg_lbd = rbind(convg_lbd,lbd_c)
        
        
        if (dist(rbind(convg_lbd[iterationnumber,], convg_lbd[(iterationnumber-1),]), method = "manhattan") < rpt_tol*tau_cand[s]) {
          print(paste('Convergence achieved at the', (iterationnumber-1), 'th iteration.'))
          break
        } else if (iterationnumber > 5) {
          print(paste('Convergence still cannot be achieved after ', (iterationnumber-1), 'iterations.'))
          break
        }
      }
      
      
      # store the current lbd vector
      conv_sln[s,] = lbd_c
      
      # compute the MSE at the current s over tune
      result = c(lbd_c, mse.tune(yTune))
      result
    }
    message('Coordinate descent algo used time', (Sys.time()-CD_st_time), 'minutes.')
    
    conv_sln = opt_res[,1:(P+Q)]
    
    MSE = opt_res[,(P+Q+1)]
    
    # save the tuning and optimized lambda
    
    # the selected tau is
    opted_tau[repetationnumber, kiter - (min.K.used-1)] = tau_cand[which.min(MSE)]
    
    tune_MSE_all[repetationnumber, kiter - (min.K.used-1),] = MSE
    
    ### The optimized lambda 
    opted_lbd[repetationnumber, kiter - (min.K.used-1), ] =  conv_sln[which.min(MSE),]
    
    lbd_c_all[repetationnumber, kiter - (min.K.used-1),,] = conv_sln
    
    
    #### Do the prediction over the testing set ####
    
    final_eval[repetationnumber, kiter - (min.K.used-1)] = mse.test(y=yTest)
    
    message('Current simulation of index:', repetationnumber, ', with K = ', K.in.use, ' is done using time', (Sys.time()-rep_st_time))
    message('Current time is ', Sys.time())
    
    
  } # end of kiter loop
  
  
}# end of repetationnumber loop



for (repetationnumber in 1:repeatlimit){
  
  for (kiter in 1:number.of.K.used){
    for (diter in 1:(P+Q)){
      if (opted_lbd[repetationnumber,kiter,diter] < 1e-2 * opted_tau[repetationnumber,kiter]){
        opted_lbd[repetationnumber,kiter,diter] = 0
      }
    } # end for diter loop  
    
    # check if all relavent predictors are selected
    
    condition1 = identical(which(opted_lbd[repetationnumber,kiter,] > 0), relavent.indecies)
    
    # check if all irrelavent predictors are excluded
    
    condition2 = identical(which(opted_lbd[repetationnumber,kiter,] == 0), irrelavent.indecies)
    
    # if both conditions are met, this is a correctly selected model
    if (condition1 & condition2 ==T) {
      csm[repetationnumber, kiter] = 1
    }
  }# end for kiter loop  
  
  
  # extract the indexes with the smallest prediction error among the all the K choices
  csm.chosenone[repetationnumber] = which.min(final_eval[repetationnumber, ])
  
  # csm.all records the correctly selected models out of the total number of repetitions
  csm.all[repetationnumber] = csm[repetationnumber, csm.chosenone[repetationnumber]]
  
  
  # get the chosen stuffs
  opted_tau.chosenone[repetationnumber] = opted_tau[repetationnumber, csm.chosenone[repetationnumber]]
  opted_lbd.chosenone[repetationnumber,] = opted_lbd[repetationnumber, csm.chosenone[repetationnumber], ]
  tune_MSE_all.chosenone[repetationnumber,] = tune_MSE_all[repetationnumber, csm.chosenone[repetationnumber],]
  lbd_c_all.chosenone[repetationnumber,,] = lbd_c_all[repetationnumber, csm.chosenone[repetationnumber],,]
  final_eval.chosenone[repetationnumber] = final_eval[repetationnumber, csm.chosenone[repetationnumber]]
  K_result.chosenone[repetationnumber] = K_result[repetationnumber, csm.chosenone[repetationnumber]]
}# end for repetationnumber loop    

#### save the objects ####

if(method.for.Dist == 'D1'){
} else if (method.for.Dist == 'Dv'){
  save(opted_tau, file = '~/Functional Kernel Regression/FAM-KRCD/sparse/1,1,1,1/me10sn2/opted_tau.Rdata')
  save(opted_tau.chosenone, file = '~/Functional Kernel Regression/FAM-KRCD/sparse/1,1,1,1/me10sn2/opted_tau.chosenone.Rdata')
  save(opted_lbd, file = '~/Functional Kernel Regression/FAM-KRCD/sparse/1,1,1,1/me10sn2/opted_lbd.Rdata')
  save(opted_lbd.chosenone, file = '~/Functional Kernel Regression/FAM-KRCD/sparse/1,1,1,1/me10sn2/opted_lbd.chosenone.Rdata')
  save(tune_MSE_all, file = '~/Functional Kernel Regression/FAM-KRCD/sparse/1,1,1,1/me10sn2/tune_MSE_all.Rdata')
  save(tune_MSE_all.chosenone, file = '~/Functional Kernel Regression/FAM-KRCD/sparse/1,1,1,1/me10sn2/tune_MSE_all.chosenone.Rdata')
  save(lbd_c_all, file = '~/Functional Kernel Regression/FAM-KRCD/sparse/1,1,1,1/me10sn2/lbd_c_all.Rdata')
  save(lbd_c_all.chosenone, file = '~/Functional Kernel Regression/FAM-KRCD/sparse/1,1,1,1/me10sn2/lbd_c_all.chosenone.Rdata')
  save(final_eval, file = '~/Functional Kernel Regression/FAM-KRCD/sparse/1,1,1,1/me10sn2/final_eval.Rdata')
  save(final_eval.chosenone, file = '~/Functional Kernel Regression/FAM-KRCD/sparse/1,1,1,1/me10sn2/final_eval.chosen.Rdata')
  save(K_result, file = '~/Functional Kernel Regression/FAM-KRCD/sparse/1,1,1,1/me10sn2/K_result.Rdata')
  save(K_result.chosenone, file = '~/Functional Kernel Regression/FAM-KRCD/sparse/1,1,1,1/me10sn2/K_result.chosenone.Rdata')
  save(csm, file = '~/Functional Kernel Regression/FAM-KRCD/sparse/1,1,1,1/me10sn2/csm.Rdata')
  save(csm.chosenone, file = '~/Functional Kernel Regression/FAM-KRCD/sparse/1,1,1,1/me10sn2/csm.chosenone.Rdata')
  save(csm.all, file = '~/Functional Kernel Regression/FAM-KRCD/sparse/1,1,1,1/me10sn2/csm.all.Rdata')
  
}


# record the number of K used throughout all trials
K_selected = rep(0,number.of.K.used)

for (kiter in 1:number.of.K.used){
  K_selected[kiter] = sum(K_result.chosenone == kiter+1)
}
K_selected

# record the number of K used throughout all correctly selected trials
K_correctly_selected = rep(0,number.of.K.used)

for (kiter in 1:number.of.K.used){
  K_correctly_selected[kiter] = sum(K_result.chosenone[csm.all == 1] == kiter+1)
}
K_correctly_selected


# record the number of chosen predictors throughout all trials
num_selected = colSums(opted_lbd.chosenone >  0)
num_selected

csm


print(round(cbind(opted_lbd.chosenone,opted_tau.chosenone),3))
print(paste(round(mean(final_eval.chosenone),2) , round(sd(final_eval.chosenone)/sqrt(repeatlimit),2)))

message(cat(num_selected, sep = ' & '),' & \\textbf{', sum(csm.all), '} &', 
        round(mean(final_eval.chosenone),3), '(', round(sd(final_eval.chosenone)/sqrt(repeatlimit),3),') \\', '\\')



message(cat(K_selected, sep = ' & '), ' \\', '\\')
message(cat(K_correctly_selected, sep = ' & '), ' \\', '\\')


stopCluster(cl)
