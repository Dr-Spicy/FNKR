library(data.table)
library(fdapace)
library(doParallel)
library(foreach)
library(nloptr)
library(doSNOW)

set.seed(27)

# env
Run.env = 'PC'
# which y
response = "Ash"
# sparsify?
Sparsified.status = "Dense"
# full x?
full.spec = "No"
##################################################################################################################################

##################################################################################################################################

#### Define the simulation parameters ####

# Set the number of subjects (N) and the size of train, tune, test
N <- length(Ash); n1 <- 150; n2 <- 50; n3 <- 50;
# number of functional predictors (P) and number of non-fucntional predictors (Q)
P <- 7; Q <- ncol(Proc.pooled);
# Get the right type of fluoresence data
if (full.spec == "Yes") {
    fluo = fluo.full
} else if (full.spec == "No") {
    fluo = fluo.trunc
}
# Set the method to choose the smoothing bandwidth when estimating the mean func and Cov func during FPCA
methodforbwCov = 'GCV'; methodforbwMu = 'GCV';
# method to calculate the distance, Dv use the inverse estimated var, D1 use just 1
method.for.Dist = 'Dv'
# min and max number of eigen components per function predictor to estimate during FPCA
min.K.used <- 2;  max.K.used <- 4; number.of.K.used <- max.K.used - min.K.used + 1;
# Set the convergence condition of the optimization during the coordinate descent algo
rpt_tol <- 1e-2
# Define the tau candidates as the searching grid during coordinate descent algo
if (method.for.Dist == 'D1'){
    tau_cand <- seq(0.6, 1.9, 0.05)
} else if (method.for.Dist == 'Dv') {
    tau_cand <- c( seq(6,11, 0.1))} # tau_cand <- c(seq(7.6,8.4,0.2), seq(8.6,10, 0.1), seq(10.2, 11, .2))}


# Get the length of the searching grid of tau
S <- length(tau_cand)



#### Define the functional predictor's parameters ####
# number of eigen components per func pred to generate, K is a P-dim vector (potentially different # for each func pred)
numberofK <- 4; K <- rep(numberofK,P)
# Set the continuum of the functional curves are defined on
sc <- as.numeric(dimnames(fluo)[[2]])
# number of obs per subject curve
if (Sparsified.status == 'Sparse') {
    sparsity_lvl <- seq(20, 30, 1);
} else if (Sparsified.status == 'Dense') {
    sparsity_lvl <- dim(fluo)[2]
}



#### load the functions and modules ####
if (Run.env == 'PC') {
    
    source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/ksi.to.x.R")
    source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/ksi.to.yTrue.R")
    source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/FPCAmaxtime.R")
    source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/pred.FPCs.R")
    source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/get.C.mat.R")
    source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/est.withz.R")
    source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/opt2.R")
    source("C:/Users/Han/Desktop/Box Sync/Prof. Wu/fdapace/R/mse.prdc.R")
    load("C:/Users/Han/Desktop/Box Sync/Prof. Wu/real data for FKR/sugar_Process/Ash.Rdata")
    load("C:/Users/Han/Desktop/Box Sync/Prof. Wu/real data for FKR/sugar_Process/fluo.full.Rdata")
    load("C:/Users/Han/Desktop/Box Sync/Prof. Wu/real data for FKR/sugar_Process/fluo.trunc.Rdata")
    load("C:/Users/Han/Desktop/Box Sync/Prof. Wu/real data for FKR/sugar_Process/Proc.pooled.Rdata")
} else if (Run.env == 'Statlab') {
    
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
lbd_c_all = array(0, dim=c(number.of.K.used,S,(P+Q)));     lbd_c_all.chosenone = array(0, dim=c(S,(P+Q)));
# Record the MSE from tune set, used to pick the best K choice
tune_MSE_all = array(0, dim = c(number.of.K.used, S)); tune_MSE_all.chosenone = array(0, dim = c(S));
# Record the optimized tau value  
opted_tau = array(0, dim = c(number.of.K.used));       opted_tau.chosenone = array(0, dim = c(1));
# Record the optimized lambda vector corresponding to optimized tau
opted_lbd = array(0, dim = c(number.of.K.used, (P+Q)));    opted_lbd.chosenone = array(0, dim = c((P+Q)));
# Record the MSE from test set, used to compare performance
final_eval = array(0, dim = c(number.of.K.used));      final_eval.chosenone = array(0, dim = c(1));

K_result = array(0, dim = c(number.of.K.used));          K_result.chosenone = array(0, dim = c(1));
csm = array(0, dim = c(number.of.K.used));              csm.chosenone = array(0, dim = c(1));    csm.all = array(0, dim = c(1));

# Record the distances between all pairs of functional curves in train and all vector predictors in train
# Those distance matrices are 3-D tensors, 1st dim is length of the train, 2nd dim is length of train/tune/test, 3rd dim is the length of lambda vector
Dtrain_all = array(0, dim = c(n1,n1,P)); DzTrain = array(0, dim = c(n1,n1,Q))
Dtune_all = array(0, dim = c(n1,n2,P)); DzTune = array(0, dim = c(n1,n2,Q))
Dtest_all = array(0, dim = c(n1,n3,P)); DzTest = array(0, dim = c(n1,n3,Q))


#### Prep for parallel computing ####

# number of cores to call
ncore = P*2
# Create & Register cluster 
registerDoSEQ()
cl <- makeCluster(ncore, type = "SOCK")
registerDoParallel(cl, cores = ncore)



##################################################################################################################################
# start to VS

rep_st_time = Sys.time()

#### Define all the needed matrices to store intermediate results####
xTrue_all = list(); xSparse_all = list();
training = list(); tuning = list(); testing = list();
FPCAsparse = list(); C = list(); pred_zeta = list(); pred_eta = list();

#### Get the sparsified X ####
## Work on the all the predictors one by one

for (diter in 1:P){
    # Convert the functional X from fluo
    xTrue_all[[diter]] <- fluo[,,diter]
    # Spasify the functional X
    xSparse_all[[diter]] <- Sparsify(samp = xTrue_all[[diter]], pts = sc, sparsity = sparsity_lvl)
}


#### partition functional X(fluo) into train, tune, test ####

condition = TRUE
while (condition == TRUE){
    
    ### sample out the index of the training set
    id_train = sample(seq(1,N),n1)
    
    # find the left over set after training
    set1 = setdiff(seq(1,N),id_train)
    
    # sample out the index of the tuning set
    id_tune = sample(set1, n2)
    
    # get the index of testing set by set substraction
    id_test = setdiff(set1, id_tune)
    
    # generate the indicator of each predictor's status of range
    temp = rep(TRUE,P)
    
    for (diter in 1:P){
        ### separate the simulated data into training, tuning and testing
        training[[diter]] = list(x= xSparse_all[[diter]]$Ly[id_train], 
                                 t= xSparse_all[[diter]]$Lt[id_train])
        tuning[[diter]] = list(x= xSparse_all[[diter]]$Ly[id_tune], 
                               t= xSparse_all[[diter]]$Lt[id_tune])
        testing[[diter]] = list(x= xSparse_all[[diter]]$Ly[id_test], 
                                t= xSparse_all[[diter]]$Lt[id_test])
        
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

remove(xTrue_all); remove(xSparse_all);

#### partition Zmatrix to zTrain, zTune, zTest in the same fashion as the functional predictors ####
Zmatrix = Proc.pooled
zTrain <- Zmatrix[id_train,];   zTune <- Zmatrix[id_tune,];   zTest <- Zmatrix[id_test,];

#### partition y to yTrain, yTune, yTest in the same fashion as the functional predictors ####
if (response == "Ash"){
    y = Ash
} else if (response == "Color") {
    y = Color 
}
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
              list(plot = F, dataType = Sparsified.status, methodSelectK = max.K.used, kernel = 'gauss',
                   nRegGrid = 101, useBinnedData = 'OFF', useBinnedCov = T, maxK = K[diter], lean = T,
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

for (kiter in 1:number.of.K.used){
    # Specify the current working K
    K.in.use = kiter+min.K.used-1
    K_result[kiter] = K.in.use 
    
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
                                          "xtol_rel" = 0.5e-3, maxeval = 10e1),
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
    opted_tau[kiter] = tau_cand[which.min(MSE)]
    
    tune_MSE_all[kiter,] = MSE
    
    ### The optimized lambda 
    opted_lbd[ kiter, ] =  conv_sln[which.min(MSE),]
    
    lbd_c_all[kiter,,] = conv_sln
    
    
    #### Do the prediction over the testing set ####
    
    final_eval[kiter] = mse.test(y=yTest)
    
    message('Current simulation with K = ', K.in.use, ' is done using time', (Sys.time()-rep_st_time))
    message('Current time is ', Sys.time())
    
    
} # end of kiter loop







for (kiter in min.K.used:max.K.used-1){
    for (diter in 1:(P+Q)){
        if (opted_lbd[kiter,diter] < rpt_tol * opted_tau[kiter]){
            opted_lbd[kiter,diter] = 0
        }
    } # end for diter loop
    
    
}# end for kiter loop



# extract the indexes with the smallest prediction error among the all the K choices
csm.chosenone = which.min(final_eval)

# csm.all records the correctly selected models out of the total number of repetitions
csm.all = csm[ csm.chosenone]


# get the chosen stuffs
opted_tau.chosenone = opted_tau[csm.chosenone]
opted_lbd.chosenone = opted_lbd[csm.chosenone, ]
tune_MSE_all.chosenone = tune_MSE_all[csm.chosenone,]
lbd_c_all.chosenone = lbd_c_all[csm.chosenone,,]
final_eval.chosenone = final_eval[csm.chosenone]
K_result.chosenone = K_result[csm.chosenone]
# overall variation
sum((y-mean(y))**2)

# prediction error fraction of overall variation 
N*final_eval.chosenone/sum((y-mean(y))**2)  

#### save the objects ####

# if(method.for.Dist == 'D1'){
#     save(opted_tau, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/D1/me100sn4.bigtau/opted_tau.Rdata')
#     save(opted_tau.chosenone, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/D1/me100sn4.bigtau/opted_tau.chosenone.Rdata')
#     save(opted_lbd, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/D1/me100sn4.bigtau/opted_lbd.Rdata')
#     save(opted_lbd.chosenone, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/D1/me100sn4.bigtau/opted_lbd.chosenone.Rdata')
#     save(tune_MSE_all, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/D1/me100sn4.bigtau/tune_MSE_all.Rdata')
#     save(tune_MSE_all.chosenone, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/D1/me100sn4.bigtau/tune_MSE_all.chosenone.Rdata')
#     save(lbd_c_all, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/D1/me100sn4.bigtau/lbd_c_all.Rdata')
#     save(lbd_c_all.chosenone, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/D1/me100sn4.bigtau/lbd_c_all.chosenone.Rdata')
#     save(final_eval, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/D1/me100sn4.bigtau/final_eval.Rdata')
#     save(final_eval.chosenone, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/D1/me100sn4.bigtau/final_eval.chosen.Rdata')
#     save(K_result, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/D1/me100sn4.bigtau/K_result.Rdata')
#     save(K_result.chosenone, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/D1/me100sn4.bigtau/K_result.chosenone.Rdata')
#     save(csm, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/D1/me100sn4.bigtau/csm.Rdata')
#     save(csm.chosenone, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/D1/me100sn4.bigtau/csm.chosenone.Rdata')
# } else if (method.for.Dist == 'Dv'){
#     save(opted_tau, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/Dv/me100sn4.bigtau/opted_tau.Rdata')
#     save(opted_tau.chosenone, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/Dv/me100sn4.bigtau/opted_tau.chosenone.Rdata')
#     save(opted_lbd, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/Dv/me100sn4.bigtau/opted_lbd.Rdata')
#     save(opted_lbd.chosenone, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/Dv/me100sn4.bigtau/opted_lbd.chosenone.Rdata')
#     save(tune_MSE_all, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/Dv/me100sn4.bigtau/tune_MSE_all.Rdata')
#     save(tune_MSE_all.chosenone, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/Dv/me100sn4.bigtau/tune_MSE_all.chosenone.Rdata')
#     save(lbd_c_all, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/Dv/me100sn4.bigtau/lbd_c_all.Rdata')
#     save(lbd_c_all.chosenone, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/Dv/me100sn4.bigtau/lbd_c_all.chosenone.Rdata')
#     save(final_eval, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/Dv/me100sn4.bigtau/final_eval.Rdata')
#     save(final_eval.chosenone, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/Dv/me100sn4.bigtau/final_eval.chosen.Rdata')
#     save(K_result, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/Dv/me100sn4.bigtau/K_result.Rdata')
#     save(K_result.chosenone, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/Dv/me100sn4.bigtau/K_result.chosenone.Rdata')
#     save(csm, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/Dv/me100sn4.bigtau/csm.Rdata')
#     save(csm.chosenone, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/Dv/me100sn4.bigtau/csm.chosenone.Rdata')
#     save(csm.all, file = '~/Functional Kernel Regression/Z10X10/X1L13.X2Q1.sq2X31X42.Z1Q/400-618/Dv/me100sn4.bigtau/csm.all.Rdata')
#     
# }







# record the number of chosen predictors throughout all trials
dim(opted_lbd.chosenone) = c(1,P+Q)
num_selected = colSums(opted_lbd.chosenone >  0)
num_selected

K_result.chosenone


print(round(cbind(opted_lbd.chosenone,opted_tau.chosenone),3))
# print(paste(round(mean(final_eval.chosenone),2) , round(sd(final_eval.chosenone)/sqrt(repeatlimit),2)))

message(cat(num_selected, sep = ' & '),' & \\textbf{', sum(csm.all), '} &', 
        round(mean(final_eval.chosenone),3), '(', round(sd(final_eval.chosenone)/sqrt(1),3),') \\', '\\')



# message(cat(K_selected, sep = ' & '), ' \\', '\\')
# message(cat(K_correctly_selected, sep = ' & '), ' \\', '\\')


stopCluster(cl)

plot(tau_cand, tune_MSE_all.chosenone)
plot(tune_MSE_all.chosenone)





##################################################################################################################################

# axx = list(nticks = 10, range = c(275, 560)) 
# axy = list(nticks = 4, range = c(230, 340))
# axz = list(nticks = 15, range = c(0, 750)) 
# 
# fig <- plot_ly(x = as.numeric(data$EmAx), y = as.numeric(data$ExAx),
#                z = fluo[1,,]) %>% add_surface() %>% layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz))
# 
# fig
