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
repeatlimit <- 100
# Set the number of subjects (N) and the size of train, tune, test
N <- 400; n1 <- .25*N; n2 <- .25*N; n3 <- .5*N;
# number of functional predictors (P) and number of non-fucntional predictors (Q)
P <- 5; Q <- 5; R <-0;
# number of the grid points per f predictor to sample from (M) 
M <- 51;
# set the variance of eigen components of all the functional predictors
FPCv <- c(16,4,1,0.25);
# Decide the S/N ratio of the measurement error (as to X) as me
me.ratio <- 8;    me <- sqrt(sum(FPCv)/me.ratio);
# Decide the S/N ratio of the regression model error (as to Y) as noise
noise.ratio <- 8; 


#### Define the functional predictor's parameters ####

# number of obs per subject curve
sparsity_lvl <- seq(M,M,1);
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
  source("~/Functional Kernel Regression/FAM-comparison/fregre.gsam.vs1.R")
  source("~/Functional Kernel Regression/FAM-comparison/gam.vs.R")
  source("~/Functional Kernel Regression/FAM-comparison/select.glm.R")
}


#### Prep for parallel computing ####

# number of cores to call
ncore = P
# Create & Register cluster 
registerDoSEQ()
cl <- makeCluster(ncore, type = "SOCK")
registerDoParallel(cl, cores = ncore)

# define the weights
weights_aa=c(1,1,1,1)

matres=array(0,dim=c(repeatlimit,10))
RMSE = array(0, dim=c(repeatlimit,1))
dimnames(matres)=list(1:repeatlimit,c(paste0("X",1:5),paste0("Z",1:5)))

repetationnumber = 1

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
  X1 = xSparse_all[[1]]$xNoisy;
  X2 = xSparse_all[[2]]$xNoisy;
  X3 = xSparse_all[[3]]$xNoisy;
  X4 = xSparse_all[[4]]$xNoisy;
  X5 = xSparse_all[[5]]$xNoisy;
  X1 = matrix(unlist(X1), ncol = M, byrow =T)
  X2 = matrix(unlist(X2), ncol = M, byrow =T)
  X3 = matrix(unlist(X3), ncol = M, byrow =T)
  X4 = matrix(unlist(X4), ncol = M, byrow =T)
  X5 = matrix(unlist(X5), ncol = M, byrow =T)
  X1 = fdata(X1, argvals = sc)
  X2 = fdata(X2, argvals = sc)
  X3 = fdata(X3, argvals = sc)
  X4 = fdata(X4, argvals = sc)
  X5 = fdata(X5, argvals = sc)
  
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
  # g.var <- sd(yTrue)**2; noise <- sqrt(g.var/noise.ratio);
  noise <- sqrt(8/noise.ratio)
  y <- yTrue +  rnorm(N,0,noise)
  
  #### Put the x and y together ####
  ldatasim=list(df=data.frame(Z1=Z1,Z2=Z2,Z3=Z3,Z4=Z4,Z5=Z5),
                X1=X1,X2=X2, X3=X3,X4=X4,X5=X5)
  ldatasim$df$y = y
  #### star the k-loop to go through all K candidate during FPCA ####
  
  b.x=list(X1=create.pc.basis(ldatasim$X1,1:4),X2=create.pc.basis(ldatasim$X2,1:4),X3=create.pc.basis(ldatasim$X3,1:4),
           X4=create.pc.basis(ldatasim$X4,1:4),X5=create.pc.basis(ldatasim$X5,1:4))
  
  res=fregre.gsam.vs1("y",ldatasim,b.x)
  #res=fda.usc:::fregre.gsam.vs(data=ldatasim,y="y")  # Error with basis.x=b.x?
  
  nam=unique(substr(attr(res$model$terms,"term.labels"),1,2))
  #nam=unique(substr(attr(res$terms,"term.labels"),1,2))
  matres[repetationnumber,dimnames(matres)[[2]] %in% nam]=1
  
  RMSE[repetationnumber,1] = mean(res$model$residuals^2)
  
  message('Current simulation of index:', repetationnumber, ' is done using time', (Sys.time()-rep_st_time))
  message('Current time is ', Sys.time())
}# end of repetationnumber loop

ARMSE = t(apply(RMSE,2,mean))
selection.res = t(apply(matres,2,mean))


save(ARMSE, file = '~/Functional Kernel Regression/FAM-comparison/indpt/1,1,1,1/me8sn8/ARMSE.Rdata')
save(selection.res, file = '~/Functional Kernel Regression/FAM-comparison/indpt/1,1,1,1/me8sn8/selection.res.Rdata')

message(cat(selection.res, sep = ' & '),' & ', 
        round(ARMSE,3), '(', round(sd(RMSE)/sqrt(repeatlimit),3),') \\', '\\')