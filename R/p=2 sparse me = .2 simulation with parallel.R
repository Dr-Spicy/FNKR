library(data.table)
library(fdapace)
library(doParallel)
library(foreach)


## to unregister the cl
registerDoSEQ()
cl <- makeCluster(32)
registerDoParallel(cl)


#### generating a toy dense functional dataset from scratch####

# Set the number of subjects (N) and the
# number of predictors (P)
# number of measurements per predictor(M)
N <- 1000;  # when N=500, PACE time = 657s
P <- 2;
M <- 501;
set.seed(100)

# Define the continuum
s <- seq(0,10,length.out = M)

#### Now, work on the 1st predictor ####
# Define the number of eigencomponents
K_t <- 2
# Define the mean and 2 eigencomponents
meanFunct <- function(s) s + sin(s)
eigFunct1 <- function(s) +cos(s*pi/10) / sqrt(5)
eigFunct2 <- function(s) -sin(s*pi/10) / sqrt(5)


# Create FPC scores
lambda_sq = c(3,1)
col1 <- rnorm(N, sd = lambda_sq[1]/2)
col2 <- rnorm(N, sd = lambda_sq[2]/2)
Ksi <- as.matrix(cbind(col1, col2))
Ksi <- apply(Ksi, 2, scale)
Ksi <- Ksi %*% diag(lambda_sq)

# Create X_true based on sampled time as a list
x1True <- Ksi %*% t(matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2)) + t(matrix(rep(meanFunct(s),N), nrow=M))

# Create the Y_true thru a non-linear relationship
yTrue <- apply(exp(Ksi),1,sum)

#### Now, work on the 2nd predictor ####
# Define the number of eigencomponents
K_t <- 2
# Define the mean and 2 eigencomponents
meanFunct <- function(s) s + cos(s)
eigFunct1 <- function(s) -cos(4*s*pi/10) / sqrt(5)
eigFunct2 <- function(s) sin(4*s*pi/10) / sqrt(5)


# Create FPC scores
lambda_sq = c(2,2)
col3 <- rnorm(N, sd = lambda_sq[1]/2)
col4 <- rnorm(N, sd = lambda_sq[2]/2)
Ksi2 <- as.matrix(cbind(col3, col4))
Ksi2 <- apply(Ksi2, 2, scale)
Ksi2 <- Ksi2 %*% diag(lambda_sq)

# Create X_true based on sampled time as a list
x2True <- Ksi2 %*% t(matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2)) + t(matrix(rep(meanFunct(s),N), nrow=M))


# Create the Y_true thru a non-linear relationship
yTrue <- yTrue + apply((Ksi2)^2,1,sum)

# sample number of measurement per subjects
num_m = sample(x=matrix(rep(c(3,4,5),N),nrow= N), size = N)

# sample ti's for each subject based on the # of measurement on [0,10] as a list
sampled_t = list()
for (i in 1:N){
  sampled_t[[i]] = sort(runif(num_m[i],0,10)) # ascending order
}

# Create X_true based on sampled time as a list
xTrue <- list()
for (i in 1:N){
  xTrue[[i]] = Ksi[i,] %*% t(matrix(c(eigFunct1(sampled_t[[i]]),eigFunct2(sampled_t[[i]])), ncol=K_t))+ meanFunct(sampled_t[[i]])
}

# Add measurement errors in X on X_true
xNoisy <- xTrue
for (i in 1:N){
  xNoisy[[i]] = xNoisy[[i]] + rnorm(xTrue[[i]],sd = 0.2)
}

# Create the Y_true thru a non-linear relationship
yTrue <- apply(Ksi[1:N,]^2,1,sum)

### separate the simulated data into training and testing
training = list(x= xNoisy[c(1:(N*.9))], t= sampled_t[c(1:(N*.9))], y = yTrue[1:(N*.9)])
testing = list(x= xNoisy[c((N*.9+1):N)], t= sampled_t[c((N*.9+1):N)], y = yTrue[(N*.9+1):N])

#### running FPCA on a dense dataset with no noise ####
#L3 <- MakeFPCAInputs(IDs = rep(1:N, each=M), tVec=rep(s,N), yVec = t(xTrue))
#FPCAdense <- FPCA(L3$Ly, L3$Lt,
#                  list(plot = T, dataType ='Dense'))


####Running FPCA on training set ####

# Do FPCA on this sparse sample with 95% FVE cut
# Notice that sparse FPCA will smooth the data internally (Yao et al., 2005)
# Smoothing is the main computational cost behind sparse FPCA
FPCAsparse <- FPCA(training$x, training$t,
                   list(plot = T, FVEthreshold = 0.95, nRegGrid = 501, verbose = T))

#### Predict the first K* FPC scores according to conditional expectation on the testing ####
pred_zeta <- predict(FPCAsparse, testing$x, testing$t, K=FPCAsparse$selectK)
#pred_zeta <- apply(pred_zeta, 2, scale)

# Get the central S matrix
temp <- t(FPCAsparse$phi) %*% FPCAsparse$phi
S <- temp %*% diag(1/FPCAsparse$lambda) %*% temp

# Calculate all the pair-wise distances between all training and all testing

Dtbt <- matrix(nrow = length(training$x), ncol = length(testing$x))
start_time = Sys.time()
Dtbt = foreach (i = 1:ncol(Dtbt), .combine=cbind) %:%
  foreach (j = 1:nrow(Dtbt), .combine=rbind) %dopar%{
    temp = FPCAsparse$xiEst[j,] - pred_zeta[i,]
    Dtbt[j,i] = t(temp) %*% S %*% temp
  }
End_time = Sys.time()
End_time - start_time


Dtbt = Dtbt/((M-1)/10)^2


#### save the objects ####
save(FPCAdense1, file = '2x1_sparse_me_FPCA.Rdata')
#FPCAdense <- readRDS('2x1_dense_noerr_FPCA.Rdata')

save(FPCAdense2, file = '2x2_sparse_me_FPCA.Rdata')
#FPCAdense <- readRDS('2x2_dense_noerr_FPCA.Rdata')

save(pred_zeta_1, file = '2x1_sparse_me_pred_zeta.RData')
#load('2x1_dense_noerr_pred_zeta.RData')

save(pred_zeta_2, file = '2x2_sparse_me_pred_zeta.RData')
#load('2x1_dense_noerr_pred_zeta.RData')

save(Dtbt1, file = '2x1_sparse_me_Dtbt.RData')
#load('2x1_sparse_me_Dtbt.RData')

save(Dtbt2, file = '2x2_sparse_me_Dtbt.RData')
#load('2x2_sparse_me_Dtbt.RData')



#### Functional NW estimate ####

# define the bandwidth
h = seq(0,5,by = 0.01) + 0.01

# define the y_hat with different bandwith
y_hat = matrix(NA, ncol = length(h), nrow = 0.1*N)

y_hat <- foreach (i =1:length(h), .combine = 'cbind') %dopar% {
  # get the exp(-d^2/h) matrix
  temp = exp(-Dtbt/(2*h[i]^2))
  # get the predicted y of testing by NW
  y_hat = (t(temp) %*% yTrue[1:(.9*N)]) / apply(temp,2,sum)
}


# define the RMSE
temp <- (y_hat - matrix(rep(yTrue[451:500], length(h)), nrow = 0.1*N))^2
RMSE <- sqrt(apply(temp,2,sum))

# find out the minimum of RMSE and do a finer serach of bandwidth

which(RMSE == min(na.omit(RMSE))) # 18

# define the bandwidth
h = seq(h[17], h[19], by = 0.0001)

# define the y_hat with different bandwith
y_hat = matrix(NA, ncol = length(h), nrow = 0.1*N)

y_hat <- foreach (i =1:length(h), .combine = 'cbind') %dopar% {
  # get the exp(-d^2/h) matrix
  temp = exp(-Dtbt/(2*h[i]^2))
  # get the predicted y of testing by NW
  y_hat = (t(temp) %*% yTrue[1:(.9*N)]) / apply(temp,2,sum)
}


# define the RMSE
temp <- (y_hat - matrix(rep(yTrue[451:500], length(h)), nrow = 0.1*N))^2
RMSE <- sqrt(apply(temp,2,sum))

# find the minimum RMSE
which(RMSE == min(RMSE)) # 58

# Get the optimal h

h[58] # 0.1757

# Get the optimal y_hat
y_hat0 <- y_hat[,58]



#
yTrue[451:500]
y_hat0


