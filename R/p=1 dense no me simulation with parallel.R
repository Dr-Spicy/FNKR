library(data.table)
library(fdapace)
library(doParallel)
library(foreach)
library(ggplot2)
library(nloptr)

## to unregister the cl
registerDoSEQ()
cl <- makeCluster(2)
registerDoParallel(cl)


#### generating a toy dense functional dataset from scratch####

# Set the number of subjects (N) and the
# number of predictors (P)
# number of measurements per predictor(M)
N <- 5000;  # when N=500, PACE time = 657s
P <- 1;
M <- 25001;
set.seed(100)

# Define the continuum
s <- seq(0,2*pi,length.out = M)

#### Now, work on the 1st predictor ####
# Define the number of eigencomponents
K_t <- 4

# Define the mean and 2 eigencomponents
meanFunct <- function(s) s + sin(s)
eigFunct1 <- function(s) cos(s) / sqrt(pi)
eigFunct2 <- function(s) sin(s) / sqrt(pi)
eigFunct3 <- function(s) cos(2*s) / sqrt(pi)
eigFunct4 <- function(s) sin(2*s) / sqrt(pi)

# Create FPC scores
lambda_sq = sort(sample(c(1,2,3), K_t, replace = TRUE),decreasing = T)
Ksi <- matrix(rnorm(N*K_t,mean=0,sd=lambda_sq), N, K_t)
#col1 <- rnorm(N, sd = lambda_sq[1])
#col2 <- rnorm(N, sd = lambda_sq[2])
#Ksi <- as.matrix(cbind(col1, col2))
Ksi <- apply(Ksi, 2, scale)
#Ksi <- Ksi %*% diag(lambda_sq)

# Create X_true based on sampled time as a list
xTrue <- Ksi %*% t(matrix(c(eigFunct1(s),eigFunct2(s),eigFunct3(s),eigFunct4(s)
), ncol=K_t)) + t(matrix(rep(meanFunct(s),N), nrow=M))

#### running FPCA on a dense dataset with no noise ####
#L3 <- MakeFPCAInputs(IDs = rep(1:N, each=M), tVec=rep(s,N), yVec = t(xTrue))
#FPCAdense <- FPCA(L3$Ly, L3$Lt,
#                  list(plot = T, dataType ='Dense'))


### separate the simulated data into training, tuning and testing
training = MakeFPCAInputs(IDs = rep(1:(.6*N), each=M), tVec=rep(s,(.6*N)), yVec = t(xTrue[1:(.6*N),]))
tuning = MakeFPCAInputs(IDs = rep((.6*N+1):(.8*N), each=M), tVec=rep(s,(.2*N)), yVec = t(xTrue[(.6*N+1):(.8*N),]))
testing = MakeFPCAInputs(IDs = rep((.8*N+1):N, each=M), tVec=rep(s,(.2*N)), yVec = t(xTrue[(.8*N+1):N,]))


####Running FPCA on training set ####

# Do FPCA on this sparse sample with 99% FVE cut
# Notice that sparse FPCA will smooth the data internally (Yao et al., 2005)
# Smoothing is the main computational cost behind sparse FPCA
FPCAdense <- FPCA(training$x, training$t,
                  list(plot = T, dataType ='Dense', nRegGrid = 501, useBinnedData = 'OFF', error = F))

# Get the central S matrix
temp <- t(FPCAdense$phi) %*% FPCAdense$phi
S <- temp %*% diag(1/FPCAdense$lambda) %*% temp

#### Predict the first K* FPC scores according to conditional expectation on the tuning ####
pred_zeta <- predict(FPCAdense, tuning$Ly, tuning$Lt, K=FPCAdense$selectK)

#### Predict the first K* FPC scores according to conditional expectation on the tuning ####
pred_eta <- predict(FPCAdense, testing$Ly, testing$Lt, K=FPCAdense$selectK)

#### Calculate all the pair-wise distances between all training samples ####
Dtrain <- matrix(0, nrow = length(training$Ly), ncol = length(training$Ly))
start_time = Sys.time()
Dtrain = foreach (i = 1:ncol(Dtrain), .combine=cbind) %:%
  foreach (j = 1:nrow(Dtrain), .combine=rbind) %dopar%{
    temp = FPCAdense$xiEst[j,] - FPCAdense$xiEst[i,]
    Dtrain[j,i] = t(temp) %*% S %*% temp
  }
End_time = Sys.time()
End_time - start_time

Dtrain = Dtrain/((M-1)/10)^2


#### Calculate all the pair-wise distances between all training samples and all tuning samples####
Dtune <- matrix(nrow = length(training$Ly), ncol = length(tuning$Ly))
start_time = Sys.time()
Dtune = foreach (i = 1:ncol(Dtune), .combine=cbind) %:%
  foreach (j = 1:nrow(Dtune), .combine=rbind) %dopar%{
    temp = FPCAdense$xiEst[j,] - pred_zeta[i,]
    Dtune[j,i] = t(temp) %*% S %*% temp
}
End_time = Sys.time()
End_time - start_time

Dtune = Dtune/((M-1)/10)^2

#### Calculate all the pair-wise distances between all training samples and all testing samples####
Dtest <- matrix(nrow = length(training$Ly), ncol = length(testing$Ly))
start_time = Sys.time()
Dtest = foreach (i = 1:ncol(Dtest), .combine=cbind) %:%
  foreach (j = 1:nrow(Dtest), .combine=rbind) %dopar%{
    temp = FPCAdense$xiEst[j,] - pred_zeta[i,]
    Dtest[j,i] = t(temp) %*% S %*% temp
  }
End_time = Sys.time()
End_time - start_time

Dtest = Dtest/((M-1)/10)^2

# Create the Y_true thru a non-linear relationship
yTrue <- apply(exp(Ksi[1:N,]),1,sum)

#### Find the optimal total budget over the tuning sets ####

pi_cand = c(1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1, 5, 1e1, 5e1, 1e2, 5e2, 1e3)

S = length(pi_cand)

foreach (s = 1:nrow(Dtest), .combine=rbind) %dopar%
  


#### save the objects ####
#save(FPCAdense, file = '1x_dense_FPCA.Rdata')
load('1x_dense_FPCA.Rdata')

#save(pred_zeta, file = '1x_dense_pred_zeta.RData')
load('1x_dense_pred_zeta.RData')
#save(Dtbt, file = '1x_dense_Dtbt.RData')
load('1x_dense_Dtbt.RData')


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
temp <- (y_hat - matrix(rep(yTrue[901:1000], length(h)), nrow = 0.1*N))^2
RMSE <- sqrt(apply(temp,2,sum))

# find out the minimum of RMSE and do a finer serach of bandwidth

which(RMSE == min(na.omit(RMSE))) # 8

# define the bandwidth
h = seq(h[7], h[9], by = 0.00001)

# define the y_hat with different bandwith
y_hat = matrix(NA, ncol = length(h), nrow = 0.1*N)

y_hat <- foreach (i =1:length(h), .combine = 'cbind') %dopar% {
  # get the exp(-d^2/h) matrix
  temp = exp(-Dtbt/(2*h[i]^2))
  # get the predicted y of testing by NW
  y_hat = (t(temp) %*% yTrue[1:(.9*N)]) / apply(temp,2,sum)
}


# define the RMSE
temp <- (y_hat - matrix(rep(yTrue[901:1000], length(h)), nrow = 0.1*N))^2
RMSE <- sqrt(apply(temp,2,sum))

# find the minimum RMSE
which(RMSE == min(RMSE)) # 1186

# Get the optimal h

h[1186] # .08185

# Get the optimal y_hat
y_hat0 <- y_hat[,1186]

#
yTrue[901:1000]
y_hat0

min(RMSE) # 59.98

plot(y_hat0, col = 'red', main = 'p=1 dense no me simulation')
lines(yTrue[901:1000], col ='black')

plot(y_hat0, col = 'red', main = 'p=1 dense no me simulation', ylim = c(0, 100))
lines(yTrue[901:1000], col ='black')


## Calculate the coefficient of determinant 
SSTO = sum((yTrue[901:1000] - mean(yTrue[901:1000]))^2)
SSRes = sum((y_hat0 - yTrue[901:1000])^2)
R2 = 1 - SSRes/SSTO  # 0.997


rbenchmark::benchmark(sapply(seq_len(4), function(k) sum(a[,k]*a[,k])*(2*10/(M-1))),
                      replications = 10)


rbenchmark::benchmark(for (k in 1:4) {sum(a[,k]*a[,k])*(2*10/(M-1))},
                      replications = 10)
