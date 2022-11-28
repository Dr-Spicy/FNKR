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
N <- 1000;  # when N=500, PACE time = 657s
P <- 3;
M <- 501;
set.seed(7)

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
lambda_sq = c(2,1)
col1 <- rnorm(N, sd = lambda_sq[1])
col2 <- rnorm(N, sd = lambda_sq[2])
Ksi <- as.matrix(cbind(col1, col2))
Ksi <- apply(Ksi, 2, scale)
Ksi <- Ksi %*% diag(lambda_sq)

# Create X_true based on sampled time as a list
x1True <- Ksi %*% t(matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2)) + t(matrix(rep(meanFunct(s),N), nrow=M))

#### running FPCA on a dense dataset with no noise ####
#L3 <- MakeFPCAInputs(IDs = rep(1:N, each=M), tVec=rep(s,N), yVec = t(xTrue))
#FPCAdense <- FPCA(L3$Ly, L3$Lt,
#                  list(plot = T, dataType ='Dense'))

# Create the Y_true thru a non-linear relationship
yTrue <- apply(exp(Ksi[1:N,]),1,sum)

### separate the simulated data into training and testing
training = MakeFPCAInputs(IDs = rep(1:(.9*N), each=M), tVec=rep(s,(.9*N)), yVec = t(x1True[1:(.9*N),]))
testing = MakeFPCAInputs(IDs = rep((.9*N+1):N, each=M), tVec=rep(s,(.1*N)), yVec = t(x1True[(.9*N+1):N,]))


####Running FPCA on training set ####

# Do FPCA on this sparse sample with 99% FVE cut
# Notice that sparse FPCA will smooth the data internally (Yao et al., 2005)
# Smoothing is the main computational cost behind sparse FPCA
FPCAdense <- FPCA(training$Ly, training$Lt,
                  list(plot = T, dataType ='Dense', nRegGrid = 1001))

#### Predict the first K* FPC scores according to conditional expectation on the testing ####
pred_zeta <- predict(FPCAdense, testing$Ly, testing$Lt, K=FPCAdense$selectK)
#pred_zeta <- apply(pred_zeta, 2, scale)

# Get the central S matrix
temp <- t(FPCAdense$phi) %*% FPCAdense$phi
S <- temp %*% diag(1/FPCAdense$lambda) %*% temp

# Calculate all the pair-wise distances between all training and all testing

Dtbt <- matrix(nrow = length(training$Ly), ncol = length(testing$Ly))
start_time = Sys.time()
Dtbt = foreach (i = 1:ncol(Dtbt), .combine=cbind) %:%
  foreach (j = 1:nrow(Dtbt), .combine=rbind) %dopar%{
    temp = FPCAdense$xiEst[j,] - pred_zeta[i,]
    Dtbt[j,i] = t(temp) %*% S %*% temp
  }
End_time = Sys.time()
End_time - start_time


Dtbt = Dtbt/((M-1)/10)^2
#### Functional NW estimate ####

# define the bandwidth
h = .5

# get the exp(-d^2/h) matrix
temp = exp(-Dtbt/(2*h^2))
# get the predicted y of testing by NW
y_hat = (t(temp) %*% yTrue[1:(.9*N)]) / apply(temp,2,sum)

yTrue[901:1000]




