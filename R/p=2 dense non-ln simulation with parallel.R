library(data.table)
library(fdapace)
library(doParallel)
library(foreach)


## to unregister the cl
registerDoSEQ()
cl <- makeCluster(8)
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
col1 <- rnorm(N, sd = lambda_sq[1])
col2 <- rnorm(N, sd = lambda_sq[2])
Ksi <- as.matrix(cbind(col1, col2))
Ksi <- apply(Ksi, 2, scale)
Ksi <- Ksi %*% diag(lambda_sq)

# Create X_true based on sampled time as a list
x1True <- Ksi %*% t(matrix(c(eigFunct1(s),eigFunct2(s)), ncol=2)) + t(matrix(rep(meanFunct(s),N), nrow=M))

# Create the Y_true thru a non-linear relationship
yTrue <- apply((Ksi)^2,1,sum)

#### Now, work on the 2nd predictor ####
# Define the number of eigencomponents
K_t <- 2
# Define the mean and 2 eigencomponents
meanFunct <- function(s) s + cos(s)
eigFunct1 <- function(s) -cos(4*s*pi/10) / sqrt(5)
eigFunct2 <- function(s) sin(4*s*pi/10) / sqrt(5)
#eigFunct1 <- function(s) -cos(s*pi/2) / sqrt(2)
#eigFunct2 <- function(s) sin(s*pi/2) / sqrt(2)

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

# sample ti's for each subject based on the # of measurement on [0,10] as a list
train_t = list()
test_t = list()
for (i in 1:(.9*N)){
  train_t[[i]] = s
}
for (i in 1:(.1*N)){
  test_t[[i]] = s
}


### separate the simulated data into training and testing
training = list(x1 = lapply(seq_len(N*.9), function(i) t(x1True)[,i]),
                x2 = lapply(seq_len(N*.9), function(i) t(x2True)[,i]),
                t = train_t, y = yTrue[1:(N*.9)])
testing = list(x1 = lapply(seq_len(N*.1), function(i) t(x1True)[,i]),
               x2 = lapply(seq_len(N*.1), function(i) t(x2True)[,i]),
               t= test_t, y = yTrue[(N*.9+1):N])

####Running FPCA on training set ####

# Do FPCA on this sparse sample with 99% FVE cut
# Notice that sparse FPCA will smooth the data internally (Yao et al., 2005)
# Smoothing is the main computational cost behind sparse FPCA

FPCAdense1 <- FPCA(training$x1, training$t,
                  list(plot = T, dataType ='Dense', nRegGrid = 1001, useBinnedData = 'OFF', error = F))

FPCAdense2 <- FPCA(training$x2, training$t,
                   list(plot = T, dataType ='Dense', nRegGrid = 1001, useBinnedData = 'OFF', error = F))


#### Predict the first K* FPC scores according to conditional expectation on the testing ####
pred_zeta_1 <- predict(FPCAdense1, testing$x1, testing$t, K=FPCAdense1$selectK)
pred_zeta_2 <- predict(FPCAdense2, testing$x2, testing$t, K=FPCAdense2$selectK)
#pred_zeta_1 <- apply(pred_zeta_1, 2, scale)
#pred_zeta_2 <- apply(pred_zeta_2, 2, scale)

# Get the central S matrix
temp1 <- t(FPCAdense1$phi) %*% FPCAdense1$phi
S1 <- temp1 %*% diag(1/FPCAdense1$lambda) %*% temp1
temp2 <- t(FPCAdense2$phi) %*% FPCAdense2$phi
S2 <- temp2 %*% diag(1/FPCAdense2$lambda) %*% temp2


# Calculate all the pair-wise distances between all training and all testing

Dtbt1 <- matrix(nrow = length(training$x1), ncol = length(testing$x1))
start_time = Sys.time()
Dtbt1 = foreach (i = 1:ncol(Dtbt1), .combine=cbind) %:%
  foreach (j = 1:nrow(Dtbt1), .combine=rbind) %dopar%{
    temp = FPCAdense1$xiEst[j,] - pred_zeta_1[i,]
    Dtbt1[j,i] = t(temp) %*% S1 %*% temp
  }
End_time = Sys.time()
End_time - start_time

Dtbt1 = Dtbt1/((M-1)/10)^2

Dtbt2 <- matrix(data = 0, nrow = length(training$x2), ncol = length(testing$x2))
start_time = Sys.time()
Dtbt2 = foreach (i = 1:ncol(Dtbt1), .combine=cbind) %:%
  foreach (j = 1:nrow(Dtbt1), .combine=rbind) %dopar%{
    temp = FPCAdense2$xiEst[j,] - pred_zeta_2[i,]
    Dtbt2[j,i] = t(temp) %*% S2 %*% temp
  }
End_time = Sys.time()
End_time - start_time


Dtbt2 = Dtbt2/((M-1)/10)^2

#### save the objects ####
#save(FPCAdense1, file = '2x1_dense_noerr_FPCA.Rdata')
load('2x1_dense_noerr_FPCA.Rdata')

#save(FPCAdense2, file = '2x2_dense_noerr_FPCA.Rdata')
load('2x2_dense_noerr_FPCA.Rdata')

#save(pred_zeta_1, file = '2x1_dense_noerr_pred_zeta.RData')
load('2x1_dense_noerr_pred_zeta.RData')

#save(pred_zeta_2, file = '2x2_dense_noerr_pred_zeta.RData')
load('2x2_dense_noerr_pred_zeta.RData')

#save(Dtbt1, file = '2x1_dense_noerr_Dtbt.RData')
load('2x1_dense_noerr_Dtbt.RData')

#save(Dtbt2, file = '2x2_dense_noerr_Dtbt.RData')
load('2x2_dense_noerr_Dtbt.RData')

##### Functional NW estimate ####

# define the bandwidth grid
h = seq(0,2,by = 1e-3) + 1e-3
h_grid <- matrix(nrow = length(h)^2, ncol = 2)
h_grid[,1] <- rep(h, length(h))
h_grid[,2] <- rep(h, each = length(h))

# define the y_hat with different bandwith
y_hat = matrix(NA, ncol = length(h), nrow = 0.1*N)

y_hat <- foreach (i =1:length(h), .combine = 'cbind') %dopar% {
  # get the exp(-d^2/h) matrix
  temp = exp(-Dtbt1/(2*h_grid[i,1]^2) - Dtbt2/(2*h_grid[i,2]^2))
  # get the predicted y of testing by NW
  y_hat = (t(temp) %*% yTrue[1:(.9*N)]) / apply(temp,2,sum)
}


# define the RMSE
temp <- (y_hat - matrix(rep(yTrue[901:1000], length(h)), nrow = 0.1*N))^2
RMSE <- sqrt(apply(temp,2,sum))

# find out the minimum of RMSE and do a finer serach of bandwidth

temp1 = which(RMSE == min(RMSE))
temp1 # 19436

h_grid[temp1,] # 4.17, 0.2

## zoom in the bandwidth (4.95~5.1), (0.0001~0.1)
h1 = seq(4.1,4.3, by =0.001)+0.001
h2 = seq(0.1,0.3, by =0.001)+0.001

h_grid <- matrix(nrow = length(h1)*length(h2), ncol = 2)
h_grid[,1] <- rep(h1, length(h2))
h_grid[,2] <- rep(h2, each = length(h1))

# define the y_hat with different bandwith
y_hat = matrix(NA, nrow = 0.1*N)

y_hat <- foreach (i =1:length(h), .combine = 'cbind') %dopar% {
  # get the exp(-d^2/h) matrix
  temp = exp(-Dtbt1/(2*h_grid[i,1]^2) - Dtbt2/(2*h_grid[i,2]^2))
  # get the predicted y of testing by NW
  y_hat = (t(temp) %*% yTrue[1:(.9*N)]) / apply(temp,2,sum)
}


# define the RMSE
temp <- (y_hat - matrix(rep(yTrue[901:1000], ncol(y_hat)), nrow = 0.1*N))^2
RMSE <- sqrt(apply(temp,2,sum))


# find the minimum RMSE
which(RMSE == min(RMSE)) # 19965

# Get the optimal h

h_grid[19965,] # 4.166 0.200

# Get the optimal y_hat
y_hat0 <- y_hat[,19965]
y_hat0
#
yTrue[901:1000]


###
