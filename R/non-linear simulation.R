library(data.table)
library(fdapace)

#### generating a toy dense functional dataset from scratch####

# Set the number of subjects (N) and the
# number of predictors (P)
# number of measurements per predictor(M) 
N <- 100;  # when N=500, PACE time = 657s
P <- 3;
M <- 201;
set.seed(127)

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
lambda_sq = c(sqrt(2),1)
col1 <- rnorm(N, sd = lambda_sq[1])
col2 <- rnorm(N, sd = lambda_sq[2])
Ksi <- as.matrix(cbind(col1, col2))
Ksi <- apply(Ksi, 2, scale)
#Ksi <- Ksi %*% diag(lambda_sq)

# sample number of measurement per subjects
num_m = sample(x=matrix(rep(c(M),N),nrow= N), size = N)

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
  xNoisy[[i]] = xNoisy[[i]] + rnorm(xTrue[[i]],sd = 0)
}

# Create the Y_true thru a non-linear relationship
yTrue <- apply(Ksi[1:N,]^2,1,sum)

### separate the simulated data into training and testing 
training = list(x= xNoisy[c(1:(N*.9))], t= sampled_t[c(1:(N*.9))], y = yTrue[1:(N*.9)])
testing = list(x= xNoisy[c((N*.9+1):N)], t= sampled_t[c((N*.9+1):N)], y = yTrue[(N*.9+1):N])


####Running FPCA on training set ####

# Do FPCA on this sparse sample with 95% FVE cut 
# Notice that sparse FPCA will smooth the data internally (Yao et al., 2005)
# Smoothing is the main computational cost behind sparse FPCA
FPCAsparse <- FPCA(training$x, training$t, list(plot = F, FVEthreshold = 0.95,
                                           nRegGrid =201)) # nRegGrid only affects time for Cov quadratically

#### Predict the first K* FPC scores according to conditional expectation on the testing ####
pred_zeta <- predict(FPCAsparse, testing$x, testing$t, K=FPCAsparse$selectK)
#pred_zeta <- apply(pred_zeta, 2, scale)

# Get the central S matrix
temp <- t(FPCAsparse$phi) %*% FPCAsparse$phi
S <- temp %*% diag(1/FPCAsparse$lambda) %*% temp

# Calculate all the pair-wise distances between all training and all testing

Dtbt <- matrix(nrow = length(training$x), ncol = length(testing$x))

for (i in 1:ncol(Dtbt)){
  for (j in 1:nrow(Dtbt)){
    temp = FPCAsparse$xiEst[j,] - pred_zeta[i,]
    Dtbt[j,i] = t(temp) %*% S %*% temp
  }
}

Dtbt = Dtbt/20^2
#### Functional NW estimate ####

# define the bandwidth
h = .3

# get the exp(-d^2/h) matrix
temp = exp(-Dtbt/(2*h^2))

# get the predicted y of testing by NW 
y_hat = (t(temp) %*% training$y) / apply(temp,2,sum)
  
testing$y



##
Dtbt_a = Dtbt/max(Dtbt)
# get the exp(-d^2/h) matrix
temp = exp(-Dtbt_a/(2*h))

# get the predicted y of testing by NW 
y_hat = (t(temp) %*% training$y) / apply(temp,2,sum)






















      
# Calculate the distance matrix D among all training data
D <- diag(rep(0,0.9*N))
for (i in 2:(0.9*N)){
  j = 1
  while (j<i) {
    temp = FPCAsparse$xiEst[i,] - FPCAsparse$xiEst[j,]
    D[i,j] = t(temp) %*% S %*% temp
    j=j+1
  }
}

#### Now do NW estimates on the training data









#### Further functionality ####
# FPCA calculates the bandwidth utilized by each smoother using generalized CV or k-fold CV automatically.
# Dense data are not smoothed by default. The argument methodMuCovEst can be switched between smooth and 
# cross-sectional if one wants to utilize different estimation techniques when work with dense data. 

# The bandwidth used for estimating the smoothed mean and the smoothed covariance are available under
# ...bwMu and bwCov respectively. Users can nevertheless provide their own bandwidth estimates:
FPCAsparseMuBW5 <- FPCA(ySparse$yNoisy, ySparse$Lt, optns = list(userBwMu = 5))
# Visualising the fitted trajectories is a good way to see if the new bandwidth made any sense:
par(mfrow=c(1,2))
CreatePathPlot( FPCAsparse, subset = 1:5, main = "GCV bandwidth", pch = 16)
CreatePathPlot( FPCAsparseMuBW5, subset = 1:3, main = "User-defined bandwidth", pch = 16)
par(mfrow=c(1,1))
# Use rectangular kernel
FPCAsparseRect <- FPCA(ySparse$yNoisy, ySparse$Lt, optns = list(kernel = 'rect'))

# FPCA returns automatically the smallest number of components required to explain 99.99% of a sample’s variance.
# Using the function selectK one can determine the number of relevant components according to AIC, 
# BIC or a different Fraction-of-Variance-Explained threshold. For example
SelectK( FPCAsparse, criterion = 'FVE', FVEthreshold = 0.95) # K = 3
SelectK( FPCAsparse, criterion = 'AIC') # K = 2
SelectK( FPCAsparse, criterion = 'BIC') # K = 2

# When working with functional data (usually not very sparse) the estimation of derivatives is often of interest.
# Using fitted.FPCA one can directly obtain numerical derivatives by defining the appropriate order p;
# fdapace provides for the first two derivatives ( p =1 or 2).
# Because the numerically differentiated data are smoothed the user can define smoothing specific arguments 
#  (see ?fitted.FPCA for more information); the derivation is done by using the derivative of the linear fit.
# Similarly using the function FPCAder , one can augment an FPCA object with functional derivatives of a sample’s mean function and eigenfunctions.

# equivalent: fitted(FPCAsparse, derOptns=list(p = 0));
fittedCurvesP0 <- fitted(FPCAsparse)
# Get first order derivatives of fitted curves, smooth using Epanechnikov kernel
fittedCurcesP1 <- fitted(FPCAsparse, K=3, derOptns=list(p = 1, kernelType = 'epan'))




































# Set the # of predictor trajectories
N <- 4
# Set the grid size for X(t) and t 
Nxt <- 40;
Nt <- 40;
set.seed(127)

# Define the continuum
time <- seq(0,10,length.out = Nt)

# Create FPC scores
U1 <- runif(N, 0, 2*pi)
U2 <- runif(N, 0, 2*pi)
Ksi1 <- cos(U1)
Ksi2 <- sin(U1)
Ksi3 <- cos(U2)
Ksi4 <- sin(U2)
Ksi <- as.matrix(cbind(Ksi1, Ksi2, Ksi3, Ksi4));
Ksi <- apply(Ksi, 2, scale)


# Define the eigenfunctions for X(t)
T = 10
eigFunct1 <- function(t) sin(2*t*pi/T) 
eigFunct2 <- function(t) cos(2*t*pi/T) 
eigFunct3 <- function(t) sin(4*t*pi/T)
eigFunct4 <- function(t) cos(4*t*pi/T) 

# Get the X(t)
X_t = Ksi %*% t(matrix(c(eigFunct1(time),eigFunct2(time),eigFunct3(time),eigFunct4(time)), ncol=4))


#### running FPCA on a dense dataset with no noise ####
L3 <- MakeFPCAInputs(IDs = rep(1:N, each=Nt), tVec=rep(time,N), yVec = t(X_t))
FPCAdense <- FPCA(L3$Ly, L3$Lt, 
                  list(plot = F, dataType ='Dense'))

SelectK(FPCAdense, criterion = 'AIC') # K = 4

# Now, we just wanna set K=2 and do the sparse FPCA again
FPCAdense_K2  <- FPCA(L3$Ly, L3$Lt, 
                        list(plot = TRUE, kernel = "gauss",
                             maxK = 2, nRegGrid=Nt)) 

# Extract the fitted Cov matrix
fitted_cov <- FPCAsparse0_K2$fittedCov

# Get the inverse Cov matrix by direct inverse
Inv_fitted_cov <- solve(fitted_cov) # Singluar!!!

# Extract the fitted eighen-functions and eigen-values
eigen_f <- FPCAdense_K2$phi
eigen_v <- FPCAdense_K2$lambda

# Get the inverse Cov matrix by our way
Inv_fitted_cov_right <- eigen_f %*% diag(eigen_v) %*% t(eigen_f)

# define a new X_0 function
X0 <- function(t) sin(t)*cos(t)

# Define the distance

dist_0 <- function(input1, input2, grid = time){
  (input1-input2(grid)) %*% Inv_fitted_cov %*% t(input1-input2(grid))
}

# define the F1() to F4() functions, i.e. the non-linear relation b/t X(t) and Y

F1 <- function(input, grid = time) {
  op =  sum(cos(time - input-1))
  return(op)
}

F2 <- function(input, grid = time) {
  op =  sum(sin(time - input-1))
  return(op)
}

F3 <- function(input, grid = time) {
  op =  sum(cos(time - input-5))
  return(op)
}

F4 <- function(input, grid = time) {
  op =  sum(sin(time - input-5))
  return(op)
}

# Define true Y w/out noise

Y_0 <- F1(X_t[1,]) + abs(F2(X_t[2,])) + log(abs(F3(X_t[3,]))) + sqrt(abs(F4(X_t[4,])))

# Specify noise for y 
sigma = 1

# get the real Y

Y <- Y_0 + rnorm(1, sd=sigma)

#### running FPCA on a sparse dataset with no noise ####
# Create sparse sample  
# Each subject has one to six readings 

X_t_Sparse <- Sparsify(samp = X_t[,3:Nt-1], pts = t(time[3:Nt-1]), sparsity = c(1:6))

FPCAsparse0 <- FPCA(X_t_Sparse$Ly, X_t_Sparse$Lt,list(plot = TRUE,
                                                kernel = "gauss")) # K = 8 

SelectK(FPCAsparse0, criterion = 'AIC') # K = 2

# Now, we just wanna set K=2 and do the sparse FPCA again
FPCAsparse0_K2  <- FPCA(X_t_Sparse$Ly, X_t_Sparse$Lt, 
                        list(plot = TRUE, kernel = "gauss",
                             maxK = 2, nRegGrid=Nt)) 

# Extract the fitted Cov matrix
fitted_cov <- FPCAsparse0_K2$fittedCov

# Get the inverse Cov matrix by direct inverse
Inv_fitted_cov <- solve(fitted_cov) # Singluar!!!

# Extract the fitted eighen-functions and eigen-values
eigen_f <- FPCAsparse0_K2$phi
eigen_v <- FPCAsparse0_K2$lambda

# Get the inverse Cov matrix by our way
Inv_fitted_cov_right <- eigen_f %*% diag(eigen_v) %*% t(eigen_f)

# 