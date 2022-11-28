#### Define g function fors functional X to y ####
#' Three functions to link functional X to scalar y
#' g1. is linear
#' g2. is quadratic
#' g3 and g4 is 2-way interaction between 2 functional X

### new g functions to even the contribution of each eigen compoentns ###


#' Generate the true response yTrue linearly
#' 
#' @param x x is the 'Ksi_all[[diter]]'
#' @param var.vec The vector of the variance of FPC scores, defalt is FPCv. 
#' @param eigen.index.cntr the indices of contributing eigencompoennts.
#'
#' @export g.result The contribution to yTrue

glinear.eq <- function(x, var.vec = FPCv, eigen.index.cntr = linear.index.of.contirbuting.eigen){
   
   if (length(eigen.index.cntr) == 1) {
      g.result = x[,eigen.index.cntr] * as.vector(1/sqrt(var.vec[eigen.index.cntr]))
   }
   else {
      g.result = x[,eigen.index.cntr] %*% as.vector(1/sqrt(var.vec[eigen.index.cntr]))
      # each eigen comp add 1 to the total variance, factor it out
      g.result = g.result/sqrt(length(eigen.index.cntr))
   }
   return(g.result*sqrt(1))
}

#' Generate the true response yTrue quadratically
#' 
#' @param x x is the 'Ksi_all[[diter]]'
#' @param var.vec The vector of the variance of FPC scores, defalt is FPCv. 
#' @param eigen.index.cntr the indices of contributing eigencompoennts.
#'
#' @export g.result The contribution to yTrue

gquadratic.eq <- function(x, var.vec = FPCv, eigen.index.cntr = quadratic.index.of.contirbuting.eigen){
   
   if (length(eigen.index.cntr) == 1) {
      g.result = x[,eigen.index.cntr]**2 * (1/var.vec[eigen.index.cntr])
      # each eigen comp add 1 to the total variance, factor it out
      g.result = g.result/sqrt(2*length(eigen.index.cntr))
   }
   else {
      g.result = x[,eigen.index.cntr]**2 %*% (1/var.vec[eigen.index.cntr])
      # each eigen comp add 1 to the total variance, factor it out
      g.result = g.result/sqrt(2*length(eigen.index.cntr))
   }
   return(g.result)
}

#' Generate the true response yTrue by a 2-way interaction
#' 
#' @param diter1 is the 1th index of 'Ksi_all[[diter]]' in the interaction
#' @param diter2 is the 2th index of 'Ksi_all[[diter]]' in the interaction
#' @param var.vec The vector of the variance of FPC scores, defalt is FPCv. 
#' @param eigen.index.cntr the indices of contributing eigencompoennts.
#'
#' @export g.result The contribution to yTrue

interaction.eq <- function(diter1, diter2, var.vec = FPCv, eigen.index.cntr = interaction.index.of.contributing.eigen){
   x1 = Ksi_all[[diter1]]; x2  = Ksi_all[[diter2]];
   g.result = x1[,eigen.index.cntr[1]] * x2[,eigen.index.cntr[2]]/sqrt(var.vec[eigen.index.cntr[1]] * var.vec[eigen.index.cntr[2]])
   return(g.result)
}

#### Define g function fors functional X to y ####
#' Three functions to link functional X to scalar y
#' g1. is linear
#' g2. is quadratic
#' g3 and g4 is 2-way interaction between 2 functional X

### new g functions to even the contribution of each eigen compoentns ###


#' Generate the true response yTrue linearly
#' 
#' @param x x is the 'Ksi_all[[diter]]'
#' @param var.vec The vector of the variance of FPC scores, defalt is FPCv. 
#' @param eigen.index.cntr the indices of contributing eigencompoennts.
#'
#' @export g.result The contribution to yTrue

glinear.eq.mean <- function(x, var.vec = FPCv, eigen.index.cntr = linear.index.of.contirbuting.eigen){
   
   if (length(eigen.index.cntr) == 1) {
      g.result = mean(x[,eigen.index.cntr]) * as.vector(1/sqrt(var.vec[eigen.index.cntr]))%*% t(rep(1,N))
   }
   else {
      g.result = apply(x[,eigen.index.cntr], 2, mean) %*% as.vector(1/sqrt(var.vec[eigen.index.cntr])) %*% t(rep(1,N))
      # each eigen comp add 1 to the total variance, factor it out
      g.result = g.result/sqrt(length(eigen.index.cntr))
   }
   return(g.result*sqrt(1))
}

#' Generate the true response yTrue quadratically
#' 
#' @param x x is the 'Ksi_all[[diter]]'
#' @param var.vec The vector of the variance of FPC scores, defalt is FPCv. 
#' @param eigen.index.cntr the indices of contributing eigencompoennts.
#'
#' @export g.result The contribution to yTrue

gquadratic.eq.mean <- function(x, var.vec = FPCv, eigen.index.cntr = quadratic.index.of.contirbuting.eigen){
   
   if (length(eigen.index.cntr) == 1) {
      g.result = mean(x[,eigen.index.cntr])**2 * (1/var.vec[eigen.index.cntr])%*% t(rep(1,N))
      # each eigen comp add 1 to the total variance, factor it out
      g.result = g.result/sqrt(2*length(eigen.index.cntr))
   }
   else {
      g.result = (apply(x[,eigen.index.cntr], 2, mean))**2 %*% (1/var.vec[eigen.index.cntr])%*% t(rep(1,N))
      # each eigen comp add 1 to the total variance, factor it out
      g.result = g.result/sqrt(2*length(eigen.index.cntr))
   }
   return(g.result)
}
