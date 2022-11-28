#' Get the central C matrices
#' @param method.for.Dist method to calc the C, D1 is not using the inversed estimated variance, whereas Dv uses it. 
#' @param FPCAobj the FPCA obj 
#' @param K.in.use the number of eigencomponents used in FPCA procedure
#' 
#' @export res The retured C mat in K.in.use by K.in.use

get.C.mat <- function(method.for.Dist, FPCAobj, K.in.use){
   res <- foreach(diter = 1:P) %do% {
      
      if (method.for.Dist == 'D1'){
         zs = diag(1/1, ncol = K.in.use, nrow = K.in.use)
      } else if (method.for.Dist == 'Dv') {
         zs = diag(1/FPCAsparse[[diter]]$lambda[1:K.in.use], ncol = K.in.use, nrow = K.in.use)
      }
      zs
   }
   return(res)
}



