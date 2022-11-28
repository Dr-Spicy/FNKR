#' Generate the predicted FPC scores via conditional expectation with parallel
#' 
#' @param FPCAobject The FPCAobject (sparse/dense) generated on the training set
#' @param dataset Either the tuning or the testing based on the situation


pred.FPCs <- function(FPCAobject, dataset){
   result = foreach(diter = 1:P, .packages = 'fdapace') %dopar% {
      zs = predict(FPCAobject[[diter]], dataset[[diter]]$x, dataset[[diter]]$t, K=FPCAobject[[diter]]$selectK)
      zs$scores
   }
   return(result)
}
