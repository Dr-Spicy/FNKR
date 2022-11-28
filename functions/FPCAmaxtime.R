#' define a function to print out the maximal time FPCA used out of all P predictors

FPCAmaxtime <- function(){
   temp = rep(0, P)
   for (diter in 1:P){
      
      temp[diter] = FPCAsparse[[diter]]$timings[1]
   }
   FPCAtotaltime = max(temp)
   
   message('FPCA used', round(FPCAtotaltime/60,4), 'minutes.')
}
