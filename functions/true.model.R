#### Model list ####

# set the indicies of the contributing eigencomponents of the functional predictors
linear.index.of.contirbuting.eigen = c(1,3)
quadratic.index.of.contirbuting.eigen = c(1,2)
interaction.index.of.contributing.eigen = c(1,2)

# set the relavent index
if (model.index == 1 || model.index == 2) {
   # indecies for true model
   relavent.indecies = as.integer(c(1,2,3,4,11)); irrelavent.indecies = setdiff(seq(1,P+Q+R), relavent.indecies);
} else if (model.index == 3 || model.index == 4) {
   # indecies for true model
   relavent.indecies = as.integer(c(1,2,5,11,12)); irrelavent.indecies = setdiff(seq(1,P+Q+R), relavent.indecies);
} else if (model.index == 5 || model.index == 6) {
  # indecies for true model
  relavent.indecies = as.integer(c(1,2,11,20)); irrelavent.indecies = setdiff(seq(1,P+Q+R), relavent.indecies);
} else if (model.index %in% c('Linear and Dense 1', 'Linear and Sparse 1') ){
  relavent.indecies = as.integer(c(1,2)); irrelavent.indecies = setdiff(seq(1,P+Q+R), relavent.indecies);
  d = 0; e = 0 ; f = 0;
  #### Define function for functional predictor ####
  beta1 <- function(t){
    return(10*sin(1*t-1)/36)
  }
  
  beta2 <- function(t){
    return(2*cos(2*t))
  }
  
  beta3 <- function(t){
    return(0)
  }
  
  beta4 <- function(t){
    return(d*sin(2*t))
  }
  
  beta5 <- function(t){
    return(e*sin(pi*t))
  }
  
  beta6 <- function(t){
    return(f*t^2)
  } 
} else if (model.index == 'FAM'){
  relavent.indecies = as.integer(c(1,2, 6, 7)); irrelavent.indecies = setdiff(seq(1,P+Q+R), relavent.indecies);
  #### Define function for functional predictor ####
  beta1 <- function(t){
    return(10*sin(1*t-1)/36)
  }
  
  beta2 <- function(t){
    return(2*cos(2*t))
  }
  
  beta3 <- function(t){
    return(0)
  }
  
  beta4 <- function(t){
    return(0)
  }
  
  beta5 <- function(t){
    return(0)
  }
  
}


#' Generate the true response yTrue based on the model index
#' 
#' @param model.index 
#' @param Ksi_all The list contains FPC scores, N by number of eigencomponents matrix, for the i-th functional predictor, which is the i-th element. 
#' @param y a vector recording the scalar response values, with the added model error. 
#'
#' @export yTrue a n by 1 vector  recording the scalar response values from X and Z by the true model, free of the regression model error.



Generate.yTrue <- function(model.index, Ksi_all){
  yTrue = rep(0,N); 
 
 if (model.index %in% c('Linear and Dense 1', 'Linear and Sparse 1')){
   
   ## 1st functional predictor
   diter = 1
   # apply the defined linear fucntion on ksi_11 and ksi_13
   temp = xTrue_all[[diter]] %*% betamatrix[diter,] * (max(sc)- min(sc))/(M-1)
   # update the yTrue
   yTrue = yTrue + temp
   
   # ## 2th predictor
   diter = 2
   # apply the defined linear fucntion on ksi_11 and ksi_13
   temp = xTrue_all[[diter]] %*% betamatrix[diter,] * (max(sc)- min(sc))/(M-1)
   # update the yTrue
   yTrue = yTrue + temp
   
 } else if (model.index == 'FAM') {
   
   ## 1st functional predictor
   diter = 1
   # apply the defined linear fucntion on ksi_11 and ksi_13
   temp = weights_aa[diter] * xTrue_all[[diter]] %*% betamatrix[diter,] * (max(sc)- min(sc))/(M-1)
   # update the yTrue
   yTrue = yTrue + temp
   
   # ## 2th predictor
   diter = 2
   # apply the defined linear fucntion on ksi_11 and ksi_13
   temp = weights_aa[diter] * xTrue_all[[diter]] %*% betamatrix[diter,] * (max(sc)- min(sc))/(M-1)
   # update the yTrue
   yTrue = yTrue + temp
   
   ## 1st scalar predictor
   temp = weights_aa[3] * sqrt(2) * Z1
   # update the yTrue
   yTrue = yTrue + temp 
   
   ## 2nd scalar predictor
   temp = weights_aa[4]  * Z2^2
   # update the yTrue
   yTrue = yTrue + temp 
   
 } else if (model.index == 1) {
    
    ## 1st functional predictor
    diter = 1
    # apply the defined linear fucntion on ksi_11 and ksi_13
    temp = glinear.eq(Ksi_all[[diter]])
    # update the yTrue
    yTrue = yTrue + t(temp)
    
    # ## 2th predictor
    diter = 2
    # apply the defined quadratic fucntion on ksi_11 and ksi_13
    temp = gquadratic.eq(Ksi_all[[diter]])
    # update the yTrue
    yTrue = yTrue + t(temp)
    
    ## add the functional interaction b/t ksi_31 and ksi_42
    temp = interaction.eq(3,4)
    # update the yTrue
    yTrue = yTrue + sqrt(2)*t(temp)
    
    # add the quadratic Z_1 
    temp = Zmatrix[,1]^2/sqrt(2)
    # update the yTrue
    yTrue = yTrue + 1*t(temp)
    
    
    
  } else if (model.index == 3) {
    
    ## 1st functional predictor
    diter = 1
    # apply the defined linear fucntion on ksi_11 and ksi_13
    temp = glinear.eq(Ksi_all[[diter]])
    # update the yTrue
    yTrue = yTrue + t(temp)
    
    # ## 2th predictor
    diter = 2
    # apply the defined quadratic fucntion on ksi_11 and ksi_13
    temp = gquadratic.eq(Ksi_all[[diter]])
    # update the yTrue
    yTrue = yTrue + t(temp)
    
    # add the interaction b/t ksi_52 and Z_2
    diterp = c(5,2); xz.interaction.index.of.contributing.eigen = 2;
    temp = Ksi_all[[diterp[1]]][,xz.interaction.index.of.contributing.eigen] * Zmatrix[,diterp[2]] / sqrt(FPCv[xz.interaction.index.of.contributing.eigen])
    # update the yTrue
    yTrue = yTrue + sqrt(2)*t(temp)
    
    # add the Z_1 
    temp = Zmatrix[,1]^2/sqrt(2)
    # update the yTrue
    yTrue = yTrue + 1*t(temp)
    
  } else if (model.index == 5) {
    
    ## 1st functional predictor
    diter = 1
    # apply the defined linear fucntion on ksi_11 and ksi_13
    temp = glinear.eq(Ksi_all[[diter]])
    # update the yTrue
    yTrue = yTrue + t(temp)
    
    # ## 2th predictor
    diter = 2
    # apply the defined quadratic fucntion on ksi_11 and ksi_13
    temp = gquadratic.eq(Ksi_all[[diter]])
    # update the yTrue
    yTrue = yTrue + t(temp)
    
    # add the arctan Z1 and Zd1 and Zd2
    temp = atan(2*(Zmatrix[,1]-0)/(-1 * (Zdmatrix[,5] == 0) + (Zdmatrix[,5] == 1)+ 1 * (Zdmatrix[,5] == 2) + 4 * (Zdmatrix[,5] == 3)))
    # update the yTrue
    yTrue = yTrue + sqrt(5)*t(temp)
    
  } else if (model.index == 2) {
    ## 1st functional predictor
    diter = 1
    # apply the defined linear fucntion on ksi_11 and ksi_13
    temp = glinear.eq(Ksi_all[[diter]])
    # update the yTrue
    yTrue = yTrue + sqrt(2)*t(temp)
    
    # ## 2th predictor
    diter = 2
    # apply the defined quadratic fucntion on ksi_11 and ksi_13
    temp = gquadratic.eq(Ksi_all[[diter]])
    # update the yTrue
    yTrue = yTrue + t(temp)
    
    ## add the functional interaction b/t ksi_31 and ksi_42
    temp = interaction.eq(3,4)
    # update the yTrue
    yTrue = yTrue + sqrt(2)*t(temp)
    
    # add the quadratic Z_1 
    temp = Zmatrix[,1]^2/sqrt(2)
    # update the yTrue
    yTrue = yTrue + 1.1*t(temp)
  } else if (model.index == 4) {
    
    ## 1st functional predictor
    diter = 1
    # apply the defined linear fucntion on ksi_11 and ksi_13
    temp = glinear.eq(Ksi_all[[diter]])
    # update the yTrue
    yTrue = yTrue + sqrt(2)*t(temp)
    
    # ## 2th predictor
    diter = 2
    # apply the defined quadratic fucntion on ksi_11 and ksi_13
    temp = gquadratic.eq(Ksi_all[[diter]])
    # update the yTrue
    yTrue = yTrue + t(temp)
    
    # add the interaction b/t ksi_52 and Z_2
    diterp = c(5,2); xz.interaction.index.of.contributing.eigen = 2;
    temp = (Ksi_all[[diterp[1]]][,xz.interaction.index.of.contributing.eigen]-2) * Zmatrix[,diterp[2]] / sqrt(FPCv[xz.interaction.index.of.contributing.eigen])
    # update the yTrue
    yTrue = yTrue + sqrt(2)*t(temp)
    
    # add the Z_1 
    temp = Zmatrix[,1]^2/sqrt(2)
    # update the yTrue
    yTrue = yTrue + 1.1*t(temp)
  } else if (model.index == 6) {
    
    ## 1st functional predictor
    diter = 1
    # apply the defined linear fucntion on ksi_11 and ksi_13
    temp = glinear.eq(Ksi_all[[diter]])
    # update the yTrue
    yTrue = yTrue + t(temp)*sqrt(2)
    
    # ## 2th predictor
    diter = 2
    # apply the defined quadratic fucntion on ksi_11 and ksi_13
    temp = gquadratic.eq(Ksi_all[[diter]])
    # update the yTrue
    yTrue = yTrue + t(temp)*sqrt(2)
    
    # add the arctan Z1 and Zd1 and Zd2
    temp = atan(1*(Zmatrix[,1]-1.35)/(-2 * (Zdmatrix[,5] == 0) + (Zdmatrix[,5] == 1)+ 2 * (Zdmatrix[,5] == 2) + 0.5 * (Zdmatrix[,5] == 3)))
    # update the yTrue
    yTrue = yTrue + sqrt(6)*t(temp)
    
  }
  
  return(yTrue)
}
