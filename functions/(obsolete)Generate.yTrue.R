#' Generate the true response yTrue based on the model index
#' 
#' @param model.index 
#' @param Ksi_all The list contains FPC scores, N by number of eigencomponents matrix, for the i-th functional predictor, which is the i-th element. 
#' @param y a vector recording the scalar response values, with the added model error. 
#'
#' @export yTrue a n by 1 vector  recording the scalar response values from X and Z by the true model, free of the regression model error.



Generate.yTrue <- function(model.index, Ksi_all){
   
   if (model.index == 1) {
     
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
      
      
      
   } else if (model.index == 2) {
      
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
   }
   
   return(yTrue)
}