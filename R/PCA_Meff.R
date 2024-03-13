#' Compute Meff statistic based on PCA to determine the number of effective markers
#'
#' @param eigenV The eigenvalues obtained from the linkage disequilibrium matrix
#' @param cut.off The threshold for percentage of the sum of the variances explained by the markers
#' @return The effective number of markers identified by the SimpleM approach
#' @name PCA_Meff
#' @keywords internal

PCA_Meff <- function(eigenV, cut.off){
  totalEigenV <- sum(eigenV)
  cuttmp <- cut.off*totalEigenV
  num_Eigens <- length(eigenV)
  EigenSum <- 0
  ind_Eigen <- 0
  
  for(i in 1:num_Eigens){
    if(EigenSum <= cuttmp){
      EigenSum <- EigenSum + eigenV[i]
      ind_Eigen <- i
    }
    else{
      break
    }
  }	
  return(ind_Eigen)
}
