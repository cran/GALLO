#' Compute the eigenvalues for the linkage disequilibrium (LD) matrix and use as input for PCA_Meff function to compute the effective number of markers
#'
#' @param mat.r A matrix composed by the LD between markers
#' @param cut.off The threshold for percentage of the sum of the variances explained by the markers
#' @return The effective number of markers identified by the SimpleM approach
#' @name inferCut
#' @keywords internal

inferCut <- function(mat.r, cut.off){

  eigen_res <- eigen(mat.r)

  eigenValues_abs <- abs(eigen_res$values)
  
  Meff_PCA <- PCA_Meff(eigenValues_abs, cut.off)
  
  return(Meff_PCA)
}