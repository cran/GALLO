#' Function to perform Weighted Z-score Approach with LD Information
#'
#' @param marker_ld A data frame containing the pairwise linkage disequilibrium between markers in a chromosome
#' @param marker_pvalues A vector with the p-values for the SNPs annotated within each gene
#' @return A vector of p-values for each gene annotated within the defined coordinates
#' @name WZ_ld
#' @keywords internal

WZ_ld<-function(marker_ld, marker_pvalues){

  z_scores <- qnorm(1 - marker_pvalues)
  
  weights <- rowMeans(marker_ld)  
  
  weighted_z_score <- sum(z_scores * weights) / sqrt(sum(weights^2))
  
  gene_pvalue_weighted_z <- 1 - pnorm(weighted_z_score)
  
  return(gene_pvalue_weighted_z)
}
