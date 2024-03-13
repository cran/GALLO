#' unction to perform Meta-Analysis with LD Correlation Coefficients
#'
#' @param marker_ld A dataframe containing the pairwise linkage disequilibrium between markers in a chromosome
#' @param marker_pvalues A vector with the p-values for the SNPs annotated within each gene
#' @return A vector of p-values for each gene annotated within the defined coordinates
#' @name meta_LD
#' @keywords internal

meta_LD<-function(marker_ld, marker_pvalues){
 
  weights <- rowMeans(marker_ld) 
  
  weighted_average_pvalue <- sum(marker_pvalues * weights) / sum(weights)
  
  gene_pvalue_weighted_avg <- weighted_average_pvalue
  
  return(gene_pvalue_weighted_avg)
}