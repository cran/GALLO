#' Estimate the number of independent segments in a chromosome based on the effective population size
#'
#' @param chr.table A table containing the chromosomes and the chromosomal length (in centiMorgans).
#' @param chr_length The name of the column where the length of the chromosomes are informed.
#' @param Ne The effective population size.
#' @details This function uses a adapted version of the formula proposed by Goddard et al. (2011) to estimate the independent number of segments in a chromosome based on the effective population size.  
#' @return A data frame with the effevtive number of segments in each chromosome.
#'@references Goddard et al. (2011) Journal of animal breeding and genetics, Volume 128, Issue 6, Pages 409-421.
#'(\doi{10.1111/j.1439-0388.2011.00964.x})
#' @name Nseg_chr
#' @export


Nseg_chr<-function(chr.table,chr_length, Ne){
  
  for(i in 1:nrow(chr.table)){
    
    L_i<-chr.table[i,chr_length]
    chr.table[i,"Nseg"]<-(2*Ne*L_i)/log10(Ne*L_i)
    
  }
  
  return(chr.table)
}





