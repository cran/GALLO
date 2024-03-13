#' Compute a multi-trait test statistic for pleiotropic effects using summary statistics from association tests
#'
#' @param data A data frame with the first column containing the SNP name and the remaining columns the signed t-values obtained for each marker in the association studies individually performed for each trait.
#' @details This function tests a null hypothesis stating that each SNP does not affect any of the traits included in the input file. The method applied here is an implementation of the statistic proposed at Bolormaa et al. (2014) and is approximately distributed as a chi-squared with n degrees of freedom, where n is equal the number of traits included in the input file.   
#' @return A data frame with the multi-trait chi-squared statistics and the correspondent p-value obtained for each SNP.
#' @importFrom stats cor
#' @importFrom stats pchisq
#'@references Bolormaa et al. (2014) Plos Genetics, Volume 10, Issue 3, e1004198.
#'(\doi{10.1371/journal.pgen.1004198})
#' @name PleioChiTest
#' @export

PleioChiTest<-function(data){
  
  ns <- nrow(data)
  nt <- ncol(data)-1
  
  awm <- matrix(0, nrow = ns, ncol = nt)
  
  awm <- as.matrix(data[,-c(1)])
  
  V<-cor(awm)
  
  Vi <- solve(V)
  
  pleio <- NULL
  for (i in 1:ns) {
    rowawm <- awm[i, ]
    rowawmvi <- rowawm %*% Vi
    pleio[i] <- sum(rowawmvi * rowawm)
  }
  
  pleio.out<-data.frame(pleio_chisq=pleio)
  
  pleio.out$pleio_pval<-1-pchisq(pleio, df = nt, lower.tail = TRUE)
  
  pleio.out<-cbind(data[,1,drop=F],pleio.out)
  
  return(pleio.out)
}
