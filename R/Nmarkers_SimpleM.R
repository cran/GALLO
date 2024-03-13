#' Estimate the number of effective markers in a chromosome based on an adapted version of the simpleM methodology
#'
#' @param ld.file A data frame with the pairwise linkage disequilibrium (LD) values for a chromosome. The column names SNP_A, SNP_B, and R are mandatory, where the SNP_A and SNP_B contained the markers names and the R column the LD values between the two markers.
#' @param PCA_cutoff A cutoff for the total of the variance explained by the markers.
#' @details This function estimate the effective number of markers in a chromosome using adapted version of the simpleM methodology described in Gao et al. (2008). The function use as input a data frame composed by three mandatory columns (SNP_A, SNP_B, and R). This data frame can be obtained using PLINK or any other software to compute LD between markers. Additionally, a threshold for percentage of the sum of the variances explained by the markers must be provided. The number of effective markers identified by this approach can be used in multiple testing corrections, such as Bonferroni.
#' @return The effective number of markers identified by the SimpleM approach
#' @importFrom data.table setkeyv
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#'@references Gao et al. (2008) Genet Epidemiol, Volume 32, Issue 4, Pages 361-369.
#'(\doi{10.1002/gepi.20310})
#' @name Nmarkers_SimpleM
#' @export

Nmarkers_SimpleM<-function(ld.file, PCA_cutoff=0.995){
  
  snps<-unique(c(ld.file$SNP_A,ld.file$SNP_B))
  
  ld.file<-as.data.table(ld.file)
  
  setkeyv(ld.file, c("SNP_A", "SNP_B"))
  
  mat.r<-matrix(NA,ncol=length(snps), nrow=length(snps), dimnames = list(snps,snps))
  
  diag(mat.r)<-1
  
  for (i in snps) {
    tmp.ldz<-ld.file[which(ld.file$SNP_A==i),]
    mat.r[i, tmp.ldz$SNP_B] <- tmp.ldz$R
    
    tmp.ldz<-ld.file[which(ld.file$SNP_B==i),]
    mat.r[i, tmp.ldz$SNP_A] <- tmp.ldz$R
  }
  
  
  mat.r[lower.tri(mat.r)] <- mat.r[upper.tri(mat.r)]
  
  cut.off <- PCA_cutoff
  
  out.cut<-inferCut(mat.r, cut.off)
  
  return(out.cut)
}