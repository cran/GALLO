#' Search genes and QTLs around candidate regions
#'
#' Takes a list of candidate markers and or regions (haplotypes, CNVs, windows, etc.) and search for genes or QTLs in a determined interval
#' @param db_file The data frame obtained using the import_gff_gtf() function
#' @param marker_file The file with the SNP or haplotype positions. Detail: For SNP files, the columns “CHR” and “BP” with the chromosome and base pair position, respectively, are mandatory. For the haplotype, the following columns are mandatory: “CHR”, “BP1” and “BP2”
#' @param method “gene” or “qtl”
#' @param marker "snp" or "haplotype"
#' @param interval The interval in base pair which can be included upstream and downstream from the markers or haplotype coordinates.
#' @param nThreads Number of threads to be used
#' @param  verbose Logical value defining if messages should of not be printed during the analysis (default=TRUE)
#' @return A dataframe with the genes or QTLs mapped within the specified intervals
#' @name find_genes_qtls_around_markers
#' @importFrom utils read.delim
#' @importFrom parallel stopCluster
#' @import DT
#' @import webshot
#' @examples
#' data(QTLmarkers)
#' data(gffQTLs)
#' out.qtls<-find_genes_qtls_around_markers(db_file=gffQTLs, marker_file=QTLmarkers,
#' method = "qtl", marker = "snp",
#' interval = 500000, nThreads = 1)
#' @export

find_genes_qtls_around_markers<-function(db_file,marker_file,method=c("gene","qtl"),marker=c("snp","haplotype"),interval=0,nThreads=NULL, verbose=TRUE){
    interval=interval
    method <- match.arg(method)
    marker <- match.arg(marker)
    if(verbose==TRUE){
        cat(paste("You are using the method:", method, "with", marker))
        cat("\n")
    }

    if(marker=="snp"){
    #Creating gene data frame
        if (method=="gene"){
            chr_list<-unique(marker_file$CHR)
            output.final<-sub_genes_markers(chr_list,db_file,marker_file,nThreads=nThreads,int=interval)
        }else {
            if(verbose==TRUE){
                cat(paste("Starting QTL searching using ", interval, " bp", " as interval", sep=""))
                cat("\n")
            }
            chr_list<-unique(marker_file$CHR)
            output.final<-sub_qtl_markers(chr_list,db_file,marker_file,nThreads=nThreads,int=interval)
            if(verbose==TRUE){
                cat("Preparing output file for QTL annotation")
                cat("\n")
            }
            output.final<-splitQTL_comment(output.final=output.final)
        }
        }else{

    if (method=="gene"){
        chr_list<-unique(marker_file$CHR)
        output.final<-sub_genes_windows(chr_list,db_file,marker_file,nThreads=nThreads,int=interval)
    }else{
        if(verbose==TRUE){
            cat(paste("Starting QTL searching using ", interval, " bp", " as interval", sep=""))
            cat("\n")
        }
        chr_list<-unique(marker_file$CHR)
        output.final<-sub_qtl_windows(chr_list,db_file,marker_file,nThreads=nThreads,int=interval)
        if(verbose==TRUE){
            cat("Preparing output file for QTL annotation")
        }
        output.final<-splitQTL_comment(output.final=output.final)
    }
        }
  return(as.data.frame(output.final))
}
