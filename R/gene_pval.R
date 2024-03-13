#' Estimate a gene-level p-value using Weighted Z-score approach and Meta-analysis with LD correlation coefficients approach
#'
#' @param data A data frame with the results of the association test performed for each marker
#' @param db_file A data frame obtained from the import_gff_gtf containing the gtf information
#' @param marker_ld A data frame containing the pairwise linkage disequilibrium between markers in a chromosome
#' @param interval The interval (in base pairs) used to annotated markers downstream and upstream from the genes coordinates
#' @param p The name of the column containing the P-values for each marker
#' @details Requires a table with p-values from a association test, a gtf file file the gene coordinates in the same assembly used to map the variants used in the association study, and a data frame with pairwise linkage disequilibrium (LD) values between markers. This analysis must be performed for each chromosome individually. The data frame with the results of the association study must have three mandatory columns names as CHR, BP and SNP containing the chromosome, base pair position and marker name, respectively. The gtf file must be imported by the import_gff_gtf() function from GALLO or can be customized by the user, since it has the same columns names. The LD table must contain three mandatory columns, SNP_A, SNP_B and R. where, the first two columns must contain the marker names and the third column, the LD value between these markers. This dtaa frame can be obtained using PLINK or any other software which computes pairwise LD between markers in the same chromosome. In the absence of LD values between any two SNPs in the data frame, a LD equal zero is assumed 
#' @return A data frame with the gene level p-values obtained using the Weighted Z-score approach (P_WZ_ld) and Meta-analysis with LD correlation coefficients approach (P_meta_LD)
#' @name gene_pval
#' @export


gene_pval<-function(data,db_file,marker_ld,interval,p){
  
  chr.list<-unique(data$CHR)
  
  tmp.gtf<-db_file[which(db_file$chr==chr.list),]
    
  marker.in.gene<-find_markers_genes(db_file=tmp.gtf, marker_file =data, int=interval)
    

  gene.pval<-NULL
  for(g in unique(marker.in.gene$gene_id)){
    
    tmp.mark<-as.data.frame(marker.in.gene[which(marker.in.gene$gene_id==g),])
    
    uniq.mark<-unique(tmp.mark$SNP)
    
    gene.ld<-marker_ld[which(marker_ld$SNP_A%in%uniq.mark & marker_ld$SNP_B%in%uniq.mark),]
    
    matrix.ld<-matrix(0, ncol=length(uniq.mark), nrow=length(uniq.mark),dimnames = list(uniq.mark,uniq.mark))
    
    
    diag(matrix.ld)<-1
    
    for (i in uniq.mark) {
      
      
      if(nrow(gene.ld[which(gene.ld$SNP_A==i),])>0){
        tmp.ldz<-gene.ld[which(gene.ld$SNP_A==i),]
        tmp.ldz<-tmp.ldz[which(tmp.ldz$SNP_B%in%uniq.mark),]
        
        matrix.ld[i, tmp.ldz$SNP_B] <- tmp.ldz$R
      }
      
      if(nrow(gene.ld[which(gene.ld$SNP_B==i),])>0){
        tmp.ldz<-gene.ld[which(gene.ld$SNP_B==i),]
        tmp.ldz<-tmp.ldz[which(tmp.ldz$SNP_A%in%uniq.mark),]
        
        matrix.ld[i, tmp.ldz$SNP_A] <- tmp.ldz$R 
      }
    }
    
    matrix.ld[lower.tri(matrix.ld)] <- matrix.ld[upper.tri(matrix.ld)]
    
    p.val_WZ_ld<-WZ_ld(matrix.ld,tmp.mark[!duplicated(tmp.mark$SNP),p])
    
    p.val_meta_LD<-meta_LD(matrix.ld,tmp.mark[!duplicated(tmp.mark$SNP),p])
    
    tmp.out<-data.frame(Gene=g, P_WZ_ld=p.val_WZ_ld, P_meta_LD=p.val_meta_LD)
    
    gene.pval<-rbind(tmp.out,gene.pval)
    
  }
  
  
  marker.in.gene<-as.data.frame(marker.in.gene)
  
  gene.pval$gene_name<-marker.in.gene[match(gene.pval$Gene,marker.in.gene$gene_id),"gene_name"]
  
  gene.pval<-gene.pval[,c("Gene","gene_name","P_WZ_ld","P_meta_LD")]
  
  return(gene.pval)

}