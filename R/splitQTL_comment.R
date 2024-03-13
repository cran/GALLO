#' Sub-function to split comment column from QTL output
#'
#' Takes a list of candidate markers and search for genes a determined interval
#' @param output.final Output from QTL annotation
#' @return A data frame with the extra_info column content, from the gff file, broken in several additional columns
#' @importFrom stringr str_split_fixed
#' @keywords internal

splitQTL_comment<-function(output.final){
    
  output_qtls<-as.data.frame(output.final)
    output_qtls$QTL_type<-as.character(output_qtls$QTL_type)
    
    for (i in seq_len(nrow(output_qtls))) {
      split_content <- stringr::str_split(output_qtls$extra_info[i], ";")
      
      split_content<-unlist(split_content)
      
      split_content<-split_content[which(split_content!="")]
      
      list_qtl<-stringr::str_split(split_content, "=")
      
      max_length <- max(sapply(list_qtl, length))
      
      my_matrix <- matrix(nrow = length(list_qtl), ncol = max_length,
                          dimnames = list(NULL, c("First_Element", "Second_Element")))
      for (j in seq_along(list_qtl)) {
        sublist <- list_qtl[[j]]
        my_matrix[j, 1:length(sublist)] <- sublist
      }
      
      df <- as.data.frame(my_matrix, stringsAsFactors = FALSE)
      
      paste.val<-NULL
      if(anyDuplicated(df$First_Element) > 0){
        
        dup.id<-df$First_Element[duplicated(df$First_Element)]
        
        for(j in 1:length(dup.id)){
          
          tmp.val<-paste(df[which(df$First_Element==dup.id[j]),"Second_Element"],collapse=",")
          
          paste.val<-c(paste.val,tmp.val)
          
        }
        
        df<-df[!duplicated(df$First_Element),]
        
        df[which(df$First_Element%in%dup.id),"Second_Element"]<-paste.val
        
      }
      
      
      output_qtls[i,df$First_Element]<-df$Second_Element
      
    }
    
    
    output_qtls<-output_qtls[,-which(colnames(output_qtls)%in%"extra_info")]
    
    output_qtls$QTL_type<-gsub("_Association","",output_qtls$QTL_type)
    
    output_qtls$QTL_type<-gsub("_QTL","",output_qtls$QTL_type)
    
    return(output_qtls)
    
}
