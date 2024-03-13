#' Function to annotate markers and the respective p-values within genes
#'
#' @param db_file A data frame obtained from the import_gff_gtf containing the gtf information
#' @param marker_file A data frame with the results of the association test performed for each marker
#' @param int The interval (in base pairs) used to annotated markers downstream and upstream from the genes coordinates
#' @return A data frame containing the markers mapped within the selected interval for each gene in the annotation file
#' @name find_markers_genes
#' @keywords internal

find_markers_genes<-function (db_file, marker_file, int = 0) 
{
  tmp_gene <- data.table::as.data.table(db_file)
  tmp_markers <- data.table::as.data.table(marker_file)
  tmp_markers$tmpBP1 <- tmp_markers$BP - int
  tmp_markers$tmpBP2 <- tmp_markers$BP + int
  data.table::setkey(tmp_markers, tmpBP1, tmpBP2)
  out <- data.table::foverlaps(tmp_gene, tmp_markers, 
                               by.x = c("start_pos", "end_pos"), by.y = data.table::key(tmp_markers), nomatch = 0)
  out[, -c("tmpBP1", "tmpBP2")]
}