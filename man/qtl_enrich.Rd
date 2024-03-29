% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qtl_enrich.R
\name{qtl_enrich}
\alias{qtl_enrich}
\title{Performs a QTL enrichment analysis based on a hypergeometric test for each QTL class}
\usage{
qtl_enrich(
  qtl_db,
  qtl_file,
  qtl_type = c("QTL_type", "Name"),
  enrich_type = c("genome", "chromosome"),
  chr.subset = NULL,
  nThreads = NULL,
  padj = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
  verbose = TRUE
)
}
\arguments{
\item{qtl_db}{The object obtained using the import_gff_gtf() function}

\item{qtl_file}{The output from find_genes_qtls_around_markers function}

\item{qtl_type}{A character indicating which type of enrichment will be performed. QTL_type indicates that the enrichment processes will be performed for the QTL classes, while Name indicates that the enrichment analysis will be performed for each trait individually}

\item{enrich_type}{A character indicating if the enrichment analysis will be performed for all the chromosomes ("genome") or for a subset of chromosomes ("chromosome). If the "genome" option is selected, the results reported are the merge of all chromosomes}

\item{chr.subset}{If enrich_type equal "chromosome", it is possible to define a subset of chromosomes to be analyzed. The default is equal NULL. Therefore, all the chromosomes will be analyzed}

\item{nThreads}{The number of threads to be used.}

\item{padj}{The algorithm for multiple testing correction to be adopted ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")}

\item{verbose}{Logical value defining if messages should of not be printed during the analysis (default=TRUE)}
}
\value{
A data frame with the p-value for the enrichment result
}
\description{
Takes the output from find_genes_qtls_around_markers and run a QTL enrichment analysis
}
\details{
The simple bias of investigation for some traits (such as milk production related traits in the QTL database for cattle) may result in a larger proportion of records in the database. Consequently, the simple investigation of the proportion of each QTL type might not be totally useful. In order to reduce the impact of this bias, a QTL enrichment analysis can be performed. The QTL enrichment analysis performed by GALLO package is based in a hypergeometric test using the number of annoatted QTLs within the candidate regions and the total number of the same QTL in the QTL database.
}
\examples{
\donttest{data(QTLmarkers)
data(gffQTLs)
out.qtls<-find_genes_qtls_around_markers(
db_file=gffQTLs,marker_file=QTLmarkers,
method = "qtl",marker = "snp",
interval = 500000, nThreads = 1)

out.enrich<-qtl_enrich(qtl_db=gffQTLs,
qtl_file=out.qtls, qtl_type = "Name",
enrich_type = "chromosome",chr.subset = NULL,
padj = "fdr",nThreads = 1)}
}
