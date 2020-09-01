## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=T, results='asis'---------------------------------------------------
#Installation
#devtools::install_github("pablobio/GALLO")

#Loading the package
library(GALLO)

## ----echo=T, results='asis'---------------------------------------------------
#Importing QTL markers from example dataset
data("QTLmarkers")

DT::datatable(QTLmarkers, rownames = FALSE, extensions = 'FixedColumns',
  options = list(scrollX = TRUE))

dim(QTLmarkers)

#Importing QTL windows from example dataset
data("QTLwindows")
DT::datatable(QTLwindows, rownames = FALSE, extensions = 'FixedColumns',
  options = list(scrollX = TRUE))

dim(QTLwindows)

## ----echo=T, results='asis'---------------------------------------------------
#Importing QTL annotation database
data("gffQTLs")

#Printing the first 100 rows
DT::datatable(gffQTLs[1:100,], rownames = FALSE, extensions = 'FixedColumns',
  options = list(scrollX = TRUE))

dim(gffQTLs)

#Importing gene annotation database
data("gtfGenes")

#Printing the first 100 rows
DT::datatable(gtfGenes[1:100,], rownames = FALSE, extensions = 'FixedColumns',
  options = list(scrollX = TRUE))

dim(gtfGenes)


## ----echo=T, results='asis'---------------------------------------------------
#An example of how to import a QTL annotation file
#qtl.inp <- import_gff_gtf(db_file="QTL_db.gff",file_type="gff")

#An example of how to import a gene annotation file
#qtf.inp <- import_gff_gtf(db_file="Gene_db.gtf",file_type="gtf")

## ----echo=T, results='asis'---------------------------------------------------

#Running gene annotation
out.genes<-find_genes_qtls_around_markers(db_file=gtfGenes,
marker_file=QTLmarkers, method = "gene", 
marker = "snp", interval = 500000, nThreads = NULL)

#Checking the first rows from the output file
DT::datatable(out.genes, rownames = FALSE, 
  extensions = 'FixedColumns',
  options = list(scrollX = TRUE))

#Checking the dimensions of the output file
dim(out.genes)

## ----echo=T, results='asis'---------------------------------------------------
#Running QTL annotation
out.qtls<-find_genes_qtls_around_markers(db_file=gffQTLs,
marker_file=QTLmarkers, method = "qtl",
marker = "snp", interval = 500000, nThreads = NULL)

#Checking the first rows from the output file
DT::datatable(out.qtls, rownames = FALSE, 
  extensions = 'FixedColumns',
  options = list(scrollX = TRUE))

#Checking the dimensions of the output file
dim(out.qtls)

## ----echo=T-------------------------------------------------------------------
#Removing the duplicated records for Gene ID
out.genes.unique<-out.genes[!duplicated(out.genes[c("Reference","gene_id")]),]
#Calculating the number of shared genes between different studies
overlap.genes<-overlapping_among_groups(file=out.genes.unique, x="Reference", y="gene_id")

overlap.genes

## ----echo=T-------------------------------------------------------------------
#Creating grouping labels (names of the references)
group.labels<-unique(out.genes.unique$Reference)

## ---- fig.align="center", fig.width=8, fig.height=8, fig.cap="Gene overlapping among studies."----
#plotting the overlapping information
plot_overlapping(overlap.genes, nmatrix=2, ntext=3, group=group.labels, labelcex = 1)


## ----echo=T, fig.height = 40, fig.width = 16, fig.align = "center"------------
#plotting the percentage of each QTL class annoatted
oldpar <- par(mar=c(1,30,1,1))
plot_qtl_info(out.qtls, qtl_plot = "qtl_type", cex=2)
par(oldpar)

## ---- fig.width=8, fig.height=6, fig.cap="Percentage of each trait annotated as a Reproduction QTL in the candidate intervals."----
#Setting margin parameter to better fit the axis labels
oldpar<-par(mar=c(5,20,1,1))
#plotting the percentage of each trait annoatted as a Reproduction QTL
plot_qtl_info(out.qtls, qtl_plot = "qtl_name", qtl_class="Reproduction")
par(oldpar)

## ---- results='hide', eval=F--------------------------------------------------
#  #Plotting percentage of the top 10 most frequent traits in all QTL classes
#  #(This is just an example code, the user do not need to execute
#  #this command for this tutorial)
#  QTL_classes<-unique(out.qtls$QTL_type)
#  
#  for(c in QTL_classes){
#    tmp.file.name<-paste(c,".png",sep="")
#    png(tmp.file.name,w=1500,h=900)
#    plot_qtl_info(out.qtls, qtl_plot = "qtl_name", qtl_class=c, n=10)
#    dev.off
#  }
#  

## ----echo=T, results='asis'---------------------------------------------------
#QTL enrichment analysis 
out.enrich<-qtl_enrich(qtl_db=gffQTLs, 
                       qtl_file=out.qtls, qtl_type = "Name",
                       enrich_type = "chromosome", chr.subset = NULL, 
                       padj = "fdr",nThreads = NULL)

#Checking the enriched QTLs
DT::datatable(out.enrich[order(out.enrich$pvalue),], 
              rownames = FALSE,
              extensions = 'FixedColumns',
              options = list(scrollX = TRUE))

## ---- echo=T------------------------------------------------------------------
#Creating a new ID composed by the trait and the chromosome
out.enrich$ID<-paste(out.enrich$QTL," - ","CHR",out.enrich$CHR,sep="")

#Match the QTL classes and filtering the Reproduction related QTLs
out.enrich.filtered<-out.enrich[which(out.enrich$adj.pval<0.05),]

## ---- fig.width=15, fig.height=10, fig.cap="Top 10 enriched traits identified in the QTL enrichment analysis. The area of the bubbles represents the number of observed QTLs for that class, while the color represents the p-value scale (the darker the color, smaller the p-value). Additionally, the x-axis shows the richness factor for each QTL, representing the ratio of number of QTLs and the expected number of that QTL."----
#Plotting the enrichment results for the QTL enrichment analysis
QTLenrich_plot(out.enrich.filtered, x="ID", pval="adj.pval")

## ---- fig.width=9, fig.height=7, fig.cap="Chord plot showing the relationship between the studies (right-hand side) and the enriched QTLs (abbreviations in the left -hand side)."----

#Filtering the output from QTL annotation

#Creating a new ID to filter the top 10 enriched QTLs
out.qtls$ID<-paste(out.qtls$Name," - ","CHR",out.qtls$CHR,sep="")
out.enrich.filtered<-out.enrich.filtered[order(out.enrich.filtered$adj.pval),]
out.qtls.filtered<-out.qtls[which(out.qtls$ID%in%out.enrich.filtered$ID[1:10]),]

#Creating color scheme based on the References
out.qtls.filtered[which(out.qtls.filtered$Reference=="Feugang et al. (2010)"),
                  "color_ref"]<-"purple"

out.qtls.filtered[which(out.qtls.filtered$Reference== "Buzanskas et al. (2017)"),
                  "color_ref"]<-"pink"


#Creating a color vector filled with black for all the traits abbreviation 
#and with the respective colors for each reference
color.grid<-c(rep("black",length(unique(out.qtls.filtered$Abbrev))),
              unique(out.qtls.filtered$color_ref))

#Naming the vector
names(color.grid)<-c(unique(out.qtls.filtered$Abbrev),
                     unique(out.qtls.filtered$Reference))

#Plotting the relationship plot using the grid color created above
relationship_plot(qtl_file=out.qtls.filtered, x="Abbrev",
                  y="Reference",cex=1,gap=2.5,degree = 60,
                  canvas.xlim = c(-5, 5), 
                  canvas.ylim = c(-3, 3), grid.col = color.grid)

