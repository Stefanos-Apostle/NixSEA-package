
####The following code will load in the package from github since it is not on CRAN yet
library(devtools)
install_github("Stefanos-Apostle/NixSEA-package")

library(NixSEA)
#######





###Below is the data and code to run the DESeq2 vigenette from the following;
###https://bioc.ism.ac.jp/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf
library("DESeq2")
library( "parathyroidSE" )

data("parathyroidGenesSE")
se <- parathyroidGenesSE
colnames(se) <- se$run

ddsFull <- DESeqDataSet( se, design = ~ patient + treatment )
ddsCollapsed <- collapseReplicates( ddsFull,
                                    groupby = ddsFull$sample,
                                    run = ddsFull$run )
dds <- ddsCollapsed[ , ddsCollapsed$time == "48h" ]
dds$time <- droplevels( dds$time )
dds$treatment <- relevel( dds$treatment, "Control" )
dds <- DESeq(dds)
res <- results( dds )

plotMA( res, ylim = c(-1, 1) )
plotDispEsts( dds, ylim = c(1e-6, 1e1) )
hist( res$pvalue, breaks=20, col="grey" )

res$ensembl <- sapply( strsplit( rownames(res), split="\\+" ), "[", 1 )
library( "biomaRt" )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = res$ensembl,
                  mart = ensembl )
idx <- match( res$ensembl, genemap$ensembl_gene_id )
res$entrez <- genemap$entrezgene_id[ idx ]
res$hgnc_symbol <- genemap$hgnc_symbol[ idx ]


rld <- rlog(dds)
rld_df <- as.data.frame(assay(rld))

library(dplyr)
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = rownames(rld_df),
                  mart = ensembl)
idx <- match( rownames(rld_df), genemap$ensembl_gene_id )
rld_df$entrez <- genemap$entrezgene[ idx ]
rld_df$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
rld_df <- distinct(rld_df, rld_df$hgnc_symbol, .keep_all = TRUE)
rld_df <- rld_df[-which(is.na(rld_df$hgnc_symbol) == TRUE),]
row.names(rld_df) <- rld_df$hgnc_symbol
rld_df <- rld_df[, 1:12]

##########################################################################



###Dependencies
library(msigdbr)
library(fgsea)
library(ggplot2)
library(gplots)


###Using the msigdb package, we can load in different specied and categoris/subcategories from the database;
###For details on the capabilities of this packafe, see the following;
####https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html

msigdb <- msigdbr(species = "Homo sapiens", category = "H")



###NixSEA package

###Make a Rank file from the DESeq2 results.
RNK_file <- RNK_make(res)

###Extracts the data from the database loaded in to only select data which is from genesets containing a keyword
###Changing the the keyword will allow you to search specific genesets; i.e. "MITOCHONDRIA","INSULIN", or "" for all)
sets <- keyword_geneset(msigdb, "")

###This function simply builds a list of genesets from the previously filtered data
genesets <- geneset_list(sets)

###Preforms GSEA analysis and plots the enrichment of the significant genesets
figure1 <- enrichment_analysis(genesets, RNK_file, sets, figure_header = "Significant Hallmark Geneset Enrichment")

###Creates a heatmap of the leading edge genes from the significantly enriched gene sets found above.
figure2 <- LE_heatmap(genesets, RNK_file, res, heatmap_header = "Expression of Significant Hallmark Leading Edge GSEA")







































