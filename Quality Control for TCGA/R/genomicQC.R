#'Updated Function to get correlation matrix
#'@param  MultiAssayExperiment object subsetted by the 2 assays with tumor samples only, disjoint copy number ranges 
#' @return Matrix with correlations between expression and copy number data

library(MultiAssayExperiment)
library (RaggedExperiment)
library(TCGAutils)
library(Biobase)
library(GenomicFeatures)
library(Homo.sapiens)

genomicQC<- function(MAE){
    cnames<- colnames(MAE)[[1]]
    subsettingList<- list(rnatumoronly = TCGAsampleSelect(cnames, 01), CNVSNP = TCGAsampleSelect(colnames(ov1)[[2]], "01"))
    ov1 <- ov[, subsettingList, ]
    ov1 <- intersectColumns(ov1) #same number of samples   
    
#1) Find the corresponding feature names for the copy number ranges
    gr.all <- GenomicFeatures::genes(Homo.sapiens, column="SYMBOL")
    seqlevelsStyle(gr.all) <- "NCBI"
    geneLogic <- unname(unlist(mcols(gr.all)[[1]])) %in% rownames(MAE[[1]])
    grSubsettor <- gr.all[geneLogic]
    names(grSubsettor) <- mcols(grSubsettor)[["SYMBOL"]] # adding gene symbols as rownames
   #RaggedExperiment--> use qreduceAssay
    newSub <- assay(sub, mcolname = "Segment_Mean", ranges = grSubsettor)
    features<- rownames(MAE[[1]]) %in% unname(unlist(mcols(gr.all)[[1]]))
    MAE[[1]] <- MAE[[1]][features, ]
    MAE[[2]] <- newSub
    #Matrices
    Exmatrix<-exprs(MAE[[1]])
    CNmatrix<-assay(MAE[[2]], mcolname = "Segment_Mean")
    correlations<-as.matrix(cor(Exmatrix, CNmatrix, method= "pearson", use = "pairwise.complete.obs"))
    return(correlations)
}




