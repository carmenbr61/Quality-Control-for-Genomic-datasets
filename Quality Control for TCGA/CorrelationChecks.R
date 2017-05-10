#Correlation checks across all assays in the ovarian cancer dataset
source("https://bioconductor.org/biocLite.R")
biocLite("MultiAssayExperiment")
library(BiocGenerics)
library(TCGAutils)
library(Biobase)
library(GenomicRanges)
library(MultiAssayExperiment)
library (RaggedExperiment)
ov<-readRDS("~/Documents/R-work/ovMAEO.rds")
names(ov)
ov<-updateObject(ov)
 
#-------------------------Correlations between ExpressionSets--------------
#RNASeqGene and RNASeq2GeneNorm
ov1<-ov[,, c("RNASeqGene", "RNASeq2GeneNorm")]
#Ran with and without subsetting by primary solid tumor samples
cnames<- colnames(ov1)[[1]]
cnames2<- colnames(ov1)[[2]]
rnatumoronly<-TCGAsampleSelect(cnames, 01)
rnatumoronly2<-TCGAsampleSelect(cnames2, 01)
subsettingList <- list(RNASeqGene = rnatumoronly, RNASeq2GeneNorm = rnatumoronly2)
ov2 <- ov1[, subsettingList, ]
ov3 <- intersectColumns(ov2)#Complete cases only
#Find gene intersection between these two assays-->same features and columns 
ov3<-intersectRows(ov3)
#------Matrices-------------
ov3a <- assays(ov3)
# ov3a[[2]] <- ov3a[[2]][match(rownames(ov3a[[1]]), rownames(ov3a[[2]])), ]
identical(rownames(ov3a[[1]]), rownames(ov3a[[2]]))
# colnames(ov3a[[1]]) <- make.names(TCGAbarcode(colnames(ov3a[[1]])))
# colnames(ov3a[[2]]) <- make.names(TCGAbarcode(colnames(ov3a[[2]])))
all.equal(rownames(ov3a[[1]]), rownames(ov3a[[2]]))


# ii <- Reduce(intersect, rownames(ov3))
# identical(sort(ii), sort(i))
# rownames(ov3)[[1]]
# ov3[ii, ,]
# ov4 <- ov3[ii, , ]
# ov4a<-assays(ov4)
a<-ov3a[[1]]
b<-ov3a[[2]]
cc<-cor(a,b)
hist(cc,main = "RNASeqGene v. RNASeq2GeneNorm")


#RNASeq2GeneNorm and mRNAArray

#RNASeqGene and mRNAArray
ov1<-ov[,, c("RNASeqGene", "mRNAArray")]
cnames<- colnames(ov1)[[1]]
cnames2<- colnames(ov1)[[2]]
rnatumoronly<-TCGAsampleSelect(cnames, 01)
rnatumoronly2<-TCGAsampleSelect(cnames2, 01)
subsettingList <- list(RNASeqGene = rnatumoronly, mRNAArray = rnatumoronly2)
ov2 <- ov1[, subsettingList, ]
ov3 <- ov2[, complete.cases(ov2), ]#Complete cases only
#Find gene intersection between these two assays--same number of features and columns 
i<- intersect(rownames(ov3)[[1]], rownames(ov3)[[2]])
ov3<-ov3[i, , ]
#------Matrices-------------
ov3a <- assays(ov3)
ov3a[[2]] <- ov3a[[2]][match(rownames(ov3a[[1]]), rownames(ov3a[[2]])), ]
identical(rownames(ov3a[[1]]), rownames(ov3a[[2]]))
colnames(ov3a[[1]]) <- make.names(TCGAbarcode(colnames(ov3a[[1]])))
colnames(ov3a[[2]]) <- make.names(TCGAbarcode(colnames(ov3a[[2]])))
all.equal(colnames(ov3a[[1]]), colnames(ov3a[[2]]))

a<-ov3a[[1]]
b<-ov3a[[2]]
c<-cor(a,b)
hist(c,main = "RNASeqGene v. mRNAArray")

#--------------------Correlations between CN and expression data assays-----

#---CNVSNP/CNASNP and mRNAArray---------> NOT Finished!
ov1<-ov[,, c("mRNAArray", "CNVSNP")]
cnames<- colnames(ov1)[[1]]
rnatumoronly<-TCGAsampleSelect(cnames, 01)
subsettingList<- list(mRNAArray = rnatumoronly, CNVSNP = TCGAsampleSelect(colnames(ov1)[[2]], "01"))
ov2 <- ov1[, subsettingList, ]
ov3 <- ov2[, complete.cases(ov2), ]
#Find gene intersection between these two assays--same number of features and columns !
#Need to find the corresponding names for the copy number ranges -so that the rows are the gene names and the columns are the samples(participants)
library(GenomicFeatures)
library(Homo.sapiens)
gr.all = genes(Homo.sapiens, column="SYMBOL")
seqlevelsStyle(gr.all) <- "NCBI"
geneLogic<- unname(unlist(mcols(gr.all)[[1]])) %in% rownames(ov3[[1]])
grSubsettor<- gr.all[geneLogic]
names(grSubsettor) <- mcols(grSubsettor)[["SYMBOL"]] # adding gene symbols as rownames
#subs<- ov3[[2]][grSubsettor, ]
newSub<- qreduceAssay(ov3[[2]],grSubsettor, weightedmean, "Segment_Mean") #needs check
#Genes in common only-- intersection. Not sure how to use intersect() in this case
featuresmRNA <- rownames(ov3[[1]]) %in% unname(unlist(mcols(gr.all)[[1]]))
ov3[[1]] <- ov3[[1]][featuresmRNA, ]
ov3[[2]] <- newSub

mRNA<- exprs(ov3[[1]])
#CNV<-assay(ov3[[2]], mcolname = "Segment_Mean")
#CNA<-assay(ov3[[2]], mcolname = "Segment_Mean")
corr<-cor(mRNA,CNV,  method= "pearson", use = "pairwise.complete.obs")
hist(corr)


#---Gistica and RNASeq2GeneNorm---> intersect function does not work in this case.
ov1<-ov[,, c("RNASeq2GeneNorm", "gistica")]
i<- intersect(rownames(ov1)[[1]], rownames(ov1)[[2]])
ov3<-ov3[i, , ]
m1<-exprs(ov1[[1]])
m2<-assay(ov1[[2]])


#---Gistica and mRNAArray

