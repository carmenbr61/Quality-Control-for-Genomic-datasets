# MIXTURE MODELS QUALITY CONTROL
#Carmen R.

#---------------------------Step 1: Prepare matrices  + get Pearson's correlations----------------------------------------------
#Use bladder MultiAssayExpriment  to create two matrices- one for RNASeqGene, one for copy number. 
#Their columns and rows should correspond exactly.

#--------------Bladder Urothelial Carcinoma, also try in UCS and UVM
download.file("http://s3.amazonaws.com/multiassayexperiments/blcaMAEO.rds", destfile = "~/Documents/R-work/blcaMAEO.rds")
blca<-readRDS("~/Documents/R-work/blcaMAEO.rds")
names(blca)
#Subset by assays
blca2<-blca[,, c("RNASeqGene", "CNVSNP")]
#Starning samples:
#RNA --> 67 samples
#Copy num --> 797 samples 


#Get tumor samples and complete cases only
library(TCGAutils)
library(Biobase)
#sampleTypes 
#RNASeqGene
cnames<- colnames(blca2)[[1]]
rnatumoronly<-TCGAsampleSelect(cnames, 01) #56 primary solid tumor samples in RNASeqGene
subsettingList <- list(RNASeqGene = rnatumoronly, CNVSNP = TCGAsampleSelect(colnames(blca2)[[2]], "01"))
bba <- blca2[, subsettingList, ]
#[1] RNASeqGene: ExpressionSet with 20502 rows and 56 columns 
#[2] CNVSNP: RangedRaggedAssay with 79410 rows and 408 columns 

blca_complete <- bba[, complete.cases(bba), ]#Complete cases only


#Need to find the corresponding names for the copy number ranges -so that the rows are the gene names and the columns are the samples(participants)
library(GenomicFeatures)
library(Homo.sapiens)

gr.all = genes(Homo.sapiens, column="SYMBOL")
seqlevelsStyle(gr.all) <- "NCBI"
geneLogic <- unname(unlist(mcols(gr.all)[[1]])) %in% rownames(blca_complete[[1]])
length(geneLogic[geneLogic==TRUE]) #17,882
grSubsettor <- gr.all[geneLogic]
names(grSubsettor) <- mcols(grSubsettor)[["SYMBOL"]] # adding gene symbols as rownames
sub <- blca_complete[[2]][grSubsettor, ]
newSub <- assay(sub, mcolname = "Segment_Mean", ranges = grSubsettor)
featuresRNA <- rownames(blca_complete[[1]]) %in% unname(unlist(mcols(gr.all)[[1]]))
blca_complete[[1]] <- blca_complete[[1]][featuresRNA, ]
blca_complete[[2]] <- newSub


#------Matrices-------------
RNAmatrix<-exprs(blca_complete[[1]])
dim(RNAmatrix)  #check: 17882 rows ( genes) and 55 tumor samples.
CNVmatrix<-assay(blca_complete[[2]], mcolname = "Segment_Mean")
dim(CNVmatrix)#check: 17882 rows ( genes) and 55 tumor samples.

rnamesL<-rownames(RNAmatrix) %in% rownames(CNVmatrix) #okay!
CNVmatrixs<-scale(CNVmatrix)
RNAmatrixs<-scale(RNAmatrix)
RNAmatrixs<-RNAmatrixs[sort(rownames(RNAmatrixs)), ]
CNVmatrixs<-CNVmatrixs[sort(rownames(CNVmatrixs)), ]


RNAmatrixt<-RNAmatrix[sort(rownames(RNAmatrix)), ]
CNVmatrixt<-CNVmatrix[sort(rownames(CNVmatrix)), ]
#checks
rownames(RNAmatrixs)[1:10]
rownames(CNVmatrixs)[1:10]


#------Correlations --> OVERALL!----
correlations<-as.matrix(cor(RNAmatrix, CNVmatrix, method= "pearson", use = "pairwise.complete.obs"))


#TRIANGULAR MATRICES
# Calculate correlations between non matching columns
nmatchingcols1<-correlations
nmatchingcols1[lower.tri(nmatchingcols1, diag = FALSE)]<- NA

nmatchingcols2<-correlations
nmatchingcols2[upper.tri(nmatchingcols2, diag =  FALSE)]<- NA
    
# Calculate correlations between matching columns
matchingcols<- as.matrix(diag(correlations))
rownames(matchingcols)<-rownames(correlations)
colnames(matchingcols)<- "Cor_rna_cn"


#----Histograms---
#Overall
hist(correlations, breaks = 25, xlab = "Correlations between RNASeqGene and CNVSNP for bladderC")
qplot(correlations, geom="histogram", binwidth = 0.5,main = "Histogram of correlations")




#Non-matched
hist(nmatchingcols1, col = "red", xlab = "Correlations between RNASeqGene and CNVSNP for bladderC")
hist(nmatchingcols2, col = "red")
#Matched 
hist(matchingcols,breaks = 25, col = "green", xlab = "Correlations between RNASeqGene and CNVSNP for bladderC")
#Overlap matched and non-matched samples
hist(nmatchingcols1, breaks = 25, main = "", xlab = "Correlations between RNASeqGene and CNVSNP for bladderC")
hist(matchingcols,breaks = 25, col = "green", add = T)
box()



#-----SAME AS ABOVE BUT USING  CNASNP
blca3<-blca[,, c("RNASeqGene", "CNASNP")]
subsettingList2 <- list(RNASeqGene = rnatumoronly, CNASNP = TCGAsampleSelect(colnames(blca3)[[2]], "01"))
bba2 <- blca3[, subsettingList2, ]
blca_complete2 <- bba2[, complete.cases(bba2), ]#Complete cases only
#copy number
geneLogic2 <- unname(unlist(mcols(gr.all)[[1]])) %in% rownames(blca_complete2[[1]])
grSubsettor2 <- gr.all[geneLogic2]
names(grSubsettor2) <- mcols(grSubsettor2)[["SYMBOL"]] # adding gene symbols as rownames
sub2 <- blca_complete2[[2]][grSubsettor2, ]
newSub2 <- assay(sub2, mcolname = "Segment_Mean", ranges = grSubsettor2)
featuresRNA2 <- rownames(blca_complete2[[1]]) %in% unname(unlist(mcols(gr.all)[[1]]))
blca_complete2[[1]] <- blca_complete2[[1]][featuresRNA2, ]
blca_complete2[[2]] <- newSub2




#using gistic a and gistic t;
names(blca)
blca3<-blca[,, c("RNASeqGene", "gistica", "gistict")]
cnames<- colnames(blca2)[[1]]
rnatumoronly<-TCGAsampleSelect(cnames, 01) 5
subsettingList <- list(RNASeqGene = rnatumoronly, CNVSNP = TCGAsampleSelect(colnames(blca2)[[2]], "01"))
bba <- blca2[, subsettingList, ]

blca_complete <- bba[, complete.cases(bba), ]#Complete cases only


#Need to find the corresponding names for the copy number ranges -so that the rows are the gene names and the columns are the samples(participants)
library(GenomicFeatures)
library(Homo.sapiens)

gr.all = genes(Homo.sapiens, column="SYMBOL")
seqlevelsStyle(gr.all) <- "NCBI"
geneLogic <- unname(unlist(mcols(gr.all)[[1]])) %in% rownames(blca_complete[[1]])
length(geneLogic[geneLogic==TRUE]) #17,882
grSubsettor <- gr.all[geneLogic]
names(grSubsettor) <- mcols(grSubsettor)[["SYMBOL"]] # adding gene symbols as rownames
sub <- blca_complete[[2]][grSubsettor, ]
newSub <- assay(sub, mcolname = "Segment_Mean", ranges = grSubsettor)
featuresRNA <- rownames(blca_complete[[1]]) %in% unname(unlist(mcols(gr.all)[[1]]))
blca_complete[[1]] <- blca_complete[[1]][featuresRNA, ]
blca_complete[[2]] <- newSub

#--------------Ovarian serous cystadenocarcinoma-----
download.file("http://s3.amazonaws.com/multiassayexperiments/ovMAEO.rds", destfile = "~/Documents/R-work/ovMAEO.rds")
ov<-readRDS("~/Documents/R-work/ovMAEO.rds")
names(ov)
ov1<-ov[,, c("RNASeqGene", "CNVSNP")]

#Get tumor samples and complete cases only
cnames_ov<- colnames(ov1)[[1]]
rnatumoronly_ov<-TCGAsampleSelect(cnames_ov, 01) #295 tumor samples in RNASeqGene
subsettingList_ov <- list(RNASeqGene = rnatumoronly_ov, CNVSNP = TCGAsampleSelect(colnames(ov1)[[2]], "01"))
ov2 <- ov1[, subsettingList_ov, ]
# [1] RNASeqGene: ExpressionSet with 19990 rows and 295 columns 
#[2] CNVSNP: RangedRaggedAssay with 200363 rows and 573 columns

ov_complete <- ov2[, complete.cases(ov2), ]#Complete cases only
# [1] RNASeqGene: ExpressionSet with 19990 rows and 292 columns 
#[2] CNVSNP: RangedRaggedAssay with 101827 rows and 292 columns 

#Need to find the corresponding names for the copy number ranges -so that the rows are the gene names and the columns are the samples(participants)
geneLogic_ov <- unname(unlist(mcols(gr.all)[[1]])) %in% rownames(ov_complete[[1]])#17,736

grSubsettor_ov <- gr.all[geneLogic_ov]

names(grSubsettor_ov) <- mcols(grSubsettor_ov)[["SYMBOL"]] # adding gene symbols as rownames

subsample_ov<- ov_complete[[2]][grSubsettor_ov, ]

newSub_ov <- assay(subsample_ov, mcolname = "Segment_Mean", ranges = grSubsettor_ov)

#Genes in common only!!
featuresRNA_ov <- rownames(ov_complete[[1]]) %in% unname(unlist(mcols(gr.all)[[1]]))
ov_complete[[1]] <- ov_complete[[1]][featuresRNA_ov, ]
ov_complete[[2]] <- newSub_ov


sum(isDisjoint(subsample_ov))
disjoin(subsample_ov)
s<-disjoin(subsample_ov, mcolname = "Segment_Mean")



RNAmatrix_ov<- exprs(ov_complete[[1]])
CNVmatrix_ov<-assay(ov_complete[[2]], mcolname = "Segment_Mean")

corr<-cor(RNAmatrix_ov,CNVmatrix_ov,  method= "pearson", use = "pairwise.complete.obs")
hist(corr)
hist(upper.tri(corr))



#library("RaggedExperiment")
#re<- as(subsample_ov, "RaggedExperiment")
#cn<-qreduceAssay(re, grSubsettor_ov, simplify= mean)

# if (sum(isDisjoint(sub)) !=  length(colnames(MAE_complete)[[2]])  ) {
#     sub<- disjoin(sub,mcolname = "Segment_Mean" )
# }








