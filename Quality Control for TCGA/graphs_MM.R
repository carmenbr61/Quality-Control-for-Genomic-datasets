#Graphs
#Input: Matrix with sample wise correlations between expression data and copy number.
#Output: Histograms of correlations between matching and non-matching samples  
#from the triangular matrices of the correlation matrix obtained in the other function

#Extract triangular matrices to find the distribution of correlations of non-matching samples
plotQC<- function(corr){
#---Calculate correlations between non matching columns
nm1<-corr
nm1[lower.tri(nm1, diag = FALSE)]<- NA
nm2<-corr
nm2[upper.tri(nm2, diag =  FALSE)]<- NA
#---Calculate correlations between matching columns
mcol<- as.matrix(diag(corr))
rownames(mcol)<-rownames(corr)
colnames(mcol)<- "Cor_rna_cn"
#----Histograms---
par(mfrow=c(1,2))
hist(corr, breaks = 25, main = "Histogram of Correlations", xlab = "Correlations between expression and copy number");qqnorm(corr);qqline(corr,col = 2)
par(mfrow=c(3,1))
hist(nm1,breaks = 25, main = "Off-diagonal (upper triangular)")
hist(nm2, breaks = 25, main = "Off-diagonal (lower triangular)")
hist(mcol,breaks = 25, main = "Diagonal (matching)")   
}



