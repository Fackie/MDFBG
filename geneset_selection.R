# data loading and pre-filtering----
library(curatedBreastData)
setwd("C:/Users/Fackie/Desktop/Genesets selection/Genesets selection")
source("jackstraw.R")
source("utilities.R")
data("curatedBreastDataExprSetList")

breastdata <- curatedBreastDataExprSetList
studysamples <- sapply(X = breastdata,FUN = function(x){dim(exprs(x))[2]})
filterstudy <- names(breastdata)[studysamples>100]
breastdata <- breastdata[filterstudy] #34 samples to 12 samples

#remove samples with NA values
index <- sapply(X = breastdata,FUN = function(x){sum(is.na(exprs(x)))==0})
curateddata <- breastdata[index]
curatedstudy <- filterstudy[index]

#convert probedata to genedata
convertdata <- list()
for(i in 1:9) 
{
  data <- curateddata[[i]]
  data <- collapseDupProbes(expr=exprs(data),keys=data@featureData$gene_symbol,method = c("highestVariance"), debug = TRUE, removeNA_keys = TRUE,varMetric = c("everything"))
  convertdata$i <- as.matrix(data[[1]])
  names(convertdata) <- curatedstudy[1:i]
}

# jackstraw selecting variable genes----
### here we can find that 6 of 9 datasets have 12722 gene features,
### thus remain 3 datasets have 18260 gene features,
### so we use the first 6 datasets to do the following procedures 

choosendata <- convertdata[c(1,2,4,5,7,8)]
save(choosendata,file = "choosendata.Rdata")
a <-  1
result <- list()
system.time(for(i in choosendata) #about 20 mins
#foreach(i = curateddata)%do%
{
  #i=choosendata[[1]]
  fc <- Run_readdata(i,type='geneprofile')
  fc <- Run_scaling(fc)
  fc <- Run_jackstraw(fc,iterations = 100,prop = 0.01,PCs = 10)
  fc <- Run_empvalue(fc,PCs = 1:10)
  fc <- Run_selectgenes(object = fc,PCs = 1,pval = 0.005)
  result$i <- fc
  names(result) <- filterstudy[1:a]
  a <- a+1
  #return(fc)
})
save(result,file="Results/10PCs_100iters_0.01prop_normalization_result.Rdata")

### intrastudy-intersection selecting significant genesets
{
load(file = "Results/10PCs_100iters_0.01prop_normalization_result.Rdata")
for(i in 1:6)
  result[[i]] <- Run_selectgenes(object = result[[i]],PCs = 1:10,pval = 0.0001) 
intra <- c()
for(i in result)
  intra <- union(intra,i@gene.select) 

genesets <- read.table("Results/msigdb.v6.2.symbols.gmt",sep = '\t',row.names = 1,colClasses = 'character')
genesets <- t(t(genesets))
geneset.name <- row.names(genesets)
dataname <- rownames(choosendata[[1]])#12722 genes

ratiores <- c()
geneset.select <- c()
for(i in 1:17810)
{
  if(sum(genesets[i,]=='') !=0)
    tmpnums <- min(which(genesets[i,]=='')) - 1 else 
    tmpnums <- length(genesets[i,])
  tmp <- genesets[i,1:tmpnums]
  tmp <- intersect(tmp,dataname)
  genenums <- length(tmp)
  if (genenums ==0)
    ratio <- 0 else #intra: 4850 genes
  {intranums <- length(intersect(intra,tmp))
  ratio <- intranums/genenums}
  ratiores <- c(ratiores,ratio)
  if(ratio >=0.9)# 90% genes selected by jackstraw in study
    geneset.select <- c(geneset.select,geneset.name[i])
}

resgeneset <- genesets[geneset.select,]
#chosen genesets
write.table(resgeneset,file = 'resgeneset.txt',sep = '\t',quote = FALSE,col.names = FALSE)
}

### load and filter selected genesets by corcondance(RRHO)
{
genesets <- readLines(file("resgeneset.txt"))
featuregene <- strsplit(genesets,split = "\t")
names(featuregene)=as.character(1:length(featuregene))
nm <- c()
curatedgeneset <- list()
genesetlength <- c()
for(i in 1:length(featuregene))
{
  tmp <- featuregene[[i]]
  tmp <- tmp[!tmp=='']
  if(length(tmp)>=10 && length(tmp)<=100)
  {
    curatedgeneset$i <- intersect(tmp[3:length(tmp)],dataname)
    nm <- c(nm,tmp[1])
    genesetlength <- c(genesetlength,length(tmp))
    names(curatedgeneset) <- nm
  }
}

#Calculate PCA score of Geneset-Dataset
load("Results/choosendata.Rdata")
source('Fusionutility.R')
score <- list()
a <- 1
for(i in choosendata)
{
  raw.data <- i
  data <- t(scale(t(raw.data),scale = FALSE)) #!!!!!!!!!!!!!!!
  score$i <- list()
  b <-  1
  for(j in curatedgeneset)
  {
    score$i$j <- list()
    intergenes <- intersect(j,row.names(data))
    Data <- data[intergenes,]
    pca <- svd(x = t(Data),nv = 1,nu=1)
    embeddings <- pca$u
    loadings <- pca$v ##### loadings
    row.names(embeddings) <- colnames(Data)
    row.names(loadings) <-  row.names(Data)
    res <- new(Class = "Datafusion")
    score$i$j$cell.embeddings <-  embeddings
    score$i$j$gene.loadings <- loadings
    names(score$i) <- names(curatedgeneset)[1:b]
    b <- b+1
  }
  names(score) <- names(choosendata)[1:a]
  a <- a+1
}
#RRHO algotithm examing corcondance of each genesets
##1.Total RRHO score
# Calculate co-direction RRHO score than get mean value

}


