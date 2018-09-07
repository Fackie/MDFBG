#load and curate selected genesets
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

save(curatedgeneset,file = "curatedgeneset.Rdata")

setwd("C:/Users/Fackie/Desktop/Genesets selection/Genesets selection")
source("jackstraw.R")
source("utilities.R")
load("Results/choosendata.Rdata")
#load("Results/featuregene.Rdata")
#calculate choosen dataset scores on 424 geneset features
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
save(score,file = "6Datasets_424Genesets_FirstPC.Rdata")

