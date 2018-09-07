Datafusion <- setClass(
  Class = "Datafusion",
  slots = list(
   # raw.data = "matrix",
   # scale.data = "matrix",
    gene.loadings = "matrix",
    cell.embeddings = "matrix")
)


numericListOverlap<- function(test1, test2){
  t1 <- c(test1)
  names(t1) <- rownames(test1)
  t2 <- c(test2)
  names(t2) <- rownames(test2)
  
  t1 <- sort(t1,decreasing = TRUE)
  t2 <- sort(t2,decreasing = TRUE)
  
  i <- names(t1)
  j <- names(t2)
  
  n<- length(t1)
  
  overlap<- function(a,b) {
    count<-as.integer(sum(as.numeric(i[1:a] %in% j[1:b])))    
    log.pval<- 1-phyper(q=count-1, m=a, n=n-a, k=b)        
    return(c(counts=count, log.pval=as.numeric(log.pval)))    
  }
  
  indexes<- expand.grid(i=seq(1,n), j=seq(1,n))
  overlaps<- apply(indexes, 1, function(x) overlap(x['i'], x['j']))
  
  nrows<- sqrt(ncol(overlaps))
  matrix.counts<- matrix(overlaps['counts',], ncol=nrows)  
  matrix.log.pvals<- matrix(overlaps['log.pval',], ncol=nrows)  
  
  return(list(counts=matrix.counts, log.pval=matrix.log.pvals))
}

test1 <- score$study_2034_GPL96_all$GNF2_BUB3$gene.loadings
test2 <- score$study_17705_GPL96_JBI_Tissue_BC_Tamoxifen$GNF2_BUB3$gene.loadings

a=numericListOverlap(test1,test2)
tp1=a$log.pval
levelplot(tp1)

rrhopvalue<- function(i,j) {

  
  count<-as.integer(sum(as.numeric(a[1:a] %in% sample2[1:b])))    
  log.pval<- -phyper(q=count-1, m=a, n=n-a, k=b, lower.tail=FALSE, log.p=TRUE)         
  signs<- 1L
  
  return(c(counts=count, 
           log.pval=as.numeric(log.pval),
           signs=as.integer(signs)
  ))    
}