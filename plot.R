library(reshape2)
library(ggplot2)

res <- list()
ind <- 1
for(i in score)
{
  m1 <- length(i)
  m2 <- length(i[[1]]$cell.embeddings)
  
  res$i <- matrix(nrow = m1,ncol = m2,data = 0)
  for(j in 1:m1)
    res$i[j,]=i[[j]]$cell.embeddings
  m = melt(res$i)
  
  g = ggplot(m, aes(x=Var1, y=Var2, fill=value))
  g1=g+geom_tile(xlab('Samples'),ylab("Genesets"),title=names(score)[ind])
  print(g1) 
  names(res) <- names(score)[1:ind]
  ind <- ind + 1
}

