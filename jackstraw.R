require(foreach)
require(irlba)
require(doParallel)

##Container class## ----
Jackstraw <- setClass(
  Class = "Jackstraw",
  slots = list(
          probes.genes = "character",
          raw.data = "matrix",
          scale.data = "matrix",
          Fstat = "matrix",
          gene.loadings = "matrix",
          cell.embeddings = "matrix",
          nullpvalue = "matrix",
          gene.select = "character")
)

##Assistant function## ----
Pre_shuffle <- function(x) 
{
  x2 <- x
  x2 <- t(x = x)
  ind <- order(c(col(x = x2)), runif(n = length(x = x2)))
  x2 <- matrix(
    data = x2[ind],
    nrow = nrow(x = x),
    ncol = ncol(x = x),
    byrow = TRUE
  )
  return(x2)
}

Pre_permutation <- function(prop, PCs, Data) 
  # prop: proportion of shuffle data
  # PCs: PCs related
{
  randGenes <- sample(x = rownames(x = Data), size = nrow(x = Data) * prop)
  Data[randGenes, ] <- Pre_shuffle(x = Data[randGenes, ])
  pca <- svd(t(Data))
  embeddings <- pca$u
  loadings <- pca$v
  row.names(embeddings) <- colnames(Data)
  row.names(loadings) <-  row.names(Data)
  return(loadings[randGenes, 1:PCs])
}

Pre_pca <- function(object,PCs)
  # object: data container
  # Data: gene * cell matrix 
  # PCs:  pc numbers
{
  Data <- object@scale.data
  pca <- irlba(A = t(Data),nv = PCs)
  embeddings <- pca$u
  loadings <- pca$v
  row.names(embeddings) <- colnames(Data)
  row.names(loadings) <-  row.names(Data)
  object@cell.embeddings <-  embeddings
  object@gene.loadings <- loadings
  return(object)
}

##Command function## ----
Run_readdata <- function(data,type='RNAseq'){
  object = new(Class = "Jackstraw")
  if (class(data) == 'matrix')
    object@raw.data = data else
    #some study samples are NAs!!!
    object@raw.data <-  as.matrix(data)
    #object@probes.genes <- as.vector(data@featureData@data$gene_symbol) if probes
  # if(type =='Microaray')
  #   object@probes.genes <- infos
   # row.names(object@raw.data) = as.vector(data@featureData@data$gene_symbol)
  return(object)
}

Run_scaling <- function(object)#,normmethod = "Log")
# maybe shouldnt be normalized?
{
  # if (normmethod == "Log")
  # {x=object@raw.data
  # object@scale.data = log(t(t(x)/colSums(x))*1000000+1)
  # }
  object@scale.data <- t(scale(t(object@raw.data),scale = FALSE))
  return(object)
}

Run_jackstraw <- function(object,iterations=100, prop = 0.01, Data, PCs = 5) 
# iterations: setting permuatation times
# prop: setting permutation proportion
# PCs:  pc numbers
{
  #permutation F-stat/parellarization
  raw.F1 = foreach(i=1:iterations) %do% Pre_permutation(PCs = PCs,Data = object@scale.data,prop = prop)
  F1 <- sapply(X = 1:PCs,FUN = function(x) {
    return(as.numeric(x = unlist(x = lapply(
      X = 1:iterations,
      FUN = function(y) {
        return(raw.F1[[y]][, x])
      }
    ))))
  })
  F1 <- as.matrix(x = F1)
  object@Fstat <- F1
  object <- Pre_pca(object,PCs)
  return(object)
}

Run_empvalue <- function(object,PCs) 
{
  emp <- sapply(
    X = PCs,
    FUN = function(x) {
      return(unlist(x = sapply(
        X = abs(object@gene.loadings[,x]),
        FUN = function(x,null){
          return(sum(null> x) / length(null))
        },
        null = abs(object@Fstat[,x])
      )))
    }
  )
  object@nullpvalue <- emp
  return(object)
}

Run_selectgenes <- function(object,PCs=1,pval= 0.001)
{
  empvalue <- object@nullpvalue
  if (length(PCs) == 1) {
    pvals.min <- empvalue[, PCs]
  }
  if (length(PCs) > 1) {
    pvals.min <- apply(X = empvalue[, PCs], MARGIN = 1, FUN = min)
  }
  names(pvals.min) <- rownames(object@scale.data)
  genes.use <- names(x = pvals.min)[pvals.min < pval]
  object@gene.select <- genes.use
  return(object)
}

