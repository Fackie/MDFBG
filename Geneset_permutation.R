# CCA distance?
# load(file = "6Datasets_424Genesets_FirstPC.Rdata")
# library(CCA)
# totalgene <- dataname
# num <- length(curatedgeneset)
# totalcor <-  vector(length=num)
# ind = 1
# for(p in curatedgeneset)
# {
#   for(i in 1:5)
#   {
#     datai = choosendata[[i]][p,]
#     for(j in (i+1):6)
#     {
#       dataj = choosendata[[j]][p,]
#       totalcor[ind] = totalcor[ind] + rcc(X = datai,Y = dataj,lambda1 = 1,lambda2 = 1)$cor[1]
#     }
#   }
#   ind <- ind +1
# }
# 
# times=1000
# varcor = vector(length = 1000)
# pnum=25
# ind <- 1
# for(t in 1:times)
# {
#   choosenum <- floor(runif(pnum,max = totalnum))
#   for(i in 1:5)
#   {
#     datai = choosendata[[i]][choosenum,]
#     for(j in (i+1):6) 
#     {
#       dataj = choosendata[[j]][choosenum,]
#       totalcor[ind] = totalcor[ind] + rcc(X = datai,Y = dataj,lambda1 = 1,lambda2 = 1)$cor[1]
#     }
#   }
#   ind <- ind + 1
# }

