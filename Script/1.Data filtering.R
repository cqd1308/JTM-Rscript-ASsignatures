options(stringsAsFactors = F)
library(impute)
####read input files. The original expression matrix file exceeds the upload limit of GitHub, so we splited it into 2 files.
rt1<-read.table('asexp-1.txt',header = T,check.names = F,row.names = 1)
rt2<-read.table('asexp-2.txt',header = T,check.names = F,row.names = 1)
rt<-rbind(rt1,rt2)
rt=as.matrix(rt)
exp=rt
dimnames=list(rownames(exp),colnames(exp))
exp=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

#impute missing expression data
mat=impute.knn(exp)
####filter data according to sd and means
data=mat$data
ob<-data[rowMeans(data)>0.05,]
ob<-as.data.frame(ob)
ob$sd<-apply(ob,1,sd)
ob<-ob[ob$sd>0.1,]
ob<-ob[,-ncol(ob)]
save(ob,file='ob.Rdata')
