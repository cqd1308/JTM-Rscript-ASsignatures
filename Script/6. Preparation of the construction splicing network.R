options(stringsAsFactors = F)
library(tidyverse)
library(limma)
library(survival)
load('clini.Rdata')
rt = read.table("SFexp.txt", row.names=1 ,header=T,sep="\t",check.names=F)   #import the expression data of splicing factors
load('corob.Rdata')    ###load the expression data of the 99 survival-related splicing events

AS=t(corob)
rownames(AS)=gsub("\\|","\\-",rownames(AS))
group=sapply(strsplit(colnames(rt),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
rtT=rt[,group==0]
colnames(rtT)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(rtT))
rtT<-as.matrix(t(rtT))
rtT<-avereps(rtT)
rtT<-as.data.frame(t(rtT))
sameSample=intersect(colnames(rtT),colnames(AS))
rtT<-rtT[,sameSample]
same<-intersect(colnames(rtT),clini$Id)
rtT<-rtT[,same]
rownames(clini)<-clini$Id
clini<-clini[same,]
rtT<-as.data.frame(t(rtT))
rtT$fustat<-clini$fustat
rtT$futime<-clini$futime
rtT$futime<-rtT$futime/365
rtT$fustat<-ifelse(rtT$fustat=='Alive',0,1)
###Cox regression analysis for splicing factors
pFilter=0.01
coxTab<-data.frame()
for(i in colnames(rtT[,1:(ncol(rtT)-2)])){
  cox <- coxph(Surv(futime, fustat) ~ rtT[,i], data = rtT)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  coxTab=rbind(coxTab,
               cbind(id=i,
                     z=coxSummary$coefficients[,"z"],
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
coxTab$pvalue<-as.numeric(coxTab$pvalue)
coxTab<-coxTab[coxTab$pvalue<0.01,]
coxTab<-dplyr::filter(coxTab,rownames(coxTab)!='NA')

corFilter=0.2             
pvalueFilter=0.01        


SF = read.table("SFexp.txt", row.names=1 ,header=T,sep="\t",check.names=F)   #????SF????????
SF<-SF[coxTab$id,]
load('corob.Rdata')
AS=t(corob)
rownames(AS)=gsub("\\|","\\-",rownames(AS))
group=sapply(strsplit(colnames(SF),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
SF=SF[,group==0]
colnames(SF)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(SF))
sameSample=intersect(colnames(SF),colnames(AS))
SF1=SF[,sameSample]
AS1=AS[,sameSample]

####analyze correlation between alternative splcing events and spling factors
outTab=data.frame()
for(i in row.names(SF1)){
  if(sd(SF1[i,])>1){
    for(j in row.names(AS1)){
      x=as.numeric(SF1[i,])
      y=as.numeric(AS1[j,])
      corT=cor.test(x,y,method = 'spearman')
      cor=corT$estimate
      pvalue=corT$p.value
      if((cor>corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(SF=i,AS=j,cor,pvalue,Regulation="postive"))
      }
      if((cor< -corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(SF=i,AS=j,cor,pvalue,Regulation="negative"))
      }
    }
  }
}



load('coxgene.Rdata')  ####import Cox results in step 4 to identify risk or protective alternative splicing events
rownames(coxgene)<-coxgene$id
rownames(coxgene)=gsub("\\|","\\-",rownames(coxgene))
asUp=coxgene[coxgene$z>0,]
asDown=coxgene[coxgene$z<0,]
SFLabel=cbind(rownames(SF),"SF")
ASupLabel=cbind(rownames(asUp),"ASup")
ASdownLabel=cbind(rownames(asDown),"ASdown")
nodeLabel=rbind(c("ID","Classify"),SFLabel,ASupLabel,ASdownLabel)
write.table(nodeLabel,file="nodeType.txt",sep="\t",quote=F,col.names=F,row.names=F)

