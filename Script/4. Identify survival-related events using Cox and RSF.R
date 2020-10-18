options(stringsAsFactors = F)
options(expressions =500000)
library(tidyverse)
library(randomForestSRC)
library(survival)                                            

load('clini.Rdata')  ####load clinical data
load('ob.Rdata')
rownames(clini)<-clini$Id
rt<-ob
group<-c(rep('Normal',59),rep('Tumor',513))
ob<-rt[,-c(1:59)]
samesample<-intersect(colnames(ob),clini$Id)
ob<-ob[,samesample]
clinical<-clini[samesample,]
rf<-cbind(clinical[,c(4,7)],t(ob))
rf$fustat<-ifelse(rf$fustat=='Alive',0,1)
rf$futime<-rf$futime/365
######randomsurvivalForest analysis
vs.as <- var.select(Surv(futime,fustat) ~ .,rf,ntree=1000,importance = TRUE)
rfgene<-vs.as$varselect[vs.as$topvars,]

#####univariate Cox analysis
pFilter=0.05
outTab=data.frame()
for(i in colnames(rf[,3:ncol(rf)])){
  cox <- coxph(Surv(futime, fustat) ~ rf[,i], data = rf)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     z=coxSummary$coefficients[,"z"],
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
outTab = outTab[is.na(outTab$pvalue)==FALSE,]
outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
coxgene=outTab[as.numeric(as.vector(outTab$pvalue))<pFilter,]
sameas<-intersect(coxgene$id,rownames(rfgene))

save(sameas,file='sameas.Rdata')
save(coxgene,file='coxgene.Rdata')