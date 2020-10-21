options(stringsAsFactors = F)
options(expressions =500000)
library(tidyverse)
library(survival)                                                
library(tableone)
library(knitr)
library(caret)
library(rbsurv)
library(survminer)
library(timeROC)

load('sameas.Rdata')   ###load the overlapping splicing events identified by Cox and RSF in step 4
load('clini.Rdata')
rownames(clini)<-clini$Id
load('ob.Rdata')
group<-c(rep('Normal',59),rep('Tumor',513))
ob<-ob[,-c(1:59)]
samesample<-intersect(colnames(ob),clini$Id)
ob<-ob[,samesample]
clinical<-clini[samesample,]
rf<-cbind(clinical[,c(4,7)],t(ob))
rf$fustat<-ifelse(rf$fustat=='Alive',0,1)
rf$futime<-rf$futime/365

####generate training set and test set randomly
set.seed(2975)
random<-createDataPartition(c(rep('a',408)),times = 1,p=0.75,list=TRUE)
train<-rf[random$Resample1,]
test<-rf[-c(random$Resample1),]
clini$Gender<-as.factor(ifelse(clini$Gender=='FEMALE','Female','Male'))
clini$Smoking.history<-as.factor(ifelse(clini$Smoking.history=='smoker',1,0))
clini$Age<-as.numeric(clini$Age)
clini$fustat<-as.factor(clini$fustat)
clini$Stage<-as.factor(clini$Stage)
clini$futime<-as.numeric(clini$futime)
clini<-clini[,-1]
pdtrain<-clini[rownames(train),][,c(1:6)]
pdtrain$group<-as.factor(1)
pdtest<-clini[rownames(test),][,c(1:6)]
pdtest$group<-as.factor(2)
pd<-rbind(pdtrain,pdtest)
write.csv(pd,file='pd.csv',row.names = T)  ####write a csv file and change blank column into NA in this file, save the changed file named as "pd-nanote.csv"
pdnote<-read.csv('pd-nanote.csv',sep=',',header = T,row.names = 1,na.strings = 'NA')
#####generate baseline table
a<-CreateTableOne(vars = c('Gender','Age','fustat','Smoking.history','Stage','futime'),
                  data=pdnote,smd=T,
                  strata = 'group',
                  factorVars = c('Gender','fustat','Smoking.history','Stage'))
summary(a)
print(a,showAllLevels = TRUE)
print(a, nonnormal = c("Age","futime")) 
a_csv<- print(a, nonnormal = c("Age","futime"),
              smd=F, 
              showAllLevels = TRUE,
              quote = FALSE, 
              noSpaces = TRUE, 
              printToggle = FALSE)
kable(a_csv,  
      align = 'c', 
      caption = 'Table 1: Comparison of grouped samples')
write.csv(a_csv, file = "baseline characteristics.csv")

###run rbsurv algorithm to identify alternative splicing events for model construction
x<-as.matrix(t(train[,sameas]))
time<-train$futime
stat<-train$fustat
fit <- rbsurv(time=time,status=stat,x=x,method="efron",n.fold=4,n.iter =100)

fitmodel<-fit$model
fitmodel<-fitmodel[fitmodel$Selected=='*       ',]
fitgene<-as.numeric(fitmodel$Gene)
fitgene<-rownames(x)[fitgene]

###calculate the coffecients of each splicing events and riskscores of patients using multivariate Cox regression
fitsurvtrain<-t(x[fitgene,])
fitsurvtrain<-fitsurvtrain[rownames(pdtrain),]
fitsurvtrain<-cbind(pdtrain[,c(3,6)],fitsurvtrain)
fitsurvtrain$fustat<-ifelse(fitsurvtrain$fustat=='Alive',0,1)
fitsurvtrain$futime<-fitsurvtrain$futime/365
regr<-coxph(Surv(futime,fustat)~.,data=fitsurvtrain)
coef<-regr$coefficients
myFun=function(x){crossprod(as.numeric(x),coef)}
riskscore<-apply(fitsurvtrain[,3:ncol(fitsurvtrain)],1,myFun)
fitsurvtrain$riskscore<-riskscore
risk<-ifelse(fitsurvtrain$riskscore>median(fitsurvtrain$riskscore),'high','low')
fitsurvtrain$risk<-risk

###draw KM curve
diff=survdiff(Surv(futime, fustat) ~risk,data = fitsurvtrain)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit.survival <- survfit(Surv(futime, fustat) ~ risk, data = fitsurvtrain)
surPlot=ggsurvplot(fit.survival, 
                   data=fitsurvtrain,
                   conf.int=FALSE,
                   pval=paste0("p=",pValue),
                   pval.size=5,
                   risk.table=TRUE,
                   legend.labs=c("High risk", "Low risk"),
                   legend.title="Risk",
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("#E64B35FF", "#4DBBD5FF"),
                   risk.table.height=.25)
pdf(file='survivaltrain.pdf',onefile = FALSE,width = 6.5,height =5.5)
print(surPlot)
dev.off()

######draw ROC curve
library(timeROC)
roctrain<-timeROC(T=fitsurvtrain$futime,
                  delta=fitsurvtrain$fustat,
                  marker=fitsurvtrain$riskscore,
                  cause=1,
                  weighting="marginal",
                  times=c(1,3,5,10),
                  ROC = TRUE)
pdf('roctrain.pdf')
plot(roctrain,time=1,title=FALSE,lwd=2.5,col='#E64B35FF')
plot(roctrain,time=3,title=FALSE,add=TRUE,lwd=2.5,col='#4DBBD5FF')
plot(roctrain,time=5,title=FALSE,add=TRUE,lwd=2.5,col='#00A087FF')
plot(roctrain,time=10,title=FALSE,add=TRUE,lwd=2.5,col='#3C5488FF')
legend('bottomright',c(paste0('AUC at 1 year: ',sprintf("%.3f", roctrain$AUC[1])),
                       paste0('AUC at 3 years: ',sprintf("%.3f", roctrain$AUC[2])),
                       paste0('AUC at 5 years: ',sprintf("%.3f", roctrain$AUC[3])),
                       paste0('AUC at 10 years: ',sprintf("%.3f", roctrain$AUC[4]))),
       col=c('#E64B35FF','#4DBBD5FF','#00A087FF','#3C5488FF'),lwd = 3,bty='n')
dev.off()


#####validate the model in the test set
fitsurvtest<-test[,c('fustat','futime',fitgene)]
riskscoretest<-apply(fitsurvtest[,3:ncol(fitsurvtest)],1,myFun)

fitsurvtest$riskscore<-riskscoretest
risk<-ifelse(fitsurvtest$riskscore>median(fitsurvtrain$riskscore),'high','low')
fitsurvtest$risk<-risk
####draw KM curve for the test set
diff=survdiff(Surv(futime, fustat) ~risk,data = fitsurvtest)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit.survival <- survfit(Surv(futime, fustat) ~ risk, data = fitsurvtest)

surPlot=ggsurvplot(fit.survival, 
                   data=fitsurvtest,
                   conf.int=FALSE,
                   pval=paste0("p=",pValue),
                   pval.size=5,
                   risk.table=TRUE,
                   legend.labs=c("High risk", "Low risk"),
                   legend.title="Risk",
                   xlab="Time(years)",
                   break.time.by = 1,
                   risk.table.title="",
                   palette=c("#E64B35FF", "#4DBBD5FF"),
                   risk.table.height=.25)
pdf(file='survivaltest.pdf',onefile = FALSE,width = 6.5,height =5.5)
print(surPlot)
dev.off()

####draw ROC curve for the test set
roctest<-timeROC(T=fitsurvtest$futime,
                 delta=fitsurvtest$fustat,
                 marker=fitsurvtest$riskscore,
                 cause=1,
                 weighting="marginal",
                 times=c(1,3,5,10),
                 ROC = TRUE)
pdf('roctest.pdf')
plot(roctest,time=1,title=FALSE,lwd=2.5,col='#E64B35FF')
plot(roctest,time=3,title=FALSE,add=TRUE,lwd=2.5,col='#4DBBD5FF')
plot(roctest,time=5,title=FALSE,add=TRUE,lwd=2.5,col='#00A087FF')
plot(roctest,time=10,title=FALSE,add=TRUE,lwd=2.5,col='#3C5488FF')
legend('bottomright',c(paste0('AUC at 1 year: ',sprintf("%.3f", roctest$AUC[1])),
                       paste0('AUC at 3 years: ',sprintf("%.3f", roctest$AUC[2])),
                       paste0('AUC at 5 years: ',sprintf("%.3f", roctest$AUC[3])),
                       paste0('AUC at 10 years: ',sprintf("%.3f", roctest$AUC[4]))),
       col=c('#E64B35FF','#4DBBD5FF','#00A087FF','#3C5488FF'),lwd = 2,bty='n')
dev.off()

save(fitgene,file='fitgene.Rdata')
save(coef,file='coef.Rdata')
