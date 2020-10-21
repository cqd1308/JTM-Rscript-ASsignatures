options(stringsAsFactors = FALSE)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)

load('clini.Rdata')
load('coef.Rdata')
load('fitgene.Rdata')
load('ob.Rdata')

myFUN<-function(x){crossprod(x,coef)}
ob<-as.data.frame(t(ob[fitgene,]))
ob<-ob[-c(1:59),]
riskscore<-apply(ob,1,myFUN)
ob$riskscore<-riskscore
rownames(clini)<-clini$Id
clini<-clini[,-1]
samesample<-intersect(rownames(ob),rownames(clini))
clini<-clini[samesample,]
ob<-ob[samesample,]
plotdata<-cbind(riskscore=ob$riskscore,clini)
#############Stage
plotstage<-dplyr::select(plotdata,riskscore,Stage)
plotstage<-plotstage[!is.na(plotstage$Stage),]

plotstage$Stage<-gsub('Stage IV','3',plotstage$Stage)
plotstage$Stage<-gsub('Stage III','3',plotstage$Stage)
plotstage$Stage<-gsub('Stage II','2',plotstage$Stage)
plotstage$Stage<-gsub('Stage I','2',plotstage$Stage)
plotstage$Stage<-gsub('3','III-IV',plotstage$Stage)
plotstage$Stage<-gsub('2','I-II',plotstage$Stage)

p1<-ggplot(plotstage,aes(Stage,riskscore))+
  geom_violin(color='white',aes(fill=Stage))+
  geom_boxplot(fill='white',width=0.4,position = position_dodge(0.9),outlier.size =  0.2)+
  stat_compare_means(data=plotstage,aes(x=Stage,y=riskscore,label=paste('P = ',format(as.numeric(..p.format..),scientific = TRUE),sep='')),
                     method='wilcox')+
  theme_bw()+ylab('Risk score')+xlab('')+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))+
  theme(legend.position = 'none')+
  ggtitle('AJCC stage')+
  theme(plot.title = element_text(size = 10,hjust = 0.5))+
  scale_fill_manual(values = c("#00B5E2FF", "#FFCD00FF"))
##################N  
plotN<-dplyr::select(plotdata,riskscore,N)
plotN<-plotN[!is.na(plotN$N),]
plotN$N<-gsub('N3','N1',plotN$N)
plotN$N<-gsub('N2','N1',plotN$N)
plotN$N<-gsub('N1','N1',plotN$N)
plotN$N<-gsub('N1','N1-3',plotN$N)
p2<-ggplot(plotN,aes(N,riskscore))+
  geom_violin(color='white',aes(fill=N))+
  geom_boxplot(fill='white',width=0.4,position = position_dodge(0.9),outlier.size =  0.2)+
  stat_compare_means(data=plotN,aes(x=N,y=riskscore),
                     label='p.format',method='wilcox')+
  theme_bw()+ylab('Risk score')+xlab('')+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))+ggtitle('LNM')+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(legend.position = 'none')+scale_fill_manual(values=c("#00B5E2FF", "#FFCD00FF"))+
  ggtitle('N stage')
##################T 
plotT<-dplyr::select(plotdata,riskscore,T)
plotT<-plotT[!is.na(plotT$T),]
plotT$T<-gsub('T4','T3',plotT$T)
plotT$T<-gsub('T3','T3-4',plotT$T)
plotT$T<-gsub('T2','T1',plotT$T)
plotT$T<-gsub('T1','T1-2',plotT$T)
p3<-ggplot(plotT,aes(T,riskscore))+
  geom_violin(color='white',aes(fill=T))+
  geom_boxplot(fill='white',width=0.4,position = position_dodge(0.9),outlier.size =  0.2)+
  stat_compare_means(data=plotT,aes(x=T,y=riskscore),
                     label = 'p.format',method='wilcox')+
  theme_bw()+ylab('Risk score')+xlab('')+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))+ggtitle('T stage')+
  theme(plot.title = element_text(hjust = 0.5,size=10))+
  theme(legend.position = 'none')+scale_fill_manual(values = c("#00B5E2FF", "#FFCD00FF"))+
  ggtitle('T stage')
##################M  
plotM<-dplyr::select(plotdata,riskscore,M)
plotM<-plotM[!is.na(plotM$M),]
p4<-ggplot(plotM,aes(M,riskscore))+
  geom_violin(color='white',aes(fill=M))+
  geom_boxplot(fill='white',width=0.4,position = position_dodge(0.9),outlier.size =  0.2)+
  stat_compare_means(data=plotM,aes(x=M,y=riskscore),
                     label = 'p.format',method='wilcox')+
  theme_bw()+ylab('Risk score')+xlab('')+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))+ggtitle('M stage')+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(legend.position = 'none')+scale_fill_manual(values = c("#00B5E2FF", "#FFCD00FF"))+
  ggtitle('M stage')
##################Age 
plotage<-dplyr::select(plotdata,riskscore,Age)
plotage<-plotage[!is.na(plotage$Age),]
plotage$Age<-ifelse(plotage$Age<=67,'<=67','>67')
p5<-ggplot(plotage,aes(Age,riskscore))+
  geom_violin(color='white',aes(fill=Age))+
  geom_boxplot(fill='white',width=0.4,position = position_dodge(0.9),outlier.size =  0.2)+
  stat_compare_means(data=plotage,aes(x=Age,y=riskscore),
                     label = 'p.format',method='wilcox')+
  theme_bw()+ylab('Risk score')+xlab('')+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))+ggtitle('Age')+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(legend.position = 'none')+scale_fill_manual(values = c("#00B5E2FF", "#FFCD00FF"))+
  ggtitle('Age')
##################Status   
plotstatus<-dplyr::select(plotdata,riskscore,fustat)
plotstatus<-plotstatus[!is.na(plotstatus$fustat),]
p6<-ggplot(plotstatus,aes(fustat,riskscore))+
  geom_violin(color='white',aes(fill=fustat))+
  geom_boxplot(fill='white',width=0.4,position = position_dodge(0.9),outlier.size =  0.2)+
  stat_compare_means(data=plotstatus,aes(x=fustat,y=riskscore,label=paste('P = ',format(as.numeric(..p.format..),scientific=TRUE),sep='')),
                     method='wilcox')+
  theme_bw()+ylab('Risk score')+xlab('')+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))+ggtitle('Status')+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(legend.position = 'none')+scale_fill_manual(values = c("#00B5E2FF", "#FFCD00FF"))+
  ggtitle('Status')
##################Gender
plotgender<-dplyr::select(plotdata,riskscore,Gender)
plotgender<-plotgender[!is.na(plotgender$Gender),]
#gendercomparison<-list(c('Male','Female'))
p7<-ggplot(plotgender,aes(Gender,riskscore))+
  geom_violin(color='white',aes(fill=Gender))+
  geom_boxplot(fill='white',width=0.4,position = position_dodge(0.9),outlier.size =  0.2)+
  stat_compare_means(data=plotgender,aes(x=Gender,y=riskscore),
                     label = 'p.format',method='wilcox')+
  theme_bw()+ylab('Risk score')+xlab('')+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))+ggtitle('Gender')+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(legend.position = 'none')+scale_fill_manual(values = c("#00B5E2FF", "#FFCD00FF"))+
  ggtitle('Gender')

##################Smoking
plotsmoke<-dplyr::select(plotdata,riskscore,Smoking.history)
plotsmoke<-plotsmoke[!is.na(plotsmoke$Smoking.history),]
plotsmoke$Smoking.history<-ifelse(plotsmoke$Smoking.history=='smoker','Smoker','Non-smoker')
p8<-ggplot(plotsmoke,aes(Smoking.history,riskscore))+
  geom_violin(color='white',aes(fill=Smoking.history))+
  geom_boxplot(fill='white',width=0.4,position = position_dodge(0.9),outlier.size =  0.2)+
  stat_compare_means(data=plotsmoke,aes(x=Smoking.history,y=riskscore),
                     label = 'p.format',method='wilcox')+
  theme_bw()+ylab('Risk score')+xlab('')+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))+ggtitle('Smoking history')+
  theme(plot.title = element_text(hjust = 0.5,size = 10))+
  theme(legend.position = 'none')+scale_fill_manual(values = c("#00B5E2FF", "#FFCD00FF"))+
  ggtitle('Smoking history')


pdf('Figure7.pdf',width = 10,height = 5)
(p1|p3|p2|p4)/(p6|p8|p7|p5)+plot_annotation(tag_levels = 'A')
dev.off()
########plot nomogram
options(stringsAsFactors = FALSE)
library(rms)
nomo<-dplyr::select(plotdata,futime,fustat,Age,Gender,Stage,riskscore)
nomo$Stage<-gsub('Stage IV','Stage III',nomo$Stage)
nomo$Stage<-gsub('Stage III','Stage III-IV',nomo$Stage)
table(nomo$Stage)
nomo$futime<-nomo$futime/365
nomo$fustat<-ifelse(nomo$fustat=='Alive',0,1)
dd <- datadist(nomo)
options(datadist="dd")
f <- cph(Surv(futime, fustat) ~ Age+Gender+Stage+riskscore, x=T, y=T, surv=T, data=nomo, time.inc=1,singular.ok = T)
surv <- Survival(f)
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(2, x), function(x) surv(3, x)), 
                lp=F, funlabel=c("1-year survival", "2-year survival", "3-year survival"), 
                maxscale=100, 
                fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05))  
pdf(file='nomogram.pdf',height=7.5,width=11)
plot(nom)
dev.off()
#######forest plot
library(survival)
rt<-dplyr::select(clini,futime,fustat,Age,Gender,Smoking.history,Stage,T,N,M)
rt$Gender<-ifelse(rt$Gender=='MALE',1,0)
rt$riskscore<-plotdata$riskscore
colnames(rt)[5]<-'Smoking history'
rt$fustat<-ifelse(rt$fustat=='Alive',0,1)
rt$`Smoking history`<-ifelse(rt$`Smoking history`=='smoker',1,0)

rt$Stage<-gsub('Stage IV',4,rt$Stage)
rt$Stage<-gsub('Stage III',3,rt$Stage)
rt$Stage<-gsub('Stage II',2,rt$Stage)
rt$Stage<-gsub('Stage I',1,rt$Stage)
rt$Stage<-as.numeric(rt$Stage)
rt$T<-as.numeric(str_sub(rt$T,2,2))
rt$N<-as.numeric(str_sub(rt$N,2,2))
rt$M<-as.numeric(str_sub(rt$M,2,2))
######univariate Cox regression
uniTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
gene <- uniTab$id
hr <- sprintf("%.3f",as.numeric(uniTab$"HR"))
hrLow  <- sprintf("%.3f",as.numeric(uniTab$"HR.95L"))
hrHigh <- sprintf("%.3f",as.numeric(uniTab$"HR.95H"))
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(as.numeric(uniTab$pvalue)<0.001, "<0.001", sprintf("%.3f", as.numeric(uniTab$pvalue)))

pdf(file='unicox.pdf', width = 6.3,height = 4.5)
n <- nrow(uniTab)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, '#0099B4FF', '#0099B4FF')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()
#######multivariate Cox regression
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=as.data.frame(cbind(id=row.names(multiTab),multiTab))

gene <- multiTab$id
hr <- sprintf("%.3f",as.numeric(multiTab$"HR"))
hrLow  <- sprintf("%.3f",as.numeric(multiTab$"HR.95L"))
hrHigh <- sprintf("%.3f",as.numeric(multiTab$"HR.95H"))
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(as.numeric(multiTab$pvalue)<0.001, "<0.001", sprintf("%.3f", as.numeric(multiTab$pvalue)))

pdf(file='multicox.pdf', width = 6.3,height = 4.5)
n <- nrow(multiTab)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, '#E18727FF', '#E18727FF')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1)
dev.off()