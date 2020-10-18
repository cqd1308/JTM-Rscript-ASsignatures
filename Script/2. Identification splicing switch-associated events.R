options(stringsAsFactors = F)
options(expressions =500000)
library(tidyverse)
library(randomForest)
library(doBy)
library(Boruta)
library(ROSE)
library(pheatmap)

load('ob.Rdata')  #load the expression file generated in step 1
rt<-ob
group<-c(rep('Normal',59),rep('Tumor',513))
dat<-as.data.frame(t(rt))

dat$Group<-group
Id<-colnames(dat)[-ncol(dat)]
###When the colnames have symbols such as "|" or "-", the Boruta algorithm reports error. Thus, we removed the symbols.
colnames(dat)<-str_remove_all(colnames(dat),'\\|')
colnames(dat)<-str_remove_all(colnames(dat),'\\-')
dat$Group<-as.factor(dat$Group)

####oversample data using ROSE package
over<-ovun.sample(Group~.,data=dat,method='over',N=1026)$data

####Boruta feture selection
boruta_train<-Boruta(Group~.,data=over,maxRuns=1000)
final_boruta <- TentativeRoughFix(boruta_train)

#####importance extraction
imp<-attStats(final_boruta)
imptruename1<-imp
imptruename<-imptruename1
rownames(imptruename)<-Id
imp<-imptruename[imptruename$decision=='Confirmed',]
imp<-imp[order(imp$medianImp,decreasing = T),]

###draw heat map
overplot<-over[,c(rownames(imptruename1)[imptruename1$decision=='Confirmed'],'Group')]
overplot<-overplot[order(overplot$Group,decreasing = T),]
overgroup<-data.frame(overplot$Group)
rownames(overgroup)<-rownames(overplot)
names(overgroup)<-'Group'
overplot1<-as.data.frame(t(overplot[,1:(ncol(overplot)-1)]))
overgroup$Group<-factor(overgroup$Group,levels = c('Normal','Tumor'))
pheatmap(overplot1,color = colorRampPalette(c('navy','white','firebrick3'))(50),
         annotation = overgroup,cluster_cols = FALSE,cluster_rows = TRUE,
         show_rownames = FALSE,show_colnames = FALSE,
         width = 12,height = 6,filename = 'Heatmap1.pdf')


rt<-dat[,c(rownames(imptruename1),'Group')]
colnames(rt)<-c(rownames(imptruename),'Group')
rt<-rt[,c(rownames(imp),'Group')]
group<-rt$Group
table(group)

rt<-select(rt,-Group)
rt<-as.data.frame(t(rt))
###extract parent genes of alternative splicing events
gene<-sapply(str_split(rownames(rt),'\\|'),'[',1)
gene<-gene[duplicated(gene)]
gene<-gene[!duplicated(gene)]
rt$gene<-sapply(str_split(rownames(rt),'\\|'),'[',1)
rt<-rt[rt$gene %in% gene,]

####calculate correlations between two splicing events from the same parent gene and extract genes with more than two splicing isoforms
tmp1<-data.frame()
manyase<-as.matrix(rt[1,])
manyase<-manyase[-1,]
for (i in rt$gene){ob<-rt[rt$gene==i,]
if (nrow(ob)< 3) {tmp2<-cor.test(as.numeric(ob[1,1:(ncol(ob)-1)]),as.numeric(ob[2,1:(ncol(ob)-1)]),method = 'spearman')
tmp2<-data.frame(pval=tmp2$p.value,cor=tmp2$estimate,gene=i)
tmp1<-rbind(tmp1,tmp2)}
else {manyase<-rbind(manyase,as.matrix(ob))}}
tmp3<-tmp1[order(tmp1$cor,tmp1$gene,decreasing = FALSE),]


####calculate the correlations of splicing events from genes with more than two isoforms one by one manually
manyase<-manyase[!duplicated(rownames(manyase)),]
rownamesase<-rownames(manyase)
manyase<-apply(manyase,2,as.numeric)
rownames(manyase)<-rownamesase
manyase<-as.data.frame(manyase)
manyase$gene<-sapply(str_split(rownames(manyase),'\\|'),'[',1)
FN1<-t(manyase[1:3,(1:ncol(manyase)-1)])
corFN1<-cor(FN1,method = 'spearman')
QKI<-t(manyase[manyase$gene=='QKI',(1:ncol(manyase)-1)])
corQKI<-cor(QKI,method = 'spearman')
SULT1A3<-t(manyase[manyase$gene=='SULT1A3',(1:ncol(manyase)-1)])
corSULT1A3<-cor(SULT1A3,method = 'spearman')
VEGFA<-t(manyase[manyase$gene=='VEGFA',(1:ncol(manyase)-1)])
corVEGFA<-cor(VEGFA,method = 'spearman')
FAM72A<-t(manyase[manyase$gene=='FAM72A',(1:ncol(manyase)-1)])
corFAM72A<-cor(FAM72A,method = 'spearman')
CD44<-t(manyase[manyase$gene=='CD44',(1:ncol(manyase)-1)])
corCD44<-cor(CD44,method = 'spearman')
C20orf96<-t(manyase[manyase$gene=='C20orf96',(1:ncol(manyase)-1)])
corC20orf96<-cor(C20orf96,method = 'spearman')
TUBB3<-t(manyase[manyase$gene=='TUBB3',(1:ncol(manyase)-1)])
corTUBB3<-cor(TUBB3,method = 'spearman')
###final list of splicing events with correlations = -1
pair<-c('ARHGEF39','ARHGEF4','ARIH2','SULT2B1','UPK3B','VSTM4','JAM2','PPP3CB','RASEF','AP1S2','MTA3','AZI2','MGAT4B','PVRL2',
        'POM121','MSL1','INF2','ADAMTS2','LDB1','ITIH5','ATP8B3','SYP')
pair<-rt[rt$gene %in% pair,]
multiname<-c('QKI|78404|AT','QKI|78405|AT','FAM72A|9576|AP','FAM72A|9575|AP')
pair<-rbind(pair,rt[multiname,])

###draw heatmap
group<-as.data.frame(dat$Group,row.names = colnames(pair)[1:(ncol(pair)-1)])
names(group)<-'Group'
#group<-group[order(group$Group,decreasing = TRUE),,drop=FALSE]
group$Group<-factor(group$Group,levels = c('Normal','Tumor'))
pair<-pair[order(pair$gene),]
pair<-pair[,rownames(group)]
pheatmap(pair[,1:(ncol(pair)-1)],color = colorRampPalette(c('navy','white','firebrick3'))(50),
         annotation = group,cluster_cols = FALSE,cluster_rows = FALSE,fontsize_row = 7,
         show_rownames = TRUE,show_colnames = FALSE,gaps_row = seq(from=2,to=46,by=2),
         width = 12,height = 6,filename = 'Heatmap2.pdf')