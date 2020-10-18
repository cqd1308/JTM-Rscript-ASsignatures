options(stringsAsFactors = F)
library(tidyverse)
library(randomForest)
library(doBy)
library(Boruta)
library(pROC)
load('ob.Rdata')
load('clinilnm.Rdata')        #####clinical information concerning lymph node metastasis
clinilnm<-clinilnm[!duplicated(clinilnm$Id),]
rownames(clinilnm)<-clinilnm$Id
group<-c(rep('Normal',59),rep('Tumor',513))
ob<-as.data.frame(t(ob[,group=='Tumor']))
samesample<-intersect(rownames(ob),rownames(clinilnm))
clinical<-clinilnm[samesample,]
ob<-ob[samesample,]
ob$LNM<-clinical$LNM
Id<-colnames(ob)[-ncol(ob)]
colnames(ob)<-str_remove_all(colnames(ob),'\\|')
colnames(ob)<-str_remove_all(colnames(ob),'\\-')
ob$LNM<-as.factor(ob$LNM)
####run Boruta algorithm and extract importance
boruta_train<-Boruta(LNM~.,data=ob,maxRuns=1000)
final_boruta <- TentativeRoughFix(boruta_train)
imp<-attStats(final_boruta)
imptruename1<-imp
imptruename<-imptruename1
rownames(imptruename)<-Id
imp<-imp[imp$decision=='Confirmed',]
imp<-imp[order(imp$medianImp,decreasing = T),]
bor_ob<-ob[,-ncol(ob)]
bor_ob<-bor_ob[,rownames(imp)]

bor_ob$LNM<-ob$LNM

###this is an adaption of function rfcv of randomForest package, which sequentially reduced the number of splicing events according to the importance of Boruta analysis 
borutarfcv<-function (trainx, trainy, cv.fold = 5, scale = "log", step = 0.5, 
                      mtry = function(p) max(1, floor(sqrt(p))), recursive = FALSE, 
                      ...) 
{
  classRF <- is.factor(trainy)
  n <- nrow(trainx)
  p <- ncol(trainx)
  q<-nrow(imptruename1)
  if (scale == "log") {
    k <- floor(log(p, base = 1/step))
    n.var <- round(p * step^(0:(k - 1)))
    same <- diff(n.var) == 0
    if (any(same)) 
      n.var <- n.var[-which(same)]
    if (!1 %in% n.var) 
      n.var <- c(n.var, 1)
  }
  else {
    n.var <- seq(from = p, to = 1, by = step)
  }
  k <- length(n.var)
  cv.pred <- vector(k, mode = "list")
  for (i in 1:k) cv.pred[[i]] <- trainy
  if (classRF) {
    f <- trainy
  }
  else {
    f <- factor(rep(1:5, length = length(trainy))[order(order(trainy))])
  }
  nlvl <- table(f)
  idx <- numeric(n)
  for (i in 1:length(nlvl)) {
    idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold, 
                                                length = nlvl[i]))
  }
  for (i in 1:cv.fold) {
    all.rf <- randomForest(trainx[idx != i, , drop = FALSE], 
                           trainy[idx != i], trainx[idx == i, , drop = FALSE], 
                           trainy[idx == i], mtry = mtry(p), importance = TRUE, 
                           ...)
    cv.pred[[1]][idx == i] <- all.rf$test$predicted
    impvar<-imptruename1[rownames(imp),]
    impvar<-impvar[order(impvar$medianImp,decreasing = T),]
    for (j in 2:k) {
      imp.idx <- rownames(impvar)[1:n.var[j]]
      sub.rf <- randomForest(trainx[idx != i, imp.idx, 
                                    drop = FALSE], trainy[idx != i], trainx[idx == 
                                                                              i, imp.idx, drop = FALSE], trainy[idx == i], 
                             mtry = mtry(n.var[j]), importance = recursive, 
                             ...)
      cv.pred[[j]][idx == i] <- sub.rf$test$predicted
      if (recursive) {
        impvar <- (1:length(imp.idx))[order(importance(sub.rf, 
                                                       type = 1), decreasing = TRUE)]
      }
      NULL
    }
    NULL
  }
  if (classRF) {
    error.cv <- sapply(cv.pred, function(x) mean(trainy != 
                                                   x))
  }
  else {
    error.cv <- sapply(cv.pred, function(x) mean((trainy - 
                                                    x)^2))
  }
  names(error.cv) <- names(cv.pred) <- n.var
  list(n.var = n.var, error.cv = error.cv, predicted = cv.pred)
}
bor_cross<-replicate(5,borutarfcv(bor_ob[-ncol(bor_ob)],bor_ob$LNM,cv.fold = 5,scale='number',step=-1), simplify = FALSE)
cverror <- data.frame(sapply(bor_cross, '[[', 'error.cv'))
cverror$variables <- rownames(cverror)
# cverror <- melt(cverror, id = 'names')
cverror<-gather(cverror,key='X',value = 'error',X1:X5)
cverror$variables<-as.numeric(as.character(cverror$variables))
errorsum<-summary_by(cverror,error~variables,FUN = mean)
###draw the plot showing cross validation error and number of features
ggplot(errorsum,aes(variables,error.mean)) +geom_point(color='brown')+
  geom_line(color = 'brown') +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +  
  labs(title = '',x = 'Number of features', y = 'Cross validation error')+
  scale_x_continuous(breaks = c(seq(0,20,2)))+
  geom_vline(xintercept = 12,lty=2,color='grey31')+
  ggsave('error.pdf',width=5,height = 4)
#####five-fold cross-validation and roc curves plotting
number<-c(1:502)
group<-c(rep(c('g1','g2','g3','g4','g5'),100),'g1','g2')
set.seed(333)
groupsam<-sample(group,502,rep = F)
groupres<-as.data.frame(cbind(number,groupsam))
groupres$number<-as.numeric(groupres$number)
i<-'g1'
#for (i in c('g1','g2','g3','g4','g5')){
train<-groupres[groupres$groupsam!=i,]
train<-bor_ob[train$number,]
test<-groupres[groupres$groupsam==i,]
test<-bor_ob[test$number,]
rftrain<-randomForest(LNM~.,data =train,ntree=10000)
rttest<-as.data.frame(predict(rftrain,test[,1:(ncol(test)-1)],type='prob'))
rttest$predict <- names(rttest)[1:2][apply(rttest[,1:2], 1, which.max)]
rttest$observed<-test[,ncol(test)]
rocg5 <- roc(rttest$observed,as.numeric(rttest$`1`))

mergelist<-list(Fold1=rocg1,Fold2=rocg2,Fold3=rocg3,Fold4=rocg4,Fold5=rocg5)
ggroc(mergelist,aes='color',legacy.axes = TRUE,lwd=0.7)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey41", linetype="twodash")+
  facet_wrap(.~name,nrow=2)+theme_bw()+
  theme(panel.grid.major = element_line(linetype = 'dashed'),
        panel.grid.minor = element_blank())+scale_x_continuous(breaks = c(0.0,0.2,0.4,0.6,0.8,1.0))+scale_y_continuous(breaks = c(0.0,0.2,0.4,0.6,0.8,1.0))+
  theme(legend.position = 'none')+
  ggsave('crossvalid-lnm.pdf',width=4.8,height = 3.7)