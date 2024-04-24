####libs------------------------------------------
library(latticeExtra)
library(gridExtra)
library(heplots)
library(corrplot)
library(mvnormtest)
library(MASS)
library(clinfun)
library(caTools)
library(VGAM)





####DATA------------------------------------------
setwd("G:\\PROJECTWORK_MULTIVARIATE")
data = read.table("Egyptian-skulls.txt")
colnames(data) = c("mb","bh","bl","nh","period")
data$period = as.factor(data$period)
levels(data$period)=c("4000BC","3300BC","1850BC") # 1,2,3 respectively
attach(data)
p1 = data[period=="4000BC",-5] #period1
p2 = data[period=="3300BC",-5] #period2
p3 = data[period=="1850BC",-5] #period3
xbar1 = colMeans(p1)
xbar2 = colMeans(p2)
xbar3 = colMeans(p3)
S1=cov(p1)
S2=cov(p2)
S3=cov(p3)



####ANALYSIS---------------------------
###boxplots.........................
par(mfrow=c(2,2))
boxplot(cbind(p1[[1]],p2[[1]],p3[[1]]),names=c("4000BC","3300BC","1850BC"),horizontal=T,col=c("red","green","blue"),xlab="maximum breadth(mm.)")
  #for period3 mb is higher on average than other periods and it is -ve skewed ,for period2 there is an outlier
boxplot(cbind(p1[[2]],p2[[2]],p3[[2]]),names=c("4000BC","3300BC","1850BC"),horizontal=T,col=c("red","green","blue"),xlab="basibregmatic height(mm.)")
  #not very much difference between periods
boxplot(cbind(p1[[3]],p2[[3]],p3[[3]]),names=c("4000BC","3300BC","1850BC"),horizontal=T,col=c("red","green","blue"),xlab="basialveolar length(mm.)")
 #average length decreasing over periods
boxplot(cbind(p1[[4]],p2[[4]],p3[[4]]),names=c("4000BC","3300BC","1850BC"),horizontal=T,col=c("red","green","blue"),xlab="nasal height(mm.)")
 #not much difference between periods
par(mfrow=c(1,1))


###scatterplot matrix...................
splom(p1,main="4000BC",varnames=c("maximum\nbreadth","basibregmatic\nheight","basialveolar\nlength","nasal\nheight"))
splom(p2,main="3300BC",varnames=c("maximum\nbreadth","basibregmatic\nheight","basialveolar\nlength","nasal\nheight"))
splom(p3,main="1850BC",varnames=c("maximum\nbreadth","basibregmatic\nheight","basialveolar\nlength","nasal\nheight"))
  #the variables appear not to be much correlated

###univariate normality......................
dns<-function(x) dnorm(x)
plt1=densityplot(~scale(p1[[1]])+scale(p2[[1]])+scale(p3[[1]]),grid=T,col=c("red","green","blue"),
            key=list(space="top",columns=3,lines=list(col=c("red","green","blue")),text=list(c("4000BC","3300BC","1850BC"))),xlab="scaled maximum breadth")+
    layer(panel.curve(dns,lty=2)) #not very much evidence against normality ,for period 1850BC density curve is -ve skewed
shapiro.test(p1[[1]])
shapiro.test(p2[[1]]) #p value 0.04
shapiro.test(p3[[1]]) #p value 0.03
plt2=densityplot(~scale(p1[[2]])+scale(p2[[2]])+scale(p3[[2]]),grid=T,col=c("red","green","blue"),
            key=list(space="top",columns=3,lines=list(col=c("red","green","blue")),text=list(c("4000BC","3300BC","1850BC"))),xlab="scaled basibregmatic height")+
  layer(panel.curve(dns,lty=2)) #not very much evidence against normality
shapiro.test(p1[[2]])
shapiro.test(p2[[2]]) 
shapiro.test(p3[[2]])
plt3=densityplot(~scale(p1[[3]])+scale(p2[[3]])+scale(p3[[3]]),grid=T,col=c("red","green","blue"),
            key=list(space="top",columns=3,lines=list(col=c("red","green","blue")),text=list(c("4000BC","3300BC","1850BC"))),xlab="scaled  basialveolar length")+
  layer(panel.curve(dns,lty=2)) #not very much evidence against normality
shapiro.test(p1[[3]])
shapiro.test(p2[[3]])
shapiro.test(p3[[3]])
plt4=densityplot(~scale(p1[[4]])+scale(p2[[4]])+scale(p3[[4]]),grid=T,col=c("red","green","blue"),
            key=list(space="top",columns=3,lines=list(col=c("red","green","blue")),text=list(c("4000BC","3300BC","1850BC"))),xlab="scaled  nasal height")+
  layer(panel.curve(dns,lty=2)) #not very much evidence against normality
shapiro.test(p1[[4]])
shapiro.test(p2[[4]])
shapiro.test(p3[[4]])
grid.arrange(plt1,plt2,plt3,plt4,nrow=2)
pval=matrix(c((shapiro.test(p1[[1]]))$p.value,(shapiro.test(p2[[1]]))$p.value,(shapiro.test(p3[[1]]))$p.value,
              (shapiro.test(p1[[2]]))$p.value,(shapiro.test(p2[[2]]))$p.value,(shapiro.test(p3[[2]]))$p.value,
              (shapiro.test(p1[[3]]))$p.value,(shapiro.test(p2[[3]]))$p.value,(shapiro.test(p3[[3]]))$p.value,
              (shapiro.test(p1[[4]]))$p.value,(shapiro.test(p2[[4]]))$p.value,(shapiro.test(p3[[4]]))$p.value),nrow=3)
rownames(pval)= c("4000BC","3300BC","1850BC")
colnames(pval)= c("mb","bh","bl","nh")
corrplot(pval,is.corr = F,addCoef.col = T,method="circle",col.lim=c(0,1),col=rainbow(100))




###Checking multivariate normality......................................... 
mshapiro.test(t(p1)) #Multivariate Shapiro Wilks ,pvalue 0.05208
mshapiro.test(t(p2)) #pvalue 0.01798
mshapiro.test(t(p3)) #pvalue 0.01038


###detecting outliers - using Mahalonobis distance...........................
md_p1 = mahalanobis(p1, center = xbar1, cov = S1)
md_p2 = mahalanobis(p2, center = xbar2, cov = S2)
md_p3 = mahalanobis(p3, center = xbar3, cov = S3)

chisq4<-function(p) qchisq(p, df = 4)
dns2<-function(x) dchisq(x,df=4)
cut_off = qchisq(0.05, df = 4,lower.tail = F)
cut_off2 = qchisq(0.10, df = 4,lower.tail = F)
qqPlot(md_p1,"chisq",df=4,id=F,ylab="Mahalonobis distances for period1",envelope = list(col="red",alpha=0),col.lines = "red") 
mtext(paste("KS test : p value = ",round((ks.test(md_p1,pchisq,df=4))$p.value,4)),line=-5,col="grey50")
qqPlot(md_p2,"chisq",df=4,id=F,ylab="Mahalonobis distances for period2",envelope = list(col="green",alpha=0),col.lines = "green")
mtext(paste("KS test : p value = ",round((ks.test(md_p2,pchisq,df=4))$p.value,4)),line=-5,col="grey50")
qqPlot(md_p3,"chisq",df=4,id=F,ylab="Mahalonobis distances for period3",envelope = list(col="blue",alpha=0),col.lines = "blue")
mtext(paste("KS test : p value = ",round((ks.test(md_p3,pchisq,df=4))$p.value,4)),line=-5,col="grey50")
out1=(md_p1>cut_off)
out2=(md_p2>cut_off)
out3=(md_p3>cut_off)
densityplot(~md_p1 + md_p2 + md_p3,grid=T,col=c("red","green","blue"),
            key=list(space="top",lines=list(col=c("red","green","blue")),text=list(c("4000BC","3300BC","1850BC"))),xlab="Mahalonobis distances")+
  layer(panel.curve(dns2,lty=2))+
  layer(panel.abline(v=cut_off,lty = 3,col="grey50"))+
  layer(panel.abline(v=cut_off2,lty = 3,col="grey50"))+
  layer(panel.text(x=8,y=.12,"level 10%"))+
  layer(panel.text(x=10,y=.08,"level 5%"))+
  layer(panel.text(md_p1[out1], y=0.001, labels = (1:30)[out1], pos = 3, col = "red"))+
  layer(panel.text(md_p2[out2], y=0.001, labels = (1:30)[out2], pos = 1, col = "green"))+
  layer(panel.text(md_p3[out3], y=0.001, labels = (1:30)[out3], pos = 4, col = "blue"))
  #Comment: not much outliers
    #compared to other observations in period2 , obs.4th is very much different ; it was also indicated by boxplots of mb



###Homogeneity of Covariance Matrices.............................
boxM(cbind(mb,bh,bl,nh)~period, data = data)
  #hypothesis of homogeneity accepted
par(mfrow=c(2,2))
corrplot(S1,is.corr = F,main="4000BC",addCoef.col = T,method="color",col.lim=c(-4,40),col=rainbow(100),type="lower")
corrplot(S2,is.corr = F,main="3300BC",addCoef.col = T,method="color",col.lim=c(-4,40),col=rainbow(100),type="lower")  #cov(bl,mb) & cov(bh,nh) quite high here
corrplot(S3,is.corr = F,main="1850BC",addCoef.col = T,method="color",col.lim=c(-4,40),col=rainbow(100),type="lower")
pooledcov=(S1+S2+S3)/3
corrplot(pooledcov,is.corr = F,main="pooled",addCoef.col = T,method="color",col.lim=c(-4,40),col=rainbow(100),type="lower")
par(mfrow=c(1,1))





####MANOVA-------------------------------------------------
#already checked the assumption of normal within group distribution , homogeneity of covariance matrix
  #interested to know if there is any difference from the skulls of three periods
summary(manova(cbind(mb,bh,bl,nh)~period, data = data))  
  #F value is significant at 5% level, implies it is very likely that there is difference between the mean of the periods
par(mfrow=c(2,2))
plot(as.factor(c("4000BC","3300BC","1850BC")),c(xbar1[1],xbar2[1],xbar3[1]),xlab="period",ylab ="maximum breadth(mm.)") 
mtext(paste("ANOVA : p value = ",round(((anova(lm(mb~period)))$`Pr(>F)`)[1],4)),col="grey50")
plot(as.factor(c("4000BC","3300BC","1850BC")),c(xbar1[2],xbar2[2],xbar3[2]),xlab="period",ylab ="basibregmatic height(mm.)") 
mtext(paste("ANOVA : p value = ",round(((anova(lm(bh~period)))$`Pr(>F)`)[1],4)),col="grey50")
plot(as.factor(c("4000BC","3300BC","1850BC")),c(xbar1[3],xbar2[3],xbar3[3]),xlab="period",ylab ="basialveolar length(mm.)") 
mtext(paste("ANOVA : p value = ",round(((anova(lm(bl~period)))$`Pr(>F)`)[1],4)),col="grey50")
plot(as.factor(c("4000BC","3300BC","1850BC")),c(xbar1[4],xbar2[4],xbar3[4]),xlab="period",ylab ="nasal height(mm.)") 
mtext(paste("ANOVA : p value = ",round(((anova(lm(nh~period)))$`Pr(>F)`)[1],4)),col="grey50")
par(mfrow=c(1,1))
###can suspect monotone trend for mb and bl ? .........................
jonckheere.test(bl,rep(1:3,each=30),alternative = "decreasing") #pvalue 0.01182
jonckheere.test(mb,rep(1:3,each=30),alternative = "increasing") #pvalue 0.004274 ,quite evidence in favour of alternative



####PCA---------------------------------------------------
(PCA=prcomp(~ .-period, data = data, scale. = TRUE))
plot(PCA) #80.89% var explained by first three PC, 58.91%only by first two
biplot(PCA)



####SPLITTING -----------------------------------------------
set.seed(72899) #in case we need to reproduce calculations
training_dataset1<-as.data.frame(subset(data[1:30,], (sample.split(1:30, SplitRatio = 0.8)) == T))  
testing_dataset1<-as.data.frame(subset(data[1:30,], (sample.split(1:30, SplitRatio = 0.8)) == F))
set.seed(101) #in case we need to reproduce calculations
training_dataset2<-as.data.frame(subset(data[31:60,], (sample.split(1:30, SplitRatio = 0.8)) == T))  
testing_dataset2<-as.data.frame(subset(data[31:60,], (sample.split(1:30, SplitRatio = 0.8)) == F))
set.seed(70541) #in case we need to reproduce calculations
training_dataset3<-as.data.frame(subset(data[61:90,], (sample.split(1:30, SplitRatio = 0.8)) == T))  
testing_dataset3<-as.data.frame(subset(data[61:90,], (sample.split(1:30, SplitRatio = 0.8)) == F))

data1=rbind(training_dataset1,training_dataset2,training_dataset3)
data2=rbind(testing_dataset1,testing_dataset2,testing_dataset3)




####LDA---------------------------------------------------
  #assumption of homogeneity of cov matrix already satisfied
(LDA=lda(period ~ .-period ,data = data1))
plot(LDA,xlim=c(-4,4) ,col = c("red","green","blue")[period]) #94.03% discrimination done by LD1 only                                                  #also in LD1 , coefficients corresponding to mb and bl are higher than others , as suspected from boxplot comparisons
###based on train set..................
predicted=predict(LDA)$class
(classify=table(data1$period,predicted))
(misclass= 1- sum(diag(classify))/sum(classify)) #misclassification rate 58.33%
rownames(classify)=c("true 4000BC","true 3300BC","true 1850BC")
colnames(classify)=c("pred. 4000BC","pred. 3300BC","pred. 1850BC")
corrplot(classify,is.corr = F,xlab="l",col.lim = c(0,72),method="pie",addCoef.col = T)
###based on test set...................
predicted=predict(LDA,data2)$class
(classify=table(data2$period,predicted))
(misclass= 1- sum(diag(classify))/sum(classify)) #misclassification rate 44.44%
rownames(classify)=c("true 4000BC","true 3300BC","true 1850BC")
colnames(classify)=c("pred. 4000BC","pred. 3300BC","pred. 1850BC")
corrplot(classify,is.corr = F,xlab="l",col.lim = c(0,18),method="pie",addCoef.col = T)
  #Note : for period2 misclassification occurs most frequently 
        #most correct classifications are from period3 
###Leave-One-Out Cross-Validation............
predictedCV = sapply(1:90,
                   function(i) {
                     return(predict(lda(period ~ .-period ,data = data[-i,]),
                                    data[i,])$class)
                   })
(classify=table(data$period,predictedCV))
(misclass= 1- sum(diag(classify))/sum(classify)) #misclassification rate 61.11%
rownames(classify)=c("true 4000BC","true 3300BC","true 1850BC")
colnames(classify)=c("pred. 4000BC","pred. 3300BC","pred. 1850BC")
corrplot(classify,is.corr = F,xlab="l",col.lim = c(0,90),method="pie",addCoef.col = T)





###what if we apply LDA based on period3 only?...................
LDA2=lda(c(rep("not 1850BC",48),rep("1850BC",24)) ~ .-period ,data = data1) 
predicted=predict(LDA2)$class
(classify=table(c(rep("not 1850BC",48),rep("1850BC",24)),predicted))
(misclass= 1- sum(diag(classify))/sum(classify))
rownames(classify)=c("true 1850BC","true not 1850BC")
colnames(classify)=c("pred. 1850BC","pred. not 1850BC")
corrplot(classify,is.corr = F,xlab="l",col.lim = c(0,72),method="pie",addCoef.col = T)  #misclassification 27.78%
predicted=predict(LDA2,data2)$class
(classify=table(c(rep("not 1850BC",12),rep("1850BC",6)),predicted))
(misclass= 1- sum(diag(classify))/sum(classify))
rownames(classify)=c("true 1850BC","true not 1850BC")
colnames(classify)=c("pred. 1850BC","pred. not 1850BC")
corrplot(classify,is.corr = F,xlab="l",col.lim = c(0,18),method="pie",addCoef.col = T)  #misclassification 22.22%
###Leave-One-Out Cross-Validation............
predictedCV2 = sapply(1:90,
                     function(i) {
                       return(predict(lda((c(rep("not 1850BC",60),rep("1850BC",30)))[-i] ~ .-period ,data = data[-i,]),
                                      data[i,])$class)
                     })
(classify=table(c(rep("not 1850BC",60),rep("1850BC",30)),predictedCV2))
(misclass= 1- sum(diag(classify))/sum(classify)) #misclassification rate 28.89%
rownames(classify)=c("true 1850BC","true not 1850BC")
colnames(classify)=c("pred. 1850BC","pred. not 1850BC")
corrplot(classify,is.corr = F,xlab="l",col.lim = c(0,90),method="pie",addCoef.col = T)

 







####QDA-----------------------------------------------------
   #but,assumption of homogeneity of covariance matrix already satisfied ,then no QDA






####MULTINOMIAL LOGISTIC-----------------------------------------
  #here the response ~ multinomial(n=3,-,-,-) and predictors are continuous , so we can try this
xyplot(period ~ mb + bh + bl + nh,data, outer = TRUE, jitter.y = TRUE, xlab ="(mm.)",scales = list(x = "free"),type=c("p","smooth"),col.line="red",lty=2)
  #as found earlier , mb and bl distinguishing periods well

mult_logis <- vglm(period ~ mb + bh + bl + nh , multinomial(refLevel = 3), data = data1) #ref level 1850BC
coef(mult_logis,matrix=T)  #note : coefficients for mb and bl higher than others 
                            # what is Hauck-Donner effect ? (this term appears in summary)
                              #Ans : a Wald test statistic is not monotonely increasing as a function of increasing distance between the parameter estimate and the null value.
###based on train data.............
pred=predict(mult_logis) #give the log(ratio of likelihood) w.r.t. ref category 
pred_multi=vector()
for(i in 1:72){
  if( (pred[i,1]<=0)*(pred[i,2]<=0) ) pred_multi[i]="1850BC"
  if( (pred[i,1]*pred[i,2]<=0) ) pred_multi[i]=(c("4000BC","3300BC"))[sign(pred[i,1]-pred[i,2])]
  if( (pred[i,1]>=0)*(pred[i,2]>=0) ) pred_multi[i]=(c("4000BC","3300BC"))[sign(pred[i,1]-pred[i,2])]
}
(classify=(table(data1$period,pred_multi))[,3:1])
(misclass= 1- sum(diag(classify))/sum(classify)) #misclassification rate 56.94%
rownames(classify)=c("true 4000BC","true 3300BC","true 1850BC")
colnames(classify)=c("pred. 4000BC","pred. 3300BC","pred. 1850BC")
corrplot(classify,is.corr = F,xlab="l",col.lim = c(0,72),method="pie",addCoef.col = T) #same as earlier
###based on test data.............
pred=predict(mult_logis,data2) #give the log(ratio of likelihood) w.r.t. ref category 
pred_multi=vector()
for(i in 1:18){
  if( (pred[i,1]<=0)*(pred[i,2]<=0) ) pred_multi[i]="1850BC"
  if( (pred[i,1]*pred[i,2]<=0) ) pred_multi[i]=(c("4000BC","3300BC"))[sign(pred[i,1]-pred[i,2])]
  if( (pred[i,1]>=0)*(pred[i,2]>=0) ) pred_multi[i]=(c("4000BC","3300BC"))[sign(pred[i,1]-pred[i,2])]
}
(classify=(table(data2$period,pred_multi))[,3:1])
(misclass= 1- sum(diag(classify))/sum(classify)) #misclassification rate 44.44%
rownames(classify)=c("true 4000BC","true 3300BC","true 1850BC")
colnames(classify)=c("pred. 4000BC","pred. 3300BC","pred. 1850BC")
corrplot(classify,is.corr = F,xlab="l",col.lim = c(0,18),method="pie",addCoef.col = T) #same as earlier

###multinomial logistic based  on mb and bl only.....................................
mult_logis2 <- vglm(period ~ mb + bl , multinomial(refLevel = 3), data = data1) #ref level 1850BC
coef(mult_logis2,matrix=T)  
pred=predict(mult_logis2) #give the log(ratio of likelihood) w.r.t. ref category 
pred_multi=vector()
for(i in 1:72){
  if( (pred[i,1]<=0)*(pred[i,2]<=0) ) pred_multi[i]="1850BC"
  if( (pred[i,1]*pred[i,2]<=0) ) pred_multi[i]=(c("4000BC","3300BC"))[sign(pred[i,1]-pred[i,2])]
  if( (pred[i,1]>=0)*(pred[i,2]>=0) ) pred_multi[i]=(c("4000BC","3300BC"))[sign(pred[i,1]-pred[i,2])]
}
(classify=(table(data1$period,pred_multi))[,3:1])
(misclass= 1- sum(diag(classify))/sum(classify)) #misclassification rate 56.94%
rownames(classify)=c("true 4000BC","true 3300BC","true 1850BC")
colnames(classify)=c("pred. 4000BC","pred. 3300BC","pred. 1850BC")
corrplot(classify,is.corr = F,xlab="l",col.lim = c(0,72),method="pie",addCoef.col = T) #same as earlier
pred=predict(mult_logis2,data2) #give the log(ratio of likelihood) w.r.t. ref category 
pred_multi=vector()
for(i in 1:18){
  if( (pred[i,1]<=0)*(pred[i,2]<=0) ) pred_multi[i]="1850BC"
  if( (pred[i,1]*pred[i,2]<=0) ) pred_multi[i]=(c("4000BC","3300BC"))[sign(pred[i,1]-pred[i,2])]
  if( (pred[i,1]>=0)*(pred[i,2]>=0) ) pred_multi[i]=(c("4000BC","3300BC"))[sign(pred[i,1]-pred[i,2])]
}
(classify=(table(data2$period,pred_multi))[,3:1])
(misclass= 1- sum(diag(classify))/sum(classify)) #misclassification rate 50.00%
rownames(classify)=c("true 4000BC","true 3300BC","true 1850BC")
colnames(classify)=c("pred. 4000BC","pred. 3300BC","pred. 1850BC")
corrplot(classify,is.corr = F,xlab="l",col.lim = c(0,72),method="pie",addCoef.col = T) #same as earlier




####WHY SUCH PERFORMANCE??-------------------------------------------------------------------------------------------
plt5=densityplot(~(p1[[1]])+(p2[[1]])+(p3[[1]]),grid=T,col=c("red","green","blue"),
                 key=list(space="top",columns=3,lines=list(col=c("red","green","blue")),text=list(c("4000BC","3300BC","1850BC"))),xlab="maximum breadth(mm.)")
plt6=densityplot(~(p1[[2]])+(p2[[2]])+(p3[[2]]),grid=T,col=c("red","green","blue"),
                 key=list(space="top",columns=3,lines=list(col=c("red","green","blue")),text=list(c("4000BC","3300BC","1850BC"))),xlab="basibregmatic height(mm.)")
plt7=densityplot(~(p1[[3]])+(p2[[3]])+(p3[[3]]),grid=T,col=c("red","green","blue"),
                 key=list(space="top",columns=3,lines=list(col=c("red","green","blue")),text=list(c("4000BC","3300BC","1850BC"))),xlab="basialveolar length(mm.)")
plt8=densityplot(~(p1[[4]])+(p2[[4]])+(p3[[4]]),grid=T,col=c("red","green","blue"),
                 key=list(space="top",columns=3,lines=list(col=c("red","green","blue")),text=list(c("4000BC","3300BC","1850BC"))),xlab="nasal height(mm.)")
grid.arrange(plt5,plt6,plt7,plt8,nrow=2)
  #as we see , the densityplots are too much overlapping
    #only period3 (blue cureve) shows some separation from others , 
      #mainly based on mb and bl (the left two plots)

####HENCE WE CAN'T IMPROVE MUCH FURTHER---------------------------------------------


#### AN INTUITIVE CLASSIFICATION RULE---------------------------------------------
plot(as.factor(c("4000BC","3300BC","1850BC")),c(xbar1[1],xbar2[1],xbar3[1]),xlab="period",ylab ="mb & bl (mm.)",lty=2,ylim=c(95,135))
plot(as.factor(c("4000BC","3300BC","1850BC")),c(xbar1[3],xbar2[3],xbar3[3]),lty=3,add=T) 
mtext("bl",col="grey50",line=-1,side=1)
mtext("mb",col="grey50",line=-1)
  #The values of 'mb - bl' increases on average from period1 to period3
boxplot(cbind(data1[1:24,1]-data1[1:24,3],data1[25:48,1]-data1[25:48,3],data1[49:72,1]-data1[49:72,3]),names=c("4000BC","3300BC","1850BC"),horizontal=T,col=c("red","green","blue"),xlab="(maximum breadth - basialveolar length)(mm.)")
xyplot(jitter(rep(0,72))~jitter(mb-bl),ylim=c(-0.2,0.2),ylab=NULL,xlab="mb - bl (mm.)",col=(c("red","green","blue"))[period],grid=T)+
    layer(panel.abline(v=quantile(data1[[1]]-data1[[3]],c(0.3333,0.6667)),lty=2))+
    layer(panel.text(x=20,y=0.1,"predict 4000BC",col="grey50"))+
    layer(panel.text(x=35,y=-0.1,"predict 3300BC",col="grey50"))+
    layer(panel.text(x=45,y=0.1,"predict 1850BC",col="grey50"))
  # divide the whole range of mb-bl from train data in equal three parts using 33.33%quantile and 66.67%quantile
    # assign category 'period1' ,'period2' ,'period3' accordingly as prediction
predicted=vector()
for(i in 1:72){
  predicted[i]="3300BC"
  if( (data1[i,1]-data1[i,3])<=quantile(data1[[1]]-data1[[3]],0.3333) ) predicted[i]="4000BC"
  if( (data1[i,1]-data1[i,3])>=quantile(data1[[1]]-data1[[3]],0.6667) ) predicted[i]="1850BC"
}
classify=(table(data1$period,predicted))[,3:1]  #performance on train data
misclass= 1- sum(diag(classify))/sum(classify) 
rownames(classify)=c("true 4000BC","true 3300BC","true 1850BC")
colnames(classify)=c("pred. 4000BC","pred. 3300BC","pred. 1850BC")
corrplot(classify,is.corr = F,xlab="l",col.lim = c(0,72),method="pie",addCoef.col = T) 
predicted=vector()
for(i in 1:18){
  predicted[i]="3300BC"
  if( (data2[i,1]-data2[i,3])<=quantile(data1[[1]]-data1[[3]],0.3333) ) predicted[i]="4000BC"
  if( (data2[i,1]-data2[i,3])>=quantile(data1[[1]]-data1[[3]],0.6667) ) predicted[i]="1850BC"
}
classify=(table(data2$period,predicted))[,3:1]  #performance on test data
misclass= 1- sum(diag(classify))/sum(classify) 
rownames(classify)=c("true 4000BC","true 3300BC","true 1850BC")
colnames(classify)=c("pred. 4000BC","pred. 3300BC","pred. 1850BC")
corrplot(classify,is.corr = F,xlab="l",col.lim = c(0,18),method="pie",addCoef.col = T) 
