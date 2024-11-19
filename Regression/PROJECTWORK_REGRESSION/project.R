library(readr)
library(lattice)
library(latticeExtra)
library(corrplot)
library(car)
library(MASS)
library(caTools)
library(leaps)
library(glmnet)


#data---------------------------------------------
setwd("G:\\PROJECTWORK_DEEPAYAN")
read_csv("car data.csv")->Data
attach(Data)
dim(Data)
summary(is.na(Data))


#plots--------------------
summary(Selling_Price)
boxplot(Selling_Price,horizontal = T,xlab="Selling_Price")  #has outliers
histogram(Selling_Price,col="red")  #response is positive skewed
boxplot(cbind(sqrt(Selling_Price),Selling_Price^(1/3),log(Selling_Price)),names=c("sq. root","cube root","log"),xlab="boxplot for various transformation")  #no outliers for log

boxplot(Present_Price,horizontal = T,xlab="Present_Price")  
histogram(Present_Price)  #same as selling price

boxplot(Driven_kms,horizontal = T,xlab="Driven_kms")  #positive skew  
histogram(Driven_kms)


plot(table(Year),col="darkgrey") #very few obs available in past or most currernt year compared to intermediate years
plot(table(Fuel_Type),col="darkred")  #most obs from petrol, negligible for CNG
plot(table(Transmission),col="darkred") #very few obs for automatic compared to manual
plot(table(Selling_type),col="blue")  #not enough information for owner>0
plot(table(Owner),col="blue")  #not enough information for owner>0


xyplot(Selling_Price~Present_Price|Car_Name)  #too many categories, not enough observations in each category


xyplot(Selling_Price~Present_Price,type=c("p","r"),col.line="red",lwd=1.5)+
  layer(panel.smoother(x,y,method="loess",se=F,col.line = "black")) #most observations in a small region of plot,due to outliers 
xyplot(log(Selling_Price)~log(Present_Price),type=c("p","r"),col.line="red",lwd=1.5)+
  layer(panel.smoother(x,y,method="loess",se=F,col.line = "black")) #strong linear relationship


xyplot(log(Selling_Price)~Year,type=c("p","r"),col.line="red",lwd=1.5)+
  layer(panel.smoother(x,y,method="loess",se=F,col.line = "black")) #strong linear relationship

xyplot(log(Selling_Price)~log(Driven_kms),type=c("p","r"),col.line="red",lwd=1.5)
+layer(panel.smoother(x,y,method="loess",se=F,col.line = "black")) #little linear relationship
xyplot(log(Selling_Price)~log(Driven_kms),group=Fuel_Type,auto.key = T)+
  layer(panel.abline(lm(log(Selling_Price[Fuel_Type=="CNG"])~log(Driven_kms[Fuel_Type=="CNG"])),col="blue"))+
  layer(panel.abline(lm(log(Selling_Price[Fuel_Type=="Diesel"])~log(Driven_kms[Fuel_Type=="Diesel"])),col="pink",lwd=2))+
  layer(panel.abline(lm(log(Selling_Price[Fuel_Type=="Petrol"])~log(Driven_kms[Fuel_Type=="Petrol"])),col="darkgreen")) 


bwplot(log(Selling_Price)~factor(Fuel_Type),fill="grey") #diesel>petrol
xyplot(log(Selling_Price)~log(Present_Price),groups =Fuel_Type,auto.key = T)
bwplot(log(Selling_Price)~factor(Transmission),fill="grey") #automatic>manual
xyplot(log(Selling_Price)~log(Present_Price),groups =Transmission,auto.key = T)
bwplot(log(Selling_Price)~factor(Selling_type),fill="grey") #dealer>individual
xyplot(log(Selling_Price)~log(Present_Price),groups =Selling_type,auto.key = T)



#collinearity-------------------------
X<-cbind(logpresent=log(Present_Price),Year,logkms=log(Driven_kms),
         model.matrix(~0+Fuel_Type),
         model.matrix(~0+Transmission),
         model.matrix(~0+Selling_type)) 
logselling=log(Selling_Price)
splom(cbind(logselling,X[,1:3]))  #others are dummy variables
#k-1 dummy var. for k levels
corrplot(cor(cbind(logselling,X[,c(-4,-7,-9)])),method="color",addCoef.col = T,diag=F) #strong collinearity between two fuel type as indicator_petrol+indicator_diesel=1,indicator_cng=0 for most obs
corrplot(cor(cbind(logselling,X[,c(-5,-7,-9)]->x)),method="color",addCoef.col = T,diag=F) #change the choice of dummy variables for collinearity
x1=x[,1]
x2=x[,2]
x3=x[,3]
x4=x[,4]
x5=x[,5]
x6=x[,6]
x7=x[,7]


#regression--------------------------
fm<-lm(logselling~x1+x2+x3+x4+x5+x6+x7)
summary(fm)
anova(fm,lm(logselling~1))


xyplot(rstudent(fm)~fm$fitted.values,main="fitted vs residual",grid=T,type="p")+
  layer(panel.smoother(x,y,form=y~x,method = "loess",span=.75,se=T,col = "red",alpha.se = .1))   #little deviation from linearity?
xyplot(abs(rstudent(fm))~fm$fitted.values,main="fitted vs |residual|",grid=T,type="p")+
  layer(panel.smoother(x,y,form=y~x,method = "loess",span=.75,degree=1,se=T,col="red",alpha.se = .1)) #heteroscedastic


qqPlot(rstudent(fm),distribution = "t",df=fm$df.residual-1,envelope = list(border=F),main="Q-Q plot of error",id=F) #heavy tailed?,but symmetric distn
dns<-function(x) dnorm(x)  #df is too high to distinguish between t & normal density
densityplot(rstudent(fm),col.line="red",grid=T,main="densityplot of error")+
  layer(panel.curve(dns,col="grey")) #symmetric & unimodal


Index=1:length(logselling) #a seq of size= sample size
id<-(hatvalues(fm)>3*mean(hatvalues(fm)))
xyplot(hatvalues(fm)~Index,type = "h",col="blue")+
  layer(panel.text(x[id], y[id], labels = Index[id], pos = 3, col = "grey50"))  #two extremely high leverage point , 19th 36th


id<-(cooks.distance(fm)>.025)
xyplot(cooks.distance(fm)~Index,type = "h",col="blue")+
  layer(panel.text(x[id], y[id], labels = Index[id], pos = 3, col = "grey50"))  #one high influence point , 86th


id=c(19,36,86)
xyplot(hatvalues(fm) ~ rstudent(fm), grid = TRUE,col="red")+
  layer(panel.text(x[id], y[id], labels = Index[id], pos = 1, col = "grey50"))  


avPlots(fm,main = "Added-Variable plots",layout=c(2,2)) 
crPlots(fm,terms = ~.-x4-x5-x6-x7)


(sqrt(vif(fm))) # multicollinearity not severe


#correction for heteroscedasticity----------------
xyplot(abs(rstudent(fm))~fm$fitted.values,main="fitted vs |residual|",grid=T,type="p")+
  layer(panel.smoother(x,y,form=y~x,method = "loess",span=.75,se=F,degree=1,col="red"))+
  layer(panel.abline(lm(abs(rstudent(fm))~fm$fitted.values)->aprxsig))
wt= 1/(aprxsig$fitted.values)^2 #sigma approximated by least sq line of "|residual| vs fitted"
wm<-lm(logselling~x1+x2+x3+x4+x5+x6+x7,weights =wt)
summary(wm)
anova(wm,lm(logselling~1))

xyplot(rstudent(wm)~wm$fitted.values,main="fitted vs residual",grid=T,type="p")+
  layer(panel.smoother(x,y,form=y~x,method = "loess",span=.6,se=T,col = "red",alpha.se = .1))   #deviation from linearity?
xyplot(abs(rstudent(wm))~wm$fitted.values,main="fitted vs |residual|",grid=T,type="p")+
  layer(panel.smoother(x,y,form=y~x,method = "loess",span=.8,degree=1,se=T,col="red",alpha.se = .1)) #heteroscedasticity rectified
ncvTest(wm)


qqPlot(rstudent(wm),distribution = "t",df=wm$df.residual-1,envelope = list(border=F),main="Q-Q plot of error",id=F) #heavy tailed?,but symmetric distn
densityplot(rstudent(wm),col.line="red",grid=T,main="densityplot of error")+
  layer(panel.curve(dns,col="black",lty=2)) #symmetric & unimodal
shapiro.test(rstudent(wm))
ks.test(rstudent(wm),pnorm)

id1<-(hatvalues(wm)>0.11)
xyplot(hatvalues(wm)~Index,type = "h",col="blue")+
  layer(panel.text(x[id1], y[id1], labels = Index[id1], pos = 3, col = "grey50"))  #two extremely high leverage point , 19th 36th


id2<-(cooks.distance(wm)>.04)
xyplot(cooks.distance(wm)~Index,type = "h",col="blue")+
  layer(panel.text(x[id2], y[id2], labels = Index[id2], pos = 3, col = "grey50"))  #one high influence point , 86th


id=c(id1,id2)
xyplot(hatvalues(wm) ~ rstudent(wm), grid = TRUE,col="red")+
  layer(panel.text(x[id], y[id], labels = Index[id], pos = 1, col = "grey50"))  

(sqrt(vif(wm))) # multicollinearity not severe


avPlots(wm) 
crPlots(wm,terms = ~.-x4-x5-x6-x7) 

{#cross validation
Yhat=vector()
ST=vector()
for(i in 1:length(x1)){
   y=logselling[-i]
   z1=x1[-i]
   z2=x2[-i]
   z3=x3[-i]
   z4=x4[-i]
   z5=x5[-i]
   z6=x6[-i]
   z7=x7[-i]
   w=wt[-i]
   l=lm(y~z1+z2+z3+z4+z5+z6+z7,weights = w)
   d=data.frame(z1=x1[i],z2=x2[i],z3=x3[i],z4=x4[i],z5=x5[i],z6=x6[i],z7=x7[i])
   Yhat[i]= predict(l,newdata=d)
   ST[i]= mean(y) }
SE=sum((logselling-Yhat)^2)
ST=sum((logselling-ST)^2)
(predictedRsq= 1- SE/ST)
xyplot(Yhat~logselling,grid=T,xlab = "actual log(Selling_Price)",ylab="predicted log(Selling_Price)",main="leave-one-out cross validation")+layer(panel.abline(a=0,b=1,col="red",lty=2)) 
}



#trying to explain the most unusual observations------------- 
Data[c(19,36,86),]
Data[which(Fuel_Type=="CNG"),]
Data[which(Fuel_Type=="Petrol" & Selling_type=="Individual" & Transmission=="Automatic"),]


xyplot((covratio(wm))~Index,grid=T)+
  layer(panel.abline(a=1,b=0,col="red"))+
  layer(panel.text(x[id], y[id], labels = Index[id], pos = 3, col = "grey50")) #removing 19th & 36th will improve precision of estimation a lot






#trying to get sparser model-------------------

#splitting
{
  set.seed(101) #in case we need to reproduce calculations
  sample <- sample.split(logselling, SplitRatio = 0.8)
  training_dataset<-as.data.frame(subset(cbind(logselling,x), sample == TRUE))  
  testing_dataset<-as.data.frame(subset(cbind(logselling,x), sample == FALSE))
  colnames(training_dataset)=c("y","x1","x2","x3","x4","x5","x6","x7")
  colnames(testing_dataset)=c("y","x1","x2","x3","x4","x5","x6","x7")
  wt2=wt[sample] 
  
  
  nm0<-lm(y~x1+x2+x3+x4+x5+x6+x7,data=training_dataset,weights = wt2)
  summary(nm0)
  P=predict(nm0,newdata = testing_dataset) 
  A= testing_dataset[,1] 
  (MAPE= 100*mean(abs((P-A)/A)))  }

#LASSO
{
  nm3<-with(training_dataset,glmnet(cbind(x1, x2, x3,x4,x5,x6,x7), y, alpha = 1))
  plot(nm3, xvar = "lambda", label = TRUE) 
  plot(nm3, xvar = "dev", label = TRUE)
  nmlamda <- with(training_dataset,cv.glmnet(cbind(x1, x2, x3,x4,x5,x6,x7), y, alpha = 1,nfolds = 50))
  lambda.min = nmlamda$lambda.min
  lambda.1se = nmlamda$lambda.1se
  plot(nmlamda) 
  mtext(paste(round(log(lambda.min),2)),side=2,line=-4,cex=.9,col="darkgrey")
  mtext(paste(round(log(lambda.1se),2)),side=2,line=-10,cex=.9,col="darkgrey")
  
  
  nmlasso<-with(training_dataset,glmnet(cbind(x1, x2, x3,x4,x5,x6,x7), y, alpha = 1,lambda =lambda.1se ))
  coef(nmlasso)
  as.logical(coef(nmlasso))}

#best subset selection
{
  reg.all <-regsubsets(y ~x1+x2+x3+x4+x5+x6+x7 ,data = training_dataset,weights = wt2)
  summary(reg.all)
  
    
  xyplot(bic ~ seq_along(bic), data = summary(reg.all),ylab="BIC",xlab="number of predictor",main="best subset selection", grid = TRUE, type = "o", pch = 16)  
  with(summary(reg.all), {
    w <- which; is.na(w) <- (w == FALSE); wbic <- w * bic
    levelplot(wbic, xlim = as.character(round(bic)),
              scales = list(x = list(rot = 90)) ,col.regions=gray(50:90/100),xlab = "BIC", ylab = "terms in model",main="best subset selection")
  })
  #best subset contains 5predictor x1,x2,x3,x5,x7
  #without losing much bic, we can reduce number of predictors to 2, namely x1,x2
}

#five predictor
{
  nm1<-lm(y~x1+x2+x3+x5+x7,data=training_dataset,weights = wt2)
  summary(nm1)
  P=predict(nm1,newdata = testing_dataset) 
  (MAPE= 100*mean(abs((P-A)/A)))
  xyplot(P~A,xlab="actual",ylab="predicted",grid=T,main="Test_Data")+
    layer(panel.abline(a=0,b=1,col="red",lty=2))
  xyplot(training_dataset$y~nm1$fitted.values,grid=T,xlab="actual",ylab="fitted",main="Train_Data")+
    layer(panel.abline(a=0,b=1,col="red",lty=2))  }

#two predictor
{
  nm2<-lm(y~x1+x2,data=training_dataset,weights = wt2)
  summary(nm2)
  P=predict(nm2,newdata = testing_dataset) 
  (MAPE= 100*mean(abs((P-A)/A)))
  
  xyplot(P~A,xlab="actual",ylab="predicted",grid=T,main="Test_Data")+
    layer(panel.abline(a=0,b=1,col="red",lty=2))
  xyplot(training_dataset$y~nm2$fitted.values,grid=T,xlab="actual",ylab="fitted",main="Train_Data")+
    layer(panel.abline(a=0,b=1,col="red",lty=2))
  
  xyplot(rstudent(nm2)~nm2$fitted.values,main="fitted vs residual",grid=T,type="p")+
    layer(panel.smoother(x,y,form=y~x,method = "loess",span=.75,se=T,col = "red",alpha.se = .1))   
  xyplot(abs(rstudent(nm2))~nm2$fitted.values,main="fitted vs |residual|",grid=T,type="p")+
    layer(panel.smoother(x,y,form=y~x,method = "loess",span=.75,degree=1,se=T,col="red",alpha.se = .1)) #homoscedastic
  qqPlot(rstudent(nm2),distribution = "t",df=nm2$df.residual-1,envelope = list(border=F),main="Q-Q plot of error",id=F) #heavy tailed? negative skewed? 
  densityplot(rstudent(nm2),col.line="red",grid=T,main="densityplot of error")+
    layer(panel.curve(dns,lty=2)) #symmetric & unimodal
  
  avPlots(nm2,id=F)
  crPlots(nm2)   
  
  sum((nm2$residuals)^2)/nm2$df.residual #SSE/df
  H=model.matrix(nm2)
  round(solve(t(H)%*%H),3) # M=(X'X)^-1
}

#ROBUST REGRESSION---------------
{
  rm<-rlm(y~x1+x2,data=training_dataset,psi=psi.bisquare,weights = wt2)
  summary(rm) 
  P=predict(rm,newdata = testing_dataset) 
  (MAPE= 100*mean(abs((P-A)/A)))
  xyplot(P~A,xlab="actual",ylab="predicted",grid=T,main="Test_Data")+
    layer(panel.abline(a=0,b=1,col="red",lty=2))
  xyplot(training_dataset$y~rm$fitted.values,grid=T,xlab="actual",ylab="fitted",main="Train_Data")+
    layer(panel.abline(a=0,b=1,col="red",lty=2))  }






#APPENDIX---------------------
random=sample(10000,500)
coeff=matrix(nrow=500,ncol=7)
for(i in 1:500){ #checking which variables are selected by lasso most of the time
  set.seed(random[i])
  
  {
    sample <- sample.split(logselling, SplitRatio = 0.8)
    training_dataset<-as.data.frame(subset(cbind(logselling,x), sample == TRUE))  
    testing_dataset<-as.data.frame(subset(cbind(logselling,x), sample == FALSE))
    colnames(training_dataset)=c("y","x1","x2","x3","x4","x5","x6","x7")
    colnames(testing_dataset)=c("y","x1","x2","x3","x4","x5","x6","x7")
    wt2=wt[sample] } #splitting
  
  
  {
    lambda.1se = (with(training_dataset,cv.glmnet(cbind(x1, x2, x3,x4,x5,x6,x7), y, alpha = 1,nfolds = 50)))$lambda.1se
    nmlasso<-with(training_dataset,glmnet(cbind(x1, x2, x3,x4,x5,x6,x7), y, alpha = 1,lambda =lambda.1se ))
    coeff[i,]=coef(nmlasso)[-1]  }  #LASSO
}
boxplot(coeff,col="blue",horizontal=F,names=c("b1","b2","b3","b4","b5","b6","b7"))
abline(h=0,col="red",lty=2)





Rsq5=vector()
MAPE5=vector()
Rsq2=vector()
MAPE2=vector()
for(i in 1:500){
  set.seed(random[i])
  {
    sample <- sample.split(logselling, SplitRatio = 0.8)
    training_dataset<-as.data.frame(subset(cbind(logselling,x), sample == TRUE))  
    testing_dataset<-as.data.frame(subset(cbind(logselling,x), sample == FALSE))
    colnames(training_dataset)=c("y","x1","x2","x3","x4","x5","x6","x7")
    colnames(testing_dataset)=c("y","x1","x2","x3","x4","x5","x6","x7")
    wt2=wt[sample] } # same splitting
  
  
  {
    nm0<-lm(y~x1+x2+x3+x5+x7,data=training_dataset,weights = wt2)
    Rsq5[i]=summary(nm0)$r.sq
    P=predict(nm0,newdata = testing_dataset) 
    A= testing_dataset[,1] 
    MAPE5[i]= 100*mean(abs((P-A)/A)) } #five parameters


  {
    nm1<-lm(y~x1+x2,data=training_dataset,weights = wt2)
    Rsq2[i]=summary(nm1)$r.sq
    P=predict(nm1,newdata = testing_dataset) 
    MAPE2[i]= 100*mean(abs((P-A)/A)) } #two parameters
  
}
xyplot(Rsq2~Rsq5,grid=T,xlab="Rsq5",ylab="Rsq2",main="simulation of R square")+layer(panel.abline(a=0,b=1,col="red"))
xyplot(MAPE2~MAPE5,grid=T,xlab="Rsq5",ylab="Rsq2",main="simulation of MAPE")+layer(panel.abline(a=0,b=1,col="red"))
densityplot(~MAPE2+MAPE5,grid=T,auto.key = T,xlab="MAPE",main="simulation of MAPE")
