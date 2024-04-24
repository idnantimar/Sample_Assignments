library(carData)
library(lattice)
library(latticeExtra)
library(gridExtra)
library(splines)


#GIVEN DATA-------------------------------------------------------
data(UN,package = "carData")
D=na.omit(UN[,(colnames(UN)=="pctUrban") | (colnames(UN)=="ppgdp") ])
#remove rows with missing values in our required columns, only 14 rows removed instead of 51 in na.omit(UN)
attach(D)
lppgdp=log(ppgdp)
n=length(lppgdp)
plot(pctUrban~ppgdp,data = D,pch=16,col="blue")
plot(pctUrban~log(ppgdp),data = D,pch=16,col="blue")


#PRESS---------------------------------------------------
PRESS<-function(meth){#name of the method to be applied
  resi_sq=numeric()
  for(i in 1:n){
    y=pctUrban[-i]
    x=lppgdp[-i]
    fit=eval(meth)
    yi_hat=predict(fit,data.frame(x=c(lppgdp[i])))
    resi_sq[i]=(pctUrban[i]-yi_hat)^2 #compute square of difference between actual value and new predicted value
  }
  return(sum(resi_sq,na.rm = T))  #PRESS Statistic will be returned
}


#LOESS with varying degrees and spans-----------------------------------------
span=seq(0.1,0.85,0.15)
deg=0:2
for(d in deg){
  plt=list()
  for(j in 1:length(span)){
    u=quote(loess(y~x,span=span[j],degree=d,family="gaussian")) #loess formula on which PRESS will be calculated
    press_stat=PRESS(u)
    plt[[j]]<-xyplot(pctUrban ~ log(ppgdp),col.line = "red",lwd=1.75,type = c("p", "smooth"), degree = d, span=span[j],family = "gaussian",main=paste("span=",span[j],"| PRESS=",round(press_stat,3)))
  }
  grid.arrange(grobs=plt,nrow=2,ncol=3,top=paste("degree=",d),bottom="LOWESS")
}


#POLYNOMIAL REGRESSION with varying degrees-----------------------------------
deg=1:4
press_stat=numeric()
for(j in deg){
  u=quote(lm(y~poly(x,degree=j))) #polynomial on which PRESS will be calculated
  press_stat[j]=PRESS(u)
}
p1=xyplot(pctUrban~log(ppgdp),type="p",main=paste("Linear","|PRESS=",round(press_stat[1],3)))+layer(panel.smoother(x,y,method="lm",se=F,form=y~poly(x,1),col.line = "red",lwd=1.5))
p2=xyplot(pctUrban~log(ppgdp),type="p",main=paste("Quadratic","|PRESS=",round(press_stat[2],3)))+layer(panel.smoother(x,y,method="lm",se=F,form=y~poly(x,2),col.line = "red",lwd=1.5))
p3=xyplot(pctUrban~log(ppgdp),type="p",main=paste("Cubic","|PRESS=",round(press_stat[3],3)))+layer(panel.smoother(x,y,method="lm",se=F,form=y~poly(x,3),col.line = "red",lwd=1.5))
p4=xyplot(pctUrban~log(ppgdp),type="p",main=paste("4th degree","|PRESS=",round(press_stat[4],3)))+layer(panel.smoother(x,y,method="lm",se=F,form=y~poly(x,4),col.line = "red",lwd=1.5))
grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2,bottom="POLYNOMIAL REGRESSION")



#BASIS SPLINES(cubic)-------------------------------------------------
df=c(5,10,20,30)
press_stat=numeric()
for(j in 1:length(df)){
  u=quote(lm(y~1+bs(x,df=df[j])))
  press_stat[j]=PRESS(u)
}
p1=xyplot(pctUrban~log(ppgdp),type="p",main=paste("df=3","|PRESS=",round(press_stat[1],3)))+layer(panel.smoother(x,y,method="lm",se=F,form=y~bs(x,3),col.line = "red",lwd=1.5))
p2=xyplot(pctUrban~log(ppgdp),type="p",main=paste("df=4","|PRESS=",round(press_stat[2],3)))+layer(panel.smoother(x,y,method="lm",se=F,form=y~bs(x,4),col.line = "red",lwd=1.5))
p3=xyplot(pctUrban~log(ppgdp),type="p",main=paste("df=5","|PRESS=",round(press_stat[3],3)))+layer(panel.smoother(x,y,method="lm",se=F,form=y~bs(x,5),col.line = "red",lwd=1.5))
p4=xyplot(pctUrban~log(ppgdp),type="p",main=paste("df=6","|PRESS=",round(press_stat[4],3)))+layer(panel.smoother(x,y,method="lm",se=F,form=y~bs(x,6),col.line = "red",lwd=1.5))
grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2,bottom="BASIS SPLINES (cubic)")


#SMOOTHING SPLINES----------------------------------------------
press_stat.ss=numeric()
for(j in 1:length(df)){
  resi_sq=numeric()
  for(i in 1:n){ 
    y=pctUrban[-i]
    x=lppgdp[-i]
    ss.fit=smooth.spline(y~x,df=df[j])
    yi_hat=as.numeric(predict(ss.fit,data.frame(x=c(lppgdp[i])))$y)
    resi_sq[i]=(pctUrban[i]-yi_hat)^2
  }
  press_stat.ss[j]=sum(resi_sq,na.rm=T)
}
p1=xyplot(pctUrban~log(ppgdp),type="p",main=paste("df=3","|PRESS=",round(press_stat.ss[1],3)))+layer(panel.lines(predict(smooth.spline(y=pctUrban,x=lppgdp,df=3)),col.line = "red",lwd=1.5))
p2=xyplot(pctUrban~log(ppgdp),type="p",main=paste("df=4","|PRESS=",round(press_stat.ss[2],3)))+layer(panel.lines(predict(smooth.spline(y=pctUrban,x=lppgdp,df=4)),col.line = "red",lwd=1.5))
p3=xyplot(pctUrban~log(ppgdp),type="p",main=paste("df=5","|PRESS=",round(press_stat.ss[3],3)))+layer(panel.lines(predict(smooth.spline(y=pctUrban,x=lppgdp,df=5)),col.line = "red",lwd=1.5))
p4=xyplot(pctUrban~log(ppgdp),type="p",main=paste("df=6","|PRESS=",round(press_stat.ss[4],3)))+layer(panel.lines(predict(smooth.spline(y=pctUrban,x=lppgdp,df=6)),col.line = "red",lwd=1.5))
grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2,bottom="SMOOTHING SPLINES")

