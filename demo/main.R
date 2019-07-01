rm(list = ls())
setwd("/media/zheqng/Seagate Backup Plus Drive/zheqng@nwu/src/RJMCMC-my-C/simu1/canadaweather/demo/")
# setwd('/home/zheqng/src/RJMCMC-my-C/canadaweather/demo/')
library('matrixStats')
# library('GPFDA')
library('MASS')
library('mvtnorm')
source('predict.r')
source('readsts.R')

Curve.num=8
testdata$X = 1:73
# predict.gp(x,y,x.new,thet,k)
# plot.gp.predict(traindata,testdata,predictdata,row.no)
# z=z.ave(traindata,theta,N,1,lag = 100)
# plot.z(z)

# source('bayss analysis.R')
# source('convergence diagnose.R')

# RMSE.m = rep(0,Curve.num)
# corr.m = rep(0,Curve.num)
op <- par(mfrow=c(3,3))
for(m in 1:Curve.num){
  # i=N;
  # k = z[i,m]+1
  predict.sum<-predict.gp.ave(traindata,testdata,theta,N,m,lag = 100)
  # RMSE.m[m]<- rmse(testdata$Y[m,],predict.sum)
  # corr.m[m] <- predictdata[[m]]$correlation
  plot.gp.predict(traindata,testdata,predict.sum,m)
  title(paste(N.paper,"th iter ","predict",i=m,"curve"))
  print(m)
}
op<- par(mfrow = c(1,1))
# Nm = length(testdata$Y[m,])
# RMSE = mean(RMSE.m)
# RMSE.each = rep(0,3)
# RMSE.each[1] = mean(RMSE.m[1:3])
# RMSE.each[2] = mean(RMSE.m[4:6])
# RMSE.each[3] = mean(RMSE.m[7:9])
z.iter=z.ave(traindata,theta,N,lag = 200)
plot.z(z.iter)
plot.iter.label(tmp,z.iter)
title('curve clustering')
############################classification
total.index = 1:35
remain.index = total.index[-choose.curve]

remaindata = list("temp","X","Y","curve.num")
# remaindata$temp =read.table("CanadianWeather_temperature.dat")
remaindata$curve.num = length(remain.index)
remaindata$X=data.frame(t(matrix(rep(1:73,27),73,27)))
tmp1= tmp - t(matrix( rep(apply(tmp,2,mean),73),35,73))

remaindata$Y=data.frame(t(tmp1[,remain.index]))
z.remain=calc.z(remaindata,theta[[N]])
z.total =rep(0,35)
z.total[choose.curve]=z.iter[51,]
z.total[remain.index] = z.remain

totaldata = list("temp","X","Y","curve.num")
# remaindata$temp =read.table("CanadianWeather_temperature.dat")
totaldata$curve.num = length(total.index)
totaldata$X=data.frame(t(matrix(rep(1:73,35),73,35)))
totaldata$Y=data.frame(t(tmp[,total.index]))
z=c(rep(1,1,15),rep(2,1,9),rep(3,1,5),rep(2,1,3),rep(4,1,3))
plot.label(totaldata,z.total)
plot.label(totaldata,z)
for (ii in 6:15)
  lines(1:73,totaldata$Y[ii,],col = "red",type = "o",pch=20,cex=0.5)
title('expert classification')
source('bayss analysis.R')
source('convergence diagnose.R')
