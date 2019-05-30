rm(list = ls())
# setwd("/media/zheqng/Seagate Backup Plus Drive/zheqng@nwu/src/RJMCMC-my-C/simu1/fix-k-move/demo/")
setwd('/home/zheqng/src/RJMCMC-my-C/canadaweather/demo/')
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
  predict.sum<-predict.gp.ave(traindata,testdata,theta,N,m,lag = 10)
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
z=z.ave(traindata,theta,N,lag = 100)
plot.z(z)
color_label = c("red","blue","yellow","green")
plot(-100,-100,xlim=c(0,73),ylim=range(tmp))
choose.curve = c(1,2,16,17,25,26,33,34)
jj=1
for (ii in choose.curve){
  lines(1:73,tmp[,ii],col = color_label[z[52,jj]],type = "l")
  jj=jj+1
}
source('bayss analysis.R')
source('convergence diagnose.R')
