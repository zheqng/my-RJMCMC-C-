rm(list = ls())
# setwd("/media/zheqng/Seagate Backup Plus Drive/zheqng@nwu/src/RJMCMC-my-C/simu1/fix-k-move/demo/")
setwd('/home/zheqng/src/RJMCMC-my-C/fix-k-move/demo/')
library('matrixStats')
# library('GPFDA')
library('MASS')
library('mvtnorm')
source('predict.r')
source('readsts.R')
N.paper=20000
RMSE.m = rep(0,9)
corr.m = rep(0,9)
op <- par(mfrow=c(3,3))
for(m in 1:9){
  # i=N;
  # k = z[i,m]+1
  predict.sum<-predict.gp.ave(traindata,testdata,theta,N.paper,m,lag = 100)
  RMSE.m[m]<- rmse(testdata$Y[m,],predict.sum)
  # corr.m[m] <- predictdata[[m]]$correlation
  plot.gp.predict(traindata,testdata,predict.sum,m)
  title(paste(N.paper,"th iter ","predict",i=m,"curve"))
  print(m)
}
op<- par(mfrow = c(1,1))
Nm = length(testdata$Y[m,])
RMSE = mean(RMSE.m)
RMSE.each = rep(0,3)
RMSE.each[1] = mean(RMSE.m[1:3])
RMSE.each[2] = mean(RMSE.m[4:6])
RMSE.each[3] = mean(RMSE.m[7:9])
