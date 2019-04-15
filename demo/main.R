rm(list = ls())
# setwd("/media/zheqng/Seagate Backup Plus Drive/zheqng@nwu/src/RJMCMC-my-C/simu1/fix-k-move/demo/")
setwd('/home/zheqng/src/RJMCMC-my-C/fix-k-move/demo/')
library('matrixStats')
# library('GPFDA')
library('MASS')
library('mvtnorm')
source('predict.r')
source('readsts.R')
rmse.m = rep(0,9)
corr.m = rep(0,9)
op <- par(mfrow=c(3,3))
for(m in 1:9){
  i=N;
  k = z[i,m]+1
    predictdata<-predict.gp(traindata$X[m,],traindata$Y[m,],testdata$X[m,],testdata$Y[m,],theta[[i]],k)
    rmse.m[m]<- rmse(testdata$Y[m,],predictdata)
    corr.m[m] <- predictdata$correlation
    # plot.gp.i.predict(traindata,predictdata,m)
    plot.gp.predict(traindata,testdata,predictdata,m)
    title(paste(N,"th iter ","predict",i=m,"curve"))
    # plot.gp.predict(traindata,testdata,predictdata,1)
    # plot.gp(traindata ,testdata )
}
op<- par(mfrow = c(1,1))
Nm = length(testdata$Y[m,])
rmse = mean(rmse.m)
rmse.each = rep(0,3)
rmse.each[1] = mean(rmse.m[1:3])
rmse.each[2] = mean(rmse.m[4:6])
rmse.each[3] = mean(rmse.m[7:9])

corr.each = rep(0,3)
corr.each[1] = mean(corr.m[1:3])
corr.each[2] = mean(corr.m[4:6])
corr.each[3] = mean(corr.m[7:9])
log.P =   print.posterior(traindata,theta[[N]])

norm = 0.0
for( m in 1:traindata$curve.num){
  norm = norm + logSumExp(log.P[m,]+ log(theta[[N]]$pi))
}

rmse(testdata$Y,predictdata$pred.mean)
