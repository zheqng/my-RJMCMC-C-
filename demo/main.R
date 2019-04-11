rm(list = ls())
setwd("/media/zheqng/Seagate Backup Plus Drive/zheqng@nwu/src/RJMCMC-my-C/simu1/fix-k-move/demo/")
# setwd('/home/zheqng/src/RJMCMC-my-C/fix-k-move/demo/')
library('matrixStats')
# library('GPFDA')
library('MASS')
library('mvtnorm')
source('predict.r')
source('readsts.R')
op <- par(mfrow=c(3,3))
for(m in 1:9){
  i=N;
  k = z[i,m]+1
    predictdata<-predict.gp(traindata$X[m,],traindata$Y[m,],testdata$X[m,],theta[[i]],k)
    plot.gp.i.predict(traindata,predictdata,m)
    title(paste(N,"th iter ","predict",i=m,"curve"))
    # plot.gp.predict(traindata,testdata,predictdata,1)
    # plot.gp(traindata ,testdata )
}
op<- par(mfrow = c(1,1))
log.P =   print.posterior(traindata,theta[[N]])

norm = 0.0
for( m in 1:traindata$curve.num){
  norm = norm + logSumExp(log.P[m,]+ log(theta[[N]]$pi))
}
