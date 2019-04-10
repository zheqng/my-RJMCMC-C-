rm(list = ls())
setwd("/media/zheqng/Seagate Backup Plus Drive/zheqng@nwu/src/RJMCMC-my-C/simu1/fix-k-move/demo/")

library('matrixStats')
library('GPFDA')
library('MASS')
library('mvtnorm')
source('predict.r')
source('readsts.R')
for(i in 1:N){
  m=7;i=1;
  k = z[i,m]+1
    predictdata<-predict.gp(traindata$X[m,],traindata$Y[m,],testdata$X[m,],theta[[i]],k)
    i=m;plot.gp.i.predict(traindata,predictdata,i)
    title(paste(N,"th iter ","predict",i,"curve"))
    # plot.gp.predict(traindata,testdata,predictdata,1)
    # plot.gp(traindata ,testdata )
}

log.P =   print.posterior(traindata,theta[[N]])

norm = 0.0
for( m in 1:traindata$curve.num){
  norm = norm + logSumExp(log.P[m,]+ log(theta[[N]]$pi))
}
