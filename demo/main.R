rm(list = ls())
# setwd("/media/zheqng/Seagate Backup Plus Drive/zheqng@nwu/src/RJMCMC-my-C/simu1/fix-k-move/demo/")
setwd('/home/zheqng/src/RJMCMC-my-C/novel perimental/regular/demo/')
library('matrixStats')
# library('GPFDA')
library('MASS')
library('mvtnorm')
source('predict.r')
# load('result.RData')

source('readsts.R')
source('bayss analysis.R')
source('convergence diagnose.R')

load('../function simulation/simudata.RData')
# testdata.known$curve.num = curve.num.test
# testdata.unknown$curve.num =  curve.num.test
# traindata$curve.num = curve.num
##############################################
#analysis classification and alpha for testdata
# curve.num = traindata$curve.num
T  = t.or+1
K = datha.K[160000]
ortheta$K = K

start.time = Sys.time()
alpha.ave = matrix(0,curve.num.test,K)
z.iter = matrix(0,T,curve.num.test)
for(iter in seq(1,T,by=10)){
  alpha.iter = calc.alpha(testdata.known,ortheta[[iter]])
  alpha.ave = alpha.ave +  alpha.iter
  z.iter[iter,] = calc.z(alpha.iter)
  iter
}
dt <- difftime(Sys.time(), start.time, units="mins")
dt
alpha.ave = alpha.ave/ length(seq(1,T,by = 10))


 z.est = calc.z(alpha.ave[,perma[[best]]])
accuracy = (curve.num.test-sum(z.test==z.est))/curve.num.test*100
cat(paste("error of clustering is",accuracy,"%\n"))
plot.alpha(alpha.ave,z.est)
#############################
# RMSE predict
RMSE.m = rep(0,curve.num.test)
SD.m = rep(0,curve.num.test)
op <- par(mfrow=c(1,1))
for(m in 1:curve.num.test)
    {tmp = predict.gp.ave.m(testdata.known,testdata.unknown,ortheta,z.iter[seq(1,T,by = 10),],m,lag = 10)
    # result=c(list('pred.mean'=pred.mean,'pred.sd'=pred.sd,'x'=testdata[[m]]$x))
    # 
    testdata.unknown[[m]]$pred.mean = tmp$pred.mean
    testdata.unknown[[m]]$pred.sd = tmp$pred.sd
    RMSE.m[m]<- rmse(testdata.unknown[[m]])
    # SD.m[m]<-mean(testdata.unknown[[m]]$pred.sd)
}

RMSE = mean(RMSE.m)
# SD = mean(SD.m)


#################################
#plot the prediction figure
# predict.sum<-predict.gp.ave(traindata,testdata,theta,z,lag = 100)

op <- par(mfrow=c(3,3))
for(m in 1:9){
  plot.gp.predict(testdata.known,testdata.unknown,m)
  title(paste(T,"th iter ","predict",i=m,"curve"))
  print(m)
}
op<- par(mfrow = c(1,1))

