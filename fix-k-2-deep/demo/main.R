rm(list = ls())
# setwd("/media/zheqng/Seagate Backup Plus Drive/zheqng@nwu/src/RJMCMC-my-C/simu1/fix-k-move/demo/")
setwd('/home/zheqng/src/RJMCMC-my-C/novel perimental/electric load 2008/fix-k-2-deep/demo/')
library('matrixStats')
# library('GPFDA')
library('MASS')
library('mvtnorm')
source('predict.r')
# load('result.RData')

source('readsts.R')
source('bayss analysis.R')
source('convergence diagnose.R')
#######################################
tempdata = read.table("electric_load_2008.dat",fill=TRUE)
N = dim(z)[1]
curve.num = dim(tempdata)[1]/2
Nm = dim(tempdata)[2]
Time = tempdata[seq(1,200,by=2),]
elecload = tempdata[seq(2,200,by=2),]
plot.z(Time,elecload,z[N,]+1)
##################################
# load('../function simulation/simudata.RData')
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

#######################################
tempdata = read.table("~/src/RJMCMC-my-C/electric load/demo/traindata.dat",fill=TRUE)
curve.num = dim(tempdata)[1]
Nm = dim(tempdata)[2]
Time = seq(1,96)
elecload = tempdata
plot.z(Time,load,z[20000,])

theta= vector("list",N)
for(i in 1:N){
  k = datha$K[i]
  theta[[i]]$pi <- theta.tmp[1:k]
  theta[[i]]$w<- theta.tmp[(k+1):(2*k)]
  theta[[i]]$v<- theta.tmp[(2*k+1):(3*k)]
  theta[[i]]$sigma2<- theta.tmp[(3*k+1):(4*k)]
  theta[[i]]$K = k
  theta.tmp<-theta.tmp[-(1:(4*k))]
}

theta = list("pi","w","v","sigma2","K")
theta$K=2
theta$pi=c(0.25211,0.74789)
theta$w=c(0.381064, 0.387153)
theta$v=c(153.693,78.4465)
theta$sigma2=c(0.0756031, 0.0148021)


log.P = rep(0,2)
for(k in 1:(theta$K)){
  log.P[k] =my.dmultinorm(x =t( elecload[1,]),mu = rep(0,length(elecload[1,])),sigma=exp.cov.noise(t(Time),theta,k),log=TRUE)
}

#####################################
tempdata = read.table("testdata.dat",fill=TRUE)
RMSE = as.vector(rep(0,14*96-97))
SD = as.vector(rep(0,14*96-97))
aaa = 0.0
for(ii in 1:(14*96-97)){
  x = seq(ii,ii+95)
  y = tempdata[ii:(ii+95)]
  x.new = ii+96
 k= calc.z.2(x,y,theta[[10]])
 y.new = tempdata[x.new]
  predictdata = predict.gp(x,y,x.new,y.new,thet = theta[[10]],k)
  RMSE[ii] = y.new - predictdata$pred.mean
  SD[ii] = predictdata$pred.sd
  aaa = aaa+RMSE[[ii]]^2
  print(ii)
}
aaa/(length(RMSE))
mean(SD)
