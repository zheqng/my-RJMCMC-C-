# rm( list = ls())
# setwd("~/src/RJMCMC-my-C/fix-k-move/demo/")
# setwd("/media/zheqng/Seagate Backup Plus Drive/zheqng@nwu/src/RJMCMC-my-C/simu1/fix-k-move/demo/")
library("fda")
# library("mcsm")

datha<- read.table("electric_load_2008.sts",fill=TRUE,col.names = c("iter number","K","log lik","split/merge","acc/reject","prob","simu_acc/simu_rej","prob","death or reserve"))
datha.loglik = datha[,"log.lik"]
datha.K = datha[,"K"]


theta.tmp<-read.table("parameter.res",fill=TRUE)
theta.tmp <- as.matrix(theta.tmp)
theta.tmp <- as.vector(theta.tmp)
N = nrow(datha)
# N = 200
# tmp= datha.K[1:N]
# rm(datha.K)
# datha.K = tmp
# rm(tmp)
# tmp = datha.loglik[1:N]
# rm(datha.loglik)
# datha.loglik = tmp
# index = (c(0,cumsum(datha$K)[-N])*4+1):(cumsum(datha$K)*4)

# time <- read.table("time.txt",fill = TRUE)
# time<-as.matrix(time)
# time<- time[-1,]
# time<-matrix(t(time),ncol=1)
# time<-as.vector(time)
# N.half = N/2
# time<-time[(1:N.half)*2]
# time.diff<-time[2:N.half] - time[1:(N.half -1)]
# plot(time[2:N.half] - time[1:(N.half -1)])

z <-read.table("z.res",fill = TRUE)

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

theta.tmp = vector("list",N/100)
j=1
for(i in seq(1,N,by=100)){
  theta.tmp[[j]] = theta[[i]]
  j=j+1
}
theta.ori = theta
theta = vector("list",N/100)
theta = theta.tmp
# w = NULL
# for(i in 1:N){
#   w <-rbind(w,theta[[i]]$w)
# }

# load('../function simulation/simudata.RData')

traindata = list("temp","X","Y","curve.num")
traindata$temp =read.table("electric_load_2008.dat")
traindata$curve.num = dim(traindata$temp)[1]/2
traindata$X=(traindata$temp)[(1:traindata$curve.num)*2-1,]
traindata$Y=traindata$temp[(1:traindata$curve.num)*2,]

# testdata = list("temp","X","Y","curve.num")
# testdata$temp =read.table("testdata.dat")
# testdata$curve.num = dim(testdata$temp)[1]/2
# testdata$X=(testdata$temp)[(1:testdata$curve.num)*2-1,]
# testdata$Y=testdata$temp[(1:testdata$curve.num)*2,]
# 
# validedata = list("temp","X","Y","curve.num")
# validedata$temp =read.table("../../validedata.dat")
# validedata$curve.num = dim(validedata$temp)[1]/2
# validedata$X=(validedata$temp)[(1:validedata$curve.num)*2-1,]
# validedata$Y=validedata$temp[(1:validedata$curve.num)*2,]

###############################################################
# N.paper=20000
plot(datha.K,type='l',xlab="iteration",ylab="components number")
plot(datha.loglik,type='l',xlab="iteration",ylab="log likelihood")
#############################################################3
boxplot(log.lik~K,data = datha)
summary(as.factor(datha$K))


boxplot(prob ~ split.merge, data = datha)
summary(datha[datha$split.merge=="split",-c(1,2,3,7,8,9)])
summary(datha[datha$split.merge=="merge",-c(1,2,3,7,8,9)])
# min(which(datha.K==2))
###############################################################
# N=N*100
datha.loglik = datha.loglik[seq(1,N,by=100)]
datha.K = datha.K[seq(1,N,by=100)]
datha=datha[seq(1,N,by=100),]

