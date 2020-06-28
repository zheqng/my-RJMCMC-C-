rm(list = ls())
library('MCMCpack')
# setwd('/home/zheqng/function simulat')
setwd('/home/zheqng/src/RJMCMC-my-C/novel perimental/regular/function simulation/')
these = read.table(file='/home/zheqng/src/RJMCMC-my-C/novel perimental/regular/function simulation/theta.txt')
# setwd('/media/zheqng/Seagate Backup Plus Drive/zheqng@nwu/src/RJMCMC-my-C/simu1/2019.2.10/function simulation/')
library("MASS")

source('basicfunc.r')
source('plot.mix.r')


K=3
curve.num = 90
Nm = 50
theta = vector("list",K)
# for(k in 1:100){
PI  = rep(1/K,1,K);
for(k in 1:K){
  theta[[k]]$w=these[k,1]
  theta[[k]]$v=these[k,2]
  theta[[k]]$sigma2=these[k,3]
  theta[[k]]$pi = PI[k]
  
}

z.training = sample(1:K,curve.num,replace= TRUE,prob=PI)


traindata =vector( "list", curve.num)
traindata$curve.num=curve.num
for(m in 1:curve.num)
{
  k=z.training[m]
  x= runif(Nm,-4,4)
  x= sort(x)
  traindata[[m]]$x =x
  traindata[[m]]$y = mvrnorm(n=1,mu=rep(0,length(x)),
                       Sigma = cov(x,theta[[k]]))
  traindata[[m]]$k=k;
}

# save.image("simudata.RData")
#######################################################
# load("simudata.RData")
plot.mixgaussian(traindata,z.training,K=K)
plot.mixgaussian(traindata,z.training,K=K,make.pdf=TRUE)
# write to file

###################################3
#traindata
unlink('../demo/traindata.dat')
sink('../demo/traindata.dat',append = TRUE)
for(m in 1:curve.num)
{
  cat(traindata[[m]]$x,"\n");
  cat(round(traindata[[m]]$y,digits=4),"\n");

}
sink()
# sink()
####################################################################3
curve.num.test = 300
Nm.test.known = 40
Nm.test.unkown = 110
Nm.test = Nm.test.known + Nm.test.unkown
z.test=sample(1:K,curve.num.test,replace= TRUE,prob=PI)


valide.dat =vector( "list", curve.num.test)
valide.dat$curve.num=curve.num.test

for(m in 1:curve.num.test)
{
  k=z.test[m]
  valide.dat[[m]]$k=k;
  x=runif(Nm.test,-4,4)
  x = sort(x)
  valide.dat[[m]]$x = x;
  valide.dat[[m]]$y = mvrnorm(n=1,mu=rep(0,Nm.test),
                       Sigma = cov(x,theta[[k]]))
}

testdata.known = vector("list",curve.num.test)
testdata.unknown = vector("list",curve.num.test)
testdata.known$curve.num = curve.num.test
testdata.unknown$curve.num = curve.num.test
for(m in 1:curve.num.test)
{
  xtoltrain = 1:Nm.test
  xindtrain = sample(xtoltrain,Nm.test.known,replace = FALSE)
  xindtrain = sort(xindtrain)
  xresttrain = xtoltrain[-xindtrain]
testdata.known[[m]]$x = valide.dat[[m]]$x[xindtrain]
testdata.known[[m]]$y = valide.dat[[m]]$y[xindtrain]
testdata.unknown[[m]]$x = valide.dat[[m]]$x[xresttrain]
testdata.unknown[[m]]$y = valide.dat[[m]]$y[xresttrain]
}


save.image("simudata.RData")

