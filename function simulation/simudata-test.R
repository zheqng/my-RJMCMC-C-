rm(list = ls())
library('MCMCpack')
# # setwd('/home/zheqng/function simulat')
setwd('/home/zheqng/src/RJMCMC-my-C/large K/function simulation/')
# these = read.table(file='/home/zheqng/src/RJMCMC-my-C/large K/function simulation/theta.txt')
# # setwd('/media/zheqng/Seagate Backup Plus Drive/zheqng@nwu/src/RJMCMC-my-C/simu1/2019.2.10/function simulation/')
# library("MASS")
# 


load("simudata.RData")


source('basicfunc.r')
source('plot.mix.r')

for(m in 7:8)
{
  
  dat[[m]]$y = mvrnorm(n=1,mu=rep(0,length(x)),
                       Sigma = cov(x,theta[[4]]))
  dat[[m]]$k=4;
  dat[[m]]$x = (dat[[m]]$x)
  # /8*0.01
}
for(m in 13:14)
{
  
  dat[[m]]$y = mvrnorm(n=1,mu=rep(0,length(x)),
                       Sigma = cov(x,theta[[7]]))
  dat[[m]]$k=7;
  dat[[m]]$x = (dat[[m]]$x)
  # /8*0.01
}
save.image("simudatanew.RData")
# load("simudata.RData")
plot.mixgaussian(dat,step=2,K=K)
plot.mixgaussian(dat,step=2,K=K,make.pdf=TRUE)
# write to file

stepsize = length(x)
tsize = floor(stepsize/2)

# for(k in 1:3)
for(m in 7:8)
{
  xtoltrain = 1:stepsize
  xindtrain = sample(1:stepsize,size = tsize,replace = FALSE)
  xindtrain = sort(xindtrain)
  xresttrain = xtoltrain[-xindtrain]
  sink('../demo/traindata.dat',append = TRUE)
  cat(dat[[m]]$x[xindtrain],"\n")
  cat(round(dat[[m]]$y[xindtrain],digits=4),"\n")
  sink()
  sink('../demo/testdata.dat',append = TRUE)
  cat(dat[[m]]$x[xresttrain],"\n")
  cat(round(dat[[m]]$y[xresttrain],digits=4),"\n")
  sink()
}
for(m in 13:14)
{
  xtoltrain = 1:stepsize
  xindtrain = sample(1:stepsize,size = tsize,replace = FALSE)
  xindtrain = sort(xindtrain)
  xresttrain = xtoltrain[-xindtrain]
  sink('../demo/traindata.dat',append = TRUE)
  cat(dat[[m]]$x[xindtrain],"\n")
  cat(round(dat[[m]]$y[xindtrain],digits=4),"\n")
  sink()
  sink('../demo/testdata.dat',append = TRUE)
  cat(dat[[m]]$x[xresttrain],"\n")
  cat(round(dat[[m]]$y[xresttrain],digits=4),"\n")
  sink()
}
# sink()
# save.image("simudata.RData")
# load("simudata.RData")
# plot.mixgaussian(dat,step=2)
# plot.mixgaussian(dat,step=2,make.pdf=TRUE)

# for(i in 1:19) cat(i,'th',these[[i]],'\n')
