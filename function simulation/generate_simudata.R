rm(list = ls())

setwd('/media/zheqng/Seagate Backup Plus Drive/zheqng@nwu/src/RJMCMC-my-C/simu1/2019.2.10/function simulation/')
library("MASS")

source('basicfunc.r')
source('plot.mix.r')

M = 9

dat =vector( "list", M)
x=seq(from=-4,to=4,by=0.08)
for(i in 1:M)
dat[[i]]$x = x;
dat$M=M

theta = vector("list",3)
theta[[1]]$pi = 1/3;theta[[1]]$w = 1.0;theta[[1]]$v = 0.2; theta[[1]]$sigma2 = 0.0025
theta[[2]]$pi = 1/3;theta[[2]]$w = 0.5;theta[[2]]$v = 1; theta[[2]]$sigma2 = 0.001
theta[[3]]$pi = 1/3;theta[[3]]$w = 10;theta[[3]]$v = 0.2; theta[[3]]$sigma2 = 0.0005

for(k in 1:3)
  for(m in 1:3)
    {
      dat[[(k-1)*3 + m]]$y = mvrnorm(n=1,mu=rep(0,length(x)),
                                 Sigma = cov(X = dat[[(k-1)*3 + m]]$x,theta = theta[[k]]))
      dat[[(k-1)*3 + m]]$k=k;
  }

# write to file

stepsize = length(x)
tsize = floor(stepsize/2)

for(k in 1:3)
  for(m in 1:3)
  {
    xtoltrain = 1:stepsize
    xindtrain = sample(1:stepsize,size = tsize,replace = FALSE)
    xindtrain = sort(xindtrain)
    xresttrain = xtoltrain[-xindtrain]
    sink('../train/demo/traindata.dat',append = TRUE)
     cat(dat[[(k-1)*3 + m]]$x[xindtrain],"\n")
    cat(dat[[(k-1)*3 + m]]$y[xindtrain],"\n")
    sink()
    sink('../train/demo/testdata.dat',append = TRUE)
    cat(dat[[(k-1)*3 + m]]$x[xresttrain],"\n")
    cat(dat[[(k-1)*3 + m]]$y[xresttrain],"\n")
    sink()
    }
  # sink()
save.image("simudata.RData")
# load("simudata.RData")
plot.mixgaussian(dat)
plot.mixgaussian(dat,make.pdf=TRUE)


