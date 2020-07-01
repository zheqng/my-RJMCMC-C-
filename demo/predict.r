# theta[i]$w theta[i]$v theta[i]$pi theta[i]$sigma2 theta[i]$k
# theta = matrix(rep(0,4*10),ncol = 4)
#################################
# likelihood and posterior
xixj<-function(x,x.star){
  n = length(x);m = length(x.star)
  tmp = matrix(rep(0,m*n),nrow = n)
  for(i in 1:n){
    for(j in 1:m){
      tmp[i,j] = (x[i] - x.star[j])^2
      # tmp[j,i] = tmp[i,j]
    }
  }
  return(tmp)
}
exp.cov<-function(x,x.star,thet,k){
  n = length(x)
  tmp<-(thet$v[k])*exp(-xixj(t(x),t(x.star))*(thet$w[k])/2) 
  return(tmp)
}
exp.cov.noise<-function(x,thet,k){
  n = length(x)
  tmp<-(thet$v[k])*exp(-xixj(t(x),t(x))*(thet$w[k])/2) + (thet$sigma2[k])*diag(n)
  return(tmp)
}
# dist(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

my.dmultinorm<-function(x,mu = rep(0,length(x)),sigma,log=TRUE){
   Nm = length(x)
   # x = t(x)
  log.d =  -t(x)%*%solve(sigma,x)/2.0 - Nm/2.0 * log(2.0 * 3.14159) - 1/2 * log(abs(det(sigma)))
  if(log==TRUE) return(log.d)
  else return(exp(log.d))
}
log.likelihood.micro <- function(traindata,thet,m,k){
  norm = dmvnorm(traindata$Y[m,],sigma=exp.cov.noise(traindata$X[m,],thet,k),log=TRUE)
   return(norm)
}

print.posterior<-function(traindata,thet){
  K = 
  norm =0.0;log.P = matrix(rep(0,(traindata$curve.num * thet$K)),traindata$curve.num,thet$K)
  for(m in 1:(traindata$curve.num)){
    # log.P = rep(0,3)
    for(k in 1:(thet$K)){
      log.P[m,k] = my.dmultinorm(x = traindata[[m]]$y,mu = rep(0,length(traindata[[m]]$y)),sigma=exp.cov.noise(traindata[[m]]$x,thet,k),log=TRUE)
    }
  }
  log.P
  return(log.P)
}
log.likelihood <-function(traindat,thet){
  norm =0.0;
  for(m in 1:(traindat$curve.num)){
    log.P = rep(0,3)
  for(k in 1:(thet$K)){
    log.P[k] =my.dmultinorm(x = traindat[[m]]$y,mu = rep(0,length(traindat[[m]]$y)),sigma=exp.cov.noise(traindat[[m]]$x,thet,k),log=FALSE)
      # dmvnorm(traindata$Y[m,],mean=rep(0,length(traindata$X[m,])),sigma=exp.cov.noise(traindata$X[m,],thet,k),log=TRUE)
  }
    norm = norm + sum(log.P)
  }
  return(norm)
}
##############################################################
#pridict gp curves
predict.gp<-function(x,y,x.new,y.new,thet,k){
  K.star<-exp.cov(x,x.new,thet,k) 
  K<-exp.cov.noise(x,thet,k)
  K.inv <-ginv(K)
  mu<- t(K.star)%*% K.inv %*% y
  sigma2<-exp.cov(x.new,x.new,thet,k) - t(K.star) %*% K.inv %*% (K.star)
  sigma2<-diag(sigma2)
  correlation<-cor(mu,y.new)
result=c(list('pred.mean'=mu,'pred.sd'=sqrt(sigma2),'newdata'=x.new,'correlation'=correlation))

return(result)
}



predict.gp.ave.m<-function(traindata,testdata,theta,z.iter,m,lag = 10){
  T = dim(z.iter)[1]
  # T = length(theta)
  warm = floor(T/2)
  batch = seq(warm,T,by = lag)
  # mu = matrix(rep(0,51*9),traindata$curve.num,51)
  # result = vector("list",traindata$curve.num)
  Test.size = length(testdata[[m]]$x)
  # for(m in 1:traindata$curve.num)
  # {
    pred.mean = matrix(rep(0,Test.size),Test.size,1)
    pred.var = matrix(rep(0,Test.size),Test.size,1)

    
    for(i in batch){
      k = z.iter[i,m]
      predictdata<-predict.gp(traindata[[m]]$x,traindata[[m]]$y,testdata[[m]]$x,testdata[[m]]$y,theta[[i]],k)
      pred.mean<-pred.mean + predictdata$pred.mean
      pred.var <- pred.var + (predictdata$pred.sd)^2 + (predictdata$pred.mean)^2 
      # print(i)
    }
    pred.mean = pred.mean/length(batch)
    pred.var = pred.var/length(batch) - (pred.mean)^2
    pred.sd = sqrt( pred.var)
    # newdata = testdata$X[m,]
    # print(m)
    # }
    result=c(list('pred.mean'=pred.mean,'pred.sd'=pred.sd,'x'=testdata[[m]]$x))
    
  return(result)
}

#####
#plot gp predict
plot.gp.predict<- function(traindata,testdata,m){
  
  upper=testdata[[m]]$pred.mean+1.96*(testdata[[m]]$pred.sd);
  lower=testdata[[m]]$pred.mean-1.96*(testdata[[m]]$pred.sd);
  plot(-100,-100,col=0,xlim=range(traindata[[m]]$x,testdata[[m]]$x),ylim=range(upper,lower,traindata[[m]]$y),type='l', xlab="input ",ylab="response")
  #
  polygon(c(testdata[[m]]$x, rev(testdata[[m]]$x)), c(upper, rev(lower)),col = rgb(127,127,127,120, maxColorValue = 255), border = NA)
  #
  lines(t(testdata[[m]]$x),t(testdata[[m]]$y),pch=8,col="blue",cex=1.5)
  lines(testdata[[m]]$x,testdata[[m]]$pred.mean,col="red",lwd=1)  
  points(t(testdata[[m]]$x),t(testdata[[m]]$y),pch=8,col="blue",cex=0.5)
  points(testdata[[m]]$x,testdata[[m]]$pred.mean,col="red",lwd=1,cex = 0.5)  
  # legend("topright",legend = c("1st","2nd ","3rd"),col = c("sienna3","gold4","steelblue"),lty = 1)
  legend("topright",legend =c("true value","prediction"),col=c("blue","red"),cex=1,lty=1,bty="n")
}

plot.gp.test<-function(traindata,testdata,i){
  # i=1
  plot(t(traindata[[i]]$x),t(traindata[[i]]$y),type='l',ylim=range(traindata[[i]]$y))
  # for(i in 2:traindata$curve.num){
    # lines(t(traindata[[i]]$x),t(traindata[[i]]$y),type='l',ylim=range(traindata[[i]]$y))
    # points(testdata$X[i,],testdata$Y[i,])
    # title(main=expression(paste(i,"curve")))
  
  # for(i in 1:testdata$curve.num){
    points(testdata[[i]]$x,testdata[[i]]$y)
  
  
}

plot.gp.i.predict<-function(traindata,predictdata,i){
  plot(t(traindata$X[i,]),t(traindata$Y[i,]),type='l',ylim=range(traindata$Y))
  points(predictdata$newdata,predictdata$pred.mean,lwd=1)  
  
}
# rmse(mu,y.new)
rmse <- function(testdata){
 y= testdata$y
  # y = t(as.matrix(y))
  sqrt(mean( (y - testdata$pred.mean)^2))
}
###############################################
#calculation  z and alpha, plot the final posterior
calc.alpha<-function(traindata,thet){
  log.P =  print.posterior(traindata,thet)
  tmp = (log.P + log (thet$pi)) 
  log.alpha = tmp  - apply(tmp,1,logSumExp)
  alpha=exp(log.alpha) 
  return(alpha)
}
calc.test.z<-function(traindata,thet){
  alpha.value = calc.alpha(traindata,thet)
 z<-apply(alpha.value,1,which.max)
 return(z)
}
calc.z<-function(alpha.value){
  z<-apply(alpha.value,1,which.max)
  
}

plot.alpha<-function(alpha,z){
  # N.paper=N
  color=c("red","blue","green")
  # N.tmp=dim(alpha)[1]
  curve.num = dim(alpha)[1]
  plot(-100,-100,xlim=c(1,curve.num),ylim = c(0,1),xlab="iteration",ylab="classification")
  
  for(i in 1:curve.num){
      points(i,apply(alpha,1,max)[i],col=color[z[i]])
    
  }
}

