# theta[i]$w theta[i]$v theta[i]$pi theta[i]$sigma2 theta[i]$k
# theta = matrix(rep(0,4*10),ncol = 4)
xixj<-function(x,x.star){
  n = length(x);m = length(x.star)
  tmp = matrix(rep(0,m*n),nrow = n)
  for(i in 1:n){
    for(j in i:m){
      tmp[i,j] = (x[i] - x.star[j])^2
      tmp[j,i] = tmp[i,j]
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
   x = t(x)
  log.d =  -t(x)%*%solve(sigma,x)/2.0 - Nm/2.0 * log(2.0 * pi) - 1/2 * log(abs(det(sigma)))
  if(log==TRUE) return(log.d)
  else return(exp(log.d))
}
log.likelihood.micro <- function(traindata,thet,m,k){
  norm = dmvnorm(traindata$Y[m,],sigma=exp.cov.noise(traindata$X[m,],thet,k),log=TRUE)
   return(norm)
}

print.posterior<-function(traindata,thet){
  norm =0.0;log.P = matrix(rep(0,(traindata$curve.num * thet$K)),traindata$curve.num,thet$K)
  for(m in 1:(traindata$curve.num)){
    # log.P = rep(0,3)
    for(k in 1:(thet$K)){
      log.P[m,k] = my.dmultinorm(x = traindata$Y[m,],mu = rep(0,length(traindata$Y[m,])),sigma=exp.cov.noise(traindata$X[m,],thet,k),log=TRUE)
    }
  }
  log.P
  return(log.P)
}
log.likelihood <-function(traindata,thet){
  norm =0.0;
  for(m in 1:(traindata$curve.num)){
    log.P = rep(0,3)
  for(k in 1:(thet$K)){
    log.P[k] =my.dmultinorm(x = traindata$Y[m,],mu = rep(0,length(traindata$Y[m,])),sigma=exp.cov.noise(traindata$X[m,],thet,k),log=FALSE)
      # dmvnorm(traindata$Y[m,],mean=rep(0,length(traindata$X[m,])),sigma=exp.cov.noise(traindata$X[m,],thet,k),log=TRUE)
  }
    norm = norm + sum(log.P)
  }
  return(norm)
}
predict.gp<-function(x,y,x.new,thet,k){
  K.star<-exp.cov(x,x.new,thet,k) 
  K<-exp.cov.noise(x,thet,k)
  # K.inv<-solve(exp.cov.noise(x,thet,k) )
  # mu<- t(K.star)%*% K.inv %*% t(y)
  # solve(t(A)%*%A,t(A)%*%B)
  K.inv <-ginv(K)
  mu<- t(K.star)%*% K.inv %*% t(y)
  sigma2<-exp.cov(x.new,x.new,thet,k) - t(K.star) %*% K.inv %*% (K.star)
  sigma2<-diag(sigma2)
result=c(list('pred.mean'=mu,'pred.sd'=sqrt(sigma2),'newdata'=x.new))

return(result)
}

plot.gp.predict<- function(traindata,testdata,predictdata,row.no){

  upper=predictdata$pred.mean+1.96*(predictdata$pred.sd);
  lower=predictdata$pred.mean-1.96*(predictdata$pred.sd);
  plot(-100,-100,col=0,xlim=range(traindata$X[row.no,],testdata$X[row.no,]),ylim=range(upper,lower,traindata$Y[row.no,]),type='l', xlab="input ",ylab="response")
  #
  polygon(c(predictdata$newdata, rev(predictdata$newdata)), c(upper, rev(lower)),col = rgb(127,127,127,120, maxColorValue = 255), border = NA)
  #
  points(traindata$X[row.no,],traindata$Y[row.no,],pch=8,col=2,cex=1.5)
  lines(predictdata$newdata,predictdata$pred.mean,col=4,lwd=1)  
}

plot.gp<-function(traindata,testdata){
  i=1
  plot(t(traindata$X[i,]),t(traindata$Y[i,]),type='l',ylim=range(traindata$Y))
  for(i in 2:traindata$curve.num){
    lines(t(traindata$X[i,]),t(traindata$Y[i,]),type='l',ylim=range(traindata$Y))
    # points(testdata$X[i,],testdata$Y[i,])
    # title(main=expression(paste(i,"curve")))
  }
  for(i in 1:testdata$curve.num){
    points(testdata$X[i,],testdata$Y[i,])
  }
  
}

plot.gp.i.predict<-function(traindata,predictdata,i){
  plot(t(traindata$X[i,]),t(traindata$Y[i,]),type='l',ylim=range(traindata$Y))
  points(predictdata$newdata,predictdata$pred.mean,lwd=1)  
  
}
# rmse(mu,y.new)