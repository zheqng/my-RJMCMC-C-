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
  tmp<-(thet[[k]]$v)*exp(-xixj(t(x),t(x.star))*(thet[[k]]$w)/2) 
  return(tmp)
}
exp.cov.noise<-function(x,thet,k){
  n = length(x)
  tmp<-(thet[[k]]$v)*exp(-xixj(x,x)*(thet[[k]]$w)/2) + (thet[[k]]$sigma2)*diag(n)
  return(tmp)
}
# dist(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

my.dmultinorm<-function(x,mu = rep(0,length(x)),sigma,log=TRUE){
  Nm = length(x)
  x=x-mu
  x=t(x)
  log.d =  -x%*%solve(sigma,t(x))/2.0 - Nm/2.0 * log(2.0 * 3.14159) - 1/2 * log(abs(det(sigma)))
  if(log==TRUE) return(log.d)
  else return(exp(log.d))
}


print.posterior<-function(traindata,thet){
  norm =0.0;log.P = matrix(rep(0,(traindata$M * thet$K)),traindata$M,thet$K)
  for(m in 1:(traindata$M)){
    # log.P = rep(0,3)
    for(k in 1:(thet$K)){
      log.P[m,k] = my.dmultinorm(x = traindata[[m]]$y,mu =rep(0,length(traindata[[m]]$x)),sigma=exp.cov.noise(traindata[[m]]$x,thet,k),log=TRUE)
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