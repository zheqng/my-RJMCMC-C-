SE<-function(xi,xj,w){
exp(-0.5*w*( xi - xj)^2)
}

cov <- function(X,Y=NULL,theta){
  if(is.null(Y))Y=X
  theta$v * outer(X,Y,SE,theta$w) 
}



mix.posterio <- function(dat, theta, m,k) {
  result<-theta[[k]]$pi * dmvnorm(dat[[m]]$y,
    mean = rep(0, length(dat[[m]]$y)),
    sigma = cov(X=dat[[m]]$x, theta = theta[[k]]),
    log = FALSE
  )
  result
  
}
