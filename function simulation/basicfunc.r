
cov <- function(x,theta){
  (theta$v) *exp(-0.5*(theta$w)*( (outer(x,x,"-")^2)))+diag(length(x))*(theta$sigma2)
  
}



mix.posterio <- function(dat, theta, m,k) {
  result<-theta[[k]]$pi * dmvnorm(dat[[m]]$y,
                                  mean = rep(0, length(dat[[m]]$y)),
                                  sigma = cov(X=dat[[m]]$x, theta = theta[[k]]),
                                  log = FALSE
  )
  result
  
}
