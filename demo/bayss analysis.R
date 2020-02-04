library('bayess')
library('MCMCpack')
plot.trace<-function(tt,PI){
  plot(tt[,1],ylim=range(tt),
       ylab=expression(w[1]),xlab="n",type="l",col="sienna3")
  lines(tt[,2],col="gold4")
  lines(tt[,3],col="steelblue")
  plot(tt[,2],PI[,2],col="sienna3",
       xlim=range(tt),ylim=range(PI),
       xlab=expression(w[2]),ylab=expression(p[2]))
  points(tt[,3],PI[,3],col="steelblue")
}

lpost<-function(traindata,theta)
{
  N = length(theta)
  result = as.vector(rep(0,N))
  for( i in 1:N){
    for(m in 1:traindata$curve.num){
      logl = matrix(0,1,theta[[i]]$K)
      for(k in 1:theta[[i]]$K){logl[k] = log.likelihood.micro(traindata,theta[[i]],m,k)}
      result[i] = result[i] + logSumExp(log(theta[[i]]$pi) + logl) + sum(log(dinvgamma(theta[[i]]$w,0.5,0.5)))
      + sum(log(dinvgamma(theta[[i]]$v,0.5,0.5))) +   sum(log(ddirichlet(theta[[i]]$pi,as.vector(rep(1,theta[[i]]$K))  )))
    
      }
    print(i)
  }
  return(result)
  
}
 



  # library(plotrix)
  # twoord.plot(1:N,lpost,1:N,datha.K,type=c("p","p"))
  logpost<-lpost(traindata,theta)
  state = data.frame(cbind(datha.K,logpost))
  names(state) <- c("K","post")
  # boxplot(post~K,data = state)
# summary(as.factor(state$K))

# state.iter = data.frame(cbind(datha$K,datha$log.lik))


# boxplot(K ~ split.merge, data = datha)
# hist(datha[datha$split.merge=="split","K"])
# hist(datha[datha$split.merge=="merge","K"])
# # 
# llikely<-function(traindata,theta)
# {
#   N = length(theta)
#   result = as.vector(rep(0,N))
#   for( i in 1:N){
#     for(m in 1:traindata$curve.num){
#       logl = matrix(0,1,theta[[i]]$K)
#       for(k in 1:theta[[i]]$K){logl[k] = log.likelihood.micro(traindata,theta[[i]],m,k)}
#       result[i] = result[i] +  logSumExp(logl+ log(theta[[i]]$pi) )
#     }
#     print(i)
#   }
#   return(result)
#   
# }
# 
# loglikely<-llikely(traindata,theta)
# state1 = data.frame(cbind(datha.K,loglikely))
# names(state1) <- c("K","likely")

# plot(datha.loglik,type = 'l',col = "black",ylim=range(c(logpost,datha.loglik)))
# lines(state$post,col="red")
# lines(state1$likely,col = "blue")


indimap=order(logpost,decreasing=TRUE)[1]
# indimap = N.3[indimap]
map=list(K = theta[[indimap]]$K,pi = theta[[indimap]]$pi, w = theta[[indimap]]$w, v = theta[[indimap]]$v, 
         sigma2 = theta[[indimap]]$sigma2)
        


lili=alloc=matrix(0,traindata$curve.num,map$K)
for (m in 1:(traindata$curve.num)){
  logl = as.vector(rep(0,map$K))
  for(k in 1:(map$K))
  {
    logl[k] = my.dmultinorm(x = traindata$Y[m,],mu = rep(0,length(traindata$Y[m,])),sigma=exp.cov.noise(traindata$X[m,],map,k),log=TRUE)
  }
  lili[m,]=log(map$pi) + logl
  # tmp = max(lili[m,])
  lili[m,]=exp(lili[m,] - logSumExp(lili[m,]))
    # exp(lili[m,]-tmp)/sum(exp(lili[m,]-tmp))
}
LogLoss <- function(y_pred, y_true) {
  eps <- 1e-15
  y_pred <- pmax(pmin(y_pred, 1 - eps), eps)
  LogLoss <- sum(y_true * log(y_pred) )
  return(LogLoss)
}

N.3 = which(datha.K==3)
ortheta = vector("list",length(N.3))
library(combinat)
perma=permn(3)
t.or=1
for (t in N.3){
  entropies=rep(0,factorial(3))
  for (m in 1:(traindata$curve.num)){
    log.P = as.vector(rep(0,3))
    for(k in 1:3)log.P[k] = log.likelihood.micro(traindata,theta[[t]],m,k)
     alloc[m,] = log(theta[[t]]$pi)+ log.P
    # alloc[m,]=simu$p[t,]*dnorm(datha[m],mean=simu$mu[t,],
                               # sd=sqrt(simu$sig[t,]))
     # min.alloc = min(alloc[m,])
    alloc[m,]= exp(alloc[m,] - logSumExp(alloc[m,]))
      # exp(alloc[m,] - min.alloc)/sum(exp(alloc[m,]-min.alloc))
    for (i in 1:factorial(3))
      entropies[i]=entropies[i]+LogLoss(alloc[m,perma[[i]]],lili[m,])
      # sum(lili[m,]*log(alloc[m,perma[[i]]]))
  }
  best=order(entropies,decreasing=TRUE)[1]
  ortheta[[t.or]]$pi = theta[[t]]$pi[perma[[best]]]
  ortheta[[t.or]]$w = theta[[t]]$w[perma[[best]]]
  ortheta[[t.or]]$v = theta[[t]]$v[perma[[best]]]
  ortheta[[t.or]]$sigma2 = theta[[t]]$sigma2[perma[[best]]]
  ortheta[[t.or]]$K = theta[[t]]$K
  t.or = t.or+1
  if(t.or%%100 == 0)print(t.or)
}

