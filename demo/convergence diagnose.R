#’ ␣###␣ Convergence ␣ d i a g n o s t i c s
library('ggplot2')
theme_set( theme_minimal( ) )
library('tidyr')
library('gganimate')
library ('ggforce')
library('MASS')
library ('rprojroot')
library ('rstan')
library('bayesplot')
library('coda')
# w.ori = c(1,0.5,10)
# v.ori=c(0.2,1,0.2)
# sigma2.ori = c(0.0025 ,0.001 ,0.0005)
# w.ori = w.ori[c(1,3,2)]
# v.ori = v.ori[c(1,3,2)]
# sigma2.ori = sigma2.ori[c(1,3,2)]
####################################################################################
PI = NULL
w = NULL
v = NULL
sigma2 = NULL

for(i in 1:(t.or-1)){
  PI = rbind(PI,ortheta[[i]]$pi)
  w = rbind(w,ortheta[[i]]$w)
  v= rbind(v,ortheta[[i]]$v)
  sigma2 = rbind(sigma2,ortheta[[i]]$sigma2)
}
colnames(PI)<-list("pi1","pi2","pi3")
colnames(w)<-list("w1","w2","w3")
colnames(v)<-list("v1","v2","v3")
colnames(sigma2)<-list("sigma21","sigma22","sigma23")
##########################################################################################
# plot.pair.trace<-function(tt,PI,tt.title,PI.title){
#   plot(tt[,1],PI[,1],col="sienna3",
#        xlim=range(tt),ylim=range(PI),
#        xlab=(tt.title),ylab=(PI.title))
#   points(tt[,3],PI[,3],col="steelblue")
#   points(tt[,2],PI[,2],col = "gold4")
# }

plot.trace<-function(tt,title){
  N = dim(tt)[1]*100
  plot(seq(1,N,by=100),tt[,1],ylim=range(tt),
       ylab=title,xlab="n",type="l",col="sienna3")
  lines(seq(1,N,by=100),tt[,2],col="gold4")
  lines(seq(1,N,by=100),tt[,3],col="steelblue")
  legend("right",legend = c("1st","2nd ","3rd"),col = c("sienna3","gold4","steelblue"),lty = 1)
}
plot.total.trace<-function(PI,w,v,sigma2){
  op<-par(mfrow = c(2,2))
  plot.trace(PI,paste("PI"))
  plot.trace(w,paste("w"))
  plot.trace(v,paste("v"))
  plot.trace(sigma2,paste("sigma2"))
  # op<-par(mfrow=c(3,2))
  # plot.pair.trace(PI,w,paste("PI"),paste("w"))
  # plot.pair.trace(PI,v,paste("PI"),paste("v"))
  # plot.pair.trace(PI,sigma2,paste("PI"), paste("sigma2"))
  # plot.pair.trace(w,v,paste("w"),paste("v"))
  # plot.pair.trace(w,sigma2,paste("w"),paste("sigma2"))
  # plot.pair.trace(v,sigma2,paste("v"),paste("sigma2"))
  op<-par(mfrow=c(1,1))
}

# plot.simulate.trace<-function(PI,w,v,sigma2){
#   plot(mcmc(PI))
#   plot(mcmc(w))
#   # cumuplot(mcmc(w), probs = c(0.05,0.5,0.95))
#   plot(mcmc(v))
#   plot(mcmc(sigma2))
# }


calc.rhat.each<-function(tt){
  N = length(tt)
  # N = dim(tt)[1]
  warm = floor(N/2)
  N.half = floor((N-warm)/2)
  tt.split = cbind(tt[(warm+1):(warm+N.half)],tt[(warm+N.half+1):(warm+2*N.half)])
  rhat=Rhat(tt.split)
  return (rhat)
}


calc.rhat<-function(PI,w,v,sigma2){
  tmp = c(0,0,0)
  op <- par(mfrow=c(4,3))
  cat(paste("convergence analysis of"," pi\n"))
  for(i in 1:3){tmp[i]=calc.rhat.each(PI[,i])}
  print(round(tmp,digits=3))
  cat(paste("convergence analysis of"," w\n"))
  for(i in 1:3){tmp[i]=calc.rhat.each(w[,i])}
  print(round(tmp,digits=3))
  cat(paste("convergence analysis of"," v\n"))
  for(i in 1:3) {tmp[i]=calc.rhat.each(v[,i])}
  print(round(tmp,digits=3))
  cat(paste("convergence analysis of"," sigma2\n"))
  for(i in 1:3){tmp[i]= calc.rhat.each(sigma2[,i])}
  print(round(tmp,digits=3))
  op<- par(mfrow = c(1,1))
  
}

plot.rhat.trace<-function(tt,title){
  N = dim(tt)[1]*100
  plot(seq(1,N,by=100),tt[,1],ylim=range(tt),xlab = "iterations",
       ylab=title,type="l",col="sienna3")
  lines(seq(1,N,by=100),tt[,2],col="gold4")
  lines(seq(1,N,by=100),tt[,3],col="steelblue")
  legend("topright",legend = c("1st","2nd ","3rd"),col = c("sienna3","gold4","steelblue"),lty = 1)
}
plot.rhat<-function(PI,w,v,sigma2){
  op<-par(mfrow = c(2,2))
  NN = dim(w)[1]
  tmp = matrix(0,NN,3)
  
  for(j in 1:3)
  for(i in  1:NN){tmp[i,j]= calc.rhat.each(PI[1:i,j])}
  plot.rhat.trace(tmp[7:NN,],paste("rhat of PI"))
  lines((1:NN)*100,rep(1.1,NN),col="red")
  
  for(j in 1:3)
    for(i in  1:NN){tmp[i,j]= calc.rhat.each(w[1:i,j])}
  plot.rhat.trace(tmp[7:NN,],paste("rhat of w"))
  lines((1:NN)*100,rep(1.1,NN),col="red")
  
  for(j in 1:3)
    for(i in  1:NN){tmp[i,j]= calc.rhat.each(v[1:i,j])}
  plot.rhat.trace(tmp[7:NN,],paste("rhat of v"))
  lines((1:NN)*100,rep(1.1,NN),col="red")
  
  for(j in 1:3)
    for(i in  1:NN){tmp[i,j]= calc.rhat.each(sigma2[1:i,j])}
  plot.rhat.trace(tmp[7:NN,],paste("rhat of sigma2"))
  lines((1:NN)*100,rep(1.1,NN),col="red")
  

  op<-par(mfrow=c(1,1))
}



###############################################################3
N.trunc =  which(datha.K[1:1600]==3)
t.or=length(N.trunc)
t.or = t.or-1;
warm = floor(t.or/2)
pi.mean = apply(PI[warm:t.or,],2,mean)
w.mean=apply(w[warm:t.or,],2,mean)
v.mean=apply(v[warm:t.or,],2,mean)
sigma2.mean = apply(sigma2[warm:t.or,],2,mean)

w.ori = c(1 ,0.5, 10)
v.ori = c(0.2,1,0.2)
sigma2.ori = c(0.0025,0.001,0.0005)

library(combinat)
perma=permn(3)
ACC=rep(0,factorial(3))
for (i in 1:factorial(3))
  ACC[i]=sum(abs(w.ori-w.mean[perma[[i]]]))
best=order(ACC,decreasing=FALSE)[1]

pi.mean = pi.mean[perma[[best]]]
w.mean = w.mean[perma[[best]]]
abs(w.ori - w.mean)/w.ori*100
v.mean = v.mean[perma[[best]]]
abs(v.ori - v.mean)/v.ori*100
sigma2.mean = sigma2.mean[perma[[best]]]
abs(sigma2.ori - sigma2.mean)/sigma2.ori

round(pi.mean,digits=2)
round(w.mean,digits=2)
round(v.mean,digits=2)
round(sigma2.mean,digits=4)
theta.mean = list( )
theta.mean$K = 3
theta.mean$w = w.mean
theta.mean$v = v.mean
theta.mean$sigma2 = sigma2.mean
theta.mean$pi = pi.mean
####################################################################
#plot trace
# plot.simulate.trace(PI[,perma[[best]]],w[,perma[[best]]],v[,perma[[best]]],sigma2[,perma[[best]]])
# calc.convergence(PI,w,v,sigma2)
plot.total.trace(PI[,perma[[best]]],w[,perma[[best]]],v[,perma[[best]]],sigma2[,perma[[best]]])
calc.rhat(PI[,perma[[best]]],w[,perma[[best]]],v[,perma[[best]]],sigma2[,perma[[best]]])
plot.rhat(PI[,perma[[best]]],w[,perma[[best]]],v[,perma[[best]]],sigma2[,perma[[best]]])



# calc.rhat(PI[1:N.3.half,c(1,3,2)],w[1:N.3.half,c(1,3,2)],v[1:N.3.half,c(1,3,2)],sigma2[1:N.3.half,c(1,3,2)])

##########################################3

