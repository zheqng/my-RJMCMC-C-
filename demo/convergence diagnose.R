#’ ␣###␣ Convergence ␣ d i a g n o s t i c s
library('ggplot2')
theme_set( theme_minimal( ) )
library('tidyr')
library('gganimate')
library ('ggforce')
library('MASS')
library ('rprojroot')
library ('rstan')

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
plot.pair.trace<-function(tt,PI,tt.title,PI.title){
  plot(tt[,1],PI[,1],col="sienna3",
       xlim=range(tt),ylim=range(PI),
       xlab=(tt.title),ylab=(PI.title))
  points(tt[,3],PI[,3],col="steelblue")
  points(tt[,2],PI[,2],col = "gold4")
}

plot.trace<-function(tt,title){
  N = dim(tt)[1]*100
  plot(seq(1,N,by=100),tt[,1],ylim=range(tt),
       ylab=title,xlab="n",type="l",col="sienna3")
  lines(seq(1,N,by=100),tt[,2],col="gold4")
  lines(seq(1,N,by=100),tt[,3],col="steelblue")
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

plot.simulate.trace<-function(PI,w,v,sigma2){
  plot(mcmc(PI))
  plot(mcmc(w))
  # cumuplot(mcmc(w), probs = c(0.05,0.5,0.95))
  plot(mcmc(v))
  plot(mcmc(sigma2))
}


plot.simulate.panel <-function(tt,N=100,warm = 50){
  dft <- data.frame(tt[1:N, ])
df <- data.frame(id=rep(1,N),
                    iter=1:N, 
                    th1 = tt[1:N, 1],
                    th2 = tt[1:N, 2],
                    th1l = c(tt[1, 1], tt[1:(N-1), 1]),
                    th2l = c(tt[1, 2], tt[1:(N-1), 2]))
# warm <- 50
labs1 <- c('Draws', 'Steps of the sampler', '90% HPD')
ind1 <- (1:warm)*2-1
dfs <- df
dfs[ind1+1,3:4]=dfs[ind1,3:4]
p1 <- ggplot() +
  geom_point(data = dfs,
             aes(th1, th2, color ='1')) +
  geom_segment(data = df, aes(x = th1, xend = th1l, color = '2',
                                 y = th2, yend = th2l)) +
  stat_ellipse(data = dft, aes(x = X1, y = X2, color = '3'), level = 0.9) +
  coord_cartesian(xlim = range(dft$X1), ylim = c(min(dft$X2)-0.1,max(dft$X2)+0.1 )) +
  labs(x = 'theta1', y = 'theta2') +
  scale_color_manual(values = c('red', 'forestgreen','blue'), labels = labs1) +
  guides(color = guide_legend(override.aes = list(
    shape = c(16, NA, NA), linetype = c(0, 1, 1)))) +
  theme(legend.position = 'bottom', legend.title = element_blank())

p1
# animate(p1 +
# transition_reveal(id=id, along=iter) +
# shadow_trail())
}

calc.convergence.each<-function(tt,title.component){
  N = length(tt)
  # N = dim(tt)[1]
  N.half = floor(N/2)
  tt.split = cbind(tt[1:N.half],tt[(N.half+1):(2*N.half)])
# dim ( samp ) <- c ( dim ( datha.loglik ) , 1 )
# samp <- aperm ( samp , c ( 1 , 3 , 2 ) )
res <- monitor ( tt.split , probs = c ( 0.25,0.5,0.75) , digits_summary = 2 )
neff <- res [ , 'n_eff' ]
reff <-  mean ( neff / ( N  ) )
print(reff)
# aaa<-acf(tt.split,plot = FALSE)
# plot(aaa,main = title.component)
}

calc.convergence<-function(PI,w,v,sigma2){
  # op <- par(mfrow=c(4,3))
  for(i in 1:9){cat(paste("convergence analysis of",i,"th pi\n"));calc.convergence.each(PI[,i],paste("pi",i))}
  for(i in 1:9){cat(paste("convergence analysis of",i,"th w\n"));calc.convergence.each(w[,i],paste("w",i))}
  for(i in 1:9) {cat(paste("convergence analysis of",i,"th v\n"));calc.convergence.each(v[,i],paste("v",i))}
  for(i in 1:9){cat(paste("convergence analysis of",i,"th sigma2\n")); calc.convergence.each(sigma2[,i],paste("sigma2",i))}
  # op<- par(mfrow = c(1,1))
  
}

getmode <- function(tt){
  uniqv <- unique(tt)
uniqv[which.max(tabulate(match(tt,uniqv)))]
}
####################################################################
plot.simulate.trace(PI,w,v,sigma2)
calc.convergence(PI,w,v,sigma2)
plot.total.trace(PI,w,v,sigma2)

t.or = t.or-1;
warm = floor(t.or/2)
pi.mean = apply(PI[warm:t.or,],2,mean)
w.mean=apply(w[warm:t.or,],2,mean)
v.mean=apply(v[warm:t.or,],2,mean)
sigma2.mean = apply(sigma2[warm:t.or,],2,mean)

pi.mean = pi.mean[c(2,3,1)]
w.ori = c(1 ,0.5, 10)
w.mean = w.mean[c(2,3,1)]
abs(w.ori - w.mean)/w.ori*100
v.ori = c(0.2,1,0.2)
v.mean = v.mean[c(2,3,1)]
abs(v.ori - v.mean)/v.ori*100
sigma2.ori = c(0.0025,0.001,0.0005)
sigma2.mean = sigma2.mean[c(2,3,1)]
abs(sigma2.ori - sigma2.mean)/sigma2.ori

theta.mean = list( )
theta.mean$K = 3
theta.mean$w = w.mean
theta.mean$v = v.mean
theta.mean$sigma2 = sigma2.mean
theta.mean$pi = pi.mean

calc.convergence(PI[1:200,],w[1:200,],v[1:200,],sigma2[1:200,])
