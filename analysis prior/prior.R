library('MCMCpack')
x = seq(0,10,length.out = 100)
density<-dinvgamma(x,0.5,0.5)
op<-par(mfrow=c(3,2))
testhist=truehist(w,xlim=c(0.95,1.02))
lines(density(t(w),adjust=2))
 curve(dinvgamma(x,1,1.5),xlim=c(0,4))
 # op<-par(mfrow=c(1,1))
 title("w")
 ##################3
 # op<-par(mfrow=c(1,2))
 testhist=truehist(v,xlim=c(0.95,1.02))
 lines(density(t(v),adjust=2))
 curve(dinvgamma(x,1,1.5),xlim=c(0,4))
 # op<-par(mfrow=c(1,1))
 title("v")
 #################
 # op<-par(mfrow=c(1,2))
 testhist=truehist(vv,xlim=c(0,7))
 lines(density(t(vv),adjust=2))
 curve(dinvgamma(x,2,6),xlim=c(0,7))
 op<-par(mfrow=c(1,1))
 title("sigma2")
 
 

x = seq(0,1,length.out = 100)
density<-dbeta(x,2,2)
# hist(v)
density1<-dbeta(x,1.1,1.1)
plot(x,density,type='l')
plot(x,density1,type='l',col="red")



