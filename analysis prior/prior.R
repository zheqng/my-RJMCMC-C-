library('MCMCpack')
x = seq(0,10,length.out = 100)
density<-dinvgamma(x,0.5,0.5)
op<-par(mfrow=c(3,2))
w.3=w[w<3]
testhist=truehist(w.3)
lines(density(t(w.3),adjust=2))
 curve(dinvgamma(x,2,0.014),xlim=c(0,0.015))
 
 testhist=truehist(w)
 lines(density(t(w),adjust=2))
 curve(dinvgamma(x,2,0.014),xlim=c(0,6))
 # op<-par(mfrow=c(1,1))
 title("w")
 ##################3
 # op<-par(mfrow=c(1,2))
 # x.v = seq(0,4000,length.out = 100)
 testhist=truehist(v)
 lines(density(t(v),adjust=2))
 curve(dinvgamma(x,3,10000),xlim=c(0,10000))
 # op<-par(mfrow=c(1,1))
 title("v")
 #################
 # op<-par(mfrow=c(1,2))
 vv.1000=vv[vv<1000]
 testhist=truehist(vv.1000)
 lines(density(t(vv.1000),adjust=2))
 curve(dinvgamma(x,3,40),xlim=c(0,25))
 op<-par(mfrow=c(1,1))
 title("sigma2")
 
 testhist=truehist(vv)
 lines(density(t(vv),adjust=2))
 curve(dinvgamma(x,3,40),xlim=c(0,8000))
 op<-par(mfrow=c(1,1))
 title("sigma2")
 
 

x = seq(0,1,length.out = 100)
density<-dbeta(x,2,2)
# hist(v)
density1<-dbeta(x,1.1,1.1)
plot(x,density,type='l')
plot(x,density1,type='l',col="red")



