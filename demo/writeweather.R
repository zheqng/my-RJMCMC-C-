rm(list=ls())
# setwd('/home/zheqng/src/RJMCMC-my-C/canadaweather/demo/')
setwd("/media/zheqng/Seagate Backup Plus Drive/zheqng@nwu/src/RJMCMC-my-C/simu1/canadaweather/demo")
library("fda")
library("GPFDA")
library("ggplot2")

temperature = matrix(rep(0,73*35*2),73,35*2)
tmp = matrix(rep(0,73*35),73,35)
for(m in 1:73){
  # temperature[ii+1,] 
  tmp[m,]= apply(CanadianWeather$dailyAv[((m-1)*5+1):((m-1)*5+5),,1],2,sum)
}
index = 1:73
train.index=sample(index,floor(73/2))
tmp= tmp - t(matrix( rep(apply(tmp,2,mean),73),35,73))
  temperature[,2*(1:35)-1] = 1:73
  temperature[,2*(1:35)]=tmp

time = matrix(1:73,73,35)

tt = rep(CanadianWeather$place,each=73)
canadian.temp=data.frame(time = time,
                         temperature = tmp,tt=tt)
ggplot(data= canadian.temp,aes(x=as.vector(time),
                               y=as.vector(tmp),colour = tt)) +  
  geom_line()
choose.curve = c(1,2,16,17,25,26,33,34)
choose.curves=sort(c(choose.curve*2,choose.curve*2-1))
write.table(t(temperature[,choose.curves]),file="CanadianWeather_temperature.dat",quote=F,col.name=F,row.names=F)
####
Diff.temperature = diff(temperature)
Diff.time = matrix(1:72,72,35)
Diff.tt = rep(CanadianWeather$place,each=72)
canadian.diff.temp=data.frame(Diff.time = Diff.time,
                         Diff.temperature = Diff.temperature,Diff.tt=Diff.tt)
ggplot(data= canadian.diff.temp,aes(x=as.vector(Diff.time),
                               y=as.vector(Diff.temperature),colour = Diff.tt)) +  
  geom_point()
write.table(Diff.temperature,file="CanadianWeather_diff_temperature.dat",quote=F,col.name=F,row.names=F)
####
temperature.fd= smooth.basis(time,temperature,tempfdPar)
veltemp.fd = deriv.fd(temperature.fd,1)


##plot class
color_label = c("red","blue","yellow","green")
z=c(rep(1,1,15),rep(2,1,9),rep(3,1,5),rep(2,1,3),rep(4,1,3))
# plot(1:73,temperature[,16],col = color_label[z[16]],type = "l",ylim = range(temperature))
plot(-100,-100,xlim=c(0,73),ylim=range(temperature))
for (ii in 1:(35))
  lines(1:73,tmp[,ii],col = color_label[z[ii]],type = "l")
plot(-100,-100,xlim=c(0,365),ylim=range(CanadianWeather$dailyAv[,,1]))
for (ii in 1:(35))
  lines(1:365,CanadianWeather$dailyAv[,ii,1],col = color_label[z[ii]],type = "l")
sample.index = sample(1:365,365)
sample.index = sort(sample.index)
plot(-100,-100,xlim=c(0,365),ylim=range(CanadianWeather$dailyAv[,,1]))
for (ii in 1:(35))
  lines(sample.index,CanadianWeather$dailyAv[sample.index,ii,1],col = color_label[z[ii]],type = "l")
#write data
write.table(t(temperature[,16:35]),file="CanadianWeather_temperature.dat",quote=F,col.name=F,row.names=F)

##################
theta=vector("list",35)
v=as.vector(rep(0,35))
w=as.vector(rep(0,35))
vv=as.vector(rep(0,35))

time.new = seq(1,73,length.out = 200)
m=1
plot(-100,-100,col=0,xlim=c(1,73),ylim=c(-150,120),type='l', xlab="input ",ylab="response")

for(m in 1:35){
  a<-gpr(time[,m],tmp[,m],Cov="pow.ex",trace = 0)
  pred.fit<-gppredict(a,time.new)
  plot.gp.predict(time[,m],tmp[,m],time.new,pred.fit)
  # a<-gpr(1:365,CanadianWeather$dailyAv[,m,1],Cov="pow.ex",trace = 100)
  theta[[m]]=sapply(a$hyper,FUN=exp) 
  v[m]=theta[[m]][1]
  w[m]=theta[[m]][2]
  vv[m]=theta[[m]][3]
  print(m)
}

plot.gp.predict<- function(x,y,x.new,pred.fit){
  upper=pred.fit$pred.mean+1.96*(pred.fit$pred.sd);
  lower=pred.fit$pred.mean-1.96*(pred.fit$pred.sd);
  # polygon(c(predictdata$newdata, rev(predictdata$newdata)), c(upper, rev(lower)),col = rgb(127,127,127,120, maxColorValue = 255), border = NA)
  lines(x.new,pred.fit$pred.mean,col="red",lwd=1)  
  points(x,y,col="red",lwd=1,cex = 0.5)  
}


