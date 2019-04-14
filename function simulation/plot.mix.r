y.range<-function(dat){
  y.range=c(0,0)
  for(m in 1:dat$M)
  {
    range.elem<-range(dat[[m]]$y)
    y.range[1]<- min(y.range,range.elem)
    y.range[2]<-max(y.range,range.elem)
  }
  y.range
}




plot.mixgaussian<-function(dat,make.pdf = FALSE){
  # make pdf
  if(make.pdf){
    pdf('simudata.pdf', 
        width=15/2.54, height=15/2.54,
        family='GB1')
    opar <- par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0))
    on.exit(dev.off())
  } else {
    opar <- par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0))
    on.exit(par(opar))
  }
  # plot data
  mix.colors = c("red","green","blue")
  
  plot(dat[[1]]$x,dat[[1]]$y,col=mix.colors[[dat[[1]]$k]],
       type="l",ylim=y.range(dat),xlab = "",ylab="")
  for(m in 2:dat$M)
    {
      lines(dat[[m]]$x
       ,dat[[m]]$y,col=mix.colors[[dat[[m]]$k]],type="l",xlab = "",ylab="")
    }
}

y.traindata.range<-function(dat){
  y.range=c(0,0)
  for(m in 1:(dat$curve.num))
  {
    range.elem<-range(traindata$Y[m,])
    y.range[1]<- min(y.range,range.elem)
    y.range[2]<-max(y.range,range.elem)
  }
  y.range
}
plot.mixgaussian.traindata<-function(dat,make.pdf = FALSE){
  # make pdf
  if(make.pdf){
    pdf('simudata.pdf', 
        width=15/2.54, height=15/2.54,
        family='GB1')
    opar <- par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0))
    on.exit(dev.off())
  } else {
    opar <- par(mar=c(3,3,3,1), mgp=c(1.5,0.5,0))
    on.exit(par(opar))
  }
  # plot data
  mix.colors = c("red","green","blue")
  z.ori = c(1,1,1,2,2,2,3,3,3)
  op <- par(bty = 'o',oma = c(1, 1, 1, 1),mfrow=c(1,3))
  plot(t(dat$X[1,]),t(dat$Y[1,]),col=mix.colors[z.ori[1]],
       type='l',ylim=c(-2,2),xlab = "",ylab="")
  for(m in 2:3)
  {
    lines(t(dat$X[m,])
          ,t(dat$Y[m,]),col=mix.colors[z.ori[m]],type="l",xlab = "",ylab="")
  }
  plot(t(dat$X[4,]),t(dat$Y[4,]),col=mix.colors[z.ori[1]],
       type='l',ylim=c(-2,2),xlab = "",ylab="")
  for(m in 5:6)
  {
    lines(t(dat$X[m,])
          ,t(dat$Y[m,]),col=mix.colors[z.ori[1]],type="l",xlab = "",ylab="")
  }
  plot(t(dat$X[7,]),t(dat$Y[7,]),col=mix.colors[z.ori[1]],
       type='l',ylim=c(-2,2),xlab = "",ylab="")
  for(m in 8:9)
  {
    lines(t(dat$X[m,])
          ,t(dat$Y[m,]),col=mix.colors[z.ori[1]],type="l",xlab = "",ylab="")
  }
  op <- par(mfrow=c(1,1))
  
}
# plot.mixgaussian2<-function(dat){
#   df<-data.frame(x = dat[[1]]$x, y = dat[[1]]$y,type = as.character(1))
#  for(m in 2:dat$M){
#    df.new <- data.frame(x = dat[[m]]$x, y = dat[[m]]$y,type = as.character(m))
#    df <-rbind(df,df.new)
#  }
# 
#   
#   library(ggplot2)
#   ggplot(df)+geom_line(aes(x,y,colour=type))
# }
