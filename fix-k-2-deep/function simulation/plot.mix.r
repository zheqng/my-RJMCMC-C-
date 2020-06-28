
y.range<-function(dat){
  y.range=c(0,0)
  for(m in 1:dat$curve.num)
  {
    range.elem<-range(dat[[m]]$y)
    y.range[1]<- min(y.range,range.elem)
    y.range[2]<-max(y.range,range.elem)
  }
  y.range
}




plot.mixgaussian<-function(dat,z,K,make.pdf = FALSE){
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
  mix.colors =rainbow(K)
  op<-par(mfrow=c(1,K))
  # plot(dat[[1]]$x,dat[[1]]$y,col=mix.colors[[dat[[1]]$k]],
  #      type="l",ylim=y.range(dat),xlab = "",ylab="")
  # title(paste(dat[[m]]$k))
  for(k in (1:K))
  {
    index = which(z==k)
    m = index[1]
    plot(dat[[m]]$x,dat[[m]]$y,col=mix.colors[k],xlim = c(-4,4),ylim= y.range(dat),type="l",xlab = "",ylab="")
    for(j in 2:length(index)){
      m = index[j]
      lines(dat[[m]]$x,dat[[m]]$y,col=mix.colors[k],type="l",xlab = "",ylab="")
    }
    title(paste("the",k,"th component"))
  }
  op<-par(mfrow=c(1,1))
}


