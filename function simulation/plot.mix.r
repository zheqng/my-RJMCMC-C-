
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




plot.mixgaussian<-function(dat,step,K,make.pdf = FALSE){
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
  op<-par(mfrow=c(3,4))
  # plot(dat[[1]]$x,dat[[1]]$y,col=mix.colors[[dat[[1]]$k]],
  #      type="l",ylim=y.range(dat),xlab = "",ylab="")
  # title(paste(dat[[m]]$k))
  for(k in (1:K))
  {
    m=(k-1)*step+1;
    plot(dat[[m]]$x,dat[[m]]$y,col=mix.colors[[dat[[m]]$k]],ylim= y.range(dat),type="l",xlab = "",ylab="")
    for(j in 2:step){
      m = (k-1)*step +j;
      lines(dat[[m]]$x,dat[[m]]$y,col=mix.colors[[dat[[m]]$k]],type="l",xlab = "",ylab="")
    }
    title(paste(k))
  }
  op<-par(mfrow=c(1,1))
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
