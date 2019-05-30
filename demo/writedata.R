rm(list=ls())
setwd('/home/zheqng/src/RJMCMC-my-C/Epileptic/demo/')
# setwd("/media/zheqng/Seagate Backup Plus Drive/zheqng@nwu/src/RJMCMC-my-C/simu1/data/")
aaa<-read.csv('data.csv')
data=aaa[,-c(1,180)]
data = data -matrix(rep(apply(data,1,mean),178),11500,178)
label = aaa[,180]
color_label = rainbow(5)
#############################################
Mean = apply(data,1,mean)
boxplot(Mean~label)
################################
# z=c(rep(1,1,15),rep(2,1,9),rep(3,1,5),rep(2,1,3),rep(4,1,3))
# plot(1:73,temperature[,16],col = color_label[z[16]],type = "l",ylim = range(temperature))
plot(-100,-100,xlim=c(0,178),ylim=range(data))
for (ii in 1:11500)
  lines(1:178,data[ii,],col = color_label[label[ii]],type = "l")
index=1:11500
train.index = sort(sample(1:178,50))
# train.index = colnames(data[,train.index])
# traindata=matrix(rep(0,35*50),35,50)
traindata=NULL
train.index = t(as.matrix(train.index))
colnames(train.index) = colnames(data[1,train.index])
# train.index = t(train.index)
# op<-par(mfrow=c(3,2))
plot(-100,-100,xlim=c(0,178),ylim=range(data))
for(k in 1:5){
  index.k = index[label[index]==k]
  # plot(-100,-100,xlim=c(0,178),ylim=range(data[index.k[1:7],]))
  for (ii in index.k[1:2]){
    lines(train.index,data[ii,train.index],col = color_label[label[ii]],type = "l")
    traindata = rbind(traindata,train.index)
    traindata=rbind(traindata,data[ii,train.index])
  }
}
# op<-par(mfrow=c(1,1))
write.table(round(traindata,digits=4),file="data.dat",quote=F,col.name=F,row.names=F)
