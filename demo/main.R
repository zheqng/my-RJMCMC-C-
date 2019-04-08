rm(list = ls())
library('GPFDA')
library('MASS')
source('predict.r')
source('readsts.R')
for(i in 1:N){
  k = z[i,4]+1
    predictdata<-predict.gp(traindata$X[4,],traindata$Y[4,],testdata$X[4,],theta[[i]],k)
    i=9;plot.gp.i.predict(traindata,predictdata,i)
    title(paste(N,"th iter ","predict",i,"curve"))
    # plot.gp.predict(traindata,testdata,predictdata,1)
    # plot.gp(traindata ,testdata )
}
