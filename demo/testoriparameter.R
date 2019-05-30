rmse2.m = rep(0,9)
corr.m = rep(0,9)
pi.ori = c(1/3,1/3,1/3)
z.ori = c(0,0,0,1,1,1,2,2,2)
w.ori = c(1,0.5,10)
v.ori=c(0.2,1,0.2)
sigma2.ori = c(0.0025 ,0.001 ,0.0005)
theta <-list("K","pi","w","v","sigma2")
theta$K = 3
theta$pi = pi.ori
theta$w = w.ori
theta$v = v.ori
theta$sigma2 = sigma2.ori
op <- par(mfrow=c(3,3))

for(m in 1:9){
  k = z.ori[m]+1
  predictdata<-predict.gp(traindata$X[m,],traindata$Y[m,],testdata$X[m,],testdata$Y[m,],theta,k)
  rmse2.m[m]<- rmse2(testdata$Y[m,],predictdata)
  # corr.m[m] <- predictdata$correlation
  # plot.gp.i.predict(traindata,predictdata,m)
  plot.gp.predict(traindata,testdata,predictdata,m)
  title(paste("original parameter ","predict",i=m,"curve"))
  # plot.gp.predict(traindata,testdata,predictdata,1)
  # plot.gp(traindata ,testdata )
}
op<- par(mfrow = c(1,1))
Nm = length(testdata$Y[m,])
rmse = sqrt(mean(rmse2.m)/(Nm))
rmse.each = rep(0,3)
rmse.each[1] = sqrt(mean(rmse2.m[1:3])/(Nm))
rmse.each[2] = sqrt(mean(rmse2.m[4:6])/(Nm))
rmse.each[3] = sqrt(mean(rmse2.m[7:9])/(Nm))

corr.each = rep(0,3)
corr.each[1] = mean(corr.m[1:3])
corr.each[2] = mean(corr.m[4:6])
corr.each[3] = mean(corr.m[7:9])
log.P =   print.posterior(traindata,theta[[N]])
