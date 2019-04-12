#’ ␣###␣ Convergence ␣ d i a g n o s t i c s
library('ggplot2')
theme_set( theme_minimal( ) )
library('tidyr')
library('gganimate')
library ('ggforce')
library('MASS')
library ('rprojroot')
library ('rstan')

pi = NULL
w = NULL
v = NULL
sigma2 = NULL
N.start = which(datha.K==2)[1]
for(i in N.start:N){
  pi = rbind(pi,theta[[i]]$pi)
  w = rbind(w,theta[[i]]$w)
  v= rbind(v,theta[[i]]$v)
  sigma2 = rbind(sigma2,theta[[i]]$sigma2)
}
samp <- v
# dim ( samp ) <- c ( dim ( datha.loglik ) , 1 )
# samp <- aperm ( samp , c ( 1 , 3 , 2 ) )
res <- monitor ( samp , probs = c ( 0.25,0.5,0.75) , digits_summary = 2 )
neff <- res [ , 'n_eff' ]
reff <-  mean ( neff / ( N  ) )

acf(datha.loglik)
