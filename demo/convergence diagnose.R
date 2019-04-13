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
N.start = which(datha.K==3)[1]
for(i in N.start:N){
  pi = rbind(pi,theta[[i]]$pi)
  w = rbind(w,theta[[i]]$w)
  v= rbind(v,theta[[i]]$v)
  sigma2 = rbind(sigma2,theta[[i]]$sigma2)
}

plot.simulate.panel <-function(tt,N=100,warm = 50){
  dft <- data.frame(tt[1:N, ])
df <- data.frame(id=rep(1,N),
                    iter=1:N, 
                    th1 = tt[1:N, 1],
                    th2 = tt[1:N, 2],
                    th1l = c(tt[1, 1], tt[1:(N-1), 1]),
                    th2l = c(tt[1, 2], tt[1:(N-1), 2]))
# warm <- 50
labs1 <- c('Draws', 'Steps of the sampler', '90% HPD')
ind1 <- (1:warm)*2-1
dfs <- df
dfs[ind1+1,3:4]=dfs[ind1,3:4]
p1 <- ggplot() +
  geom_point(data = dfs,
             aes(th1, th2, color ='1')) +
  geom_segment(data = df, aes(x = th1, xend = th1l, color = '2',
                                 y = th2, yend = th2l)) +
  stat_ellipse(data = dft, aes(x = X1, y = X2, color = '3'), level = 0.9) +
  coord_cartesian(xlim = c(1.52,1.6), ylim = c(1.1,1.3)) +
  labs(x = 'theta1', y = 'theta2') +
  scale_color_manual(values = c('red', 'forestgreen','blue'), labels = labs1) +
  guides(color = guide_legend(override.aes = list(
    shape = c(16, NA, NA), linetype = c(0, 1, 1)))) +
  theme(legend.position = 'bottom', legend.title = element_blank())

p1
}

samp <- v
# dim ( samp ) <- c ( dim ( datha.loglik ) , 1 )
# samp <- aperm ( samp , c ( 1 , 3 , 2 ) )
res <- monitor ( samp , probs = c ( 0.25,0.5,0.75) , digits_summary = 2 )
neff <- res [ , 'n_eff' ]
reff <-  mean ( neff / ( N  ) )

acf(datha.loglik)
