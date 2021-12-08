###########################################
# Validation of analytic methods
############################################
#Author Yibai Zhao
#Date November 15, 2021

#Run some example code to test the general analytic functions
#States in the  model
#1=no cancer
#2=preclinical
#3=clinical

library("here")
library("ggplot2")
library(reshape2)
source("empirical_sensitivity_analytic_general_functions.R")
source("simulation_functions_empirical_sensitivity_expanded.R")
# source("empirical_sensitivity_analytic_timevarying_sens_general_functions.R")

rateexample=matrix(c(-.002,.002,0,0,-.16,.16,0,0,0),byrow=T,nrow=3)

##################
# screen-by-screen 
##################
# k=1
## simulation method
empirical_out_simulation_V1=lapply(seq(.1,1,by=.1),empirical.sensitivity.simulation,rate.matrix=rateexample,
                                start.dist=c(1,0,0),
                                screen.times=c(3),
                                post.screen.lookout=1,
                                clinical.cancer.state=3,
                                pre.clinical.cancer.state=2,
                                nreps=10000)
plot(seq(.1,1,by=.1),unlist(empirical_out_simulation_V1),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="red",xlim=c(0,1),
     ylim=c(0,1), 
     main = "Empirical versus true sensitivity\nScreen by screen at time=3")
lines(seq(0,1),seq(0,1))


## analytic method
empirical_analytic_out <-
  lapply(seq(.1,1,by=.01), 
         empirical.sensitivity.general,
         screen.start.time = 1, k = 3, #start.time=NULL, end.time=NULL,
         rate.matrix = rateexample,
         post.screen.lookout = 1,
         clinical.cancer.state = 3, 
         pre.clinical.cancer.state = 2,
         method='Single')
lines(seq(.1,1,by=.01),unlist(empirical_analytic_out),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="blue")
legend("bottomright",legend=c("analytic", "simulation","y=x"),col=c("blue","red","black"),lty=1,lwd=1)

##################
# Series of screen 
##################
# k=4
## simulation method
empirical_out_simulation_V2=lapply(seq(.1,1,by=.1),empirical.sensitivity.simulation,rate.matrix=rateexample,
                                start.dist=c(1,0,0),
                                screen.times=c(1,2,3,4),
                                post.screen.lookout=1,
                                clinical.cancer.state=3,
                                pre.clinical.cancer.state=2,
                                nreps=10000)
plot(seq(.1,1,by=.1),unlist(empirical_out_simulation_V2),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="red",xlim=c(0,1),
     ylim=c(0,1), 
     main = "Empirical versus true sensitivity\n Across screen time=1,2,3,4")
lines(seq(0,1),seq(0,1))


## analytic method
empirical_analytic_out <-
  lapply(seq(.1,1,by=.01), 
         empirical.sensitivity.general,
         screen.start.time = 2, k = 5, #start.time=NULL, end.time=NULL,
         rate.matrix = rateexample,
         post.screen.lookout = 1,
         clinical.cancer.state = 3, 
         pre.clinical.cancer.state = 2,
         method='All')
lines(seq(.1,1,by=.01),unlist(empirical_analytic_out),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="blue")
legend("bottomright",legend=c("analytic", "simulation","y=x"),col=c("blue","red","black"),lty=1,lwd=1)

###########
# Plots of empirical sensitity(t_k) for k=1,…,n
###########
mst <- 1.5
post.screen.lookout <- 1
# screen-by-screen 
rateexample=matrix(c(-.002,.002,0,0,-1/mst,1/mst,0,0,0),byrow=T,nrow=3)
empirical_analytic_out_k1 <-
  sapply(seq(1, 10, by = 2), function(k){
    sapply(seq(.1,1,by=.01), 
           empirical.sensitivity.general,
           screen.start.time = 1, k = k, #start.time=NULL, end.time=NULL,
           rate.matrix = rateexample,
           post.screen.lookout = post.screen.lookout,
           clinical.cancer.state = 3, 
           pre.clinical.cancer.state = 2,
           method='Single')
  })
empirical_analytic_out_k1_long <- as.data.frame(melt(empirical_analytic_out_k1))
names(empirical_analytic_out_k1_long) <- c('No.', 'k', 'empirical_sensitivity')
empirical_analytic_out_k1_long$true_sensitivity <- rep(seq(.1,1,by=.01), 5)
empirical_analytic_out_k1_long$type <- 'single screen'

# Series of screen 
empirical_analytic_out_k2 <-
  sapply(seq(1, 10, by = 2), function(k){
    sapply(seq(.1,1,by=.01), 
           empirical.sensitivity.general,
           screen.start.time = 1, k = k, #start.time=NULL, end.time=NULL,
           rate.matrix = rateexample,
           post.screen.lookout = post.screen.lookout,
           clinical.cancer.state = 3, 
           pre.clinical.cancer.state = 2,
           method='All')
  })
empirical_analytic_out_k2_long <- as.data.frame(melt(empirical_analytic_out_k2))
names(empirical_analytic_out_k2_long) <- c('No.', 'k', 'empirical_sensitivity')
empirical_analytic_out_k2_long$true_sensitivity <- rep(seq(.1,1,by=.01), 5)
empirical_analytic_out_k2_long$type <- 'series of screen'

# merge two scenarior
empirical_analytic_out_k <- rbind(empirical_analytic_out_k1_long, empirical_analytic_out_k2_long)
empirical_analytic_out_k %>% 
  ggplot(aes(x = true_sensitivity, y = empirical_sensitivity, color = k, group = k)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, lty = 'dashed') +
  xlim(0, 1) + ylim(0,1) +
  facet_grid(.~type) +
  scale_color_continuous(name = "kth screen") +
  ggtitle(paste0('MST=', mst, '; lookout time=',
                 post.screen.lookout))

diff <- (empirical_analytic_out_k1_long$empirical_sensitivity-empirical_analytic_out_k2_long$empirical_sensitivity)
dt_diff <- data.frame(diff = diff, 
                      true_sensitivity = empirical_analytic_out_k1_long$true_sensitivity, 
                      k = empirical_analytic_out_k1_long$k)
ggplot(data = dt_diff, aes(x = true_sensitivity, y = diff, color = k, group = k)) + 
  geom_line() +
  ylab('single screen - series of screen') +
  scale_color_continuous(name = "kth screen")
  
###########
# sojourn time follows Weibull distribution
###########
mst = 6
shape_incs = 2
scale_incs <- optim(par = 0.1, fn = find_paramsF, method = 'L-BFGS-B', lower = 0, upper = Inf, shape = shape_incs, mst = mst)$par
shape_decs = 0.5
scale_decs <- optim(par = 0.1, fn = find_paramsF, method = 'L-BFGS-B', lower = 0, upper = Inf, shape = shape_decs, mst = mst)$par
# shape = st_params[1]
# scale = st_params[2]
# mean(rweibull(100, shape = shape, scale = scale))
# scale*gamma(1+1/shape)

# screen-by-screen 
# k=1
## simulation method
rateexample=matrix(c(-.002,.002,0,0,-1/mst,1/mst,0,0,0),byrow=T,nrow=3)
empirical_out_simulation=lapply(seq(.1,1,by=.1),empirical.sensitivity.simulation,rate.matrix=rateexample,
                                   start.dist=c(1,0,0),
                                   screen.times=c(3),
                                   post.screen.lookout=1,
                                   clinical.cancer.state=3,
                                   pre.clinical.cancer.state=2,
                                   nreps=10000)
plot(seq(.1,1,by=.1),unlist(empirical_out_simulation),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="red",xlim=c(0,1),
     ylim=c(0,1), 
     main = "Empirical versus true sensitivity\nScreen by screen at time=3")
lines(seq(0,1),seq(0,1))

empirical_analytic_exp_out <-
  lapply(seq(.1,1,by=.01), 
         empirical.sensitivity.general,
         screen.start.time = 3, k = 1, #start.time=NULL, end.time=NULL,
         pre_onset_params = 0.002, st_params = 1/mst, funs = 'exp',
         post.screen.lookout = 1,
         method='Single')
lines(seq(.1,1,by=.01),unlist(empirical_analytic_exp_out),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="blue")

empirical_analytic_weibull_out <-
  lapply(seq(.1,1,by=.01), 
         empirical.sensitivity.general,
         screen.start.time = 3, k = 1, #start.time=NULL, end.time=NULL,
         pre_onset_params = 0.002, st_params = c(shape_incs, scale_incs), funs = 'weibull',
         post.screen.lookout = 1,
         method='Single')
lines(seq(.1,1,by=.01),unlist(empirical_analytic_weibull_out),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="green")

empirical_analytic_weibull_out2 <-
  lapply(seq(.1,1,by=.01), 
         empirical.sensitivity.general,
         screen.start.time = 3, k = 1, #start.time=NULL, end.time=NULL,
         pre_onset_params = 0.002, st_params = c(shape_decs, scale_decs), funs = 'weibull',
         post.screen.lookout = 1,
         method='Single')
lines(seq(.1,1,by=.01),unlist(empirical_analytic_weibull_out2),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="orange")
legend("bottomright",legend=c("simulation_exp", "analytic_exp", "analytic_weibull_incs", "analytic_weibull_decs", "y=x"),col=c("red","blue",'green',"orange", "black"),lty=1,lwd=1)


# Series of screen 
# k=4
## simulation method
empirical_out_simulation_V2=lapply(seq(.1,1,by=.1),empirical.sensitivity.simulation,rate.matrix=rateexample,
                                   start.dist=c(1,0,0),
                                   screen.times=c(1,2,3,4),
                                   post.screen.lookout=1,
                                   clinical.cancer.state=3,
                                   pre.clinical.cancer.state=2,
                                   nreps=10000)
plot(seq(.1,1,by=.1),unlist(empirical_out_simulation_V2),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="red",xlim=c(0,1),
     ylim=c(0,1), 
     main = "Empirical versus true sensitivity\n Across screen time=1,2,3,4")
lines(seq(0,1),seq(0,1))


## analytic method
empirical_analytic_exp_out <-
  lapply(seq(.1,1,by=.01), 
         empirical.sensitivity.general,
         screen.start.time = 1, k = 4, #start.time=NULL, end.time=NULL,
         pre_onset_params = 0.002, st_params = 1/mst, funs = 'exp',
         post.screen.lookout = 1,
         method='All')
lines(seq(.1,1,by=.01),unlist(empirical_analytic_exp_out),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="blue")


empirical_analytic_weibull_out <-
  lapply(seq(.1,1,by=.01), 
         empirical.sensitivity.general,
         screen.start.time = 3, k = 1, #start.time=NULL, end.time=NULL,
         pre_onset_params = 0.002, st_params = c(shape_incs, scale_incs), funs = 'weibull',
         post.screen.lookout = 1,
         method='All')
lines(seq(.1,1,by=.01),unlist(empirical_analytic_weibull_out),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="green")

empirical_analytic_weibull_out2 <-
  lapply(seq(.1,1,by=.01), 
         empirical.sensitivity.general,
         screen.start.time = 3, k = 1, #start.time=NULL, end.time=NULL,
         pre_onset_params = 0.002, st_params = c(shape_decs, scale_decs), funs = 'weibull',
         post.screen.lookout = 1,
         method='All')
lines(seq(.1,1,by=.01),unlist(empirical_analytic_weibull_out2),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="orange")
legend("bottomright",legend=c("simulation_exp", "analytic_exp", "analytic_weibull_incs", "analytic_weibull_decs", "y=x"),col=c("red","blue",'green',"orange", "black"),lty=1,lwd=1)

###########
# varing sensitivity over sojourn time
###########
# screen-by-screen 
# k=1
## simulation method
empirical_out_simulation_V1=lapply(seq(.1,1,by=.1),empirical.sensitivity.simulation,rate.matrix=rateexample,
                                   start.dist=c(1,0,0),
                                   screen.times=c(3),
                                   post.screen.lookout=1,
                                   clinical.cancer.state=3,
                                   pre.clinical.cancer.state=2,
                                   nreps=10000)
plot(seq(.1,1,by=.1),unlist(empirical_out_simulation_V1),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="red",xlim=c(0,1),
     ylim=c(0,1), 
     main = "Empirical versus true sensitivity\nScreen by screen at time=3")
lines(seq(0,1),seq(0,1))


## analytic method
empirical_analytic_out <-
  lapply(seq(.1,1,by=.01), 
         empirical.sensitivity.general.V2,
         screen.start.time = 3, k = 1, #start.time=NULL, end.time=NULL,
         sens.time.varying = FALSE,
         pre_onset_params = 0.002, st_params = 0.16, funs = 'exp',
         post.screen.lookout = 1,
         method='Single')
lines(seq(.1,1,by=.01),unlist(empirical_analytic_out),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="blue")
legend("bottomright",legend=c("analytic", "simulation","y=x"),col=c("blue","red","black"),lty=1,lwd=1)


## analytic method and time-varying
true_sens_time_varying <- sensF(seq(0,10,length.out=100), sensitivity=NULL, sens.time.varying=TRUE, sens_params=c(-0.5,0.5), dist = 'linear')
plot(seq(0,10,length.out=100),true_sens_time_varying,xlab=c("time"),ylab=c("sensitvity"),type="l",col="red",xlim=c(0,10),
     ylim=c(0,1), 
     main = "Empirical versus true sensitivity varying over time")
empirical_time_varying_analytic_out <-
  sapply(seq(1,10,1), function(k)
         empirical.sensitivity.general.V2(
         screen.start.time = 1, k=k, #start.time=NULL, end.time=NULL,
         sensitivity = NULL,
         sens.time.varying = TRUE,
         pre_onset_params = 0.002, st_params = 0.16, sens_params = c(-0.5, 0.5), funs = 'exp',
         post.screen.lookout = 1,
         method='Single'))
lines(seq(1,10,1),unlist(empirical_time_varying_analytic_out),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="blue")
legend("bottomright",legend=c("empirical sens with time varying", "true sensitivity"),col=c("blue","red"),lty=1,lwd=1)

