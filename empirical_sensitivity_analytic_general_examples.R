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
source("empirical_sensitivity_analytic_general_functions.R")
source("simulation_functions_empirical_sensitivity_expanded.R")

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
         screen.start.time = 3, k = 1, #start.time=NULL, end.time=NULL,
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
                                screen.times=c(2,3,4,5,6),
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


