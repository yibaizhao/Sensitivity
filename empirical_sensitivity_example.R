
###########################################
#Sample code for empirical sensitivity
############################################
#Author JL
#Date October 15, 2021

#Run some example code to test the empirical sensitvity functions
#States in the  model
#1=no cancer
#2=preclinical
#3=clinical

library("here")
source("simulation_functions_empirical_sensitivity.R")

#use the breast cancer example, fit to BCSC data
rateexample=matrix(c(-.002,.002,0,0,-.135,.135,0,0,0),byrow=T,nrow=3)

#Empirical sensitivity at 80% underlying sensitvity
empirical.sensitivity(sensitivity=.8,rate.matrix=rateexample,
                      start.dist=c(1- 0.003984064, 0.003984064,0),screen.times=c(0,1,2,3),
                      post.screen.lookout=1,
                      clinical.cancer.state=3,
                      pre.clinical.cancer.state=2,
                                nreps=10000)

#get empirical sensitvity values for a range of true sensitivities, simulation method
empirical_out=lapply(seq(.1,1,by=.1),empirical.sensitivity,rate.matrix=rateexample,start.dist=c(1- 0.003984064, 0.003984064,0),
                     screen.times=c(0,1,2,3),
       post.screen.lookout=1,
       clinical.cancer.state=3,
       pre.clinical.cancer.state=2,
       nreps=10000)
plot(seq(.1,1,by=.1),unlist(empirical_out),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="red",xlim=c(0,1),ylim=c(0,1))
lines(seq(0,1),seq(0,1))

#get empirical sensitivities for a range of true sensitivites, 1st screen, analytic method
empirical_out_first_screen=lapply(seq(.1,1,by=.01),empirical.sensitivity.first.screen,rate.matrix=rateexample,start.dist=c(1- 0.003984064, 0.003984064,0),
                     post.screen.lookout=1,
                     clinical.cancer.state=3,
                     pre.clinical.cancer.state=2)
lines(seq(.1,1,by=.01),unlist(empirical_out_first_screen),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="blue")
# legend("bottomright",legend=c("first screen, analytic","four screens, simulation","y=x"),col=c("blue","red","black"),lty=1,lwd=1)


#get empirical sensitivities for a range of true sensitivites, 1st screen, analytic method from Overleaf
empirical_out_first_screen_V2=lapply(seq(.1,1,by=.01),empirical.sensitivity.first.screen.V2,
                                     rate.matrix=rateexample,start.dist=c(1- 0.003984064, 0.003984064,0),
                                     post.screen.lookout=1,
                                     clinical.cancer.state=3,
                                     pre.clinical.cancer.state=2)
lines(seq(.1,1,by=.01),unlist(empirical_out_first_screen_V2),xlab=c("true sensitvity"),ylab=c("empirical sensitvity"),type="l",col="green")
legend("bottomright",legend=c("first screen, analytic V2", "first screen, analytic", "first screens, simulation","y=x"),col=c("green", "blue","red","black"),lty=1,lwd=1)

