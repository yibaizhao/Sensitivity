library(msm)

########################################################################################
sim.ctmc<-function(start.state,rate.matrix, end.time,start.time=0,absorbing.state=0){
  ########################################################################################
  #Author: JL
  #This function simulates from a homogeneous CTMC characterized by rate.matrix
  #INPUTS: rate.matrix=rate matrix, start.state=starting state for CTMC, end.time=time
  #         to stop data simulations; start.time= time to start data simulations
  #        absorbing.state=possible multiple absorbing states in the model
  #OUTPUTS: a list with two objects: "times"=transition times and "states"=transition states
  #
  #
  ##########################################################################################
  #browser()
  state.space<-seq(1:dim(rate.matrix)[1])
  size<-dim(rate.matrix)[1]
  cur.state<-start.state
  times<-vector()
  states<-vector()
  
  times[1]<-start.time
  states[1]<-cur.state
  cur.time<-start.time
  k<-2
  while((cur.time<end.time)&(!(cur.state%in% absorbing.state))){
    exp.rate<-(-1)*rate.matrix[cur.state,cur.state]
    if(exp.rate==0){
      cur.time<-end.time
    } else{
      cur.time<-cur.time+rexp(n=1,rate=exp.rate)
      if(cur.time<end.time){
        times[k]<-cur.time
        if(size==2){
          cur.state=as.numeric(state.space[-cur.state])
        }else{
          cur.state<-sample(state.space[-cur.state],size=1,prob=rate.matrix[cur.state,-cur.state])
        }
        states[k]<-cur.state
        k<-k+1
        
      }
    }
  }
  return.list<-list(times,states)
  names(return.list)<-c("times","states")
  return(return.list)
}

discrete.ctmc<-function(ctmc.times,ctmc.states,obs.times){
  ###################################################################################################
  #Author: 
  # This function gets the state of a CTMC at different discrete observation times
  #INPUTS: ctmc.times=the transition times for the CTMC
  #        ctmc.states=the states at each of the transition times
  #        obs.times = the discrete observation times
  #
  #OUTPUTS: a dataframe with two columns: obs.times= observation times, states=value of state at obs.times
  #WARNING: if the observation times are outside of max and min transition time, then
  #         the state is assumed to be unchanged from the closest recorded transition time
  #################################################################################################
  out<- data.frame(approx(x=ctmc.times,y=ctmc.states,xout=obs.times,rule=2,f=0,method="constant"))
  colnames(out)<-c("obs.times","states")
  return(out)
}

get.observed.datapoint<-function(underlying.state,emission.matrix){
  ###################################################################################################
  #Author: 
  # This function gets an observed data point in a HMM based on an underlying state an emission matrix
  #INPUTS: underlying.state = unobserved underlying state in HMM,
  #        emmision.matrix=a matrix with the emission probablities.
  #        the ith row corresponds to the hidden value X(t)=i, and the kth column to O(t)=k|X(t)=i
  #        thus the rows sum to 1, and k columns correspond to the k possible observed states
  #OUTPUTS: the observed data point
  #################################################################################################
  states<-seq(1:dim(emission.matrix)[2])
  probs<-emission.matrix[underlying.state,]
  # browser()
  sample(x=states,size=1,prob=probs)
}


observed.data.hmm<-function(obs.times,underlying.states,emission.matrix){
  ###################################################################################################
  #Author JL
  # This function the observed states in an HMM for multiple observation times
  # INPUTS: obs.times=a vector of observatio11n times; underlying states=corresponding underlying states
  #         emission.matrix=emission matrix for observed data
  #These functins get the observed data point based on an underlying state and a emission matrix
  #
  #
  #################################################################################################
  obs.data<-sapply(underlying.states,FUN="get.observed.datapoint",emission.matrix)
  # browser()
  out<-data.frame(obs.times, obs.data)
  colnames(out)<-c("obs.times","obs.data")
  return(out)
}


get.dx=function(sensitivity,rate.matrix,start.dist,screen.times,post.screen.lookout,clinical.cancer.state,pre.clinical.cancer.state){
  ###################################################################################################
  #Author JL
  # This function ascertains whether a person is screen or interval detected during a series of screens
  # INPUTS: rate.matrix=transition matrix of cancer model
  #         start.dist=starting probability distribution for cancer model
  #         sensitivity of the test for detecting pre-clinical cancer
  #         screen.times=times of screening
  #         post.screen.lookout = duration of time to look out for interval cancer after last screen
  #         clinical.cancer.state= state corresponding to clinical cancer
  #         pre.clinical.cancer.state=state corresponding to pre-clinical cancer
  # OUTPUTS: result=1 if no cancer discovered, 2 if screen detected, 3=interval detected
  #
  #################################################################################################
  
    nstates=dim(rate.matrix)[1]
    thetimes=c(screen.times,max(screen.times)+post.screen.lookout)
    
    emission.mat=diag(x=1,nrow=nstates,ncol=nstates)
    emission.mat[pre.clinical.cancer.state,pre.clinical.cancer.state]=sensitivity
    emission.mat[pre.clinical.cancer.state,1]=1-sensitivity
    the.start.state=sample(1:nstates,size=1,prob=start.dist)
    outctmc=sim.ctmc(start.state=the.start.state,rate.matrix=rate.matrix, 
                   end.time=(max(screen.times)+1),start.time=0,absorbing.state=clinical.cancer.state)
    discreteout=discrete.ctmc(ctmc.times = outctmc$times,ctmc.states=outctmc$states,obs.times=thetimes)
  
    obsout=observed.data.hmm(obs.times=discreteout$obs.times,underlying.states=discreteout$states,emission.matrix = emission.mat)
    #obtain 1=no cancer detected, 2=screen, 3=interval
    result=1
    result[clinical.cancer.state%in%obsout$obs.data]=3
    result[pre.clinical.cancer.state%in%head(obsout$obs.data,-1)]=2
   # browser()
    return(result)
}


###################################################################################################
#Author JL
# This function ascertains the empirical.sensitvity after a series of screens
#empirical sensitvity is defined as the fraction of all cancers that are screeen detected=screen/(screen +interval)
# INPUTS: rate.matrix=transition matrix of cancer model
#         start.dist=starting probability distribution for cancer model
#         sensitivity of the test for detecting pre-clinical cancer
#         screen.times=times of screening
#         post.screen.lookout = duration of time to look out for interval cancer after last screen
#         clinical.cancer.state= state corresponding to clinical cancer
#         pre.clinical.cancer.state=state corresponding to pre-clinical cancer
# OUTPUTS: result=1 if no cancer discovered, 2 if screen detected, 3=interval detected
#
#################################################################################################

empirical.sensitivity<-function(sensitivity,rate.matrix,start.dist,screen.times,post.screen.lookout,clinical.cancer.state,pre.clinical.cancer.state,
                                nreps){
  
  out=replicate(n=nreps,get.dx(sensitivity=sensitivity,rate.matrix=rate.matrix,start.dist=start.dist,screen.times=screen.times,
                 clinical.cancer.state=clinical.cancer.state,
                pre.clinical.cancer.state=pre.clinical.cancer.state,post.screen.lookout),simplify=T)
  
  empirical_sensitivity=sum(I(out==2))/(sum(I(out==2))+sum(I(out==3)))
  return(empirical_sensitivity)
}


###################################################################################################
#Author JL
# This function ascertains the empirical.sensitvity after the first screen, using an analytic approach
#empirical sensitvity is defined as the fraction of all cancers that are screeen detected=screen/(screen +interval)
# INPUTS: rate.matrix=transition matrix of cancer model
#         start.dist=starting probability distribution for cancer model
#         sensitivity of the test for detecting pre-clinical cancer
#         screen.times=times of screening
#         post.screen.lookout = duration of time to look out for interval cancer after last screen
#         clinical.cancer.state= state corresponding to clinical cancer
#         pre.clinical.cancer.state=state corresponding to pre-clinical cancer
# OUTPUTS: result=1 if no cancer discovered, 2 if screen detected, 3=interval detected
#
#################################################################################################
empirical.sensitivity.first.screen<-function(sensitivity,rate.matrix,start.dist,post.screen.lookout,
                                   clinical.cancer.state,pre.clinical.cancer.state){
  
  prob.mat=MatrixExp(mat=rate.matrix,t=post.screen.lookout)
  num=start.dist[pre.clinical.cancer.state]*sensitivity
  denom=start.dist[pre.clinical.cancer.state]*sensitivity+start.dist[1]*prob.mat[1,clinical.cancer.state]+
            start.dist[pre.clinical.cancer.state]*(1-sensitivity)*prob.mat[pre.clinical.cancer.state,clinical.cancer.state]
  return(num/denom)
  
  
}

empirical.sensitivity.first.screen.YZ <- function(sensitivity,rate.matrix,start.dist,post.screen.lookout,
                                             clinical.cancer.state,pre.clinical.cancer.state){
  
  screen.detect = start.dist[pre.clinical.cancer.state]*sensitivity
  lambda1 = rate.matrix[1, pre.clinical.cancer.state]
  lambda2 = rate.matrix[2, clinical.cancer.state]
  interval.cancer = 
    start.dist[1]*integrate(fxn.f.g, 
                            1, 1+post.screen.lookout, 
                            lambda1, lambda2, 
                            t=1+post.screen.lookout,
                            fun = "lambda1*exp(-lambda1*x)*lambda2*exp(-lambda2*(t-x))")$value +
    start.dist[pre.clinical.cancer.state]*(1-sensitivity)*
    integrate(fxn.f.g, 
              0, 1, 
              lambda1, lambda2, 
              t=1+post.screen.lookout,
              fun = "lambda1*exp(-lambda1*x)*lambda2*exp(-lambda2*(t-x))")$value
  
  return(screen.detect / (screen.detect + interval.cancer))
  
}

fxn.f.g <- function(x, lambda1, lambda2, t, fun = "exp(-lambda1*x)*exp(-lambda2*(t-x))"){
  fun <- parse(text=fun)
  eval(fun, envir=list(x=x, lambda1=lambda1, lambda2=lambda2))
}

empirical.sensitivity.first.screen.V2 <- function(sensitivity,rate.matrix,start.dist,post.screen.lookout,
                                                  clinical.cancer.state,pre.clinical.cancer.state){
  
  lambda1 = rate.matrix[1, pre.clinical.cancer.state]
  lambda2 = rate.matrix[2, clinical.cancer.state]
  
  screen.detect = sensitivity * integrate(fxn.f.g, 
                                          0, 1, 
                                          lambda1=lambda1, lambda2=lambda2, 
                                          t=1,
                                          fun = "lambda1*exp(-lambda1*x)*exp(-lambda2*(t-x))")$value
  interval.cancer = 
    (1-sensitivity) * integrate(fxn.f.g, 
                                0, 1, 
                                lambda1=lambda1, lambda2=lambda2, 
                                t=1,
                                fun = "lambda1*exp(-lambda1*x)*( exp(-lambda2*(t-x)) - exp(-lambda2*(t+1-x)) )")$value +
    integrate(fxn.f.g, 
              1, 2, 
              lambda1=lambda1, lambda2=lambda2, 
              t=1+post.screen.lookout,
              fun = "lambda1*exp(-lambda1*x)*(1 - exp(-lambda2*(t-x)))")$value
  
  return(screen.detect / (screen.detect + interval.cancer))
  
}

# start.dist <- exp(-rate.matrix[1,pre.clinical.cancer.state] * 1)
# integrate(fxn.f.g, 
#           1, 2, 
#           lambda1=0.1, lambda2=0.2, 
#           t=2,
#           fun = "lambda1*exp(-lambda1*x)*(1 - exp(-lambda2*(t-x)))")$value
# 
# InnerIntegral = Vectorize(function(t) { integrate(fxn.f.g, 1, t,
#                                                   lambda1=0.1, lambda2=0.2, 
#                                                   t=t,
#                                                   fun = "lambda1*exp(-lambda1*x)*lambda2*exp(-lambda2*(t-x))")$value})
# integrate(InnerIntegral , 1, 2)$value

