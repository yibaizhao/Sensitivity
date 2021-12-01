source("empirical_sensitivity_analytic_general_functions.R")

#######################################
# Given true sensitivity and mean sojourn time, solve for post-screen lookout time, so that empirical=true sensitivity.
#######################################
empirical.sensitivity.diff.V1 <- function(params, t1, rate.matrix, post.screen.lookout=NULL, sensitivity=NULL,
                                          clinical.cancer.state,pre.clinical.cancer.state, type){
  if(is.null(sensitivity)){
    sensitivity = params
  }
  
  if(is.null(post.screen.lookout)){
    post.screen.lookout = params
  }
  
  em_sens <- empirical.sensitivity.general(screen.start.time = t1, k = 1, 
                                           sensitivity = sensitivity,
                                           rate.matrix = rate.matrix,
                                           post.screen.lookout = post.screen.lookout,
                                           clinical.cancer.state = clinical.cancer.state, 
                                           pre.clinical.cancer.state = pre.clinical.cancer.state,
                                           method = type)

  return(abs(em_sens - sensitivity))
  
}


#######################################
# Given empirical sensitivity, true sensitivity, and screening interval, solve for sojourn time, so that empirical=true sensitivity.
#######################################
empirical.sensitivity.diff <- function(params, 
                                       k = 1,
                                       screen.start.time, 
                                       mean.sojourn.time=NULL, 
                                       post.screen.lookout=NULL, 
                                       sensitivity=NULL,
                                       clinical.cancer.state,pre.clinical.cancer.state,
                                       type){
  if(is.null(sensitivity)){
    sensitivity = params
  }
  
  if(is.null(post.screen.lookout)){
    post.screen.lookout = params
  }
  if(is.null(mean.sojourn.time)){
    mean.sojourn.time = params
  }
  rate.matrix = matrix(c(-.002,.002,0,0,-1/mean.sojourn.time,1/mean.sojourn.time,0,0,0),byrow=T,nrow=3)
  
  em_sens <- empirical.sensitivity.general(screen.start.time = screen.start.time, k = k, 
                                           sensitivity = sensitivity,
                                           rate.matrix = rate.matrix,
                                           post.screen.lookout = post.screen.lookout,
                                           clinical.cancer.state = 3, 
                                           pre.clinical.cancer.state = 2,
                                           method = type)
  
  c(abs(em_sens - sensitivity))
  
}


