source("empirical_sensitivity_analytic_general_functions.R")


#######################################
# Find parameters for mst follows weibull distribution
#######################################
find_paramsF <- function(params, mst){
  shape = params[1]
  scale = params[2]
  abs(scale*gamma(1+1/shape) - mst)
}

#######################################
# Given true sensitivity and mean sojourn time, solve for post-screen lookout time, so that empirical=true sensitivity.
#######################################
empirical.sensitivity.diff.V1 <- function(params, t1, pre_onset_params, st_params, funs, 
                                          post.screen.lookout=NULL, sensitivity=NULL,
                                          type){
  if(is.null(sensitivity)){
    sensitivity = params
  }
  
  if(is.null(post.screen.lookout)){
    post.screen.lookout = params
  }
  
  em_sens <- empirical.sensitivity.general(screen.start.time = t1, k = 1, 
                                           sensitivity = sensitivity,
                                           pre_onset_params, st_params, funs,                                           post.screen.lookout = post.screen.lookout,
                                           method = type)

  return(abs(em_sens - sensitivity))
  
}


#######################################
# Given empirical sensitivity, true sensitivity, and screening interval, solve for sojourn time, so that empirical=true sensitivity.
#######################################
empirical.sensitivity.diff <- function(params, 
                                       k = 1,
                                       screen.start.time, 
                                       pre_onset_params, st_params, pi = 0, funs, 
                                       post.screen.lookout=NULL, 
                                       sensitivity=NULL,
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

  em_sens <- empirical.sensitivity.general(screen.start.time = screen.start.time, k = k, 
                                           sensitivity = sensitivity,
                                           pre_onset_params, st_params, pi, funs, 
                                           post.screen.lookout = post.screen.lookout,
                                           method = type)
  
  c(abs(em_sens - sensitivity))
  
}


