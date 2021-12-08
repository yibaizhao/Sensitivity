# integral
integrand.screen.detect <- function(u, screen.times, k, lambda, st_params, pi, funs){
  ########################################################################################
  #Author: Yibai Zhao
  #This function calculate the probability of screen detection at time u
  # assumes exponential distributions for onset and sojourn times 
  #INPUTS: u: pre-clinical onset time
  #        screen.times: screen times
  #        k: k^{th} screen time, t_k
  #        lambda: mean pre-clinical onset time
  #        st_params: parameters for sojourn time function
  #        pi: probability of being indolent, such that F*=pi+(1-pi)*F --> S*=(1-pi)S
  #        funs: function of sojourn time
  #OUTPUTS: numeric value between 0 and 1
  ##########################################################################################
  
  if(funs == 'exp'){
    return(dexp(u, rate=lambda)*(1-pi)*pexp(screen.times[k]-u, rate=st_params,lower.tail=F))
  }
  if(funs == 'weibull'){
    return(dexp(u, rate=lambda)*(1-pi)*pweibull(screen.times[k]-u, shape=st_params[1], scale=st_params[2],lower.tail=F))
  }
  
}


integrand.interval.cancer.false.negative <- function(u, screen.times, k, post.screen.lookout, lambda, st_params, pi, funs){
  ########################################################################################
  #Author: Yibai Zhao
  #This function calculate the probability of interval cancer when onset at time u and k^{th} test is false negative
  # assumes exponential distributions for onset and sojourn times 
  #INPUTS: u: pre-clinical onset time
  #        screen.times: screen times
  #        k: k^{th} screen time, t_k
  #        post.screen.lookout: post screen lookout time
  #        lambda: mean pre-clinical onset time from exp(lambda)
  #        st_params: parameters for sojourn time function
  #        pi: probability of being indolent, such that F*=pi+(1-pi)*F --> S*=(1-pi)S
  #        funs: function of sojourn time
  #OUTPUTS: numeric value between 0 and 1
  ##########################################################################################
  
  screen.times_k_1 <- screen.times[k] + post.screen.lookout
  if(funs == 'exp'){
    return(dexp(u, rate=lambda)*( (1-pi)*pexp(screen.times_k_1-u, rate=st_params) - (1-pi)*pexp(screen.times[k]-u, rate=st_params) ))
  }
  if(funs == 'weibull'){
    return(dexp(u, rate=lambda)*( 
      (1-pi)*pweibull(screen.times_k_1-u, shape=st_params[1], scale=st_params[2]) - 
        (1-pi)*pweibull(screen.times[k]-u, shape=st_params[1], scale=st_params[2]) 
      ))
  }
  
}

integrand.interval.cancer.true<-function(u, screen.times, k, post.screen.lookout, lambda, st_params, pi, funs){
  ########################################################################################
  #Author: Yibai Zhao
  #This function calculate the probability of interval cancer when onset at time u and k^{th} screen test is true negative
  # assumes exponential distributions for onset and sojourn times 
  #INPUTS: u: pre-clinical onset time
  #        screen.times: screen times
  #        k: k^{th} screen time, t_k
  #        post.screen.lookout: post screen lookout time
  #        lambda: mean pre-clinical onset time from exp(lambda)
  #        st_params: parameters for sojourn time function
  #        pi: probability of being indolent, such that F*=pi+(1-pi)*F --> S*=(1-pi)S
  #        funs: function of sojourn time
  #OUTPUTS: numeric value between 0 and 1
  ##########################################################################################
  
  screen.times_k_1 <- screen.times[k] + post.screen.lookout
  if(funs == 'exp'){
    return(dexp(u, rate=lambda)*(pi + (1-pi)*pexp(screen.times_k_1-u, rate=st_params, lower.tail=T)))
  }
  if(funs == 'weibull'){
    return(dexp(u, rate=lambda)*(pi + (1-pi)*pweibull(screen.times_k_1-u, shape=st_params[1], scale=st_params[2], lower.tail=T)))
  }
  
}

screen.detect.time.k <- function(screen.times, k, lambda, st_params, pi, funs, 
                                 sensitivity){
  ########################################################################################
  #Author: Yibai Zhao
  #This function calculate the probability of screen detected at k^{th} screen time
  # assumes exponential distributions for onset and sojourn times 
  #INPUTS: screen.times: screen times
  #        k: k^{th} screen time, t_k
  #        lambda: mean pre-clinical onset time from exp(lambda)
  #        st_params: parameters for sojourn time function
  #        funs: function of sojourn time
  #        pi: probability of being indolent, such that F*=pi+(1-pi)*F --> S*=(1-pi)S
  #        sensitivity: true sensitivity
  #OUTPUTS: numeric value between 0 and 1
  ##########################################################################################
  
  screen.detect <- 0
  # k=0 #####
  if(k == 1){
    screen.detect <- screen.detect + integrate(integrand.screen.detect,
                                               lower=0, upper=screen.times[k], 
                                               screen.times, k, lambda, st_params, pi, funs)$value*sensitivity
    return(screen.detect)
  }
  # k>0 ####
  # case 1: false negative before k
  for (i in c(1:(k-1))) {
    screen.times_i_1 <- ifelse(i == 1, 0, screen.times[i-1])
    screen.detect <- screen.detect + integrate(integrand.screen.detect,
                                               lower=screen.times_i_1, upper=screen.times[i], 
                                               screen.times, k, lambda, st_params, pi, funs)$value*(1-sensitivity)^(k-i)*sensitivity
  }
  # case 2: true negative before k
  screen.detect <- screen.detect + integrate(integrand.screen.detect,
                                             lower=screen.times[k-1], upper=screen.times[k], 
                                             screen.times, k, lambda, st_params, pi, funs)$value*sensitivity
  return(screen.detect)
}

interval.cancer.time.k <- function(screen.times, k, 
                                   lambda, st_params, pi, funs,
                                   sensitivity,
                                   post.screen.lookout){
  ########################################################################################
  #Author: Yibai Zhao
  #This function calculate the probability of interval cancer at k^{th} screen time
  # assumes exponential distributions for onset and sojourn times 
  #INPUTS: screen.times: screen times
  #        k: k^{th} screen time, t_k
  #        lambda: mean pre-clinical onset time from exp(lambda)
  #        st_params: parameters for sojourn time function
  #        pi: probability of being indolent, such that F*=pi+(1-pi)*F --> S*=(1-pi)S
  #        funs: function of sojourn time
  #        sensitivity: true sensitivity
  #        post.screen.lookout: post screen lookout time
  #OUTPUTS: numeric value between 0 and 1
  ##########################################################################################
  
  interval.cancer <- 0
  
  # case 1: false negative at/before k
  for (i in c(1:k)) {
    screen.times_i_1 <- ifelse(i == 1, 0, screen.times[i-1])
    interval.cancer <- interval.cancer + integrate(integrand.interval.cancer.false.negative,
                                                   lower=screen.times_i_1, upper=screen.times[i], 
                                                   screen.times, k, post.screen.lookout, lambda, st_params, pi, funs)$value * (1-sensitivity)^{(k-i+1)}
  }
  # case 2: interval detected in (t_k, t_{k+1})
  screen.times_k_1 <- screen.times[k] + post.screen.lookout
  # if(screen.times[k] < max(screen.times)){
  #   screen.times_k_1 <- screen.times[k+1]
  # }
  interval.cancer <- interval.cancer + integrate(integrand.interval.cancer.true,
                                                 lower=screen.times[k], upper=screen.times_k_1, 
                                                 screen.times, k, post.screen.lookout, lambda, st_params, pi, funs)$value
  
  return(interval.cancer)
}

empirical.sensitivity.general <- function(screen.start.time, k, #start.time=NULL, end.time=NULL,
                                          sensitivity,
                                          pre_onset_params, st_params, pi = 0, funs,
                                          post.screen.lookout,
                                          method='Single'){
  ########################################################################################
  #Author: Yibai Zhao
  #This function calculate the empirical sensitivity using analytic method
  # assumes exponential distributions for onset and sojourn times 
  #INPUTS: screen.start.time: first screen time
  #        k: k^{th} screen time, t_k
  #        sensitivity: true sensitivity
  #        rate.matrix=transition matrix of cancer model
  #        pi: probability of being indolent, such that F*=pi+(1-pi)*F --> S*=(1-pi)S
  #        funs: function of sojourn time
  #        post.screen.lookout: post screen lookout time
  #        method: (1) 'Single': screen by screen at k^{th} screen time;
  #                (2) 'All': across k serious of screen times
  #OUTPUTS: numeric value between 0 and 1
  ##########################################################################################

  if( is.null(screen.start.time) |  (screen.start.time<=0) ) stop('screen start time must be greater than zero')

  post.screen.lookout <- as.numeric(post.screen.lookout)
  screen.times <- as.numeric(screen.start.time + cumsum(rep(post.screen.lookout, k)) - post.screen.lookout)
  k <- as.numeric(k)
  sensitivity <- as.numeric(sensitivity)
  lambda = pre_onset_params

  screen.detect = screen.detect.time.k(screen.times, k, lambda, st_params, pi, funs,
                                       sensitivity)
  interval.cancer = interval.cancer.time.k(screen.times, k, lambda, st_params, pi, funs,
                                           sensitivity,
                                           post.screen.lookout)
  
  if((method=='All') & (k>1)){ #& (1<start.time)){
    screen.detect <- interval.cancer <- 0
    for (tk in c(1:k)) {
      screen.detect <- screen.detect + screen.detect.time.k(screen.times, k=tk, lambda, st_params, pi, funs,
                                           sensitivity)
      interval.cancer <- interval.cancer + interval.cancer.time.k(screen.times, k=tk, lambda, st_params, pi, funs,
                                               sensitivity,
                                               post.screen.lookout)
      
    }
  }
  
  return(screen.detect / (screen.detect + interval.cancer))
  
}


#######################################
# Find parameters for mst follows weibull distribution
#######################################
find_paramsF <- function(shape=NULL, scale=NULL, mst){
  # shape = params[1]
  # scale = params[2]
  abs(scale*gamma(1+1/shape) - mst)
}


