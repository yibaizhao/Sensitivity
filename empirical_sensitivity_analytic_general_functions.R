# integral
integrand.screen.detect <- function(u, screen.times, k, lambda, gamma){
  ########################################################################################
  #Author: Yibai Zhao
  #This function calculate the probability of screen detection at time u
  # assumes exponential distributions for onset and sojourn times 
  #INPUTS: u: pre-clinical onset time
  #        screen.times: screen times
  #        k: k^{th} screen time, t_k
  #        lambda: mean pre-clinical onset time
  #        gamma: mean sojourn time
  #OUTPUTS: numeric value between 0 and 1
  ##########################################################################################
  
  dexp(u, rate=lambda)*pexp(screen.times[k]-u, rate=gamma,lower.tail=F)
  
}


integrand.interval.cancer.false.negative <- function(u, screen.times, k, post.screen.lookout, lambda, gamma){
  ########################################################################################
  #Author: Yibai Zhao
  #This function calculate the probability of interval cancer when onset at time u and k^{th} test is false negative
  # assumes exponential distributions for onset and sojourn times 
  #INPUTS: u: pre-clinical onset time
  #        screen.times: screen times
  #        k: k^{th} screen time, t_k
  #        post.screen.lookout: post screen lookout time
  #        lambda: mean pre-clinical onset time from exp(lambda)
  #        gamma: mean sojourn time fromo exp(gamma)
  #OUTPUTS: numeric value between 0 and 1
  ##########################################################################################
  
  screen.times_k_1 <- screen.times[k] + post.screen.lookout
  dexp(u, rate=lambda)*( pexp(screen.times_k_1-u, rate=gamma) - pexp(screen.times[k]-u, rate=gamma) )
  
}

integrand.interval.cancer.true<-function(u, screen.times, k, post.screen.lookout, lambda, gamma){
  ########################################################################################
  #Author: Yibai Zhao
  #This function calculate the probability of interval cancer when onset at time u and k^{th} screen test is true negative
  # assumes exponential distributions for onset and sojourn times 
  #INPUTS: u: pre-clinical onset time
  #        screen.times: screen times
  #        k: k^{th} screen time, t_k
  #        post.screen.lookout: post screen lookout time
  #        lambda: mean pre-clinical onset time from exp(lambda)
  #        gamma: mean sojourn time fromo exp(gamma)
  #OUTPUTS: numeric value between 0 and 1
  ##########################################################################################
  
  screen.times_k_1 <- screen.times[k] + post.screen.lookout
  dexp(u, rate=lambda)*pexp(screen.times_k_1-u, rate=gamma, lower.tail=T)
}

screen.detect.time.k <- function(screen.times, k, lambda, gamma,
                                 sensitivity,
                                 clinical.cancer.state, pre.clinical.cancer.state){
  ########################################################################################
  #Author: Yibai Zhao
  #This function calculate the probability of screen detected at k^{th} screen time
  # assumes exponential distributions for onset and sojourn times 
  #INPUTS: screen.times: screen times
  #        k: k^{th} screen time, t_k
  #        lambda: mean pre-clinical onset time from exp(lambda)
  #        gamma: mean sojourn time fromo exp(gamma)
  #        sensitivity: true sensitivity
  #OUTPUTS: numeric value between 0 and 1
  ##########################################################################################
  
  screen.detect <- 0
  # k=0 #####
  if(k == 1){
    screen.detect <- screen.detect + integrate(integrand.screen.detect,
                                               lower=0, upper=screen.times[k], screen.times, k, lambda, gamma)$value*sensitivity
    return(screen.detect)
  }
  # k>0 ####
  # case 1: false negative before k
  for (i in c(1:(k-1))) {
    screen.times_i_1 <- ifelse(i == 1, 0, screen.times[i-1])
    screen.detect <- screen.detect + integrate(integrand.screen.detect,
                                               lower=screen.times_i_1, upper=screen.times[i], screen.times, k, lambda, gamma)$value*(1-sensitivity)^(k-i)*sensitivity
  }
  # case 2: true negative before k
  screen.detect <- screen.detect + integrate(integrand.screen.detect,
                                             lower=screen.times[k-1], upper=screen.times[k], screen.times, k, lambda, gamma)$value*sensitivity
  return(screen.detect)
}

interval.cancer.time.k <- function(screen.times, k, 
                                   lambda, gamma,
                                   sensitivity,
                                   post.screen.lookout,
                                   clinical.cancer.state, pre.clinical.cancer.state){
  ########################################################################################
  #Author: Yibai Zhao
  #This function calculate the probability of interval cancer at k^{th} screen time
  # assumes exponential distributions for onset and sojourn times 
  #INPUTS: screen.times: screen times
  #        k: k^{th} screen time, t_k
  #        lambda: mean pre-clinical onset time from exp(lambda)
  #        gamma: mean sojourn time fromo exp(gamma)
  #        sensitivity: true sensitivity
  #        post.screen.lookout: post screen lookout time
  #OUTPUTS: numeric value between 0 and 1
  ##########################################################################################
  
  interval.cancer <- 0
  
  # case 1: false negative at/before k
  for (i in c(1:k)) {
    screen.times_i_1 <- ifelse(i == 1, 0, screen.times[i-1])
    interval.cancer <- interval.cancer + integrate(integrand.interval.cancer.false.negative,
                                                   lower=screen.times_i_1, upper=screen.times[i], screen.times, k, post.screen.lookout, lambda, gamma)$value * (1-sensitivity)^{k-i+1}
  }
  # case 2: interval detected in (t_k, t_{k+1})
  screen.times_k_1 <- screen.times[k] + post.screen.lookout
  if(screen.times[k] < max(screen.times)){
    screen.times_k_1 <- screen.times[k+1]
  }
  interval.cancer <- interval.cancer + integrate(integrand.interval.cancer.true,
                                                 lower=screen.times[k], upper=screen.times_k_1, screen.times, k, post.screen.lookout, lambda, gamma)$value
  
  return(interval.cancer)
}

empirical.sensitivity.general <- function(screen.start.time, k, #start.time=NULL, end.time=NULL,
                                          sensitivity,
                                          rate.matrix,
                                          post.screen.lookout,
                                          clinical.cancer.state, pre.clinical.cancer.state,
                                          method='Single'){
  ########################################################################################
  #Author: Yibai Zhao
  #This function calculate the empirical sensitivity using analytic method
  # assumes exponential distributions for onset and sojourn times 
  #INPUTS: screen.start.time: first screen time
  #        k: k^{th} screen time, t_k
  #        sensitivity: true sensitivity
  #        rate.matrix=transition matrix of cancer model
  #        post.screen.lookout: post screen lookout time
  #        method: (1) 'Single': screen by screen at k^{th} screen time;
  #                (2) 'All': across k serious of screen times
  #OUTPUTS: numeric value between 0 and 1
  ##########################################################################################

  if( is.null(screen.start.time) |  (screen.start.time<=0) ) stop('screen start time must be greater than zero')

  screen.times <- screen.start.time + cumsum(rep(post.screen.lookout, k)) - post.screen.lookout
  
  lambda = rate.matrix[1, pre.clinical.cancer.state]
  gamma = rate.matrix[pre.clinical.cancer.state, clinical.cancer.state]
  
  screen.detect = screen.detect.time.k(screen.times, k, lambda, gamma,
                                       sensitivity,
                                       clinical.cancer.state, pre.clinical.cancer.state)
  interval.cancer = interval.cancer.time.k(screen.times, k, lambda, gamma,
                                           sensitivity,
                                           post.screen.lookout,
                                           clinical.cancer.state, pre.clinical.cancer.state)
  
  if((method=='All') & (k>1)){ #& (1<start.time)){
    # screen.detect <- interval.cancer <- 0
    for (tk in c(2:k)) {
    # for (tk in c(start.time:end.time)) {
      screen.detect <- screen.detect + screen.detect.time.k(screen.times, k=tk, lambda, gamma,
                                           sensitivity,
                                           clinical.cancer.state, pre.clinical.cancer.state)
      interval.cancer <- interval.cancer + interval.cancer.time.k(screen.times, k=tk, lambda, gamma,
                                               sensitivity,
                                               post.screen.lookout,
                                               clinical.cancer.state, pre.clinical.cancer.state)
      
    }
  }
  
  return(screen.detect / (screen.detect + interval.cancer))
  
}




