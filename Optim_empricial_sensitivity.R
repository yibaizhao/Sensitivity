library("here")
library("ggplot2")
source("empirical_sensitivity_analytic_general_functions.R")

#######################################
# Given true sensitivity and mean sojourn time, solve for post-screen lookout time, so that empirical=true sensitivity.
#######################################
empirical.sensitivity.diff.V1 <- function(params, t1, rate.matrix, post.screen.lookout=NULL, sensitivity=NULL,
                                          clinical.cancer.state,pre.clinical.cancer.state){
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
                                           clinical.cancer.state = 3, 
                                           pre.clinical.cancer.state = 2,
                                           method='Single')

  return(abs(em_sens - sensitivity))
  
}

sens_seq <- seq(0.9, 0.5, by = -0.1)
mean.sojourn.time = c(2, 4, 6)
post.screen.lookout = seq(0, 10, length.out = 100)
sens_mst_lookout_seq <- expand.grid(sens_seq, mean.sojourn.time, post.screen.lookout)

out5 <- apply(sens_mst_lookout_seq, 1, function(row){
  rateexample=matrix(c(-.002,.002,0,0,-1/row[2],1/row[2],0,0,0),byrow=T,nrow=3)
  sens_empirical <- empirical.sensitivity.general(screen.start.time = 1, k = 1, 
                                                  sensitivity = row[1],
                                                  rate.matrix = rateexample,
                                                  post.screen.lookout = row[3],
                                                  clinical.cancer.state = 3, 
                                                  pre.clinical.cancer.state = 2,
                                                  method='Single')
  c(row, sens_empirical)
})
out5 <- as.data.frame(t(out5))
names(out5) <- c('true_senstivity', 'mean_sojourn_time', 'post_screen_lookout', 'empirical_senstivity')
out5$`mean_sojourn_time` <- with(out5, factor(paste0('Mean sojourn time=', mean_sojourn_time), levels = paste0('Mean sojourn time=', mean.sojourn.time)))
out5$`true_senstivity` <- with(out5, paste0('True sensitivity=', sens_seq))

# optimization
sens_seq <- seq(0.9, 0.5, by = -0.1)
mean.sojourn.time = c(2, 4, 6)
sens_mst_seq <- expand.grid(sens_seq, mean.sojourn.time)

out5_ests <- apply(sens_mst_seq, 1, function(row){
  params_init = 1
  rateexample=matrix(c(-.02,.02,0,0,-1/row[2],1/row[2],0,0,0),byrow=T,nrow=3)
  ests <- optim(par=params_init, fn=empirical.sensitivity.diff.V1,
                # method = 'Brent', lower = 1e-6, upper = 1e+6,
                method='L-BFGS-B', lower=0, upper=20,
                t1 = 1, rate.matrix = rateexample, post.screen.lookout = NULL, sensitivity = row[1],
                clinical.cancer.state = 3, pre.clinical.cancer.state = 2)$par
  sens_empirical <- empirical.sensitivity.general(screen.start.time = 1, k = 1, 
                                                  sensitivity = row[1],
                                                  rate.matrix = rateexample,
                                                  post.screen.lookout = ests,
                                                  clinical.cancer.state = 3, 
                                                  pre.clinical.cancer.state = 2,
                                                  method='Single')
  return(c(row, ests, sens_empirical))
})
out5_ests <- as.data.frame(t(out5_ests))
names(out5_ests) <- c('true_senstivity', 'mean_sojourn_time', 'post_screen_lookout', 'empirical_senstivity')
out5_ests$empirical_senstivity <- with(out5_ests, round(empirical_senstivity, 1))
out5_ests$`mean_sojourn_time` <- with(out5_ests, factor(paste0('Mean sojourn time=', mean_sojourn_time), levels = paste0('Mean sojourn time=', mean.sojourn.time)))
out5_ests$`true_senstivity` <- with(out5_ests, paste0('True sensitivity=', sens_seq))


# plot
ggplot(data = out5, aes(x = post_screen_lookout, y = empirical_senstivity, color = factor(true_senstivity))) +
  geom_line() +
  geom_hline(data = out5_ests, aes(yintercept = empirical_senstivity, color = factor(true_senstivity)), lty = 'dashed') +
  geom_vline(data = out5_ests, aes(xintercept = post_screen_lookout, color = factor(true_senstivity)), lty = 'dashed') +
  geom_point(data = out5_ests, aes(x = post_screen_lookout, y = empirical_senstivity, color = factor(true_senstivity)), size = 3) +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  facet_grid(mean_sojourn_time~.) +
  theme(legend.title = element_blank())

ggplot(data = out5, aes(x = post_screen_lookout, y = empirical_senstivity, color = factor(mean_sojourn_time))) +
  geom_line() +
  geom_hline(data = out5_ests, aes(yintercept = empirical_senstivity, color = factor(mean_sojourn_time)), lty = 'dashed') +
  geom_vline(data = out5_ests, aes(xintercept = post_screen_lookout, color = factor(mean_sojourn_time)), lty = 'dashed') +
  geom_point(data = out5_ests, aes(x = post_screen_lookout, y = empirical_senstivity, color = factor(mean_sojourn_time)), size = 3) +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  facet_grid(true_senstivity~.) +
  theme(legend.title = element_blank())

#######################################
# Given empirical sensitivity, true sensitivity, and screening interval, solve for sojourn time, so that empirical=true sensitivity.
#######################################
empirical.sensitivity.diff.V2 <- function(params, 
                                          k = 1,
                                          screen.start.time, 
                                          mean.sojourn.time=NULL, 
                                          post.screen.lookout=NULL, 
                                          sensitivity=NULL,
                                          clinical.cancer.state,pre.clinical.cancer.state){
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
                                           method='Single')
  
  c(abs(em_sens - sensitivity))
  
}

sens_seq <- seq(0.9, 0.2, by = -0.2)
mean.sojourn.time = seq(0.1, 30, length.out = 100)
post.screen.lookout = seq(1, 5, by = 1)
sens_mst_lookout_seq <- expand.grid(sens_seq, mean.sojourn.time, post.screen.lookout)

out_mst <- apply(sens_mst_lookout_seq, 1, function(row){
  rateexample=matrix(c(-.002,.002,0,0,-1/row[2],1/row[2],0,0,0),byrow=T,nrow=3)
  sens_empirical <- empirical.sensitivity.general(screen.start.time = 1, k = 1, 
                                                  sensitivity = row[1],
                                                  rate.matrix = rateexample,
                                                  post.screen.lookout = row[3],
                                                  clinical.cancer.state = 3, 
                                                  pre.clinical.cancer.state = 2,
                                                  method='Single')
  c(row, sens_empirical)
})
out_mst <- as.data.frame(t(out_mst))
names(out_mst) <- c('true_senstivity', 'mean_sojourn_time', 'post_screen_lookout', 'empirical_senstivity')
out_mst$post_screen_lookout <- with(out_mst, factor(paste0('post screen lookout time=', post_screen_lookout), levels = paste0('post screen lookout time=', post.screen.lookout)))
out_mst$`true_senstivity` <- with(out_mst, paste0('True sensitivity=', sens_seq))

# optimization
sens_lookout_seq <- expand.grid(sens_seq, post.screen.lookout)

library(rootSolve)# find all roots
out_mst_ests <- apply(sens_lookout_seq, 1, function(row){
  # ests <- uniroot(empirical.sensitivity.diff.V2, lower=0.1, upper=10, extendInt = "yes",
  #         k=1,
  #         screen.start.time=1, 
  #         mean.sojourn.time=NULL, 
  #         post.screen.lookout=as.numeric(row[2]), 
  #         sensitivity=as.numeric(row[1]),
  #         clinical.cancer.state = 3, pre.clinical.cancer.state = 2)$root
  # print(ests)
  ests <- optim(par=1, fn=empirical.sensitivity.diff.V2,
                method='L-BFGS-B', lower=0.1, upper=30,
                k=1,
                screen.start.time=1,
                mean.sojourn.time=NULL,
                post.screen.lookout=row[2],
                sensitivity=row[1],
                clinical.cancer.state = 3, pre.clinical.cancer.state = 2)$par
  rateexample=matrix(c(-.002,.002,0,0,-1/ests,1/ests,0,0,0),byrow=T,nrow=3)
  sens_empirical <- empirical.sensitivity.general(screen.start.time = 1, k = 1, 
                                                  sensitivity = row[1],
                                                  rate.matrix = rateexample,
                                                  post.screen.lookout = row[2],
                                                  clinical.cancer.state = 3, 
                                                  pre.clinical.cancer.state = 2,
                                                  method='Single')
  return(c(row, ests, sens_empirical))
})
out_mst_ests <- as.data.frame(t(out_mst_ests))
names(out_mst_ests) <- c('true_senstivity', 'post_screen_lookout', 'mean_sojourn_time', 'empirical_senstivity')
out_mst_ests$empirical_senstivity <- with(out_mst_ests, round(empirical_senstivity, 1))
out_mst_ests$post_screen_lookout <- with(out_mst_ests, factor(paste0('post screen lookout time=', post_screen_lookout), levels = paste0('post screen lookout time=', post.screen.lookout)))
out_mst_ests$`true_senstivity` <- with(out_mst_ests, paste0('True sensitivity=', sens_seq))

ggplot(data = out_mst, aes(x = mean_sojourn_time, y = empirical_senstivity, color = factor(true_senstivity))) +
  geom_line() +
  geom_segment(data = out_mst_ests, aes(x = mean_sojourn_time , y = empirical_senstivity, 
                                        xend = 0, yend = empirical_senstivity, 
                                        color = factor(true_senstivity)), lty = 'dashed') +
  geom_segment(data = out_mst_ests, aes(x = mean_sojourn_time , y = empirical_senstivity, 
                                        xend = mean_sojourn_time, yend = 0, 
                                        color = factor(true_senstivity)), lty = 'dashed') +
  geom_point(data = out_mst_ests, aes(x = mean_sojourn_time, y = empirical_senstivity, color = factor(true_senstivity)), size = 3) +
  scale_x_continuous(limits=c(0, 31), breaks = seq(0, 30, by = 1), expand=c(0, 0)) +
  scale_y_continuous(limits=c(0, 1.01), breaks = seq(0, 1, by = 0.1), expand=c(0, 0)) +
  facet_grid(post_screen_lookout~.) +
  theme(legend.title = element_blank())



