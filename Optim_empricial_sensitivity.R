library("here")
library("ggplot2")
source("empirical_sensitivity_analytic_general_functions.R")
source("Optim_empricial_sensitivity_functions.R")

#######################################
# Given true sensitivity and mean sojourn time, solve for post-screen lookout time, so that empirical=true sensitivity.
#######################################
sens_seq <- seq(0.9, 0.2, by = -0.2)
mean.sojourn.time = c(2, 4, 6)
post.screen.lookout = seq(0, 3, length.out = 100)
sens_mst_lookout_seq <- expand.grid(sens_seq, mean.sojourn.time, post.screen.lookout)

out_interval <- apply(sens_mst_lookout_seq, 1, function(row){
  rateexample=matrix(c(-.002,.002,0,0,-1/row[2],1/row[2],0,0,0),byrow=T,nrow=3)
  sens_empirical <- empirical.sensitivity.general(screen.start.time = 1, k = 1, 
                                                  sensitivity = row[1],
                                                  rate.matrix = rateexample,
                                                  post.screen.lookout = row[3],
                                                  clinical.cancer.state = 3, 
                                                  pre.clinical.cancer.state = 2,
                                                  method='All')
  c(row, sens_empirical)
})
out_interval <- as.data.frame(t(out_interval))
names(out_interval) <- c('true_senstivity', 'mean_sojourn_time', 'post_screen_lookout', 'empirical_senstivity')
out_interval$`mean_sojourn_time` <- with(out_interval, factor(paste0('Mean sojourn time=', mean_sojourn_time), levels = paste0('Mean sojourn time=', mean.sojourn.time)))
out_interval$`true_senstivity` <- with(out_interval, paste0('True sensitivity=', sens_seq))

# optimization
sens_mst_seq <- expand.grid(sens_seq, mean.sojourn.time)

out_interval_ests <- apply(sens_mst_seq, 1, function(row){
  params_init = 1
  rateexample=matrix(c(-.002,.002,0,0,-1/row[2],1/row[2],0,0,0),byrow=T,nrow=3)
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
                                                  method='All')
  return(c(row, ests, sens_empirical))
})
out_interval_ests <- as.data.frame(t(out_interval_ests))
names(out_interval_ests) <- c('true_senstivity', 'mean_sojourn_time', 'post_screen_lookout', 'empirical_senstivity')
out_interval_ests$empirical_senstivity <- with(out_interval_ests, round(empirical_senstivity, 1))
out_interval_ests$`mean_sojourn_time` <- with(out_interval_ests, factor(paste0('Mean sojourn time=', mean_sojourn_time), levels = paste0('Mean sojourn time=', mean.sojourn.time)))
out_interval_ests$`true_senstivity` <- with(out_interval_ests, paste0('True sensitivity=', sens_seq))


# plot
# g_optim_interval <-
#   ggplot(data = out_interval, aes(x = post_screen_lookout, y = empirical_senstivity, color = factor(true_senstivity))) +
#   geom_line() +
#   geom_segment(data = out_interval_ests, aes(x = post_screen_lookout , y = empirical_senstivity,
#                                         xend = 0, yend = empirical_senstivity,
#                                         color = factor(true_senstivity)), lty = 'dashed') +
#   geom_segment(data = out_interval_ests, aes(x = post_screen_lookout , y = empirical_senstivity,
#                                         xend = post_screen_lookout, yend = 0,
#                                         color = factor(true_senstivity)), lty = 'dashed') +
#   geom_point(data = out_interval_ests, aes(x = post_screen_lookout, y = empirical_senstivity, color = factor(true_senstivity)), size = 3) +
#   scale_x_continuous(breaks = seq(0, 10, by = 1), expand=c(0, 0)) +
#   scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand=c(0, 0)) +
#   facet_grid(mean_sojourn_time~.) +
#   theme_bw() +
#   labs(x = 'Post Screen Lookout Time', y = 'Empirical Sensitivity') +
#   theme(legend.title = element_blank())
# g_optim_interval

g_optim_interval <-
  ggplot(data = out_interval, aes(x = empirical_senstivity, y = post_screen_lookout, color = factor(true_senstivity))) +
  geom_line() +
  geom_segment(data = out_interval_ests, aes(y = post_screen_lookout , x = empirical_senstivity, 
                                             yend = 0, xend = empirical_senstivity, 
                                             color = factor(true_senstivity)), lty = 'dashed') +
  geom_segment(data = out_interval_ests, aes(y = post_screen_lookout , x = empirical_senstivity, 
                                             yend = post_screen_lookout, xend = 0, 
                                             color = factor(true_senstivity)), lty = 'dashed') +
  geom_point(data = out_interval_ests, aes(x = empirical_senstivity, y = post_screen_lookout, color = factor(true_senstivity)), size = 3) +
  scale_y_continuous(breaks = seq(0, 10, by = 1), expand=c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), expand=c(0, 0)) +
  facet_grid(.~mean_sojourn_time) +
  theme_bw() +
  labs(y = 'Post Screen Lookout Time', x = 'Empirical Sensitivity') +
  theme(legend.title = element_blank(),
        panel.spacing = unit(1, "lines"))
g_optim_interval

ggsave("./Results/plot_optim_interval_V2.jpg", width = 13, height = 7)

ggplot(data = out_interval, aes(x = post_screen_lookout, y = empirical_senstivity, color = factor(mean_sojourn_time))) +
  geom_line() +
  geom_hline(data = out_interval_ests, aes(yintercept = empirical_senstivity, color = factor(mean_sojourn_time)), lty = 'dashed') +
  geom_vline(data = out_interval_ests, aes(xintercept = post_screen_lookout, color = factor(mean_sojourn_time)), lty = 'dashed') +
  geom_point(data = out_interval_ests, aes(x = post_screen_lookout, y = empirical_senstivity, color = factor(mean_sojourn_time)), size = 3) +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  facet_grid(true_senstivity~.) +
  theme(legend.title = element_blank())

#######################################
# Given empirical sensitivity, true sensitivity, and screening interval, solve for sojourn time, so that empirical=true sensitivity.
#######################################
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

# library(rootSolve)# find all roots
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
                clinical.cancer.state = 3, pre.clinical.cancer.state = 2, type = 'All')$par
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

g_optim_mst <-
  ggplot(data = out_mst, aes(y = mean_sojourn_time, x = empirical_senstivity, color = factor(true_senstivity))) +
  geom_line() +
  geom_segment(data = out_mst_ests, aes(y = mean_sojourn_time , x = empirical_senstivity, 
                                        yend = 0, xend = empirical_senstivity, 
                                        color = factor(true_senstivity)), lty = 'dashed') +
  geom_segment(data = out_mst_ests, aes(y = mean_sojourn_time , x = empirical_senstivity, 
                                        yend = mean_sojourn_time, xend = 0, 
                                        color = factor(true_senstivity)), lty = 'dashed') +
  geom_point(data = out_mst_ests, aes(y = mean_sojourn_time, x = empirical_senstivity, color = factor(true_senstivity)), size = 3) +
  scale_y_continuous(limits=c(0, 31), breaks = seq(0, 30, by = 1), expand=c(0, 0)) +
  scale_x_continuous(limits=c(0, 1.01), breaks = seq(0, 1, by = 0.1), expand=c(0, 0)) +
  facet_grid(.~post_screen_lookout) +
  theme_bw() +
  labs(y = 'Mean Sojourn Time', x = 'Empirical Sensitivity') +
  theme(legend.title = element_blank(),
        panel.spacing = unit(1, "lines"))
g_optim_mst
ggsave("./Results/plot_optim_mst_V2.jpg", width = 13, height = 7)


#######################################
# Given empirical sensitivity, true sensitivity, and mst, solve for screen interval, so that empirical=true sensitivity.
# Find optimal screen interval for multiple screens
#######################################
sens_seq <- seq(0.9, 0.2, by = -0.2)
mean.sojourn.time = c(1, 2, 4, 6, 8)
post.screen.lookout = seq(0, 7, length.out = 100)
k_seq <- seq(1, 5, by=1)
sens_mst_lookout_k_seq <- expand.grid(sens_seq, mean.sojourn.time, post.screen.lookout, k_seq)

out_interval <- apply(sens_mst_lookout_k_seq, 1, function(row){
  # st_params <- optim(par = c(0.1, 0.1), fn = find_paramsF, method = 'L-BFGS-B', lower = 0, upper = Inf, mst = row[2])$par
  sens_empirical <- empirical.sensitivity.general(screen.start.time = 1, k = row[4], 
                                                  sensitivity = row[1],
                                                  pre_onset_params = 0.002, st_params = 1/row[2], funs = 'exp',
                                                  post.screen.lookout = row[3],
                                                  method='All')
  c(row, sens_empirical)
})
out_interval <- as.data.frame(t(out_interval))
names(out_interval) <- c('true_senstivity', 'mean_sojourn_time', 'post_screen_lookout', 'k', 'empirical_senstivity')
out_interval$`mean_sojourn_time` <- with(out_interval, factor(paste0('Mean sojourn time=', mean_sojourn_time), levels = paste0('Mean sojourn time=', mean.sojourn.time)))
out_interval$`true_senstivity` <- with(out_interval, paste0('True sensitivity=', true_senstivity))
out_interval$`k` <- with(out_interval, paste0('no. of screens=', k))

# optimization
sens_mst_k_seq <- expand.grid(sens_seq, mean.sojourn.time, k_seq)

out_interval_ests <- apply(sens_mst_k_seq, 1, function(row){
  ests <- optim(par=1, fn=empirical.sensitivity.diff,
                method='L-BFGS-B', lower=0, upper=10,
                k=as.numeric(row[3]),
                screen.start.time=1,
                pre_onset_params = 0.002, st_params = 1/row[2], funs = 'exp',
                sensitivity=row[1],
                type = 'All')$par
  sens_empirical <- empirical.sensitivity.general(screen.start.time = 1, k = row[3], 
                                                  sensitivity = row[1],
                                                  pre_onset_params = 0.002, st_params = 1/row[2], funs = 'exp',
                                                  post.screen.lookout = ests,
                                                  method='All')
  return(c(row, ests, sens_empirical))
})
out_interval_ests <- as.data.frame(t(out_interval_ests))
names(out_interval_ests) <- c('true_senstivity', 'mean_sojourn_time', 'k', 'post_screen_lookout', 'empirical_senstivity')
out_interval_ests$empirical_senstivity <- with(out_interval_ests, round(empirical_senstivity, 1))
out_interval_ests$`mean_sojourn_time` <- with(out_interval_ests, factor(paste0('Mean sojourn time=', mean_sojourn_time), levels = paste0('Mean sojourn time=', mean.sojourn.time)))
out_interval_ests$`true_senstivity` <- with(out_interval_ests, paste0('True sensitivity=', true_senstivity))
out_interval_ests$`k` <- with(out_interval_ests, paste0('no. of screens=', k))

g_optim_interval <-
  ggplot(data = out_interval, aes(x = empirical_senstivity, y = post_screen_lookout, color = factor(true_senstivity))) +
  geom_line() +
  geom_segment(data = out_interval_ests, aes(y = post_screen_lookout , x = empirical_senstivity, 
                                             yend = 0, xend = empirical_senstivity, 
                                             color = factor(true_senstivity)), lty = 'dashed') +
  geom_segment(data = out_interval_ests, aes(y = post_screen_lookout , x = empirical_senstivity, 
                                             yend = post_screen_lookout, xend = 0, 
                                             color = factor(true_senstivity)), lty = 'dashed') +
  geom_point(data = out_interval_ests, aes(x = empirical_senstivity, y = post_screen_lookout, color = factor(true_senstivity)), size = 3) +
  scale_y_continuous(breaks = seq(0, 6, by = 1), expand=c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), expand=c(0, 0)) +
  facet_grid(k~mean_sojourn_time) +
  theme_bw() +
  labs(y = 'Post Screen Lookout Time', x = 'Empirical Sensitivity') +
  theme(legend.title = element_blank(),
        panel.spacing = unit(1, "lines"))

g_optim_interval <-
  ggplot(data = out_interval, aes(x = post_screen_lookout, y = empirical_senstivity, color = factor(true_senstivity))) +
  geom_line() +
  geom_segment(data = out_interval_ests, aes(x = post_screen_lookout , y = empirical_senstivity, 
                                             xend = 0, yend = empirical_senstivity, 
                                             color = factor(true_senstivity)), lty = 'dashed') +
  geom_segment(data = out_interval_ests, aes(x = post_screen_lookout , y = empirical_senstivity, 
                                             xend = post_screen_lookout, yend = 0, 
                                             color = factor(true_senstivity)), lty = 'dashed') +
  geom_point(data = out_interval_ests, aes(x = post_screen_lookout, y = empirical_senstivity, color = factor(true_senstivity)), size = 3) +
  scale_x_continuous(breaks = seq(0, 7, by = 1), expand=c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand=c(0, 0)) +
  facet_grid(k~mean_sojourn_time) +
  theme_bw() +
  labs(x = 'Post Screen Lookout Time', y = 'Empirical Sensitivity') +
  theme(legend.title = element_blank(),
        panel.spacing = unit(1, "lines"))
g_optim_interval

ggsave("./Results/plot_optim_interval_across_k_V2.jpg", width = 13, height = 7)

# One single plot
## plot optimal lookout (y axis) versus true sensitivity (x axis). 
## Colors by Mean sojourn time.
ggplot(data = out_interval_ests, aes(x = true_senstivity, y = post_screen_lookout, color = factor(mean_sojourn_time))) +
  geom_point() +
  facet_grid(k~.) +
  theme_bw() +
  labs(x = 'True Sensitivity', y = 'Post Screen Lookout Time') +
  theme(legend.title = element_blank()) +
  ggtitle('Start screen time = 1')
ggsave("./Results/plot_optim_interval_across_k_dots_start1.jpg", width = 8, height = 7)

# check empirical across series of screens
out_interval %>% 
  ggplot(aes(x = k, y = empirical_senstivity,  color = factor(mean_sojourn_time))) +
  geom_line() + geom_point() +
  geom_hline(aes(yintercept = true_senstivity), lty = 'dashed') +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), expand=c(0, 0)) +
  facet_grid(paste0('lookout=', post_screen_lookout)~paste0('true sensitivity=', true_senstivity)) +
  theme(legend.title = element_blank()) +
  xlab('Across k screens') + ggtitle('Start screen time = 1')
ggsave("./Results/plot_sens_across_k_start1.jpg", width = 12, height = 8)

#######################################
# Given empirical sensitivity, true sensitivity, and mst, pi, solve for screen interval, so that empirical=true sensitivity.
# Find optimal screen interval for multiple screens
#######################################
pi <- c(0.1, 0.3, 0.5, 0.7)
sens_seq <- seq(0.9, 0.2, by = -0.2)
# mean.sojourn.time = c(1, 2, 4, 6, 8)
mean.sojourn.time <- 2
post.screen.lookout = seq(0, 7, length.out = 100)
k_seq <- seq(1, 10, by=2)
sens_mst_lookout_k_seq <- expand.grid(pi, sens_seq, mean.sojourn.time, post.screen.lookout, k_seq)

out_interval <- apply(sens_mst_lookout_k_seq, 1, function(row){
  # st_params <- optim(par = c(0.1, 0.1), fn = find_paramsF, method = 'L-BFGS-B', lower = 0, upper = Inf, mst = row[2])$par
  sens_empirical <- empirical.sensitivity.general(screen.start.time = 1, k = row[5], 
                                                  sensitivity = row[2],
                                                  pre_onset_params = 0.002, st_params = 1/row[3], pi = row[1], funs = 'exp',
                                                  post.screen.lookout = row[4],
                                                  method='All')
  c(row, sens_empirical)
})
out_interval <- as.data.frame(t(out_interval))
names(out_interval) <- c('pi', 'true_senstivity', 'mean_sojourn_time', 'post_screen_lookout', 'k', 'empirical_senstivity')
out_interval$pi <- with(out_interval, (paste0('pi=', pi)))
out_interval$`true_senstivity` <- with(out_interval, paste0('True sensitivity=', true_senstivity))
out_interval$`k` <- with(out_interval, paste0('no. of screens=', k))

# optimization
sens_mst_k_seq <- expand.grid(pi, sens_seq, mean.sojourn.time, k_seq)

out_interval_ests <- apply(sens_mst_k_seq, 1, function(row){
  ests <- optim(par=1, fn=empirical.sensitivity.diff,
                method='L-BFGS-B', lower=0, upper=10,
                k=as.numeric(row[4]),
                screen.start.time=1,
                pre_onset_params = 0.002, st_params = 1/row[3], pi = row[1], funs = 'exp',
                sensitivity=row[2],
                type = 'All')$par
  sens_empirical <- empirical.sensitivity.general(screen.start.time = 1, k = row[4], 
                                                  sensitivity = row[2],
                                                  pre_onset_params = 0.002, st_params = 1/row[3], pi = row[1], funs = 'exp',
                                                  post.screen.lookout = ests,
                                                  method='All')
  return(c(row, ests, sens_empirical))
})
out_interval_ests <- as.data.frame(t(out_interval_ests))
names(out_interval_ests) <- c('pi', 'true_senstivity', 'mean_sojourn_time', 'k', 'post_screen_lookout', 'empirical_senstivity')
out_interval_ests$empirical_senstivity <- with(out_interval_ests, round(empirical_senstivity, 1))
out_interval_ests$pi <- with(out_interval_ests, (paste0('pi=', pi)))
out_interval_ests$`true_senstivity` <- with(out_interval_ests, paste0('True sensitivity=', true_senstivity))
out_interval_ests$`k` <- with(out_interval_ests, paste0('no. of screens=', k))

g_optim_interval <-
  ggplot(data = out_interval, aes(x = empirical_senstivity, y = post_screen_lookout, color = factor(true_senstivity))) +
  geom_line() +
  geom_segment(data = out_interval_ests, aes(y = post_screen_lookout , x = empirical_senstivity, 
                                             yend = 0, xend = empirical_senstivity, 
                                             color = factor(true_senstivity)), lty = 'dashed') +
  geom_segment(data = out_interval_ests, aes(y = post_screen_lookout , x = empirical_senstivity, 
                                             yend = post_screen_lookout, xend = 0, 
                                             color = factor(true_senstivity)), lty = 'dashed') +
  geom_point(data = out_interval_ests, aes(x = empirical_senstivity, y = post_screen_lookout, color = factor(true_senstivity)), size = 3) +
  scale_y_continuous(breaks = seq(0, 6, by = 1), expand=c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), expand=c(0, 0)) +
  facet_grid(k~pi) +
  theme_bw() +
  labs(y = 'Post Screen Lookout Time', x = 'Empirical Sensitivity') +
  theme(legend.title = element_blank(),
        panel.spacing = unit(1, "lines"))

g_optim_interval <-
  ggplot(data = out_interval, aes(x = post_screen_lookout, y = empirical_senstivity, color = factor(true_senstivity))) +
  geom_line() +
  geom_segment(data = out_interval_ests, aes(x = post_screen_lookout , y = empirical_senstivity, 
                                             xend = 0, yend = empirical_senstivity, 
                                             color = factor(true_senstivity)), lty = 'dashed') +
  geom_segment(data = out_interval_ests, aes(x = post_screen_lookout , y = empirical_senstivity, 
                                             xend = post_screen_lookout, yend = 0, 
                                             color = factor(true_senstivity)), lty = 'dashed') +
  geom_point(data = out_interval_ests, aes(x = post_screen_lookout, y = empirical_senstivity, color = factor(true_senstivity)), size = 3) +
  scale_x_continuous(breaks = seq(0, 7, by = 1), expand=c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), expand=c(0, 0)) +
  facet_grid(k~pi) +
  theme_bw() +
  labs(x = 'Post Screen Lookout Time', y = 'Empirical Sensitivity') +
  theme(legend.title = element_blank(),
        panel.spacing = unit(1, "lines"))
g_optim_interval

ggsave("./Results/plot_optim_interval_across_k_V2.jpg", width = 13, height = 7)

