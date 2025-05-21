# Compare performance of drcmd to tmle, drtmle, and AIPW when outcome is MAR
# (All 4 packages can handle this setting)
#-------------------------------------------------------------------------------
# set up main params
rm(list=ls())
n_grid <- c(500,1000,2500,5000) # sample sizes
nsim <- 100 # sims per grid point
#-------------------------------------------------------------------------------
# Load in relevant stuff from drcmd

source('../R/utils.R')
source('../R/drcmd.R')
source('../R/nuis.R')
source('../R/methods.R')
library(tictoc) # for timing packages
library(SuperLearner)
#-------------------------------------------------------------------------------
# Params for drcmd

eem_ind <- TRUE # estimate with empirical efficiency maximization
sl_learners <- c('SL.glm','SL.glm.interaction','SL.earth','SL.gam')
k <- 1
#-------------------------------------------------------------------------------

drcmdests <- rep(NA,100) ; drcmdtimes <- rep(NA,100)
tmleests <- rep(NA,100) ; tmletimes <- rep(NA,100)
aipwests <- rep(NA,100) ; aipwtimes <- rep(NA,100)
drtmleests <- rep(NA,100) ; drtmletimes <- rep(NA,100)

res_df <- data.frame()
for (n in n_grid) {
for (ss in 1:nsim) {
  print(ss)
  # Simulate data
  X <- rnorm(n) ; A <- rbinom(n,1,plogis(X)) ; Y <- rnorm(n) + A + X
  Ystar <- Y + rnorm(n)/2 ; R <- rbinom(n,1,plogis(X)) ; X <- as.data.frame(X)
  YY <- Y ; Y[R==0] <- NA

  # AIPW
  tic()
  tempaipw <- AIPW::AIPW$new(Y=Y,A=A,W=X,
                                 Q.SL.library = sl_learners,
                                 g.SL.library = sl_learners,
                            k_split = 1)$fit()$summary()
  aipwests <- tempaipw$estimates$RD[1]
  temp <- toc()
  aipwtimes <- temp$toc - temp$tic

  # drcmd
  tic()
  drcmdests <- drcmd(Y,A,X,
                        default_learners = sl_learners,
                        eem_ind=eem_ind, k=8) # $results$estimates$psi_hat_ate
  summary(drcmdests)
  temp <- toc()
  drcmdtimes <- temp$toc - temp$tic

  # tmle
  tic()
  tmleests <- tmle::tmle(Y,A,X, Delta=R,
                             Q.SL.library = sl_learners,
                             g.SL.library = sl_learners,
                             g.Delta.SL.library = sl_learners)$estimates$ATE$psi
  temp <- toc()
  tmletimes <- temp$toc - temp$tic

  # drtmle
  tic()
  tempdrtmle <- drtmle::drtmle(W=X,A=A,Y=Y,
                               SL_Q = sl_learners,
                               SL_g = sl_learners,
                               SL_Qr = sl_learners,
                               SL_gr = sl_learners,
                               a_0=c(1,0))
  drtmleests <- tempdrtmle$drtmle$est[1] - tempdrtmle$drtmle$est[2]
  temp <- toc()
  drtmletimes <- temp$toc - temp$tic

  # update results
  res_df <- rbind(res_df,data.frame(n=n,ss=ss,
                                    drcmdests=drcmdests,tmleests=tmleests,
                                    aipwests=aipwests,drtmleests=drtmleests,
                                    drcmdtimes=drcmdtimes,tmletimes=tmletimes,
                                    aipwtimes=aipwtimes,drtmletimes=drtmletimes))
}
}
#-------------------------------------------------------------------------------
# Plots for estimation

# reshape to wide format, removng est from each name
results <- data.frame(res_df) %>% select(-ends_with('times'))  %>%
  pivot_longer(cols=ends_with('ests'),names_to='method',values_to='est') %>%
  mutate(method = case_when(method=='drcmdests' ~ 'drcmd',
                            method=='tmleests' ~ 'tmle',
                            method=='aipwests' ~ 'AIPW',
                            method=='drtmleests' ~ 'drtmle'))

results %>% ggplot(aes(x=est,y=method,fill=method)) +
  facet_wrap(~as.factor(paste0('n=',n)), ncol=2) +
  geom_density_ridges(alpha = 0.5,
                      scale = 1,
                      rel_min_height = 0.015,
                      quantile_lines = T,
                      quantiles = 0.5,
                      panel_scaling = F) +
  theme_bw() +
  labs(title='Comparison of DR causal inference packages: estimation',
       subtitle='Outcome missing at random conditional on measured covariates',
       fill='Package',
       x='Estimate of E[Y(1) - Y(0)]',
       y='Package') +
  theme(legend.position = 'bottom')
ggsave('figures/packages_ests.pdf',width=8,height=4,units='in')

# Also do a version where RMSE is calculated
results %>% group_by(n,method) %>%
  summarise(rmse = sqrt(mean((est - 1)^2))) %>% ungroup() %>%
  filter(method != 'AIPW') %>%
  ggplot(aes(x=n,y=rmse,color=method, group=method)) +
  geom_point() + geom_line(size=0.7) +
  theme_bw() +
  labs(title='Comparison of DR causal inference packages: RMSE',
       subtitle='Outcome missing at random conditional on measured covariates',
       fill='Package',
       x='Sample size, n',
       y='RMSE') +
  theme(legend.position = 'bottom')
#-------------------------------------------------------------------------------
# Plots for timing

results <- data.frame(res_df) %>% select(-ends_with('ests'))  %>%
  pivot_longer(cols=ends_with('times'),names_to='method',values_to='times') %>%
  mutate(method = case_when(method=='drcmdtimes' ~ 'drcmd',
                            method=='tmletimes' ~ 'tmle',
                            method=='aipwtimes' ~ 'AIPW',
                            method=='drtmletimes' ~ 'drtmle'))

results %>% ggplot(aes(x=times,y=method,fill=method)) +
  facet_wrap(~as.factor(paste0('n=',n)), ncol=2) +
  geom_density_ridges(alpha = 0.5,
                      scale = 1,
                      rel_min_height = 0.015,
                      quantile_lines = T,
                      quantiles = 0.5,
                      panel_scaling = F) +
  theme_bw() +
  labs(title='Comparison of DR causal inference packages: time',
       subtitle='Outcome missing at random conditional on measured covariates',
       fill='Package',
       x='Time to run (s)',
       y='Package') +
  theme(legend.position = 'bottom')
ggsave('figures/packages_times.pdf',width=8,height=4,units='in')


results %>% group_by(n,method) %>%
  summarise(times = mean(times)) %>% ungroup() %>%
  ggplot(aes(x=n,y=times,color=method, group=method)) +
  geom_point() + geom_line(size=0.7) +
  theme_bw() +
  labs(title='Comparison of DR causal inference packages: time',
       subtitle='Outcome missing at random conditional on measured covariates',
       color='Package',
       x='Sample size, n',
       y='Time to run (s)') +
  theme(legend.position = 'bottom')
#-------------------------------------------------------------------------------
results <- drcmd(Y,A,X, hal_ind=hal_ind,
      sl_learners = sl_learners,
      eem_ind=eem_ind)
