# Compare performance of drcmd to tmle, drtmle, and AIPW when outcome is MAR
#-------------------------------------------------------------------------------
# set up main params

n <- 1e3
#-------------------------------------------------------------------------------
# Params for dtcmd

hal_ind <- FALSE
eem_ind <- TRUE
sl_learners <- c('SL.glm','SL.glm.interaction','SL.earth')
k <- 1
#-------------------------------------------------------------------------------

drcmdests <- rep(NA,100) ; drcmdtimes <- rep(NA,100)
tmleests <- rep(NA,100) ; tmletimes <- rep(NA,100)
aipwests <- rep(NA,100) ; aipwtimes <- rep(NA,100)
drtmleests <- rep(NA,100) ; drtmletimes <- rep(NA,100)
for (ss in 1:100) {
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
                            k_split = k)$fit()$summary()
  aipwests[ss] <- tempaipw$estimates$RD[1]
  temp <- toc()
  aipwtimes[ss] <- temp$toc - temp$tic

  # drcmd
  tic()
  drcmdests[ss] <- drcmd(Y,A,X, hal_ind=hal_ind,
                        sl_learners = sl_learners,
                        eem_ind=eem_ind)$results$estimates$psi_hat_ate
  temp <- toc()
  drcmdtimes[ss] <- temp$toc - temp$tic

  # tmle
  tic()
  tmleests[ss] <- tmle::tmle(Y,A,X, Delta=R,
                             Q.SL.library = sl_learners,
                             g.SL.library = sl_learners,
                             g.Delta.SL.library = sl_learners)$estimates$ATE$psi
  temp <- toc()
  tmletimes[ss] <- temp$toc - temp$tic

  # drtmle
  tic()
  tempdrtmle <- drtmle::drtmle(W=X,A=A,Y=Y,
                               SL_Q = sl_learners,
                               SL_g = sl_learners,
                               SL_Qr = sl_learners,
                               SL_gr = sl_learners,
                               a_0=c(1,0))
  drtmleests[ss] <- tempdrtmle$drtmle$est[1] - tempdrtmle$drtmle$est[2]
  temp <- toc()
  drtmletimes[ss] <- temp$toc - temp$tic
}
#-------------------------------------------------------------------------------
# Plots for estimation

# reshape to wide format, removng est from each name
results <- data.frame(drcmdests,tmleests,aipwests,drtmleests) %>%
  pivot_longer(cols=everything(),names_to='method',values_to='est') %>%
  mutate(method = case_when(method=='drcmdests' ~ 'drcmd',
                            method=='tmleests' ~ 'tmle',
                            method=='aipwests' ~ 'AIPW',
                            method=='drtmleests' ~ 'drtmle'))

results %>% ggplot(aes(x=est,y=method,fill=method)) +
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
#-------------------------------------------------------------------------------


