rm(list=ls())
source('../R/utils.R')
source('../R/drcmd.R')
source('../R/nuis.R')
source('../R/methods.R')

#-------------------------------------------------------------------------------
# version 1: with sbp
# Data cleaning
nhefs <- causaldata::nhefs
# Get set of covariates for treatment + outcome regressions
covs <- c('sex','age','race','education','smokeyrs','smokeintensity',
          'active','exercise','wt71','sbp')
# Set outcome and treatment
outcome <- 'wt82_71' ; treatment <- 'qsmk'
# Set factor variables
factor_vars <- c('sex','race','education',
                 'active','exercise')
# Filter out the vars we want and only keep complete cases
analysis_data <- nhefs %>%
  select(covs, outcome,treatment) %>%
  mutate_at(all_of(factor_vars), as.factor) 

X <- analysis_data[,covs] %>% 
  mutate(across(everything(), ~ replace(., is.na(.), Inf)))
X <- as.data.frame(model.matrix(~ . -1, data=X)) %>%
  mutate(sbp = ifelse(sbp==Inf,NA,sbp))

A <- as.integer(analysis_data[[treatment]])
Y <- as.double(analysis_data[[outcome]])

default_learners <- c('SL.glmnet')

drcmd_res = drcmd(Y,A,X,
      default_learners=default_learners)
#-------------------------------------------------------------------------------
# version 2: without sbp

nhefs <- causaldata::nhefs
# Get set of covariates for treatment + outcome regressions
covs <- c('sex','age','race','education','smokeyrs','smokeintensity',
          'active','exercise','wt71')
# Set outcome and treatment
outcome <- 'wt82_71' ; treatment <- 'qsmk'
# Set factor variables
factor_vars <- c('sex','race','education',
                 'active','exercise')
# Filter out the vars we want and only keep complete cases
analysis_data <- nhefs %>%
  select(covs, outcome,treatment) %>%
  mutate_at(all_of(factor_vars), as.factor) 

X <- analysis_data[,covs] %>% 
  mutate(across(everything(), ~ replace(., is.na(.), Inf)))
X <- as.data.frame(model.matrix(~ . -1, data=X)) 

A <- as.integer(analysis_data[[treatment]])
Y <- as.double(analysis_data[[outcome]])

default_learners <- c('SL.glmnet')

drcmd_res = drcmd(Y,A,X,
                  default_learners=default_learners)
#-------------------------------------------------------------------------------

# Filter out the vars we want and only keep complete cases
analysis_data <- nhefs %>%
  select(covs, outcome,treatment) %>%
  mutate_at(all_of(factor_vars), as.factor) %>%
  filter(!is.na(wt82_71))

X <- analysis_data[,covs] 
X <- as.data.frame(model.matrix(~ . -1, data=X)) 

A <- as.integer(analysis_data[[treatment]])
Y <- as.double(analysis_data[[outcome]])

default_learners <- c('SL.gam','SL.glm')

drcmd_res = drcmd(Y,A,X,
                  default_learners=default_learners)
