#-------------------------------------------------------------------------------
# Script simulating from a random intercept model for the investigation
# of the power of Q-Q plots for LMEs.
#
# Adam Loy
# November 2014
#
# RADON DATA
#
# References: - Gelman, A & Pardoe, I (2006). Bayesian measures of explained
#  			variance and pooling in multilevel (hierarchical) models.
#				Technometrics, 48(2), 241--251
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Preliminaries
#-------------------------------------------------------------------------------

# required to TS test bands
library(mvtnorm)    # for the multivariate normal distribution
library(robustbase) # for robust estimates for the mean and sd


library(ggplot2)
library(lme4)     # for modeling
library(HLMdiag)  # for residuals
library(plyr)     # for formating data for lineups
library(reshape2) # for formating data for lineups
library(stringr) # for formating data for lineups


source("functions.R") 

#-------------------------------------------------------------------------------
# Simulation
#-------------------------------------------------------------------------------

# Base model
fm <- lmer(log.radon ~ basement + uranium + (1 | county), data = radon)

# Trial simulations
sim1 <- sim_random_intercept_model(.mod = fm, nsim = 1, e.dsn = "exp", 
                                   b0.dsn = "t", sigma.err = 2, sigma.b0 = 1)

mod1 <- refit(sim1, object = fm)
re1  <- std_ranef(mod1)

# Everything seems to be working based on my tests, and we don't hit the issue
# with the standardized random effects using this set up.

#-------------------------------------------------------------------------------
# Lineups with TS test bands
#-------------------------------------------------------------------------------

# this will create lineups for all of the simulation settings
# TODO: fully integrate with your old ShowNext function to record
#       our attempts.
simLineup(.mod = fm, e.dsn = "exp", alt.b0.dsn = "exp", null.b0.dsn = "norm", 
          sigma.err = 2, sigma.b0 = 1)