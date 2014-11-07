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

library(lme4)     # for modeling
library(HLMdiag)  # for residuals

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

sim_random_intercept_model <- function(.mod, nsim, e.dsn, b0.dsn, sigma.err, sigma.b0){
  vc <- VarCorr( .mod )
  
  n <- nobs(.mod)                   # No. obs.
  nre <- unname( getME(.mod, "q") ) # No. ranefs
  
  ## Simulating error terms
  if(e.dsn == "norm") {
    e  <- rnorm(n = nsim * n, mean = 0, sd = sigma.err)
  } 
  if(e.dsn == "t") {
    e  <- (sigma.err / sqrt(3)) * rt(n = nsim * n, df = 3)
  }
  if(e.dsn == "exp") {
    e  <- sigma.err * ( rexp(n = nsim * n) - 1 )
  }
  e <- matrix(e, nc = nsim)
  
  ## Simulating random intercept
  if(b0.dsn == "norm") {
    b0  <- rnorm(n = nsim * nre, mean = 0, sd = sigma.b0)
  } 
  if(b0.dsn == "t") {
    b0  <- (sigma.b0 / sqrt(3)) * rt(n = nsim * nre, df = 3)
  }
  if(b0.dsn == "exp") {
    b0  <- sigma.b0 * ( rexp(n = nsim * nre) - 1 )
  }
  b0 <- matrix(b0, nc = nsim)
  
  ## Generating y
  y <- getME(.mod, "X") %*% fixef(.mod) + getME(.mod, "Z") %*% b0 + e
  
  y_df <- as.data.frame( as.matrix(y) )
  colnames(y_df) <- paste("sim_", 1:ncol(y_df), sep = "")
  
  return( y_df )
}


#-------------------------------------------------------------------------------
# Simulation
#-------------------------------------------------------------------------------

# Base model
fm <- lmer(log.radon ~ basement + uranium + (1 | county), data = radon)

# Trial simulations
sim_random_intercept_model(.mod = fm, nsim = 1, e.dsn = "t", b0.dsn = "norm", 
                           sigma.err = 2, sigma.b0 = 1)

sim_random_intercept_model(.mod = fm, nsim = 1, e.dsn = "norm", b0.dsn = "norm", 
                           sigma.err = 2, sigma.b0 = 1)

sim_random_intercept_model(.mod = fm, nsim = 1, e.dsn = "exp", b0.dsn = "norm", 
                           sigma.err = 2, sigma.b0 = 1)

