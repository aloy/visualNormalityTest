#-------------------------------------------------------------------------------
# R functions for the paper assessing the power of lineups to assess normality
# of LMEs.
#
# Adam Loy and Heike Hofmann
# November 2014
#-------------------------------------------------------------------------------

# Calculate marginal variance of the ranefs
lev2.marginal.var <- function(.model) {
  y <- getME(.model, "y")
  X <- getME(.model, "X")
  Z <- HLMdiag:::BlockZ(.model)
  n <- nrow(X)
  ngrps   <- unname(sapply(.model@flist, function(x) length(levels(x))))
  n_rtrms <- getME(.model, "n_rtrms")
  
  # Constructing V = Cov(Y)
  sig0 <- sigma(.model)
  
  ZDZt <- crossprod( getME(.model, "A") ) # sig0^2 * crossprod( getME(.model, "A") )
  R    <- Diagonal(n) # Diagonal( n = n, x = sig0^2 )
  D    <- kronecker( Diagonal(ngrps), bdiag(VarCorr(.model)) )
  V    <- R + ZDZt
  
  # Inverting V
  V.chol <- chol( V )
  Vinv   <- chol2inv( V.chol )
  
  bse <- sig0^2 * crossprod( chol(Vinv) %*% Z %*% D ) # Marginal COV. used by Lange and Ryan
  bse.diag <- diag(bse)
  
  semat <- matrix(sqrt(bse.diag), ncol = n_rtrms, byrow = TRUE)
  
  return(semat)
}

# Calculate standardized ranefs
std_ranef <- function(.model) {
  res <- ranef(.model)[[1]]
  semat <- lev2.marginal.var(.model)
  
  RVAL <- res / semat ## ISSUE: we can get NaNs if SEs are 0
  return(RVAL)
}

# Naive simulation envelope
sim_env <- function(x, conf = .95){
  n <- length(x)
  P <- ppoints(x)
  z <- qnorm(P)
  a <- as.numeric(HLMdiag:::qqlineInfo(x)[1])
  b <- as.numeric(HLMdiag:::qqlineInfo(x)[2])
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (b/dnorm(z)) * sqrt(P * (1 - P)/n)
  fit.value <- a + b * z
  upper <- fit.value + zz * SE
  lower <- fit.value - zz * SE
  return(data.frame(lower, upper))
}

# Self test function for lineups
showPlot <- function(dframe) {
  print(qplot(sample = basement, data = dframe, stat = "qq")  + 
          facet_wrap(~ sample, ncol = 5) + 
          geom_ribbon(aes(x = x, ymin = band.lower, ymax = band.upper), alpha = .25) + 
          #	xlab("Normal Quantiles") + ylab("Sample Quantiles") +
          ylab(NULL) + xlab(NULL) + 
          theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
                axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
  )
  dframe$sample[nrow(dframe)]
}

# Simulate a random intercept model fit using lmer() with specified distributions 
# and variance components for the error terms and random effects
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