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



showPlot <- function(dframe) {
  print(qplot(x, y, data = dframe)  + 
          facet_wrap(~ sample, ncol = 5) + 
          geom_ribbon(aes(x = x, ymin = lower, ymax = upper), alpha = .25) + 
          #  xlab("Normal Quantiles") + ylab("Sample Quantiles") +
          ylab(NULL) + xlab(NULL) + 
          theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
                axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
  )
  dframe$sample[nrow(dframe)]
}

simLineup <- function(.mod, e.dsn, alt.b0.dsn, null.b0.dsn, sigma.err, sigma.b0) {
  sim_step <- sim_random_intercept_model(.mod, nsim = 1, e.dsn = e.dsn, b0.dsn = alt.b0.dsn, 
                                         sigma.err = sigma.err, sigma.b0 = sigma.b0)
  refit.b.alt <-  refit(.mod, sim_step)
  b.alt <- std_ranef(refit.b.alt)[[1]]
  
  sim.y   <- sim_random_intercept_model(.mod, nsim = 19, e.dsn = e.dsn, b0.dsn = null.b0.dsn,
                                        sigma.err = sigma.err, sigma.b0 = sigma.b0)                        
  sim.mod <- apply(sim.y, 2, refit, object = fm)            ## a list of models
  
  
  sim.b1 <- llply(sim.mod, function(x) std_ranef(x)[[1]])   ## a list of random intercepts
  sim.b1 <- melt( do.call("rbind", sim.b1) )[,-2]           ## changing to a data frame
  names(sim.b1) <- c("sample", "county")                    ## setting colnames for faceting
  sim.b1        <- arrange(sim.b1, sample)                  ## ordering by simulation
  
  sim.b1$.n <- as.numeric( str_extract(sim.b1$sample, "\\d+") )
  sim.b1 <- ddply(sim.b1, .(.n), function(df) QQ.UNcb(df$county, plot = FALSE))
  #   sim.b1 <- ddply(sim.b1, .(.n), transform, band = do.call("cbind", QQ.UNcb(county, plot = FALSE)), 
  #                   x = sort(qqnorm(county, plot.it=FALSE)$x),
  #                   y = sort(qqnorm(county, plot.it=FALSE)$y))
  #   sim.b1 <- sim.b1[,-2]
  
  b1 <- data.frame(.n = 20, QQ.UNcb(b.alt, plot = FALSE))
  #   b1 <- data.frame(.n=20, band = do.call("cbind", QQ.UNcb(b.alt, plot = FALSE)), 
  #                    x = sort(qqnorm(b.alt, plot.it=FALSE)$x), 
  #                    y = sort(qqnorm(b.alt, plot.it=FALSE)$y))
  #   
  lineup_data <- rbind(sim.b1, b1)
  
  lineup_data$sample <- sample(20,20, replace=FALSE)[lineup_data$.n]
  location <- lineup_data$sample[nrow(lineup_data)]
  #write.csv(tranef, file=sprintf("radontranef-test/radontranef-%s-%s.csv", null, location))
  showPlot(lineup_data)
}

#-------------------------------------------------------------------------------
# Functions for TS test bands, I made a slight change in one
#-------------------------------------------------------------------------------


###################################
# a few functions that are useful #
###################################

Qn.scale<-function(x){
  Qn(x,finite.corr=FALSE)
}

Qn.location<-function(x){
  s_Qn(x,mu.too=TRUE)[[1]]
}


out<-function(curve,upper,lower){
  k<-length(which(curve>upper | curve<lower))
  return(k>0)
}



###############################################################
# The confidence band function                                #
# The function returns a list with the following:             #
# 1. lower,upper      the normal scale confidence bands       #
# 2. individual.alpha the final individual significance level #
###############################################################

QQ.cb<-function(x.sample,mu=0,sigma=1,M.sim=1000,alpha=0.05,plot=TRUE){
  
  n<-length(x.sample)
  upper.ci<-rep(NA,n)
  lower.ci<-rep(NA,n)
  p.value <-matrix(NA,nrow=n,ncol=M.sim)
  sim=NULL
  
  # simulate data
  for(i in 1:M.sim)
    sim<-cbind(sim,sort(runif(n)))
  
  # widen the CI to get a simultanoues 1-alpha CI
  for(i in 1:n){
    tmp<-pbeta(sim[i,],shape1=i,shape2=n+1-i)
    p.value[i,]<-apply(cbind(tmp,1-tmp),1,min)
  }
  
  critical.values<-apply(p.value,2,min)
  C.crit<-quantile(critical.values,prob=alpha)
  
  upper.ci<-qbeta(1-C.crit,shape1=1:n,shape2=n+1-(1:n))
  lower.ci<-qbeta(C.crit,shape1=1:n,shape2=n+1-(1:n))
  
  
  # now translate back to normal
  norm.upper<-qnorm(upper.ci)
  norm.lower<-qnorm(lower.ci)
  
  if(plot==TRUE){
    q.prob<-qnorm((1:n)/(n+1))
    plot(q.prob,norm.upper,type="l",col="red",ylim=c(-3,3),xlim=c(-3,3),ylab="Sample Quantile",xlab="Sample Quantile")
    lines(q.prob,norm.lower,col="red")
    z.sample<-(x.sample-mu)/sigma
    points(q.prob,sort(z.sample),pch=19,cex=0.6)
  }
  
  return(list(lower=norm.lower,upper=norm.upper))
}


###############################################################
# The confidence band function with unknown parameters        #
# The function returns a list with the following:             #
# 1. lower,upper      the normal scale confidence bands       #
# 2. individual.alpha the final individual significance level #
###############################################################


QQ.UNcb<-function(x.sample,M.sim=1000,alpha=0.05,plot=TRUE,center.func=Qn.location,scale.func=Qn.scale){
  
  
  n<-length(x.sample)
  upper.ci<-rep(NA,n)
  lower.ci<-rep(NA,n)
  p.value <-matrix(NA,nrow=n,ncol=M.sim)
  
  sim=NULL
  for(i in 1:M.sim)
    sim<-cbind(sim,sort(rnorm(n)))
  
  
  # here is the difference... scale it -- so you are estimating the values
  center<-apply(sim,2,center.func)
  scale<-apply(sim,2,scale.func)
  sim<-sweep(sweep(sim,2,center,FUN="-"),2,scale,FUN="/")
  
  
  # convert the norm to beta
  sim<-t(apply(sim,1,pnorm))
  
  # widen the CI to get a simultanoues 1-alpha CI
  for(i in 1:n){
    tmp<-pbeta(sim[i,],shape1=i,shape2=n+1-i)
    p.value[i,]<-apply(cbind(tmp,1-tmp),1,min)
  }
  
  critical.values<-apply(p.value,2,min)
  C.crit<-quantile(critical.values,prob=alpha)
  
  upper.ci<-qbeta(1-C.crit,shape1=1:n,shape2=n+1-(1:n))
  lower.ci<-qbeta(C.crit,shape1=1:n,shape2=n+1-(1:n))
  
  
  # now translate back to normal
  norm.upper<-qnorm(upper.ci)
  norm.lower<-qnorm(lower.ci)
  
  q.prob<-qnorm((1:n)/(n+1))
  z.sample<-(x.sample-center.func(x.sample))/scale.func(x.sample)
  
  if(plot==TRUE){
#     q.prob<-qnorm((1:n)/(n+1))
    plot(q.prob,norm.upper,type="l",col="red",ylim=c(-3,3),xlim=c(-3,3),ylab="Sample Quantile",xlab="Sample Quantile")
    lines(q.prob,norm.lower,col="red")
#     z.sample<-(x.sample-center.func(x.sample))/scale.func(x.sample)
    points(q.prob,sort(z.sample),cex=0.6,pch=19)
  }
  
  # I changed this for our purposes (Adam Loy)
  return(data.frame(x=q.prob, y=sort(z.sample), lower=norm.lower, upper=norm.upper))
}