#-------------------------------------------------------------------------------
# Script constructing lineups of the Q-Q plots.
#
# Adam Loy
# May 2013
#
# RADON DATA
# 
#
# References: - Gelman, A & Pardoe, I (2006). Bayesian measures of explained
#				variance and pooling in multilevel (hierarchical) models.
#				Technometrics, 48(2), 241--251
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Preliminaries
#-------------------------------------------------------------------------------

setwd("~/Documents/Thesis/Dissertation/sociology_chapter/")
library(ggplot2)
library(lme4)     # for modeling
library(HLMdiag)  # for residuals
library(nullabor) # for lineups
library(plyr)
library(reshape2)
library(stringr)
library(mvtnorm)

BlockZ <- function(object) {
  Z <- getME(object, "Z")
  
  grp.size <- table(object@flist)
  ngrps <- length(grp.size)
  nranef <- dim(ranef(object)[[1]])[2]
  
  base.ord <- seq(from = 1, by = ngrps, length.out = nranef)
  ord <- base.ord + rep(0:(ngrps - 1), each = nranef)
  
  perm.mat <- t(as(ord, "pMatrix"))
  
  return(Z %*% perm.mat)
}


lev2.marginal.var <- function(.model) {
  y <- .model@y
  X <- getME(.model, "X")
  Z <- BlockZ(.model)
  n <- nrow(X)
  ngrps <- unname(sapply(.model@flist, function(x) length(levels(x))))
  
  # Constructing V = Cov(Y)
  sig0 <- attr(VarCorr(.model), "sc") # sigma(.model)
  
  ZDZt <- sig0^2 * crossprod( .model@A )
  R    <- Diagonal( n = n, x = sig0^2 )
  D    <- kronecker( Diagonal(ngrps), bdiag(VarCorr(.model)) )
  V    <- Diagonal(n) + ZDZt
  
  # Inverting V
  V.chol <- chol( V )
  Vinv   <- chol2inv( V.chol )

  bse <- crossprod( chol(Vinv) %*% Z %*% D ) # Marginal COV. used by Lange and Ryan
  bse.diag <- diag(bse)

  semat <- matrix(sqrt(bse.diag), ncol = 2, byrow = TRUE)

  return(semat)
}

std_ranef <- function(.model) {
	res <- ranef(.model)[[1]]
	semat <- lev2.marginal.var(.model)
	
	RVAL <- res / semat
	return(RVAL)
}

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

sim_t_hlm <- function(.mod) {
	vc <- VarCorr( .mod )
	D  <- as.matrix( bdiag(vc) )
	sig.e <- sigma(.mod)
	
	dims <- .mod@dims
	n <- dims[["n"]]
	m <- dims[["q"]] / nrow(D)

	## normal errors
	e  <- rnorm(n = n, mean = 0, sd = sig.e)

	## mutlivariate t random effects
	b <- rmvt(n = m, sigma = D, df = 3)

#	print(c(ad.test(b[,1])$p.value,ad.test(b[,2])$p.value) )

	## Generating y
	bvec <- c(b[,1], b[,2])
	y <- getME(.mod, "X") %*% fixef(.mod) + getME(.mod, "Z") %*% bvec + e
	
	return( list(y=as.numeric(y), b=b) )
}

sim_indep_ranef_hlm <- function(.mod, nsim, e.dsn, b0.dsn, b1.dsn, sigma.err, sigma.b0, sigma.b1){
  vc <- VarCorr( .mod )
	
	dims <- .mod@dims
	n <- dims[["n"]]
	m <- dims[["q"]] / dims[["nt"]]
	
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
		b0  <- rnorm(n = nsim * m, mean = 0, sd = sigma.b0)
	} 
	if(b0.dsn == "t") {
		b0  <- (sigma.b0 / sqrt(3)) * rt(n = nsim * m, df = 3)
	}
	if(b0.dsn == "exp") {
		b0  <- sigma.b0 * ( rexp(n = nsim * m) - 1 )
	}
	b0 <- matrix(b0, nc = nsim)

	## Simulating random slope
	if(b1.dsn == "norm") {
		b1  <- rnorm(n = nsim * m, mean = 0, sd = sigma.b1)
	} 
	if(b1.dsn == "t") {
		b1  <- (sigma.b1 / sqrt(3)) * rt(n = nsim * m, df = 3)
	}
	if(b1.dsn == "exp") {
		b1  <- sigma.b1 * ( rexp(n = nsim * m) - 1 )
	}
	b1 <- matrix(b1, nc = nsim)
	
	## Generating y
	b <- rbind(b0, b1)
	y <- getME(.mod, "X") %*% fixef(.mod) + getME(.mod, "Z") %*% b + e
	
	y.df <- as.data.frame( as.matrix( y) )
	colnames(y.df) <- paste("sim_", 1:ncol(y.df), sep = "")
	
	return( y.df )

}

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

#load("~/Documents/Thesis/Dissertation/sociology_chapter/data/radon.RData")
load(file.choose())
fm <- lmer(log.radon ~ basement + uranium + (1 | county) + (basement - 1 | county), data = radon)


showNext <- function(whoami) {
sim_step <- sim_t_hlm(fm)
y.b.t <- sim_step$y
refit.b.t <-  refit(fm, y.b.t)
b.t <- ranef(refit.b.t)[[1]]

sim.y   <- simulate(fm, nsim = 19)                        
sim.mod <- apply(sim.y, 2, refit, object = fm)            ## a list of models



sim.b1 <- llply(sim.mod, function(x) ranef(x)[[1]][,2])   ## a list of random slopes
sim.b1 <- melt( do.call("rbind", sim.b1) )[,-2]           ## changing to a data frame
names(sim.b1) <- c("sample", "basement")                  ## setting colnames for faceting
sim.b1        <- arrange(sim.b1, sample)                  ## ordering by simulation

sim.b1$.n <- as.numeric( str_extract(sim.b1$sample, "\\d+") )
sim.b1 <- ddply(sim.b1, .(.n), transform, band = sim_env(basement), 
	x = sort(qqnorm(basement, plot.it=FALSE)$x))

b1 <- data.frame(.n=20, transform(b.t, band = sim_env(basement), 
	x = sort(qqnorm(basement, plot.it=FALSE)$x)))

tranef <- rbind(sim.b1[,c(".n", "basement", "x", "band.lower", "band.upper")],
				b1[,c(".n", "basement", "x", "band.lower", "band.upper")])
tranef$sample <- sample(20,20, replace=FALSE)[tranef$.n]
location <- tranef$sample[nrow(tranef)]
#write.csv(tranef, file=sprintf("radontranef-test/radontranef-%s-%s.csv", null, location))
correct <- showPlot(tranef)
choice <- scan()
response <- paste(choice, collapse=",")
if (file.exists("tranef-results.csv")) {
	df <-read.csv("tranef-results.csv")
	k <- max(df$step)+1
} else k <- 1
res <- data.frame(id=whoami, step=k, b0=sim_step$b[,1],b1=sim_step$b[,2], response=response, location_no=location)
write.table(res, file="tranef-results.csv", sep=",", append=file.exists("tranef-results.csv"), row.names=FALSE, col.names=!file.exists("tranef-results.csv"))
}

for (i in 1:10)
	showNext("heike")

#df <- read.csv("tranef-results.csv")
df <- read.csv(file.choose(), stringsAsFactors=FALSE)
unique(df[,c("response", "location_no", "step")])
library(plyr)
library(nortest)
res <- ddply(df, .(step, id), summarize, response=response[1], 
	location_no=location_no[1],
	num_answers=unlist(llply(strsplit(as.character(response), ","), length)),
	correct=location_no %in% as.vector(strsplit(as.character(response), ",")[[1]]),
	ad0 = ad.test(b0)$p.value,
	ad1 = ad.test(b1)$p.value
)
mean(res$correct/res$num_answers)

qplot(ad1, 1-correct/num_answers, data=res, geom="jitter")
qplot(ad0, 1-correct/num_answers, data=res, geom="jitter")
