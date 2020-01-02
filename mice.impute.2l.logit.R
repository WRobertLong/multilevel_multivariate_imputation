mice.impute.2l.logit <- function(y, ry, x, type, ...) {

###########################################################
####
#### 		Method for logistic multilevel imputation in MICE
#### 		Author: W Robert Long 
#### 		Date:	 03 July 2016
####
###########################################################

# y  		Response vector to be imputed
# ry 		Vector of missing data pattern FALSE=missing, TRUE=observed
# x  		Matrix n x p of complete covariates.
# type	Predictor matrix


	require(lme4)

	# data for call to glmer
	dt.glmer <- cbind(y[ry],x[ry, ])

	# build formula for call to glmer:
	
	# all the variables in x
	all.x.vars <- names(x)

	# the grouping variable
	group.var <- names(type[type==-2])
	# need to check there is one and only one ? Currently not allowing crossed random effects.

	# random effects variable(s)
	re.vars <- names(type[type==2])

	# fixed effects part of formula
	response.str <- "y[ry] "
	fe.varstr <- "1 "
	
	for (var in all.x.vars[!(all.x.vars %in% group.var)]) {
		fe.varstr <- paste(fe.varstr,var,sep=" + ")
	}

	# random effects part of formula
	re.varstr <- "1"
	for (var in re.vars) {
		re.varstr <- paste(re.varstr,var,sep=" + ")
	}

	re.str <- paste("+ (", re.varstr, "|" ,group.var, ")", sep="")

	varstr <- paste(response.str, "~ ",fe.varstr, re.str, sep="")
	fm <- as.formula(varstr)

	fit.glmer <- glmer(fm, data=dt.glmer, family=binomial(link =logit))

	# from this we can extract the estimated fixed effects 
	beta <- fixef(fit.glmer)

	# Now find the transpose of Cholesky decomposition of the fixed effects 
	# variance-covariance matrix
	rv <- t(chol(vcov(fit.glmer)))

	# Now we can draw beta
	beta.star <- beta + rv %*% rnorm(ncol(rv))

	# The model matrix for fixed effects should be just those rows 
	# corresponding to the missing values in the response.
	X <- as.matrix(x[!ry, ])   	

	# now we compute the vector of log odds for the 
	# missing responses contributed by the fixed effects:
	p.fe <- X %*% beta.star

	# Now for contribution of the random effects to the log odds.
	# First we need the transpose of model matrix for the random effects.
	# This requires care because not all clusters may have observations
	# in the missing data. So we should use the entire x and y objects
	# to build the model matrix using lme4 (but without fitting a model)
	# 
	dt.for.re <- cbind(y,x)

	varstr1 <- paste("y ~ ",fe.varstr, re.str, sep="")
	fm1 <-  as.formula(varstr1)

	# get the design matrix for the random effects without fitting a model
	Zt <- lmer(fm1, data = dt.for.re,doFit = FALSE)$FL$trms[[1]]$Zt 
	# lmer rather than glmer was used here because glmer produces warnings 
	# about 0 and 1 probabilities even when doFit=FALSE is used.

	# the random effects
	z <- unlist(ranef(fit.glmer))

	# Then the vector of log odds contributed by the random effects is
	p.re <- (t(Zt) %*% z)[!ry]
	# where we use [!ry] to select just those values corresponding to the missings

	# So the probabilities are 
	p <- plogis(as.vector(p.fe + p.re))

	vec <- (runif(length(p)) <= p)

	# cast to binary and convert to factor if needed
	vec[vec] <- 1
	if (is.factor(y)) {
		vec <- factor(vec, c(0, 1), levels(y))
	}
	return(vec)

}