pkgname <- "XenoCat"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('XenoCat')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("XenoCat-package")
### * XenoCat-package

flush(stderr()); flush(stdout())

### Name: XenoCat-package
### Title: Categorizing mixed-effects models for the analysis of
###   xenograft/tumor experiments.
### Aliases: XenoCat-package XenoCat
### Keywords: package xenograft mixed-effects models

### ** Examples

data(mcf_low)
library(lme4)

# MCF-7 low dosage example dataset
xeno.summary(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + Timepoint:Growth 
+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Visualizing the categorized fit
par(mfrow=c(2,3))
xeno.draw.data(mcf_low)
xeno.draw.fixed(mcf_low_fit,colorgroups=FALSE)
xeno.draw.fit(mcf_low_fit)
xeno.draw.res(mcf_low_fit)
xeno.draw.autocor(mcf_low_fit)
xeno.draw.propres(mcf_low_fit)

# Non-categorizing fit
mcf_low_fit2 = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint + 
Treatment:Timepoint + (1|Tumor_id) + (0+Timepoint|Tumor_id))
# Visual comparison of the two fits
par(mfrow=c(2,2))
xeno.draw.fixed(mcf_low_fit, colorgroups=FALSE, maintitle="Categories fixed effects")
xeno.draw.fixed(mcf_low_fit2, colorgroups=FALSE, maintitle="No categories fixed effects")
xeno.draw.fit(mcf_low_fit, maintitle="Categories full fit")
xeno.draw.fit(mcf_low_fit2, maintitle="No categories full fit", legendposition=NULL)





graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("mcf_low")
### * mcf_low

flush(stderr()); flush(stdout())

### Name: mcf_low
### Title: MCF-7 LAR lower dosage dataset
### Aliases: mcf_low
### Keywords: datasets MCF-7 xenograft

### ** Examples

data(mcf_low)

head(mcf_low)
tail(mcf_low)

xeno.summary(mcf_low)



cleanEx()
nameEx("xeno.EM")
### * xeno.EM

flush(stderr()); flush(stdout())

### Name: xeno.EM
### Title: Linear categorizing EM-algorithm
### Aliases: xeno.EM
### Keywords: EM categories

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
# Fitted object with lme4
mcf_low_fit
# The data.frame with found categories (repeated for each observation for a tumor)
xeno.summary(mcf_low_EM)




cleanEx()
nameEx("xeno.cat")
### * xeno.cat

flush(stderr()); flush(stdout())

### Name: xeno.cat
### Title: Function for extracting the individual Growth-covariate values.
### Aliases: xeno.cat
### Keywords: categories

### ** Examples


data(mcf_low)
library(lme4)

# MCF-7 low dosage example dataset
xeno.summary(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + Timepoint:Growth 
+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Obtaining categories for the tumor labels:
xeno.cat(mcf_low_fit)




cleanEx()
nameEx("xeno.cat.fuzzy")
### * xeno.cat.fuzzy

flush(stderr()); flush(stdout())

### Name: xeno.cat.fuzzy
### Title: Mixture Gaussian model fit for distinguishing subgroups of the
###   probabilistic categorization
### Aliases: xeno.cat.fuzzy
### Keywords: categories

### ** Examples


data(mcf_low)
library(lme4)

# MCF-7 low dosage example dataset
xeno.summary(mcf_low)

# Probabilistic categorizing fit
mcf_low_EM_prob = xeno.EM(data=mcf_low, formula=Response ~ 1 + 
Treatment + Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + 
(0+Timepoint|Tumor_id), discriminate=FALSE)

mcf_low_fit_prob = lmer(data=mcf_low_EM_prob, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Obtaining probabilistic categories for the tumor labels
xeno.cat(mcf_low_fit_prob)

# Using mixture gaussian fits to the probabilistic categories
library(mixtools)
par(mfrow=c(1,2))

set.seed(123)
mcf_low_EM_prob2 = xeno.cat.fuzzy(mcf_low_fit_prob)
mcf_low_3groups = lmer(data=mcf_low_EM_prob2, Response ~ 1 + Treatment + Timepoint:Growth 
+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id)) 
xeno.draw.cat(mcf_low_3groups)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("xeno.diag")
### * xeno.diag

flush(stderr()); flush(stdout())

### Name: xeno.diag
### Title: Diagnostics function for categories and fixed effect estimates
###   in XenoCat-package
### Aliases: xeno.diag
### Keywords: ~kwd1 ~kwd2

### ** Examples

data(mcf_low)
library(lme4)


# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + Timepoint:Growth 
+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Running the diagnostics-function
xeno.diag(mcf_low_fit,
model = Response ~ 1 + Treatment + Timepoint:Growth 
+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))





cleanEx()
nameEx("xeno.draw.EM")
### * xeno.draw.EM

flush(stderr()); flush(stdout())

### Name: xeno.draw.EM
### Title: Function for drawing categorizing EM-algorithm iterations
### Aliases: xeno.draw.EM
### Keywords: visualization EM

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Drawing EM-iterations
xeno.draw.EM(
	Response ~ 1 + Treatment + Timepoint:Growth + Treatment:Timepoint:Growth + 
	(1|Tumor_id) + (0+Timepoint|Tumor_id),
	data=mcf_low,
	mfrow=c(2,2))





cleanEx()
nameEx("xeno.draw.HPDfixef")
### * xeno.draw.HPDfixef

flush(stderr()); flush(stdout())

### Name: xeno.draw.HPDfixef
### Title: Function for drawing density plots with HPD-interval,
###   MCMC-sample mean and coefficient estimate
### Aliases: xeno.draw.HPDfixef
### Keywords: visualization fixed effects model validation significance

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Significance testing with MCMC-samples
a = xeno.fixef.pvals(mcf_low_fit,digits=3,MCMCsim=100000)
a[,1:3]
a[,4:6]

# Visualizing posterior MCMC-samples, HPD 95%, mean of samples and coefficient 
# estimate
xeno.draw.HPDfixef(mcf_low_fit)




cleanEx()
nameEx("xeno.draw.autocor")
### * xeno.draw.autocor

flush(stderr()); flush(stdout())

### Name: xeno.draw.autocor
### Title: Function for drawing within-tumor autocorrelation from a fitted
###   model
### Aliases: xeno.draw.autocor
### Keywords: visualization model validation autocorrelation

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Visual validation
par(mfrow=c(1,3))
xeno.draw.res(mcf_low_fit)
xeno.draw.autocor(mcf_low_fit)
xeno.draw.propres(mcf_low_fit)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("xeno.draw.cat")
### * xeno.draw.cat

flush(stderr()); flush(stdout())

### Name: xeno.draw.cat
### Title: Function for drawing found categories from a fitted model's
###   frame
### Aliases: xeno.draw.cat
### Keywords: visualization categories

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

xeno.draw.data(mcf_low, maintitle="MCF-7 LAR Low dosage")

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

par(mfrow=c(1,2))
# Visualizing data
xeno.draw.data(mcf_low)
# Visualizing categories
xeno.draw.cat(mcf_low_fit)





graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("xeno.draw.data")
### * xeno.draw.data

flush(stderr()); flush(stdout())

### Name: xeno.draw.data
### Title: Function for drawing data from data.frame or fitted model's
###   frame
### Aliases: xeno.draw.data
### Keywords: visualization

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Drawing the data
xeno.draw.data(mcf_low, maintitle="MCF-7 LAR Low dosage")




cleanEx()
nameEx("xeno.draw.fit")
### * xeno.draw.fit

flush(stderr()); flush(stdout())

### Name: xeno.draw.fit
### Title: Function for drawing model fit
### Aliases: xeno.draw.fit
### Keywords: visualization

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Visualizing data
xeno.draw.data(mcf_low)
# Non-categorizing fit
mcf_low_fit2 = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint + 
Treatment:Timepoint + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Visual comparison of the two fits
par(mfrow=c(2,2))
xeno.draw.fixed(mcf_low_fit, colorgroups=FALSE, maintitle="Categories fixed effects")
xeno.draw.fixed(mcf_low_fit2, colorgroups=FALSE, maintitle="No categories fixed effects")
xeno.draw.fit(mcf_low_fit, maintitle="Categories full fit")
xeno.draw.fit(mcf_low_fit2, maintitle="No categories full fit", legendposition=NULL)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("xeno.draw.fixed")
### * xeno.draw.fixed

flush(stderr()); flush(stdout())

### Name: xeno.draw.fixed
### Title: Function for drawing fitted model's fixed effects
### Aliases: xeno.draw.fixed
### Keywords: visualization fixed effects

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Visualizing data
xeno.draw.data(mcf_low)
# Non-categorizing fit
mcf_low_fit2 = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint + 
Treatment:Timepoint + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Visual comparison of the two fits
par(mfrow=c(2,2))
xeno.draw.fixed(mcf_low_fit, colorgroups=FALSE, maintitle="Categories fixed effects")
xeno.draw.fixed(mcf_low_fit2, colorgroups=FALSE, maintitle="No categories fixed effects")
xeno.draw.fit(mcf_low_fit, maintitle="Categories full fit")
xeno.draw.fit(mcf_low_fit2, maintitle="No categories full fit", legendposition=NULL)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("xeno.draw.power")
### * xeno.draw.power

flush(stderr()); flush(stdout())

### Name: xeno.draw.power
### Title: Function for drawing power curves
### Aliases: xeno.draw.power
### Keywords: visualization power

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Power simulation example with too few simulations to be actually reliable;
# Actual simulations would be suggested to be done with e.g. nsim=10000 and 
# MCMCsim=10000, more if able
power_test = xeno.sim.power(mcf_low_fit, nsim=100, MCMCsim=100, Ns=c(10,25,50), T=5)
# Testing a longer treatment period than in the original dataset
power_test_longer = xeno.sim.power(mcf_low_fit, nsim=100, MCMCsim=100, 
Ns=c(10,25,50), T=10)
# Drawing power curve
xeno.draw.power(
	powermats=list(power_test, power_test_longer), 
	powerlabs=c("Example power (too low sims)", "Longer tp example (too low sims)"),
	minN = 5,
	maxN = 50,
	legendposition="right")




cleanEx()
nameEx("xeno.draw.precision")
### * xeno.draw.precision

flush(stderr()); flush(stdout())

### Name: xeno.draw.precision
### Title: Function for drawing precision curves
### Aliases: xeno.draw.precision
### Keywords: visualization precision

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
# Non-categorizing fit
mcf_low_fit2 = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint + 
Treatment:Timepoint + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Testing precision of combining offset and slope treatment effect terms for 
# the hypothesis
combined_prec = xeno.test.precision(mcf_low_fit, K=c(0,1,0,1))
combined_prec

# Testing precision by omitting time points from the end
# Precision of the slope hypothesis
slope_prec = xeno.test.precision.fit(mcf_low_fit, K=c(0,0,0,1))
slope_prec

# Precision of the offset hypothesis
offset_prec = xeno.test.precision.fit(mcf_low_fit, K=c(0,1,0,0))
offset_prec

# Treatment-term split to be specific to each time point
# K-vectors tested separately for each time point term
tps_prec = xeno.test.precision.tps(mcf_low_fit)
tps_prec


# Drawing precision curves
xeno.draw.precision(
	fits = list(mcf_low_fit, mcf_low_fit2),
	testedK = list(c(0,1,0,0), c(0,0,0,1), c(0,1,0,1)),
	Klabels = c("Offset effect", "Slope effect", "Combined effect"),
	fitlabels = c("MCF-7 LAR low dosage categorized", c("MCF-7 LAR low 
	dosage no categories")))





cleanEx()
nameEx("xeno.draw.propres")
### * xeno.draw.propres

flush(stderr()); flush(stdout())

### Name: xeno.draw.propres
### Title: Function for drawing fitted model's proportional residuals
### Aliases: xeno.draw.propres
### Keywords: visualization model validation residuals

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Visual validation
par(mfrow=c(1,3))
xeno.draw.res(mcf_low_fit)
xeno.draw.autocor(mcf_low_fit)
xeno.draw.propres(mcf_low_fit)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("xeno.draw.ranef")
### * xeno.draw.ranef

flush(stderr()); flush(stdout())

### Name: xeno.draw.ranef
### Title: Function for drawing model random effects
### Aliases: xeno.draw.ranef
### Keywords: visualization random effects

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Visualizing the random effects
xeno.draw.ranef(mcf_low_fit)




cleanEx()
nameEx("xeno.draw.res")
### * xeno.draw.res

flush(stderr()); flush(stdout())

### Name: xeno.draw.res
### Title: Function for drawing fitted model's residuals
### Aliases: xeno.draw.res
### Keywords: visualization model validation residuals

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Visual validation
par(mfrow=c(1,3))
xeno.draw.res(mcf_low_fit)
xeno.draw.autocor(mcf_low_fit)
xeno.draw.propres(mcf_low_fit)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("xeno.fixef.LRT")
### * xeno.fixef.LRT

flush(stderr()); flush(stdout())

### Name: xeno.fixef.LRT
### Title: Function for computing fixed effects' Likelihood Ratio Test
###   statistical significance
### Aliases: xeno.fixef.LRT
### Keywords: fixed effects significance model validation

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Testing significance with LRT
lrts = xeno.fixef.LRT(mcf_low_fit, digits=3)
names(lrts) = names(fixef(mcf_low_fit))
lrts




cleanEx()
nameEx("xeno.fixef.MCMC")
### * xeno.fixef.MCMC

flush(stderr()); flush(stdout())

### Name: xeno.fixef.MCMC
### Title: Fixed effects statistical significance using MCMC
### Aliases: xeno.fixef.MCMC
### Keywords: ~kwd1 ~kwd2

### ** Examples


data(mcf_low)
library(lme4)


# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + Timepoint:Growth 
+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Running the default MCMC inference
xeno.fixef.MCMC(mcf_low_fit)




cleanEx()
nameEx("xeno.fixef.perm")
### * xeno.fixef.perm

flush(stderr()); flush(stdout())

### Name: xeno.fixef.perm
### Title: Permutation of the treatment labels for statistical significance
### Aliases: xeno.fixef.perm
### Keywords: fixed effects significance model validation

### ** Examples


data(mcf_low)
library(lme4)

# MCF-7 low dosage example dataset
xeno.summary(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + Timepoint:Growth 
	+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
	Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# MCMC p-values
xeno.fixef.pvals(mcf_low_fit)
# p-values according to permutation
xeno.fixef.perm(mcf_low_fit)




cleanEx()
nameEx("xeno.fixef.pvals")
### * xeno.fixef.pvals

flush(stderr()); flush(stdout())

### Name: xeno.fixef.pvals
### Title: Function for estimating statistical significance of model's
###   fixed effects according to MCMC and/or LRT
### Aliases: xeno.fixef.pvals
### Keywords: fixed effects significance model validation

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Significance testing with MCMC-samples
a = xeno.fixef.pvals(mcf_low_fit,digits=3,MCMCsim=100000)
a[,1:3]
a[,4:6]

# Visualizing posterior MCMC-samples, HPD 95%, mean of samples and coefficient 
# estimate
xeno.draw.HPDfixef(mcf_low_fit)





cleanEx()
nameEx("xeno.formula")
### * xeno.formula

flush(stderr()); flush(stdout())

### Name: xeno.formula
### Title: Easy XenoCat model formulation with several handy options
### Aliases: xeno.formula
### Keywords: ~kwd1 ~kwd2

### ** Examples


# Default model formulation returned
xeno.formula()
# Omit b1,b2 - no target size included
xeno.formula(target=FALSE)
# Default non-categorizing model
xeno.formula(cat=FALSE)
# Default non-categorizing model for an experiment without target size
xeno.formula(cat=FALSE, target=FALSE)


# MCF-7 LAR low dosage example dataset
data(mcf_low)
# Model
frml = xeno.formula()

mcf_low_EM = xeno.EM(data=mcf_low, formula=frml)
mcf_low_fit = lmer(data=mcf_low_EM, frml)
# Fitted object with lme4
mcf_low_fit






cleanEx()
nameEx("xeno.nlmer.EM")
### * xeno.nlmer.EM

flush(stderr()); flush(stdout())

### Name: xeno.nlmer.EM
### Title: Non-linear categorizing EM-algorithm
### Aliases: xeno.nlmer.EM
### Keywords: non-linear EM

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Concept of fitting a non-linear model with the presented EM-algorithm procedure
# Defining a model with Michaelis-Menten kinetics. 
# Allowing random parameters of the M-M with the addition of an intercept term 
# that ought to catch the starting criteria
Model = function(Intercept, Offset, Treatment, Growth, VM, K, Timepoint) { 
	Intercept + Offset*Treatment + (Growth*VM*Timepoint/(K+Timepoint))
}
ModelGradient = deriv(body(Model)[[2]], namevec = c("Intercept", "Offset", "VM","K"), 
function.arg=Model)
starting_conditions = c(Intercept=20, Offset=-1, VM=100, K = 1)



mcf_nlmer_EM = xeno.nlmer.EM(
	formula = Response ~ ModelGradient(Intercept, Offset, Treatment, Growth, 
	VM, K, Timepoint) ~ (VM|Tumor_id) + (K|Tumor_id) + (Intercept|Tumor_id),
	data=mcf_low, 
	Model=Model, 
	ModelGradient=ModelGradient,
	start=starting_conditions, 
	verbose=TRUE,
	max.iter=10)

mcf_nlmer_fit = nlmer(Response ~ ModelGradient(Intercept, Offset, Treatment, Growth, 
VM, K, Timepoint) ~ (VM|Tumor_id) + (K|Tumor_id) + (Intercept|Tumor_id), data = mcf_nlmer_EM, 
start=starting_conditions)

par(mfrow=c(2,2))
xeno.draw.data(mcf_low)
xeno.nlmer.draw.fixed(mcf_nlmer_fit, mcf_low, Model=Model)
xeno.nlmer.draw.fit(mcf_nlmer_fit, mcf_low)
xeno.draw.res(mcf_nlmer_fit)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("xeno.nlmer.draw.fit")
### * xeno.nlmer.draw.fit

flush(stderr()); flush(stdout())

### Name: xeno.nlmer.draw.fit
### Title: Non-linear (nlmer) function for drawing full fit
### Aliases: xeno.nlmer.draw.fit
### Keywords: non-linear visualization

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Concept of fitting a non-linear model with the presented EM-algorithm procedure
# Defining a model with Michaelis-Menten kinetics. 
# Allowing random parameters of the M-M with the addition of an intercept 
#term that ought to catch the starting criteria
Model = function(Intercept, Offset, Treatment, Growth, VM, K, Timepoint) { 
	Intercept + Offset*Treatment + (Growth*VM*Timepoint/(K+Timepoint))
}
ModelGradient = deriv(body(Model)[[2]], namevec = c("Intercept", "Offset", "VM","K"), 
function.arg=Model)
starting_conditions = c(Intercept=20, Offset=-1, VM=100, K = 1)



mcf_nlmer_EM = xeno.nlmer.EM(
	formula = Response ~ ModelGradient(Intercept, Offset, Treatment, Growth, 
	VM, K, Timepoint) ~ (VM|Tumor_id) + (K|Tumor_id) + (Intercept|Tumor_id),
	data=mcf_low, 
	Model=Model, 
	ModelGradient=ModelGradient,
	start=starting_conditions, 
	verbose=TRUE,
	max.iter=10)

mcf_nlmer_fit = nlmer(Response ~ ModelGradient(Intercept, Offset, Treatment, Growth, 
VM, K, Timepoint) ~ (VM|Tumor_id) + (K|Tumor_id) + (Intercept|Tumor_id), data = mcf_nlmer_EM, 
start=starting_conditions)

par(mfrow=c(2,2))
xeno.draw.data(mcf_low)
xeno.nlmer.draw.fixed(mcf_nlmer_fit, mcf_low, Model=Model)
xeno.nlmer.draw.fit(mcf_nlmer_fit, mcf_low)
xeno.draw.res(mcf_nlmer_fit)





graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("xeno.nlmer.draw.fixed")
### * xeno.nlmer.draw.fixed

flush(stderr()); flush(stdout())

### Name: xeno.nlmer.draw.fixed
### Title: Non-linear (nlmer) fixed effect draw function
### Aliases: xeno.nlmer.draw.fixed
### Keywords: non-linear fixed effects visualization

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Concept of fitting a non-linear model with the presented EM-algorithm procedure
# Defining a model with Michaelis-Menten kinetics. 
# Allowing random parameters of the M-M with the addition of an intercept 
#term that ought to catch the starting criteria
Model = function(Intercept, Offset, Treatment, Growth, VM, K, Timepoint) { 
	Intercept + Offset*Treatment + (Growth*VM*Timepoint/(K+Timepoint))
}
ModelGradient = deriv(body(Model)[[2]], namevec = c("Intercept", "Offset", "VM","K"), 
function.arg=Model)
starting_conditions = c(Intercept=20, Offset=-1, VM=100, K = 1)



mcf_nlmer_EM = xeno.nlmer.EM(
	formula = Response ~ ModelGradient(Intercept, Offset, Treatment, Growth, 
	VM, K, Timepoint) ~ (VM|Tumor_id) + (K|Tumor_id) + (Intercept|Tumor_id),
	data=mcf_low, 
	Model=Model, 
	ModelGradient=ModelGradient,
	start=starting_conditions, 
	verbose=TRUE,
	max.iter=10)

mcf_nlmer_fit = nlmer(Response ~ ModelGradient(Intercept, Offset, Treatment, Growth, 
VM, K, Timepoint) ~ (VM|Tumor_id) + (K|Tumor_id) + (Intercept|Tumor_id), data = mcf_nlmer_EM, 
start=starting_conditions)

par(mfrow=c(2,2))
xeno.draw.data(mcf_low)
xeno.nlmer.draw.fixed(mcf_nlmer_fit, mcf_low, Model=Model)
xeno.nlmer.draw.fit(mcf_nlmer_fit, mcf_low)
xeno.draw.res(mcf_nlmer_fit)





graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("xeno.norm.start")
### * xeno.norm.start

flush(stderr()); flush(stdout())

### Name: xeno.norm.start
### Title: Function for response transformations and normalization
### Aliases: xeno.norm.start
### Keywords: normalization transform

### ** Examples


data(mcf_low)
library(lme4)

# MCF-7 low dosage example dataset
xeno.summary(mcf_low)

# Log-transform of response values, wich are normalized by the first measured value
mcf_low_new = xeno.norm.start(mcf_low)

# Dropping out the first time point because it is no longer informative
mcf_low_new = subset(mcf_low_new, mcf_low_new["Timepoint"]>0)
mcf_low_new["Timepoint"] = mcf_low_new["Timepoint"] - 1

# Categorizing fit for the log-transformed, normalized response
mcf_low_EM_lognorm = xeno.EM(data=mcf_low_new, formula=LogResponse ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id),
responsename="LogResponse")

mcf_low_fit_lognorm = lmer(data=mcf_low_EM_lognorm, LogResponse ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

#  Drawing the fit to the transformed responses
xeno.draw.fit(mcf_low_fit_lognorm, responsename="LogResponse")

# Normalization did not change the trend of the identified categories;
# Less tumors identified as non-growing in the control group
xeno.test.cat(mcf_low_fit_lognorm)




cleanEx()
nameEx("xeno.sim.power")
### * xeno.sim.power

flush(stderr()); flush(stdout())

### Name: xeno.sim.power
### Title: Function for power simulation for chosen tumor N and time scale
### Aliases: xeno.sim.power
### Keywords: simulation power

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Power simulation example with too few simulations to be actually reliable;
# Actual simulations would be suggested to be done with nsim=10000 and 
# MCMCsim=10000, more if able
power_test = xeno.sim.power(mcf_low_fit, nsim=100, MCMCsim=100, Ns=c(10,25,50), T=5)
# Testing a longer treatment period than in the original dataset
power_test_longer = xeno.sim.power(mcf_low_fit, nsim=100, MCMCsim=100, 
Ns=c(10,25,50), T=10)
# Drawing power curve
xeno.draw.power(
	powermats=list(power_test, power_test_longer), 
	powerlabs=c("Example power (too low sims)", "Longer tp example (too low sims)"),
	minN = 5,
	maxN = 50,
	legendposition="right")




cleanEx()
nameEx("xeno.sim.tumor")
### * xeno.sim.tumor

flush(stderr()); flush(stdout())

### Name: xeno.sim.tumor
### Title: Function for simulating a tumor according to a lme4 model fit
### Aliases: xeno.sim.tumor
### Keywords: simulation power

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

par(mfrow=c(2,2))
# Testing simulated datas of different lengths (original length T=5)
Ts = c(4:7)
for(i in 1:4){
	simulated_data = xeno.sim.tumor(mcf_low_fit, N=10, T=Ts[i], 
	responsename="Response")
	simulated_fit = lmer(data=simulated_data, Response ~ 1 + Treatment + 
	Timepoint:Growth + 	Treatment:Timepoint:Growth + (1|Tumor_id) + 
	(0+Timepoint|Tumor_id))
	xeno.draw.data(simulated_fit, maintitle="Simulated Data")
}




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("xeno.summary")
### * xeno.summary

flush(stderr()); flush(stdout())

### Name: xeno.summary
### Title: Summary function for a data.frame of a tumor growth experiment
### Aliases: xeno.summary
### Keywords: summary

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Check top and bottom of the data.frame
head(mcf_low)
tail(mcf_low)

# MCF-7 LAR low dosage example dataset
xeno.summary(mcf_low)




cleanEx()
nameEx("xeno.test.autocor")
### * xeno.test.autocor

flush(stderr()); flush(stdout())

### Name: xeno.test.autocor
### Title: Function for computing within-tumor autocorrelation
### Aliases: xeno.test.autocor
### Keywords: autocorrelation model validation

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Autocorrelation values of within-tumor residuals
autocors = xeno.test.autocor(mcf_low_fit)
autocors

# Visual validation
par(mfrow=c(1,3))
xeno.draw.res(mcf_low_fit)
xeno.draw.autocor(mcf_low_fit)
xeno.draw.propres(mcf_low_fit)





graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("xeno.test.cat")
### * xeno.test.cat

flush(stderr()); flush(stdout())

### Name: xeno.test.cat
### Title: Function for extracting 2x2 contingency table of Growth
###   categories vs. Treatment groups for Fisher's Exact Test
### Aliases: xeno.test.cat
### Keywords: categories fisher's exact test

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# 2x2 contingency table for Fisher's exact test
groups = xeno.test.cat(mcf_low_fit)
groups




cleanEx()
nameEx("xeno.test.endcor")
### * xeno.test.endcor

flush(stderr()); flush(stdout())

### Name: xeno.test.endcor
### Title: Correlation of end-point measurements
### Aliases: xeno.test.endcor
### Keywords: correlation end-point measurements

### ** Examples


data(mcf_low)
library(lme4)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + Timepoint:Growth 
	+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
	Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Both produce the same plot and result
xeno.test.endcor(mcf_low, legendposition="topleft")
xeno.test.endcor(mcf_low_fit, legendposition="topleft")





cleanEx()
nameEx("xeno.test.marker.data")
### * xeno.test.marker.data

flush(stderr()); flush(stdout())

### Name: xeno.test.marker.data
### Title: Function for extracting marker data from a data.frame
### Aliases: xeno.test.marker.data
### Keywords: biomarker

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + Timepoint:Growth 
+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Extracting KI-67 markers according to the Growth subcategories
KI67_growth = xeno.test.marker.data(mcf_low_EM, componentname="Growth", 
other="KI67", verbose=FALSE)
KI67_growth

# Could be extracted also by taking Growth from the model frame and KI67 from the 
# data.frame
KI67_growth = xeno.test.marker.fit(mcf_low_fit, mcf_low_EM, componentname="Growth", 
other="KI67", verbose=FALSE)
KI67_growth




cleanEx()
nameEx("xeno.test.marker.fit")
### * xeno.test.marker.fit

flush(stderr()); flush(stdout())

### Name: xeno.test.marker.fit
### Title: Function for extracting marker data from a fitted model
### Aliases: xeno.test.marker.fit
### Keywords: biomarker

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Extracting KI-67 markers according to the Growth subcategories
KI67_growth = xeno.test.marker.data(mcf_low_EM, componentname="Growth", other="KI67", 
verbose=FALSE)
KI67_growth

# Could be extracted also by taking Growth from the model frame and KI67 from the 
# data.frame
KI67_growth = xeno.test.marker.fit(mcf_low_fit, mcf_low_EM, componentname="Growth", 
other="KI67", verbose=FALSE)
KI67_growth




cleanEx()
nameEx("xeno.test.precision")
### * xeno.test.precision

flush(stderr()); flush(stdout())

### Name: xeno.test.precision
### Title: Function for extracting precision of fitted model for a
###   hypothesis of choice
### Aliases: xeno.test.precision
### Keywords: precision

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + Timepoint:Growth 
+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
# Non-categorizing fit
mcf_low_fit2 = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint + 
Treatment:Timepoint + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Testing precision of combining offset and slope treatment effect terms for the 
# hypothesis
combined_prec = xeno.test.precision(mcf_low_fit, K=c(0,1,0,1))
combined_prec

# Testing precision by omitting time points from the end
# Precision of the slope hypothesis
slope_prec = xeno.test.precision.fit(mcf_low_fit, K=c(0,0,0,1))
slope_prec

# Precision of the offset hypothesis
offset_prec = xeno.test.precision.fit(mcf_low_fit, K=c(0,1,0,0))
offset_prec

# Treatment-term split to be specific to each time point
# K-vectors tested separately for each time point term
tps_prec = xeno.test.precision.tps(mcf_low_fit)
tps_prec


# Drawing precision curves
xeno.draw.precision(
	fits = list(mcf_low_fit, mcf_low_fit2),
	testedK = list(c(0,1,0,0), c(0,0,0,1), c(0,1,0,1)),
	Klabels = c("Offset effect", "Slope effect", "Combined effect"),
	fitlabels = c("MCF-7 LAR low dosage categorized", c("MCF-7 LAR low dosage no 
	categories")))




cleanEx()
nameEx("xeno.test.precision.fit")
### * xeno.test.precision.fit

flush(stderr()); flush(stdout())

### Name: xeno.test.precision.fit
### Title: Function for extracting precision for different time periods for
###   vector K from a fitted lme4 object (mer)
### Aliases: xeno.test.precision.fit
### Keywords: precision

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + Timepoint:Growth 
+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
# Non-categorizing fit
mcf_low_fit2 = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint + 
Treatment:Timepoint + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Testing precision of combining offset and slope treatment effect terms for the 
# hypothesis
combined_prec = xeno.test.precision(mcf_low_fit, K=c(0,1,0,1))
combined_prec

# Testing precision by omitting time points from the end
# Precision of the slope hypothesis
slope_prec = xeno.test.precision.fit(mcf_low_fit, K=c(0,0,0,1))
slope_prec

# Precision of the offset hypothesis
offset_prec = xeno.test.precision.fit(mcf_low_fit, K=c(0,1,0,0))
offset_prec

# Treatment-term split to be specific to each time point
# K-vectors tested separately for each time point term
tps_prec = xeno.test.precision.tps(mcf_low_fit)
tps_prec


# Drawing precision curves
xeno.draw.precision(
	fits = list(mcf_low_fit, mcf_low_fit2),
	testedK = list(c(0,1,0,0), c(0,0,0,1), c(0,1,0,1)),
	Klabels = c("Offset effect", "Slope effect", "Combined effect"),
	fitlabels = c("MCF-7 LAR low dosage categorized", c("MCF-7 LAR low dosage no 
	categories")))




cleanEx()
nameEx("xeno.test.precision.tps")
### * xeno.test.precision.tps

flush(stderr()); flush(stdout())

### Name: xeno.test.precision.tps
### Title: Precision for time point specific treatment effect
### Aliases: xeno.test.precision.tps
### Keywords: precision

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + Timepoint:Growth 
+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
# Non-categorizing fit
mcf_low_fit2 = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint + 
Treatment:Timepoint + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Testing precision of combining offset and slope treatment effect terms for the 
# hypothesis
combined_prec = xeno.test.precision(mcf_low_fit, K=c(0,1,0,1))
combined_prec

# Testing precision by omitting time points from the end
# Precision of the slope hypothesis
slope_prec = xeno.test.precision.fit(mcf_low_fit, K=c(0,0,0,1))
slope_prec

# Precision of the offset hypothesis
offset_prec = xeno.test.precision.fit(mcf_low_fit, K=c(0,1,0,0))
offset_prec

# Treatment-term split to be specific to each time point
# K-vectors tested separately for each time point term
tps_prec = xeno.test.precision.tps(mcf_low_fit)
tps_prec


# Drawing precision curves
xeno.draw.precision(
	fits = list(mcf_low_fit, mcf_low_fit2),
	testedK = list(c(0,1,0,0), c(0,0,0,1), c(0,1,0,1)),
	Klabels = c("Offset effect", "Slope effect", "Combined effect"),
	fitlabels = c("MCF-7 LAR low dosage categorized", c("MCF-7 LAR low dosage no 
	categories")))




cleanEx()
nameEx("xeno.test.ran")
### * xeno.test.ran

flush(stderr()); flush(stdout())

### Name: xeno.test.ran
### Title: Function for extracting random effect components of specific
###   tumor fits and test them for correlation with defined marker values
### Aliases: xeno.test.ran
### Keywords: random effects biomarker correlation

### ** Examples

# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Testing random intercept values to the corresponding KI-67 measurements
xeno.test.ran(mcf_low_fit, orig_data=mcf_low, rand_component="(Intercept)", 
other="KI67")
# Pearson correlation by default, Spearman rank correlation with cormethod:
xeno.test.ran(mcf_low_fit, orig_data=mcf_low, rand_component="(Intercept)", 
other="KI67", 
cormethod="spearman")

# Testing random slope values to the corresponding KI-67 measurements
xeno.test.ran(mcf_low_fit, orig_data=mcf_low, rand_component="Timepoint", 
other="KI67")
# Pearson correlation by default, Spearman rank correlation with cormethod:
xeno.test.ran(mcf_low_fit, orig_data=mcf_low, rand_component="Timepoint", 
other="KI67", 
cormethod="spearman")




### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
