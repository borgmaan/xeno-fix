\name{xeno.fixef.pvals}
\alias{xeno.fixef.pvals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for estimating statistical significance of model's fixed effects according 
to MCMC and/or LRT
}
\description{
Function that wraps the fixed effect estimates, standard deviation, t statistic, 
Highest Posterior Density interval, Markov-Chain Monte Carlo (MCMC) 
and alternatively Likelihood Ratio Test (LRT) p-values. The MCMC p-values are 
estimated as in Bayeen (package languageR).
}
\usage{
xeno.fixef.pvals(fit, LRT = FALSE, MCMC = TRUE, MCMCsim = 10000, 
digits = 3, id = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
Mixed-effects model object fit with lme4 (mer).
}
  \item{LRT}{
Should Likelihood-Ratio Testing (LRT) used for fixed effect statistical significance. 
Defaults to FALSE.
}
  \item{MCMC}{
Should MCMC-samples be used for fixed effect statistical significance. 
Defaults to TRUE.
}
  \item{MCMCsim}{
How many mcmcsamp-function (lme4-package) samples should be used for estimating 
the statistical significance. Defaults to 10,000.
}
  \item{digits}{
How many digits ought to be reported of the p-value. Defaults to 3.
}
  \item{id}{
R session specific id number if progression of p-value estimation should be written to a separate text file.
}
  \item{verbose}{
Should the samples generated with mcmcsamp-function be returned. Defaults to FALSE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Return value is a matrix with columns indicating fixed effect coefficient estimate,
standard error, corresponding t value and possibly LRT p-value and/or HPD 95 percent 
interval with MCMC p-value. Rows are the fixed effects of the model.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Baayen RH. Modeling data with fixed and random effects. 
In: Analyzing linguistic data, a practical introduction to statistics using R, 
Cambridge University Press; 2008. p. 242-58.

Pinheiro JC, Bates DM. Hypothesis tests for fixed-effects terms. 
In: Chambers J, Eddy W, Hardle W, Sheather S, Tierney L (editors). 
Mixed-effects models in S and S-PLUS.  Springer-Verlag; 2000. p. 87-92.

Douglas Bates, R-help, "lmer, p-values and all that", May 19th 2006,
https://stat.ethz.ch/pipermail/r-help/2006-May/094765.html

}
\author{
Teemu D Laajala <tlaajala@cc.hut.fi>
}
\note{
The MCMC diagnostics include evaluating the coherence of the estimated model 
fixed effect coefficients and the mean of the MCMC samples. The function
xeno.draw.HPDfixef can be used for this purpose. MCMC-sample generation for 
assessing statistical significance was suggested by the lme4-author Douglas Bates,
see R-help reference for further discussion.

As MCMC sample generation is a stochastic process, the p-values and HPD-intervals
have slight variation.

Additional statistical significance validation can be done with the function 
xeno.fixef.perm and the boundary suggestions with the function pamer.fnc 
provided in the LMERConvenienceFunctions-package.
}
\section{Warning }{
Likelihood-Ratio Testing can be anti-conservative for datasets with only a few 
observations. See Pinheiro & Bates for further discussion.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{xeno.draw.HPDfixef}},
\code{\link{xeno.fixef.perm}},
\code{\link{xeno.fixef.LRT}},
\code{\link{xeno.fixef.MCMC}}
}
\examples{
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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fixed effects }
\keyword{ significance }
\keyword{ model validation }