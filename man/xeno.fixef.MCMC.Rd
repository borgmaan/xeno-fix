\name{xeno.fixef.MCMC}
\alias{xeno.fixef.MCMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fixed effects statistical significance using MCMC
}
\description{
Function that wraps the fixed effect estimates, standard deviation, t statistic, 
Highest Posterior Density interval, Markov-Chain Monte Carlo (MCMC). 
The MCMC p-values are estimated as in Bayeen (package languageR).
}
\usage{
xeno.fixef.MCMC(fit, MCMCsim = 10000, digits = 6, 
draw = TRUE, burn.draw = FALSE, burn.cut = 0, id = NULL, 
verbose = FALSE, time = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
Model fit (lme4 mer-object)
}
  \item{MCMCsim}{
The amount of MCMC samples generated. Defaults to 10,000.
}
  \item{digits}{
The number of digits that are reported for the p-values. Defaults to 6.
}
  \item{draw}{
Should a visualization of the MCMC samples be drawn (suggested for validation). Defaults to TRUE.
}
  \item{burn.draw}{
Should a visualization of the MCMC samples be drawn according to their sampling index. Defaults to FALSE.
}
  \item{burn.cut}{
Positive integer indicating if a number of initial samples should be left out. Defaults to 0 (no samples left out).
}
  \item{id}{
R session specific id number if progression of p-value estimation should be written to a separate text file.
}
  \item{verbose}{
Should output be printed to the user. Defaults to FALSE.
}
  \item{time}{
Should a brief test permutation be run instead of running the actual permutation, so that a time estimate can be given to the user. Defaults to FALSE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{xeno.fixef.pvals}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

data(mcf_low)
library(lme4)


# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + Timepoint:Growth 
+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Running the default MCMC inference
xeno.fixef.MCMC(mcf_low_fit)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
