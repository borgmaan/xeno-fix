\name{xeno.fixef.LRT}
\alias{xeno.fixef.LRT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for computing fixed effects' Likelihood Ratio Test statistical significance
}
\description{
Function for computing Likelihood Ratio Test p-values for fixed effects. The 
model is re-fit by omitting one fixed effect at a time and then statistical 
significance is estimated by comparing the likelihoods. The model fit is done 
using Maximum Likelihood (ML) in order for the likelihoods to be comparable, 
while by default the mixed-effects models are fitted using Restricted Maximum 
Likelihood (REML). Due to issues with reliability for datasets with only a 
small amount of observations, LRT is not proposed to be the default choice.
}
\usage{
xeno.fixef.LRT(fit, digits = 6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
Mixed-effects model object fit with lme4 (mer).
}
  \item{digits}{
How many digits should be reported for the p-value.
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
Pinheiro JC, Bates DM. Hypothesis tests for fixed-effects terms. 
In: Chambers J, Eddy W, Hardle W, Sheather S, Tierney L (editors). 
Mixed-effects models in S and S-PLUS.  Springer-Verlag; 2000. p. 87-92.
}
\author{
Teemu D Laajala <tlaajala@cc.hut.fi>
}
\note{
May be unreliable for datasets with only a few observations, see Pinheiro & Bates.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fixed effects }
\keyword{ significance }
\keyword{ model validation }