\name{xeno.draw.HPDfixef}
\alias{xeno.draw.HPDfixef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for drawing density plots with HPD-interval, MCMC-sample mean and 
coefficient estimate
}
\description{
Visualization of the MCMC-samples used for assessing fixed effect terms' 
statistical significance. The function draws the density plots using 
lattice-package, with the addition of the HPD-interval, the mean of the 
MCMC-samples and the model coefficient estimate indicated by vertical lines. 
The concordance of the mean and the estimate ought to be done along with 
verifying the trivial shape of the MCMC-sample distribution, to avoid 
multimodality problems in the esimation of the HPD interval.
}
\usage{
xeno.draw.HPDfixef(fit, MCMCsim = 10000, samples=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
Mixed-effects model object fit with lme4 (mer).
}
  \item{MCMCsim}{
How many mcmcsamp-function (lme4-package) simulations should be used.
}
 \item{samples}{
If samples generated with the mcmcsamp-function are provided are provided, density plots and intervals will be visualized according to them.
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
%% ~put references to the literature/web site here ~
}
\author{
Teemu D Laajala <tlaajala@cc.hut.fi>
}
\note{
Requires Lattice-package to draw the density plot.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{xeno.fixef.pvals}}
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
\keyword{ visualization }
\keyword{ fixed effects }
\keyword{ model validation }
\keyword{ significance }
