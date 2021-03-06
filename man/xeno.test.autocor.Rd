\name{xeno.test.autocor}
\alias{xeno.test.autocor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for computing within-tumor autocorrelation
}
\description{
Function that computes autocorrelation for lags of {0, 1, ... , T-1} for within-tumor 
residual autocorrelation. Residual pairs are considered only within tumors, and all 
the pairs are combined in the end over all the tumors in order to compute Pearson 
correlation.
}
\usage{
xeno.test.autocor(fit, idname = "Tumor_id", tpname = "Timepoint", 
draw = FALSE, maintitle = "Autocorrelation")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
Mixed-effects model object fit with lme4 (mer).
}
  \item{idname}{
Column name with the individual tumor or animal labels in the data.frame. 
Defaults to "Tumor_id".
}
  \item{tpname}{
Column name with the time points in the data.frame. Defaults to "Timepoint".
}
  \item{draw}{
Should the autocorrelations at different lags be drawn (same as xeno.draw.autocor).
Defaults to TRUE.
}
  \item{maintitle}{
Text title for the figure.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Vector containing the correlations at different lags.
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
%%  ~~further notes~~
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

# Autocorrelation values of within-tumor residuals
autocors = xeno.test.autocor(mcf_low_fit)
autocors

# Visual validation
par(mfrow=c(1,3))
xeno.draw.res(mcf_low_fit)
xeno.draw.autocor(mcf_low_fit)
xeno.draw.propres(mcf_low_fit)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ autocorrelation }
\keyword{ model validation }
