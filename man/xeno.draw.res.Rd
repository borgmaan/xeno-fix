\name{xeno.draw.res}
\alias{xeno.draw.res}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for drawing fitted model's residuals
}
\description{
Function for drawing the residuals of a fitted model. Can be used e.g. for finding 
outlier tumors and adjusting the model for visible structure in the residuals. 
Ideally the residuals ought to form a uniform pattern around the 0 mean.
}
\usage{
xeno.draw.res(fit, responsename = "Response", maintitle = "Residuals", 
drawline = FALSE, drawquad = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
Mixed-effects model object fit with lme4 (mer).
}
  \item{responsename}{
Column name with the response values in the data.frame. Defaults to "Response".
}
  \item{maintitle}{
Text title for the figure.
}
  \item{drawline}{
Should a linear fit of form y = a*x + b be fit and drawn with the residuals to visualize the trend. Defaults to FALSE.
}
  \item{drawquad}{
Should a quadratic trend in residuals be illustrated in a similar manner to "drawline". Defaults to FALSE.
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

# Visual validation
par(mfrow=c(1,3))
xeno.draw.res(mcf_low_fit)
xeno.draw.autocor(mcf_low_fit)
xeno.draw.propres(mcf_low_fit)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ visualization }
\keyword{ model validation }
\keyword{ residuals }
