\name{xeno.draw.autocor}
\alias{xeno.draw.autocor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for drawing within-tumor autocorrelation from a fitted model
}
\description{
Drawing function for different lags of within-tumor residual autocorrelation. 
The autocorrelation plots are shown as lines for quick visualization of the 
direction and magnitude of the autocorrelation at each lag. High absolute 
correlation at low lags ought to be avoided, though lag 0 will always 
naturally produce correlation of 1.
}
\usage{
xeno.draw.autocor(fit, orig_data, responsename = "Response", 
treatmentname = "Treatment", tpname = "Timepoint", 
idname = "Tumor_id", maintitle = "Autocorrelation")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
Mixed-effects model object fit with lme4 (mer).
}
  \item{orig_data}{
The original data.frame used in the fitting of the mixed-effects model.
}
  \item{responsename}{
Column name with the response values in the data.frame. Defaults to "Response".
}
  \item{treatmentname}{
Column name with the binary treatment group indicators in the data.frame. 
Defaults to "Treatment".
}
  \item{tpname}{
Column name with the time points in the data.frame. Defaults to "Timepoint".
}
  \item{idname}{
Column name with the individual tumor or animal labels in the data.frame. 
Defaults to "Tumor_id".
}
  \item{maintitle}{
Text title for the figure.
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
\keyword{ autocorrelation }
