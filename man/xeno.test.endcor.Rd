\name{xeno.test.endcor}
\alias{xeno.test.endcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Correlation of end-point measurements
}
\description{
Function for testing for correlation of the first and last measured response value. The
pairs of observations are tested over individual units by default with Pearson correlation.
}
\usage{
xeno.test.endcor(x, responsename = "Response", tpname = "Timepoint", 
idname = "Tumor_id", method = "pearson", subgroup = "Treatment", 
draw = TRUE, legendposition = "topright", 
maintitle = "End-point measurements")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Mixed-effects model object fit with lme4 (mer) or a data.frame object.
}
  \item{responsename}{
Column name with the response values in the data.frame. Defaults to "Response".
}
  \item{tpname}{
Column name with the time points in the data.frame. Defaults to "Timepoint".
}
  \item{idname}{
Column name with the individual tumor or animal labels in the data.frame. 
Defaults to "Tumor_id".
}
  \item{method}{
Method for computing correlation: should be either "pearson" or "spearman".
}
  \item{subgroup}{
Subgroup for indicating different colours in the plot. Defaults to distinguishing the
treatment and control groups.
}
  \item{draw}{
Should a scatter plot figure of the end-point measurements be drawn. Defaults to TRUE.
}
  \item{legendposition}{
Text string indicating where the legend should be placed; must be one of "bottomright",
"bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and
"center". NULL or FALSE value will omit legend.
}
  \item{maintitle}{
Text title for the figure.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix with 2 columns is returned where the rows indicate different tumor ids and
the columns indicate the first and the last measured response value.
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Teemu D. Laajala <tlaajala@cc.hut.fi>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{xeno.norm.start}}
}
\examples{

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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ correlation }
\keyword{ end-point measurements }
