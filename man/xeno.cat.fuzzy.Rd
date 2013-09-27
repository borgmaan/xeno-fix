\name{xeno.cat.fuzzy}
\alias{xeno.cat.fuzzy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mixture Gaussian model fit for distinguishing subgroups of the probabilistic 
categorization
}
\description{
The mixtools-package is used for mixture Gaussian model fit to the probabilistic
categories identified by the EM-algorithm (xeno.EM with parameter discriminate=FALSE).
A chosen number of Gaussian distributions are fit and the distribution with the highest
likelihood is assigned to be the corresponding discrete category. The number of 
identified categories in this fashion is not neccessarily equivalent to the number
of fitted Gaussian distributions.
}
\usage{
xeno.cat.fuzzy(fit, tpname = "Timepoint", idname = "Tumor_id", 
responsename = "Response", mu = c(0.1, 0.5, 0.9), 
sigma = c(0.1, 0.1, 0.1), draw = TRUE, legendposition = "top", 
rand=FALSE, randeff = "Timepoint", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
Mixed-effects model object fit with lme4 (mer).
}
  \item{tpname}{
Column name with the time points in the data.frame. Defaults to "Timepoint".
}
  \item{idname}{
Column name with the individual tumor or animal labels in the data.frame. 
Defaults to "Tumor_id".
}
  \item{responsename}{
Column name with the response values in the data.frame. Defaults to "Response".
}
  \item{mu}{
Initial mu-values for the Gaussian distributions.
}
  \item{sigma}{
Initial sigma-values for the Gaussian distributions.
}
  \item{draw}{
Should a figure be drawn with the probabilistic categorization values on the x-axis and y-axis 
showing the density of the estimated Gaussian distributions. Notice that the category values
are jittered slightly on the y-axis.
}
  \item{legendposition}{
Text string indicating where the legend should be placed; must be one of "bottomright",
"bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and
"center". NULL or FALSE value will omit legend.
}
 \item{rand}{
Should a random effect be used for identifying subcategories. If TRUE, the random effect 
defined by parameter randeff will be explored for subcategories. Defaults to FALSE.
}
 \item{randeff}{
Specified random effect to be explored for subcategories (if rand==TRUE), 
defaults to the random slope.
}
  \item{\dots}{
Additional parameters passed to the normalmixEM-function of the mixtools-package.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data.frame is returned where the Growth-column includes the discrete categories
according to the mixture Gaussian fit.
}
\references{
Benaglia T, Chauveau D, Hunter DR, Young DS. mixtools: An R Package for Analyzing Mixture Models. J Stat Soft 2009;32:1-29.
}
\author{
Teemu D. Laajala <tlaajala@cc.hut.fi>
}
\note{
A seed number ought to be used for reproducibility of the mixture Gaussian fit.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{xeno.cat}}
}
\examples{

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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ categories }

