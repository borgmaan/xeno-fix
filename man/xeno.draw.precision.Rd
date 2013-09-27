\name{xeno.draw.precision}
\alias{xeno.draw.precision}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for drawing precision curves
}
\description{
Function for drawing precision curves according to the pre-defined vector K and 
the model. The curves can be used for evaluation of most informative time periods, 
according to the model fits with omitting of time periods from the end of the 
longitudinal study or by fitting time point specific treatment terms.
}
\usage{
xeno.draw.precision(fits = list(), testedK = list(), Klabels = c(), 
fitlabels = c(), responsename = "Response", tpname = "Timepoint", 
treatmentname = "Treatment", idname = "Tumor_id", verbose = TRUE, 
legendposition = "topright")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fits}{
List containing mixed-effects model fits with lme4 (mer).
}
  \item{testedK}{
List of vectors containing different formulations of K for the hypothesis tests.
}
  \item{Klabels}{
Vector of strings for labeling the different K vectors (testedK-parameter).
}
  \item{fitlabels}{
Vector of strings for labeling the different fits (fits-parameter).
}
  \item{responsename}{
Column name with the response values in the data.frame. Defaults to "Response".
}
  \item{tpname}{
Column name with the time points in the data.frame. Defaults to "Timepoint".
}
  \item{treatmentname}{
Column name with the binary treatment group indicators in the data.frame. 
Defaults to "Treatment".
}
  \item{idname}{
Column name with the individual tumor or animal labels in the data.frame. 
Defaults to "Tumor_id".
}
  \item{verbose}{
Should output be printed to the user. Defaults to TRUE.
}
  \item{legendposition}{
Text string indicating where the legend should be placed; must be one of "bottomright",
"bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and 
"center". NULL or FALSE value will omit legend.
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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ visualization }
\keyword{ precision }
