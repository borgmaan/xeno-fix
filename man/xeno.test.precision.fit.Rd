\name{xeno.test.precision.fit}
\alias{xeno.test.precision.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for extracting precision for different time periods for vector K from a 
fitted lme4 object (mer)
}
\description{
Function for testing omitting of time points from the end of the experiment and 
observing how this affects the precision of a hypothesis of form H0: t(K) b = 0. 
The function will re-fit the provided model object with shorter than original time 
periods and return the desired precisions.
}
\usage{
xeno.test.precision.fit(fit, K = c(0, 0, 0, 1), 
responsename = "Response", tpname = "Timepoint", 
treatmentname = "Treatment", idname = "Tumor_id")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
Mixed-effects model object fit with lme4 (mer).
}
  \item{K}{
Vector K indicating the type of hypothesis to be testedin H0: t(K) b = 0
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
}
\details{
Returns list of vectors with amount of data points omitted, the non-centrality 
lambda parameter, critical F value, power according to the F statistic and 
precision for the provided K and fit.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Stroup WW.  Mixed model procedures to assess power, precision and sample 
size in the design of experiments. ASA Proceedings of the Biopharmaceutical 
Section. Alexandria, American Statistical Association; 1999. p. 15-24.
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
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + Timepoint:Growth 
+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
# Non-categorizing fit
mcf_low_fit2 = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint + 
Treatment:Timepoint + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Testing precision of combining offset and slope treatment effect terms for the 
# hypothesis
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
	fitlabels = c("MCF-7 LAR low dosage categorized", c("MCF-7 LAR low dosage no 
	categories")))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ precision }

