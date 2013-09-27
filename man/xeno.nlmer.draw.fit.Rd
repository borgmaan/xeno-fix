\name{xeno.nlmer.draw.fit}
\alias{xeno.nlmer.draw.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Non-linear (nlmer) function for drawing full fit
}
\description{
Function implementation for visualizing whole model fit for the non-linear mixed-effects models fit 
with the lme4-function nlmer.
}
\usage{
xeno.nlmer.draw.fit(fit, orig_data, responsename = "Response", 
treatmentname = "Treatment", idname = "Tumor_id", tpname = "Timepoint", 
draw_orig = TRUE, maintitle = "Full fit", ymax = NULL, 
per_individual = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
Non-linear mixed-effects model object fit with lme4 (mer).
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
  \item{idname}{
Column name with the individual tumor or animal labels in the data.frame. Defaults to "Tumor_id".
}
  \item{tpname}{
Column name with the time points in the data.frame. Defaults to "Timepoint".
}
  \item{draw_orig}{
Should the original data curves be drawn on the background of the figure, defaults to TRUE.
}
  \item{maintitle}{
Text title for the figure.
}
  \item{ymax}{
A cutoff point for y-axis not derived from the data.
}
  \item{per_individual}{
Should individual plots be drawn. This will open a new figure panel for each of the individuals,
so e.g. a PostScript file with many pages ought to be used to capture the individual observed
values and fits.
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

# Concept of fitting a non-linear model with the presented EM-algorithm procedure
# Defining a model with Michaelis-Menten kinetics. 
# Allowing random parameters of the M-M with the addition of an intercept 
#term that ought to catch the starting criteria
Model = function(Intercept, Offset, Treatment, Growth, VM, K, Timepoint) { 
	Intercept + Offset*Treatment + (Growth*VM*Timepoint/(K+Timepoint))
}
ModelGradient = deriv(body(Model)[[2]], namevec = c("Intercept", "Offset", "VM","K"), 
function.arg=Model)
starting_conditions = c(Intercept=20, Offset=-1, VM=100, K = 1)



mcf_nlmer_EM = xeno.nlmer.EM(
	formula = Response ~ ModelGradient(Intercept, Offset, Treatment, Growth, 
	VM, K, Timepoint) ~ (VM|Tumor_id) + (K|Tumor_id) + (Intercept|Tumor_id),
	data=mcf_low, 
	Model=Model, 
	ModelGradient=ModelGradient,
	start=starting_conditions, 
	verbose=TRUE,
	max.iter=10)

mcf_nlmer_fit = nlmer(Response ~ ModelGradient(Intercept, Offset, Treatment, Growth, 
VM, K, Timepoint) ~ (VM|Tumor_id) + (K|Tumor_id) + (Intercept|Tumor_id), data = mcf_nlmer_EM, 
start=starting_conditions)

par(mfrow=c(2,2))
xeno.draw.data(mcf_low)
xeno.nlmer.draw.fixed(mcf_nlmer_fit, mcf_low, Model=Model)
xeno.nlmer.draw.fit(mcf_nlmer_fit, mcf_low)
xeno.draw.res(mcf_nlmer_fit)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ non-linear }
\keyword{ visualization }

