\name{xeno.nlmer.EM}
\alias{xeno.nlmer.EM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Non-linear categorizing EM-algorithm
}
\description{
Implementation for the categorizing EM-algorithm for non-linear mixed-effects models. 
In comparison to the linear version, this function can be much more time consuming 
and convergence issues can occur due to problematic starting values or model 
formulation. 
}
\usage{
xeno.nlmer.EM(formula, data, Model, ModelGradient, start, 
max.iter = 100, loglh.difference = 0.01, responsename = "Response", 
treatmentname = "Treatment", tpname = "Timepoint", 
idname = "Tumor_id", verbose = FALSE, discriminate = TRUE, 
randomstart = FALSE, return.iterations = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Model formula used for fitting with lme4.
}
  \item{data}{
The original data.frame of the experiment in R long-data format. Categories will 
be added to the 
data.frame as an additional column.
}
  \item{Model}{
Model function, see the example.
}
  \item{ModelGradient}{
Model gradient, see the example.
}
  \item{start}{
Starting conditions as a vector for the fitting function.
}
  \item{max.iter}{
Maximum number of EM-iteration steps. Defaults to 100.
}
  \item{loglh.difference}{
Difference threshold between log-likelihoods of consecutive steps that is used 
as criteria for stopping the EM-algorithm. Defaults to 0.01.
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
  \item{verbose}{
Should output be printed to the user. Defaults to FALSE.
}
  \item{discriminate}{
Should the categories be discrete {0,1}. Defaults to TRUE. If false, categories 
will be probabilistic [0,1].
}
  \item{randomstart}{
Should random binary categories be assigned to the tumor units in the beginning, 
defaults to FALSE. Default start gives all tumors parameter value 1 in the beginning.
}
  \item{return.iterations}{
Should EM-iterations be returned by the function, defaults to FALSE. Used mainly 
for plotting.
}
  \item{\dots}{
Further parameters passed for the lme4 fitting function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Default function return is by default a data.frame containing the original data with
the additional column Growth indicating the identified categories. If return.iterations
parameter is TRUE, return will be a list of length 3:
\item{comp1 }{The final fit (mer-object)}
\item{comp2 }{The final data (data.frame)}
\item{comp3 }{The model fits at each step (list of mer-objects)}
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
# Allowing random parameters of the M-M with the addition of an intercept term 
# that ought to catch the starting criteria
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
\keyword{ EM }
