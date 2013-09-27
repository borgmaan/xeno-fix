\name{xeno.fixef.perm}
\alias{xeno.fixef.perm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Permutation of the treatment labels for statistical significance
}
\description{
This function can be used for validating fixed effect significance with treatment specific
effects. The permutation is done by randomly distributing the treatment-labels in the
data over the individual tumor units and re-fitting the model. The p-value of a fixed 
effect term is the percentage of its permutated t statistics that have absolute value higher 
than that of the original fitted model. For categorizing mixed-effects model, the identified 
categories are permutated while conserving the connection to their original treatment labels.
}
\usage{
xeno.fixef.perm(fit, perms = 10000, idname = "Tumor_id", 
tpname = "Timepoint", treatmentname = "Treatment", 
responsename = "Response", id = NULL, verbose = FALSE,
time = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
Mixed-effects model object fit with lme4 (mer).
}
  \item{perms}{
Number of label permutations.
}
  \item{idname}{
Column name with the individual tumor or animal labels in the data.frame. 
Defaults to "Tumor_id".
}
  \item{tpname}{
Column name with the time points in the data.frame. Defaults to "Timepoint".
}
  \item{treatmentname}{
Column name with the binary treatment group indicators in the data.frame. 
Defaults to "Treatment".
}
  \item{responsename}{
Column name with the response values in the data.frame. Defaults to "Response".
}
  \item{id}{
R session specific id number if progression of p-value estimation should be written to a separate text file.
}
  \item{verbose}{
Should the function print a line every 100 permutations. Defaults to FALSE.
}
  \item{time}{
Should a brief test permutation be run instead of running the actual permutation, so that a time estimate can be given to the user. Defaults to FALSE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector for the fixed effects is returned with the corresponding p-values.
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
Teemu D. Laajala <tlaajala@cc.hut.fi>
}
\note{
The significance of the treatment and growth category specific terms can be 
validated through the permutation. However, this method will not provide feasible
p-values for the significance of model terms that do not depend on the treatment 
and/or identified categories, such as the intercept-term.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{xeno.fixef.pvals}}
}
\examples{

data(mcf_low)
library(lme4)

# MCF-7 low dosage example dataset
xeno.summary(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + Timepoint:Growth 
	+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
	Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# MCMC p-values
xeno.fixef.pvals(mcf_low_fit)
# p-values according to permutation
xeno.fixef.perm(mcf_low_fit)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ fixed effects }
\keyword{ significance }
\keyword{ model validation }