\name{xeno.formula}
\alias{xeno.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Easy XenoCat model formulation with several handy options
}
\description{
This function will produce an R formula-object according to the provided 
parameters.
}
\usage{
xeno.formula(target = TRUE, cat = TRUE, responsename = "Response", 
treatmentname = "Treatment", tpname = "Timepoint", idname = "Tumor_id", 
terms = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{target}{
Should the fixed effects terms be included that are most suitable for
experiments with a designated target size. These terms are b1 (base level)
and b2 (offset), that require possibility for y-axis shift in the response
profiles. Defaults to TRUE.
}
  \item{cat}{
Should the categorizing Growth-covariate be included in the model. If FALSE,
the model will reduce to the conventional, non-categorizing model. Defaults to
TRUE.
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
  \item{terms}{
Should several fixed and random effects be included. The terms in the
vector with 6 elements are: b1, b2, b3, b4 (fixed effects) and u1, u2 
(random effects).
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

# Default model formulation returned
xeno.formula()
# Omit b1,b2 - no target size included
xeno.formula(target=FALSE)
# Default non-categorizing model
xeno.formula(cat=FALSE)
# Default non-categorizing model for an experiment without target size
xeno.formula(cat=FALSE, target=FALSE)


# MCF-7 LAR low dosage example dataset
data(mcf_low)
# Model
frml = xeno.formula()

mcf_low_EM = xeno.EM(data=mcf_low, formula=frml)
mcf_low_fit = lmer(data=mcf_low_EM, frml)
# Fitted object with lme4
mcf_low_fit



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
