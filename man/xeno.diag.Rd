\name{xeno.diag}
\alias{xeno.diag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diagnostics function for categories and fixed effect estimates in XenoCat-package
}
\description{
This function will perform diagnostics for detecting local or global 
optimum of the identified growth categories. The a prior solution 
(initial Growth-parameter set to unity) is compared to random start 
solutions and the user is informed whether the a prior solution 
resulted in global optimum. Additionally, the user will be given
automated feedback for the feasibility of the fixed effects 
coefficients, such as the base level b1-estimate (should be coherent
with the target size in experiment if such was used) and if the
growth inhibition ratio, |b4/b3|, has resulted in an unfeasible solution,
which hints of model assumption violations.
}
\usage{
xeno.diag(x, model = Response ~ 1 + Treatment + Timepoint:Growth 
+ Treatment:Timepoint:Growth + (1 | Tumor_id) + (0 + Timepoint | Tumor_id)
, starts = 1000, discriminate = TRUE, idname = "Tumor_id", 
tpname = "Timepoint", treatmentname = "Treatment", 
responsename = "Response", verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data.frame or a lme4-mer object (latter preferred).
}
  \item{model}{
Model formulation. Defaults to the categorizing model for an experiment
with a target size included in the experiment design.
}
  \item{starts}{
How many random starts should be initialized. Defaults to 1000.
}
  \item{discriminate}{
Should the categories be discrete {0,1}. Defaults to TRUE. If FALSE, categories 
will be probabilistic [0,1].
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
  \item{verbose}{
Should output be printed to the user. Defaults to FALSE.
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
\code{\link{xeno.EM}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(mcf_low)
library(lme4)


# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + Timepoint:Growth 
+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Running the diagnostics-function
xeno.diag(mcf_low_fit,
model = Response ~ 1 + Treatment + Timepoint:Growth 
+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
