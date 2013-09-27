\name{xeno.cat}
\alias{xeno.cat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for extracting the individual Growth-covariate values.
}
\description{
Function for extracting the Growth-covariate values for the identified categorization, 
where each id label holds a single categorization value. A vector is returned with the 
category values connected to the id labels.
}
\usage{
xeno.cat(fit, tpname = "Timepoint", idname = "Tumor_id")
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Vector with the category values. The vector elements have the id labels as a name vector.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
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

# Obtaining categories for the tumor labels:
xeno.cat(mcf_low_fit)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ categories }

