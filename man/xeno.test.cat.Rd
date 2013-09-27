\name{xeno.test.cat}
\alias{xeno.test.cat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for extracting 2x2 contingency table of Growth categories vs. Treatment 
groups for Fisher's Exact Test
}
\description{
Function that builds a 2x2 matrix with rows corresponding to the binary categories 
Growth and Non-growing and columns for the binary treatment groups assumed to be 
Control and Treatment.
}
\usage{
xeno.test.cat(x, treatmentname = "Treatment", tpname = "Timepoint", 
idname = "Tumor_id")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Mixed-effects model object fit with lme4 (mer) or a data.frame object.
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
For binary categories a 2x2 matrix with rows for Growing and Non-growing categories
and columns for the Control and Treatment groups.
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
This function will only work with the binary categorization.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fisher.test}}, 
\code{\link{xeno.test.marker.fit}}, 
\code{\link{xeno.test.marker.data}}, 
\code{\link{xeno.test.ran}}
}
\examples{
# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# 2x2 contingency table for Fisher's exact test
groups = xeno.test.cat(mcf_low_fit)
groups

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ categories }
\keyword{ fisher's exact test }
