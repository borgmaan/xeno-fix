\name{xeno.test.marker.fit}
\alias{xeno.test.marker.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for extracting marker data from a fitted model
}
\description{
Function that will extract a field indicated in the fitted model's frame and compare 
it to a marker provided in the original data.frame used in the model fitting. The 
connected values are returned to the user for further use.
}
\usage{
xeno.test.marker.fit(fit, orig_data, componentname = "Growth", 
other = "", idname = "Tumor_id", tpname = "Timepoint", valueat = 1, 
verbose = TRUE, rm.na = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
Mixed-effects model object fit with lme4 (mer).
}
  \item{orig_data}{
The original data.frame used in the fitting of the mixed-effects model.
}
  \item{componentname}{
The column name of the first component to test, by default the category covariate "Growth".
}
  \item{other}{
The column name of the second component to test, must be included in the model fitting
thus present in the model frame.
}
  \item{idname}{
Column name with the individual tumor or animal labels in the data.frame. 
Defaults to "Tumor_id".
}
  \item{tpname}{
Column name with the time points in the data.frame. Defaults to "Timepoint".
}
  \item{valueat}{
The t:th time point in the long-format data.frame that holds the corresponding 
marker value, defaults to 1. This means that the marker value at the row with 
the first measured time point is used for representing the marker for the tumor unit.
}
  \item{verbose}{
Should output be printed to the user. Defaults to TRUE.
}
  \item{rm.na}{
Should NA-values be removed. Defaults to TRUE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a vector list with length equal to the different Growth category values found.
Each list element will include the marker values found for the corresponding
sorted Growth element (sorted from smallest to the largest). Thus for binary
categorization the returned list is of length 2, with first vector holding marker
values for Growth = 0 and second vector values for Growth = 1.
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
Differs from xeno.test.marker.data in that the first component is extracted from the model 
frame instead of data.frame holding the original data.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{xeno.test.marker.data}}, \code{\link{xeno.test.cat}}, 
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

# Extracting KI-67 markers according to the Growth subcategories
KI67_growth = xeno.test.marker.data(mcf_low_EM, componentname="Growth", other="KI67", 
verbose=FALSE)
KI67_growth

# Could be extracted also by taking Growth from the model frame and KI67 from the 
# data.frame
KI67_growth = xeno.test.marker.fit(mcf_low_fit, mcf_low_EM, componentname="Growth", 
other="KI67", verbose=FALSE)
KI67_growth

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ biomarker }

