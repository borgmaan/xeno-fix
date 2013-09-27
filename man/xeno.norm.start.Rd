\name{xeno.norm.start}
\alias{xeno.norm.start}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for response transformations and normalization
}
\description{
A simple function for transforming response values in a data.frame. The default
approach includes logarithmic transform with the possilibity to normalize the 
response values by diving the measurements by the initial tumor size.
}
\usage{
xeno.norm.start(data, newfield = "LogResponse", 
responsename = "Response", tpname = "Timepoint", 
idname = "Tumor_id", log = TRUE, start = TRUE, 
type = "divide")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The data.frame in which the values are to be transformed.
}
  \item{newfield}{
Name for the new column including the transformed response values. Defaults to "LogResponse".
}
  \item{responsename}{
Column name with the response values in the data.frame. Defaults to "Response".
}
  \item{tpname}{
Column name with the time points in the data.frame. Defaults to "Timepoint".
}
  \item{idname}{
Column name with the individual tumor or animal labels in the data.frame. 
Defaults to "Tumor_id".
}
  \item{log}{
Should logarithmic transformation after normalization be done. Defaults to TRUE.
}
  \item{start}{
Should start-point normalization be done to the response values. Defaults to TRUE.
}
  \item{type}{
Type of start-point normalization to do. Defaults to dividing the response values by the
first measurement of the tumor response.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The original data.frame is returned with the addition of the new column with the 
transformed response values.
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
\code{\link{xeno.test.endcor}}
}
\examples{

data(mcf_low)
library(lme4)

# MCF-7 low dosage example dataset
xeno.summary(mcf_low)

# Log-transform of response values, wich are normalized by the first measured value
mcf_low_new = xeno.norm.start(mcf_low)

# Dropping out the first time point because it is no longer informative
mcf_low_new = subset(mcf_low_new, mcf_low_new["Timepoint"]>0)
mcf_low_new["Timepoint"] = mcf_low_new["Timepoint"] - 1

# Categorizing fit for the log-transformed, normalized response
mcf_low_EM_lognorm = xeno.EM(data=mcf_low_new, formula=LogResponse ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id),
responsename="LogResponse")

mcf_low_fit_lognorm = lmer(data=mcf_low_EM_lognorm, LogResponse ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

#  Drawing the fit to the transformed responses
xeno.draw.fit(mcf_low_fit_lognorm, responsename="LogResponse")

# Normalization did not change the trend of the identified categories;
# Less tumors identified as non-growing in the control group
xeno.test.cat(mcf_low_fit_lognorm)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ normalization }
\keyword{ transform }
