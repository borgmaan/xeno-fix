\name{xeno.draw.cat}
\alias{xeno.draw.cat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for drawing found categories from a fitted model's frame
}
\description{
Visualization of found categories from the data using the defined model and 
the EM-algorithm. For binary categorization a two-color curve plot is done, 
and for continuous categorization a continuous color range is used instead.
For non-categorizing fit the plot is non-informative.
}
\usage{
xeno.draw.cat(fit, orig_data, responsename = "Response", 
treatmentname = "Treatment", tpname = "Timepoint", idname = "Tumor_id", 
maintitle = "Found categorization", legendposition = "topleft", 
ymax = NULL, ymin = NULL, scaled.growth = FALSE, 
discrete.colours = c("blue", "red", "green", "purple"), 
continuous.colours = rainbow(1000, start = 0.7, end = 0.95))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
Mixed-effects model object fit with lme4 (mer).
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
  \item{tpname}{
Column name with the time points in the data.frame. Defaults to "Timepoint".
}
  \item{idname}{
Column name with the individual tumor or animal labels in the data.frame. 
Defaults to "Tumor_id".
}
  \item{maintitle}{
Text title for the figure.
}
  \item{legendposition}{
Text string indicating where the legend should be placed; must be one of "bottomright",
"bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and
"center". NULL or FALSE value will omit legend.
}
  \item{ymax}{
A maximum cutoff point for y-axis not derived from the data.
}
  \item{ymin}{
A minimum cutoff point for y-axis not derived from the data.
}
  \item{scaled.growth}{
Should the growth-parameter should be scaled to be within [0,1] range for the 
purposes of highlighting different tumor growth values in the probablistic 
categorization procedure. Defaults to FALSE.
}
  \item{discrete.colours}{
Vector containing discrete colours for discrete categories.
}
  \item{continuous.colours}{
Vector containing continuous colour palettes for probabilistic categories.
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

xeno.draw.data(mcf_low, maintitle="MCF-7 LAR Low dosage")

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

par(mfrow=c(1,2))
# Visualizing data
xeno.draw.data(mcf_low)
# Visualizing categories
xeno.draw.cat(mcf_low_fit)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ visualization }
\keyword{ categories }
