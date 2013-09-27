\name{xeno.draw.data}
\alias{xeno.draw.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for drawing data from data.frame or fitted model's frame
}
\description{
Visualization of the data with separate tumor/animal units as their separate 
curves and Treatment groups (as in Control or actual Treatment) distinguished.
}
\usage{
xeno.draw.data(x, responsename = "Response", treatmentname = "Treatment", 
tpname = "Timepoint", idname = "Tumor_id", maintitle = "Data", 
legendposition = "topleft")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Mixed-effects model object fit with lme4 (mer) or a data.frame used for 
the model fitting.
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
Text string indicating where the legend should be placed; must be one of "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center". NULL or FALSE value will omit legend.
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

# Drawing the data
xeno.draw.data(mcf_low, maintitle="MCF-7 LAR Low dosage")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ visualization }

