\name{xeno.draw.EM}
\alias{xeno.draw.EM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for drawing categorizing EM-algorithm iterations
}
\description{
Visualization of consecutive and final EM-algorithm step taken in the 
estimation of the binary categories. The fixed effects fit at each iterative 
step is visualized, as it is used for prediction in the Expectation-step of 
the EM-algorithm. The found categories at each step are distinguished.
}
\usage{
xeno.draw.EM(formula, data, mfrow = NULL, responsename = "Response", 
tpname = "Timepoint", idname = "Tumor_id", treatmentname = "Treatment", 
legendposition = "topleft")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Model formula used for fitting with lme4.
}
  \item{data}{
The original data.frame of the experiment in R long-data format. Categories 
will be added to the data.frame as an additional column.
}
  \item{mfrow}{
Vector for paneling in par(), first element corresponding to rows and second to columns.
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
  \item{treatmentname}{
Column name with the binary treatment group indicators in the data.frame. 
Defaults to "Treatment".
}
  \item{legendposition}{
Text string indicating where the legend should be placed; must be one of "bottomright",
"bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and
"center". NULL or FALSE value will omit legend.
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

# Drawing EM-iterations
xeno.draw.EM(
	Response ~ 1 + Treatment + Timepoint:Growth + Treatment:Timepoint:Growth + 
	(1|Tumor_id) + (0+Timepoint|Tumor_id),
	data=mcf_low,
	mfrow=c(2,2))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ visualization }
\keyword{ EM }
