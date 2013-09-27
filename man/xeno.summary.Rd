\name{xeno.summary}
\alias{xeno.summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary function for a data.frame of a tumor growth experiment
}
\description{
Function that explores a data.frame provided as function input. Reported values 
include e.g. missing observation percentages, tumor N amounts, measured time points 
and characteristics of the subgroups such as treatment groups.
}
\usage{
xeno.summary(orig_data, responsename = "Response", 
tpname = "Timepoint", treatmentname = "Treatment", 
idname = "Tumor_id")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{orig_data}{
The original data.frame used in the fitting of the mixed-effects model.
}
  \item{responsename}{
Column name with the response values in the data.frame. Defaults to "Response".
}
  \item{tpname}{
Column name with the time points in the data.frame. Defaults to "Timepoint".
}
  \item{treatmentname}{
Column name with the binary treatment group indicators in the data.frame. 
Defaults to "Treatment".
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

# Check top and bottom of the data.frame
head(mcf_low)
tail(mcf_low)

# MCF-7 LAR low dosage example dataset
xeno.summary(mcf_low)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ summary }

