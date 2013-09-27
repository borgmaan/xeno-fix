\name{xeno.EM}
\alias{xeno.EM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear categorizing EM-algorithm
}
\description{
The main categorizing mixed-effects model fitting procedure. The EM-algorithm 
assumes that there is tumor response heterogeneity present both in both treatment 
groups (Control and actual Treatment). The algorithm will try find hidden subgroups 
of profiles by performing a likelihood ratio test between binary categories 
Growing (Growth = 1) and Non-growing (Growth = 0). The prediction in 
Expectation-step for the groups is done by using the corresponding fixed effect 
fits and in Maximization-step the model is fit using lme4-package with the new 
Growth-parameter vector. Iteration is continued until a pre-determined amount of 
iterative steps has been taken or the log-likelihood difference of consecutive fits 
reduces below a threshold.
}
\usage{
xeno.EM(formula, data, max.iter = 100, loglh.difference = 0.01, 
responsename = "Response", treatmentname = "Treatment", 
tpname = "Timepoint", idname = "Tumor_id", discriminate = TRUE, 
randomstart = FALSE, return.iterations = FALSE, id = NULL,
verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
Model formula used for fitting with lme4.
}
  \item{data}{
The original data.frame of the experiment in R long-data format. Categories 
will be added to 
the data.frame as an additional column.
}
  \item{max.iter}{
Maximum number of EM-iteration steps. Defaults to 100.
}
  \item{loglh.difference}{
Difference threshold between log-likelihoods of consecutive steps that is used 
as criteria for 
stopping the EM-algorithm. Defaults to 0.01.
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
  \item{discriminate}{
Should the categories be discrete {0,1}. Defaults to TRUE. If FALSE, categories 
will be probabilistic [0,1].
}
  \item{randomstart}{
Should random binary categories be assigned to the tumor units in the beginning, 
defaults to FALSE. Default start gives all tumors parameter value 1 in the beginning.
}
  \item{return.iterations}{
Should EM-iterations be returned by the function, defaults to FALSE. Used 
mainly for plotting.
}
  \item{verbose}{
Should output be printed to the user. Defaults to TRUE.
}
  \item{id}{
R session specific id number if progression of p-value estimation should be written to a separate text file.
}
  \item{\dots}{
Further parameters passed for the lme4 fitting function.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Default function return is by default a data.frame containing the original data with
the additional column Growth indicating the identified categories. If return.iterations
parameter is TRUE, return will be a list of length 3:
\item{comp1 }{The final fit (mer-object)}
\item{comp2 }{The final data (data.frame)}
\item{comp3 }{The model fits at each step (list of mer-objects)}
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
\code{\link{xeno.nlmer.EM}}
}
\examples{
# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
# Fitted object with lme4
mcf_low_fit
# The data.frame with found categories (repeated for each observation for a tumor)
xeno.summary(mcf_low_EM)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EM }
\keyword{ categories }
