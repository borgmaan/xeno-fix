\name{xeno.test.ran}
\alias{xeno.test.ran}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for extracting random effect components of specific tumor fits and test 
them for correlation with 
defined marker values
}
\description{
Extracts such random effects such as random offset or random slope and tests for 
correlation with defined markers, e.g. KI-67 index in the example dataset. This can 
be done in respect to the whole dataset or within subsets, such as within specific 
treatment group(s) or growth categories, in search for biological subgroups.
}
\usage{
xeno.test.ran(fit, orig_data, rand_component = "Timepoint", 
other = "", treatmentname = "Treatment", idname = "Tumor_id", 
valueat = 1, cormethod = "pearson")
}
\arguments{
  \item{fit}{
Mixed-effects model object fit with lme4 (mer).
}
  \item{orig_data}{
data.frame containing the fitted data
}
  \item{rand_component}{
Random component that we're interested in, e.g. "(Intercept)" or "Timepoint". Defaults 
to "Timepoint".
}
  \item{other}{
Name of the column in orig_data: marker extracted from the data.frame that will be 
tested for correlation with the random effect.
}
  \item{treatmentname}{
Column name with treatment group. Defaults to "Treatment".
}
  \item{idname}{
Column name with separate individual (tumor) units. Defaults to "Tumor_id".
}
  \item{valueat}{
The time point value that contains the marker in the data.frame, typically the first 
row of the particular tumor should hold the value. Indexing starts at 1.
}
  \item{cormethod}{
Correlation method: defaults to "pearson", see also "spearman".
}
}
\details{
% 
}
\value{
Returns a list containing 3 elements of if not using binary categorization and 
5 if using binary categorization. First element holds two columns with whole 
data, second and third elements with subsets of data in respect to the first and 
second treatment group and the fourth and fifth elements hold subsets of data in 
respect to the first and second growth group. The elements contain 2 columns with 
first one with random effect value and second one with the marker value. The number 
of rows corresponds to the tumor units within the dataset.
}
\references{
% 
}
\author{
Teemu D Laajala <tlaajala@cc.hut.fi>
}
\note{
% 
}

\seealso{
\code{\link{xeno.test.cat}}, \code{\link{xeno.test.marker.fit}}, 
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

# Testing random intercept values to the corresponding KI-67 measurements
xeno.test.ran(mcf_low_fit, orig_data=mcf_low, rand_component="(Intercept)", 
other="KI67")
# Pearson correlation by default, Spearman rank correlation with cormethod:
xeno.test.ran(mcf_low_fit, orig_data=mcf_low, rand_component="(Intercept)", 
other="KI67", 
cormethod="spearman")

# Testing random slope values to the corresponding KI-67 measurements
xeno.test.ran(mcf_low_fit, orig_data=mcf_low, rand_component="Timepoint", 
other="KI67")
# Pearson correlation by default, Spearman rank correlation with cormethod:
xeno.test.ran(mcf_low_fit, orig_data=mcf_low, rand_component="Timepoint", 
other="KI67", 
cormethod="spearman")

}
\keyword{ random effects }
\keyword{ biomarker }
\keyword{ correlation }
