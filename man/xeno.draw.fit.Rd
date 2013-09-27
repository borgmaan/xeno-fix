\name{xeno.draw.fit}
\alias{xeno.draw.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for drawing model fit
}
\description{
Function for visualizing the whole model fit. Treatment groups (Control and 
actual Treatment group) will be visualized separately, also if categorizing 
modeling procedure has been used with binary categories, the corresponding 
Growth {0,1} values will be distinguished.
}
\usage{
xeno.draw.fit(fit, orig_data, responsename = "Response", 
treatmentname = "Treatment", tpname = "Timepoint", 
idname = "Tumor_id", maintitle = "Fixed + Random effect fit", 
legendposition = "topleft", drawgroups = TRUE, drawcat = FALSE,
colorgroups = FALSE)
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
Text string indicating where the legend should be placed; must be one of "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center". NULL or FALSE value will omit legend.
}
  \item{drawgroups}{
Should the original data curves be drawn on the background. Defaults to TRUE.
}
  \item{drawcat}{
Should the identified categories be indicated instead of the treatment groups (overrides 
drawgroups if set to TRUE). Defaults to FALSE.
}
  \item{colorgroups}{
If the original data curves are drawn on the background, should the treatment 
groups be indicated by colors. Defaults to TRUE.
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

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Visualizing data
xeno.draw.data(mcf_low)
# Non-categorizing fit
mcf_low_fit2 = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint + 
Treatment:Timepoint + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Visual comparison of the two fits
par(mfrow=c(2,2))
xeno.draw.fixed(mcf_low_fit, colorgroups=FALSE, maintitle="Categories fixed effects")
xeno.draw.fixed(mcf_low_fit2, colorgroups=FALSE, maintitle="No categories fixed effects")
xeno.draw.fit(mcf_low_fit, maintitle="Categories full fit")
xeno.draw.fit(mcf_low_fit2, maintitle="No categories full fit", legendposition=NULL)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ visualization }

