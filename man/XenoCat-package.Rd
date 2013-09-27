\name{XenoCat-package}
\alias{XenoCat-package}
\alias{XenoCat}
\docType{package}
\title{
Categorizing mixed-effects models for the analysis of xenograft/tumor experiments.
}
\description{
The modeling framework offered here is intended for the analysis of xenograft and 
more general tumor growth experiments through categorizing mixed-effects models. 
The package lme4 by Bates et al. is essential for the mixed-effects procedures. The main
framework of the package is based on the EM-algorithm implanted within the mixed models 
to produce a biologically feasible categorization of tumor growth profiles
through binary categories by default. A probabilistic approach to the categories is also
offered, which can be used for fuzzy categorization for more than two categories of 
subgroups of tumors. Basic functionality for non-linear models is provided
(indicated by the prefix xeno.nlmer), but the use of linear models is emphasized.

}
\details{
\tabular{ll}{
Package: \tab XenoCat\cr
Type: \tab Package\cr
Version: \tab 1.0.3\cr
Date: \tab 2012-05-28\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
The main functionality is divided into subcategories:
All the XenoCat-package's functions start with the xeno-prefix:
\tabular{ll}{
prefix:\tab Functionality \cr
xeno.EM:\tab Expectation-Maximization implementation for categorizing linear mixed models \cr
xeno.cat:\tab Functions for handling the identified categories \cr
xeno.draw:\tab Visualization functions \cr
xeno.fixef:\tab Helper functions for inference based on the fixed effects of the model \cr
xeno.nlmer:\tab Functions for non-linear models fit with nlmer-function of lme4-package \cr
xeno.sim:\tab Simulation functions intended for power simulations \cr
xeno.test:\tab Helper functions for testing for biomarkers, precision, etc \cr
}
}
\author{
Teemu D. Laajala <tlaajala@cc.hut.fi>
}
\references{
Laajala TD, Corander J, Saarinen NM, Makela K, Savolainen S, Suominen MI, Alhoniemi E, 
Makela S, Poutanen M, Aittokallio T. Improved statistical modeling of tumor growth 
and treatment effect in pre-clinical animal studies with highly heterogeneous 
responses in vivo. Submitted.
}
\keyword{ package }
\keyword{ xenograft }
\keyword{ mixed-effects models }
\seealso{
lme4, lattice, languageR, LMERConvenienceFunctions, mixtools
}
\examples{
data(mcf_low)
library(lme4)

# MCF-7 low dosage example dataset
xeno.summary(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + Timepoint:Growth 
+ Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Visualizing the categorized fit
par(mfrow=c(2,3))
xeno.draw.data(mcf_low)
xeno.draw.fixed(mcf_low_fit,colorgroups=FALSE)
xeno.draw.fit(mcf_low_fit)
xeno.draw.res(mcf_low_fit)
xeno.draw.autocor(mcf_low_fit)
xeno.draw.propres(mcf_low_fit)

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
