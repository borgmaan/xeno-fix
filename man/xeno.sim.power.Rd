\name{xeno.sim.power}
\alias{xeno.sim.power}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for power simulation for chosen tumor N and time scale
}
\description{
Function that wraps the simulated tumor profile procedure for different tumor N 
amounts and allows testing power for shorter or longer time periods than in the 
original experiment. Depending on the artificial dataset and MCMC simulation 
amounts the function call be rather time consuming, but the procedure can be 
split into smaller parts by running separate N vectors or by extending a parallel 
simulation/significance implementation.
}
\usage{
xeno.sim.power(fit, file = "POW_ANALYSIS.txt", nsim = 10000, 
Ns = c(10, 15, 20), pvalMCMC = TRUE, MCMCsim = 10000, alpha = 0.05, 
responsename = "Response", treatmentname = "Treatment", 
tpname = "Timepoint", idname = "Tumor_id", 
T = max(fit@frame[[tpname]]))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
Mixed-effects model object fit with lme4 (mer).
}
  \item{file}{
Filename to write temporary output to. Defaults to "POW_ANALYSIS.txt"
}
  \item{nsim}{
Number of desired amount of simulations of data, defaults to 10000.
}
  \item{Ns}{
Vector containing the N values to be tested, holds values 10, 15 and 20 by default.
}
  \item{pvalMCMC}{
Whether MCMC-samples should be used for computing statistical significance of 
fixed effect terms. If FALSE, LRT is used instead. Defaults to TRUE.
}
  \item{MCMCsim}{
Amount of MCMC simulations for estimating statistical significance, defaults to 10000.
}
  \item{alpha}{
Desired statistical significance threshold, defaults to 0.05.
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
  \item{T}{
Last measured time point. Defaults to the last original measured time point.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returned value is a matrix with rows corresponding to the indicated tumor N amounts to
test. Columns are the fixed effects of the fitted model.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Gelman A, Hill J. Multilevel power calculation using fake-data simulation. 
In: Alvarez RM, Beck NL, Wu LL (editors). Data analysis using regression 
and multilevel/hierarchical models, Cambridge University Press; 2007. p. 449-54.
}
\author{
Teemu D Laajala <tlaajala@cc.hut.fi>
}
\note{
The approach for power analysis through simulated dataset was adopted from Gelman & Hill.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{xeno.sim.tumor}}
}
\examples{
# MCF-7 LAR low dosage example dataset
data(mcf_low)

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Power simulation example with too few simulations to be actually reliable;
# Actual simulations would be suggested to be done with nsim=10000 and 
# MCMCsim=10000, more if able
power_test = xeno.sim.power(mcf_low_fit, nsim=100, MCMCsim=100, Ns=c(10,25,50), T=5)
# Testing a longer treatment period than in the original dataset
power_test_longer = xeno.sim.power(mcf_low_fit, nsim=100, MCMCsim=100, 
Ns=c(10,25,50), T=10)
# Drawing power curve
xeno.draw.power(
	powermats=list(power_test, power_test_longer), 
	powerlabs=c("Example power (too low sims)", "Longer tp example (too low sims)"),
	minN = 5,
	maxN = 50,
	legendposition="right")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ simulation }
\keyword{ power }
