\name{xeno.draw.power}
\alias{xeno.draw.power}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for drawing power curves
}
\description{
Function for drawing power curves produced with the simulated dataset power 
analysis procedure. Multiple fixed effect terms can be visualized simultaneously.
}
\usage{
xeno.draw.power(powermats = list(), powerlabs = c(), minpw = 0, 
maxpw = 1, minN = 14, maxN = 50, maintitle = "Power curves", 
drawcrit = 0.8, legendposition = "right")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{powermats}{
List of power matrices produced by the xeno.sim.power-function.
}
  \item{powerlabs}{
Labels for the power matrices provided in the powermats-parameter.
}
  \item{minpw}{
Minimum y-axis power value to be drawn. Defaults to 0.
}
  \item{maxpw}{
Maximum y-axis power value to be drawn. Defaults to 1.
}
  \item{minN}{
Minimum x-axis tumor amount N value to be drawn. Defaults to 14.
}
  \item{maxN}{
Maximum x-axis tumor amount N value to be drawn. Defaults to 50.
}
  \item{maintitle}{
Text title for the figure.
}
  \item{drawcrit}{
If non-NULL, a dashed line is drawn to indicate a certain power requirement in the graph.
Defaults to 0.8.
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

# Categorizing fit
mcf_low_EM = xeno.EM(data=mcf_low, formula=Response ~ 1 + Treatment + 
Timepoint:Growth + Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))
mcf_low_fit = lmer(data=mcf_low_EM, Response ~ 1 + Treatment + Timepoint:Growth + 
Treatment:Timepoint:Growth + (1|Tumor_id) + (0+Timepoint|Tumor_id))

# Power simulation example with too few simulations to be actually reliable;
# Actual simulations would be suggested to be done with e.g. nsim=10000 and 
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
\keyword{ visualization }
\keyword{ power }
