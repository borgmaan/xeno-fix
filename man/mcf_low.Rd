\name{mcf_low}
\alias{mcf_low}
\docType{data}
\title{
MCF-7 LAR lower dosage dataset
}
\description{
Example dataset in categorizing mixed-effects model R-package XenoCat 
}
\usage{data(mcf_low)}
\format{
  A data frame with 210 observations on the following 5 variables.
  \describe{
    \item{\code{Treatment}}{a numeric vector}
    \item{\code{Response}}{a numeric vector}
    \item{\code{Timepoint}}{a numeric vector}
    \item{\code{Tumor_id}}{a factor with levels \code{ER1} \code{ER10} \code{ER11} 
    \code{ER12} \code{ER13} \code{ER14} \code{ER15} \code{ER2} \code{ER20LAR1} 
    \code{ER20LAR10} \code{ER20LAR11} \code{ER20LAR12} \code{ER20LAR13} 
    \code{ER20LAR14} \code{ER20LAR15} \code{ER20LAR16} \code{ER20LAR17} 
    \code{ER20LAR18} \code{ER20LAR19} \code{ER20LAR2} \code{ER20LAR20} 
    \code{ER20LAR3} \code{ER20LAR4} \code{ER20LAR5} \code{ER20LAR6} \code{ER20LAR7} 
    \code{ER20LAR8} \code{ER20LAR9} \code{ER3} \code{ER4} \code{ER5} \code{ER6} 
    \code{ER7} \code{ER8} \code{ER9}}
    \item{\code{KI67}}{a numeric vector}
  }
}
\details{
Contains the tumor growth experiment of MCF-7 cell line xenograft treated with 
dietary Lariciresinol (lower dose). See the reference for more information. 
}
\source{
Saarinen NM, Warri A, Dings RP, Airio M, Smeds AI, Makela S. Dietary lariciresinol 
attenuates mammary tumor growth and reduces blood vessel density in human MCF-7 
breast cancer xenografts and carcinogen-induced mammary tumors in rats. Int J 
Cancer 2008;123:1196-204
}
\examples{
data(mcf_low)

head(mcf_low)
tail(mcf_low)

xeno.summary(mcf_low)
}
\keyword{datasets}
\keyword{MCF-7}
\keyword{xenograft}
