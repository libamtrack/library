\name{AT.CPPSC.alpha.and.beta}
\alias{AT.CPPSC.alpha.and.beta}
\title{AT.CPPSC.alpha.and.beta}
\description{Returns alpha and beta for an ion survival curve predicted by the CPPSC model. The arguments are directly passed on to the \code{\link{AT.run.CPPSC.method}} routine, some are provided default values here to simplify the use.

!! THIS ROUTINE IS LIMITED TO THE LIN-QUAD MODEL BUT SHOULD BE GENERALIZED FOR ANY GAMMA REPONSE, BOTH IN TERMS OF INPUT X RAY PARAMETERS *AND* OF RETURNED (FITTED) ION PARAMETERS!!

Cave: Even if input X ray response is LQ with cut dose, the resulting LQ will be purely linear-quadratic without high-dose linear characteristic.}
\usage{AT.CPPSC.alpha.and.beta(E.MeV.u, particle.no, fluence.cm2.or.dose.Gy, material.no, RDD.model, RDD.parameters, ER.model, gamma.model, gamma.parameters, N2 = 20, fluence.factor = 1.0, write.output = FALSE, shrink.tails = TRUE, shrink.tails.under = 1e-30, adjust.N2 = TRUE)
}
\arguments{
  \item{E.MeV.u}{ Energy}
  \item{particle.no}{ Particle index}
  \item{fluence.cm2.or.dose.Gy}{ Mean dose - will be varied for survival / response in a range 0.1x-10x.}
  \item{material.no}{ Material index}
  \item{RDD.model}{ RDD index}
  \item{RDD.parameters}{ RDD parameters}
  \item{ER.model}{ ER index}
  \item{gamma.model}{ gamma model index, must be LQ = 5}
  \item{gamma.parameters}{ gamma parameters, here alpha, beta, Dcut}
  \item{N2}{defaulted}
  \item{fluence.factor}{defaulted}
  \item{write.output}{defaulted}
  \item{shrink.tails}{defaulted}
  \item{shrink.tails.under}{defaulted}
  \item{adjust.N2}{defaulted}
}
\value{
List with follwing items
  \item{alpha}{ alpha for ion response curve}
  \item{beta}{ beta for ion response curve}
  \item{df}{data frame with complete ion response curve as computed to find alpha and beta including residuals fo that fit.}
}
\examples{
# None yet
}
