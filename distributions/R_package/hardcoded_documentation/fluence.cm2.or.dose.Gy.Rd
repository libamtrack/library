% Documentation for 'fluence.cm2.or.dose.Gy' variable
\name{fluence.cm2.or.dose.Gy}

\alias{fluence.cm2.or.dose.Gy}
\alias{Fluence or dose argument}

\title{Fluence and dose as equivalent input parameters}

\description{Many functions in libamtrack can take either the fluence or the dose as an argument. For a specific nuclid and energy they can be converted into each other (see also \code{AT.dose.Gy()} and \code{AT.fluence.cm2}). As both can only be positive, the argument \code{fluence.cm2.or.dose.Gy} will be interpreted as fluence [1/cm2] if positive and as dose [Gy] if given as negative value. 

\itemize{
\item{\code{fluence.cm2.or.dose.Gy <- 1e7}}{ means fluence = 1e7 / cm2}
\item{\code{fluence.cm2.or.dose.Gy <- -10}}{ means dose = 10 Gy}
}
}
