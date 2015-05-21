% Documentation for 'number.of.field.components' variable
\name{number.of.field.components}

\alias{number.of.field.components}
\alias{Handling of mixed particle fields}

\title{Handling of mixed particle fields}

\description{Many functions in libamtrack offer the possibility to consider mixed particle fields. These fields are composed from a finite number of specific nuclids defined by \code{\link{particle.no}}. All particles of an individual component are supposed to have exactly the same energy (given by \code{\link{particle.no}}), without spread. Thus an energy distribution has to be describe in discrete energy steps.

Below is a compilation of the arguments a mixed field is described by. Depending on the function all or a subset are necessary. If arrays of length > \code{number.of.field.components} are entered, the additional data will be ignored. If array hold too few data, the computation will not be executed.}

\arguments{
\item{number.of.field.components}{A variable describing how many components are used for the field. This also dictates the array lengths of all following variables.}
\item{E.MeV.u}{Particle energies (array of size \code{number.of.field.components}, see also \code{\link{E.MeV.u}}}
\item{particle.no}{Particle type indices (array of size \code{number.of.field.components}), see also \code{\link{particle.no}}}
\item{fluence.cm2.or.dose.Gy}{Fluences (if positive), or doses (if negative, array of size \code{number.of.field.components}), see also \code{\link{fluence.cm2.or.dose.Gy}}}
}


