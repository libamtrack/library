\name{AT.particle.no.from.particle.name.single}
\alias{AT.particle.no.from.particle.name.single}
\title{AT.particle.no.from.particle.name.single}
\description{Returns particle number for given particle index
}
\usage{AT.particle.no.from.particle.name.single(particle.name)
}
\arguments{
  \item{particle.name}{ particle name (see also \code{\link{particle.no}}).}
}
\value{
  \item{particle.no}{ particle index number (see also \code{\link{particle.name}}).}
}

\examples{
# Get particle names for some nuclides
AT.particle.no.from.particle.name.single( particle.name = "12C")
}


