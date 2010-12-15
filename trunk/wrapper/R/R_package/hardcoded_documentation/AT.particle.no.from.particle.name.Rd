\name{AT.particle.name.from.particle.no}
\alias{AT.particle.name.from.particle.no}
\title{AT.particle.name.from.particle.no}
\description{Returns particle name(s) for given particle index number(s)
}
\usage{AT.particle.name.from.particle.no(particle.no)
}
\arguments{
  \item{particle.no}{ particle index number(s) (see also \code{\link{particle.no}}).}
}
\value{
  \item{particle.name}{ particle name(s) (see also \code{\link{particle.name}}).}
}

\examples{
# Get particle names for some nuclides
AT.particle.name.from.particle.no( particle.no     = c(1001,1002,3004,6012,8016,92238))
}


