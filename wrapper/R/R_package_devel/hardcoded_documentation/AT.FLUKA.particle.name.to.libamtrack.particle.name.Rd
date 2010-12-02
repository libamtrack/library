\name{AT.FLUKA.particle.name.to.libamtrack.particle.name}
\alias{AT.FLUKA.particle.name.to.libamtrack.particle.name}
\title{AT.FLUKA.particle.name.to.libamtrack.particle.name}
\description{Coverts FLUKA style particle names to libamtrack convention.}
\usage{AT.FLUKA.particle.name.to.libamtrack.particle.name( FLUKA.particle.names)
}
\arguments{
  \item{FLUKA.particle.names}{ vector of strings containing the FLUKA-style particle names}
}
\value{
  Character string vector with libamtrack-style particle names.
}

\details{
FLUKA particle names relating to all particle with a specific Z, e.g. LI* or H*, will be translated to the 'even' isotope A = 2*Z. The only exception is H* which will be replaced by 1H. For example, LI* --> 6Li, C* --> 12C, BE* --> 8Be. This is a pure convention and does not necessarily represent the most prominent isotope. For most computations which only rely on Z it is, however, valid.
}

\examples{
AT.FLUKA.particle.name.to.libamtrack.particle.name( FLUKA.particle.names = c("H*", "B*", "B10", "C12", "BE7", "U238") )
}
