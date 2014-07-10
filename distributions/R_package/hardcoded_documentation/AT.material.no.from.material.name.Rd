\name{AT.material.no.from.material.name}
\alias{AT.material.no.from.material.name}
\title{AT.material.no.from.material.name}
\description{Returns material name(s) for given material index number(s)
}
\usage{AT.material.no.from.material.name(material.name)
}
\arguments{
  \item{material.name}{ material name(s) (see also \code{\link{material.name}}).}
}
\value{
  \item{material.no}{ material index number(s) (see also \code{\link{material.no}}).}
}

\examples{
# Get material index numbers for three materials of which one is not implemented
AT.material.no.from.material.name( material.name     = c("Water, Liquid", "PMMA", "Squirrel liver pudding"))
}


