\name{AT.material.name.from.material.no}
\alias{AT.material.name.from.material.no}
\title{AT.material.name.from.material.no}
\description{Returns material name(s) for given material index number(s)
}
\usage{AT.material.name.from.material.no(material.no)
}
\arguments{
  \item{material.no}{ material index number(s) (see also \code{\link{material.no}}).}
}
\value{
  \item{material.name}{ material name(s) (see also \code{\link{material.name}}).}
}

\examples{
# Get material names for all predefined materials
AT.material.name.from.material.no( material.no = 1:25)
}


