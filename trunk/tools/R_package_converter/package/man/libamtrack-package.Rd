\name{libamtrack-package}
\alias{libamtrack-package}
\alias{libamtrack}
\docType{package}
\title{
libamtrack package
}
\description{
This package is the R interface to the open-source, ANSI C library libamtrack (http://libamtrack.dkfz.org). libamtrack provides computational routines for the prediction of detector response and radiobiological efficiency in heavy charged particle beams. It is designed for research in proton and ion dosimetry and radiotherapy. libamtrack also provides many auxiliary physics routines for proton and ion beams.
}
\details{
\tabular{ll}{
Package: \tab libamtrack\cr
Version: \tab 0.3\cr
Date: \tab 2010-10-22\cr
Depends: \tab R (>= 2.2.0)\cr
License: \tab GPL (version 3 or later)\cr
Built: \tab R 2.10.1; i386-pc-mingw32; 2010-10-22 12:00:00 UTC; windows\cr
}

Index:
\preformatted{
AT.run.SPIFF.method     AT.run.SPIFF.method
}

}

\author{
Felix Klein <f.klein@dkfz.de>.

Maintainer: Steffen Greilich <s.greilich@dkfz.de>
}

\references{
TODO
}

\keyword{ package }

\seealso{
http://libamtrack.dkfz.org
}

\examples{

#######################################################################
############################ 1. LET ###################################
print("Compute the LET (in keV/um) of a 270 MeV/u carbon ion in Aluminum:")

AT.LET.keV.um( n           = 1,
               E.MeV.u     = 270,
			   particle.no = 6012,
			   material.no = 3)

}
