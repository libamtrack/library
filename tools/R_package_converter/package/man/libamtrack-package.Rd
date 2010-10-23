\name{libamtrack-package}
\alias{libamtrack-package}
\alias{libamtrack}
\docType{package}
\title{
libamtrack package
}
\description{
This package is the R interface to the open-source, ANSI C library libamtrack. 

libamtrack provides computational routines for the prediction of detector response and radiobiological efficiency in heavy charged particle beams. It is designed for research in proton and ion dosimetry and radiotherapy. libamtrack also provides many auxiliary physics routines for proton and ion beams.

Please note that libamtrack is still under heavy development and so is the R interface. Function can be unstable especially when arguments are pushed out of their scope. If you experience any trouble your feed back is very appreciated: \email{s.greilich@dkfz.de}.
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

\emph{Index:}

\bold{Efficiency / RBE functions:}

These function compute the relative efficiency / RBE of a mixed particle field according to a specific amorphous track model flavour and phyiscs
\preformatted{
AT.run.GSM.method
AT.run.IGK.method
AT.run.SPIFF.method
}

}

\author{
Original C to R autoconversion developed by Felix Klein \email{fklein@embl.de}
Package maintainer: Steffen Greilich \email{s.greilich@dkfz.de}
}

\references{
\cite{Greilich, Grzanka, Bassler, Andersen and Jäkel, Amorphous track models: A numerical comparison study, doi:10.1016/j.radmeas.2010.05.039}
}

\keyword{ package }

\seealso{
\url{http://libamtrack.dkfz.org}
}

\examples{

#######################################################################
############################ 1. LET ###################################
print("Compute the LET (in keV/um) of a 270 MeV/u carbon ion in Aluminum:")

AT.LET.keV.um( n           = 1,
               E.MeV.u     = 270,
			   particle.no = 6012,
			   material.no = 3)
			   
#######################################################################
############################ 2. RESPONSE###############################
print("Compute the LET (in keV/um) of a 270 MeV/u carbon ion in Aluminum:")

AT.LET.keV.um( n           = 1,
               E.MeV.u     = 270,
			   particle.no = 6012,
			   material.no = 3)
}
