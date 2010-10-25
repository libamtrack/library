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
Date: \tab 2010-10-25\cr
Depends: \tab R (>= 2.2.0)\cr
License: \tab GPL (version 3 or later)\cr
}

\bold{FUNCTION INDEX:}

\bold{Efficiency / RBE routines:}

These functions compute the relative efficiency / RBE of a mixed particle field according to a specific amorphous track model flavour and phyiscs

\tabular{ll}{
\code{\link{AT.run.GSM.method}} \tab Grid-summation ('checkerboard') method \cr
\code{\link{AT.run.IGK.method}} \tab Ion-gamma-kill ('Katz') method \cr
\code{\link{AT.run.SPIFF.method}} \tab Compund-Poison process with successive convolutions (CPP-SC, 'SPIFF') method \cr
}

\bold{Physics routines:}

These functions handle the physics of proton and ion beams needed in libamtrack:

\tabular{ll}{
\code{\link{AT.beta.from.E}} \tab Relativistic beta of an ion \cr
\code{\link{AT.LET.keV.um}} \tab Unrestricted LET / Stopping power \cr
}

}

\author{
C2R autoconversion developed by: Felix Klein \email{fklein@embl.de}

Package maintainer: Steffen Greilich \email{s.greilich@dkfz.de}
}

\references{
\cite{Greilich, Grzanka, Bassler, Andersen and Jakel, Amorphous track models: A numerical comparison study, doi:10.1016/j.radmeas.2010.05.039}
}

\keyword{ package }

\seealso{
\url{http://libamtrack.dkfz.org}
}

\examples{

###########################################################################
############################ 1. LET #######################################
print("Compute the LET (in keV/um) of a 270 MeV/u carbon ion in Aluminum:")

AT.LET.keV.um( n           = 1,
               E.MeV.u     = 270,
			   particle.no = 6012,
			   material.no = 3)
			   
#####################################################################
#################### 2. DETECTOR EFFICIENCY #########################
print("Compute the relative efficiency of Alanine in 10 MeV protons")
print("Waligorskis version of the Katz' model")

AT.run.IGK.method( number.of.field.components           = 1,                           # one component only
                   particle.no                          = 1001,                        # namely protons with
                   E.MeV.u                              = 10,                          # 10 MeV/u
                   fluence.cm2.or.dose.Gy               = c(-1.0),                     # delivering 1 Gy
                   material.no                          = 5,                           # i.e. Alanine
                   rdd.model                            = 4,                           # Katz parametrization of radial dose distribution with simplified extended targets
                   rdd.parameter                        = c(5e-8,1e-10)                # with 50 nm target size and 1e-10 dose minimum
                   er.model                             = 2,                           # Butts&Katz parametrization of track radius
                   gamma.model                          = 2,                           # General hit/target X ray response, but
                   gamma.parameters                     = c(1,500,1,1,0),              # as simple single exponential saturation (one hit, one target), saturation dose 500 Gy
                   saturation.cross.section.factor      = 1.4,                         # factor to take 'brush' around track into account
				   write.output                         = T) 
}
