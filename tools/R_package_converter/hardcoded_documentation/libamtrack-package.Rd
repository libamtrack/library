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
Version: \tab 0.4.1 (Blue Ocelot)\cr
Date: \tab 2010-10-31\cr
Depends: \tab R (>= 2.2.0)\cr
License: \tab GPL (version 3 or later)
}

\bold{FUNCTION INDEX:}

\bold{Efficiency / RBE routines:}

These functions compute the relative efficiency / RBE of a mixed particle field according to a specific amorphous track model flavour and phyiscs

\tabular{ll}{
\code{\link{AT.run.GSM.method}} \tab Grid-summation ('checkerboard') method \cr
\code{\link{AT.run.IGK.method}} \tab Ion-gamma-kill ('Katz') method \cr
\code{\link{AT.run.CPPSC.method}} \tab Compund-Poison process with successive convolutions (CPP-SC, 'SPIFF') method
}

\bold{Track-structure routines:}

These functions handle underlying physic used in amorphous track modeling:

\tabular{ll}{
\code{\link{AT.D.RDD.Gy}} \tab Dose distribution around a particle track \cr
\code{\link{AT.max.electron.ranges.m}} \tab Maximum electron range / track width
%\code{\link{AT.gamma.response}} \tab Gamma / X ray response of a system
}


\bold{Physics routines:}

These functions handle the physics of proton and ion beams needed in libamtrack.

Stopping-power routines:

\tabular{ll}{
\code{\link{AT.LET.keV.um}} \tab Unrestricted LET / Stopping power\cr
\code{\link{AT.LET.MeV.cm2.g}} \tab Unrestricted LET per mass / mass stopping power\cr
\code{\link{AT.E.MeV.from.LET}} \tab Inverse routine: energy from LET (cave!)\cr
\code{\link{AT.stopping.power.ratio}} \tab Computes stopping power ratios, also for mixed fields\cr

\code{\link{AT.CSDA.range.m}} \tab CSDA range in m\cr
\code{\link{AT.CSDA.range.g.cm2}} \tab CSDA range in mass thickness\cr
\code{\link{AT.E.MeV.from.CDSA.range}} \tab Inverse routine: energy from CSDA range\cr

\code{\link{AT.Bethe.Mass.Stopping.Power.MeV.cm2.g}} \tab Restricted (!) LET / mass stopping power from the classical Bethe formular
}

Mean LET / energy in mixed fields routines:

\tabular{ll}{
\code{\link{AT.fluence.weighted.LET.MeV.cm2.g}} \tab Computes fluence-weighted LET\cr
\code{\link{AT.dose.weighted.LET.MeV.cm2.g}} \tab Computes dose-weighted LET \cr

\code{\link{AT.fluence.weighted.E.MeV.u}} \tab Computes fluence-weighted mean energy\cr
\code{\link{AT.dose.weighted.E.MeV.u}} \tab Computes dose-weighted mean energy
}


Dose / fluence conversions:

\tabular{ll}{
\code{\link{AT.dose.Gy.from.fluence.cm2}} \tab Compute dose(s) for given fluence(s) and particle(s)\cr
\code{\link{AT.fluence.cm2.from.dose.Gy}} \tab Compute fluence(s) given dose(s) and particle(s)\cr
\code{\link{AT.total.D.Gy}} \tab Computes total dose for a mixed field\cr
\code{\link{AT.total.fluence.cm2}} \tab Computes total fluence for a mixed field
}


Beam related routines:

\tabular{ll}{
\code{\link{AT.beam.par.technical.to.physical}} \tab For double Gaussian beam, converts FWHM and particle number into fluence and sigma width\cr
\code{\link{AT.beam.par.physical.to.technical}} \tab Inverse, converts fluence and sigma width into FWHM and particle number
}


Misc physics routines:

\tabular{ll}{
\code{\link{AT.momentum.MeV.c.u.from.E.MeV.u}} \tab Momentum from kinetic energy \cr
\code{\link{AT.E.MeV.u.from.momentum.MeV.c.u}} \tab Kinetic energy from momentum \cr
\code{\link{AT.effective.charge.from.E.MeV.u}} \tab Effective charge of an ion depending on its kinetic energy \cr
\code{\link{AT.max.E.transfer.MeV}} \tab Max energy transfered from an ion to secondary electrons\cr
\code{\link{AT.mean.number.of.tracks.contrib}} \tab Mean number of tracks that desposite dose in a representative point\cr

\code{\link{AT.beta.from.E}} \tab Relativistic beta of an ion with \cr
\code{\link{AT.E.from.beta}} \tab Kinetic energy for given beta \cr
\code{\link{AT.gamma.from.E}} \tab Relativistic gamma of an ion
}


Other routines:

\tabular{ll}{
\code{\link{AT.particle.name.from.particle.no}} \tab Converts particle index numbers into particle names\cr
\code{\link{AT.particle.no.from.particle.name}} \tab Converts particle names into particle index numbers\cr
\code{\link{AT.material.name.from.material.no}} \tab Converts material index numbers into material names\cr
\code{\link{AT.material.no.from.material.name}} \tab Converts material names into material index numbers\cr
\code{\link{AT.get.material.data}} \tab Returns properties of pre-defined material
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

##############################################################################
############################ 1. LET ##########################################
cat("Compute the LET (in MeV*cm2/g) of a 270 MeV/u carbon ion in Aluminum:\n")

AT.LET.MeV.cm2.g( E.MeV.u     = 270,
			      particle.no = AT.particle.no.from.particle.name("12C"),
			      material.no = AT.material.no.from.material.name("Aluminum"))
				  
cat("... and in water:\n")

AT.LET.MeV.cm2.g( E.MeV.u     = 270,
			      particle.no = AT.particle.no.from.particle.name("12C"),
			      material.no = AT.material.no.from.material.name("Water, Liquid"))

########################################################################################
########################### 2. DOSE AROUND A TRACK #####################################
cat("Compare the Geiss parametrization for protons and Carbon at different energies:\n")

require(lattice)
df	<-	expand.grid(	E.MeV.u	       = 10^seq(0, 3, length.out = 4),       # from 1 to 1000 MeV/u in 4 steps
                        particle.no	   = c(1001,6012),                       # protons and carbons
                        r.m		       = 10^seq(-9, -2, length.out = 100),   # from 1 nm to 1 cm in 100 steps
                        material.no    = 2,                                  # Aluminium Oxide
                        rdd.model	   = 3,                                  # Geiss parametrization
                        rdd.parameter  = 5e-8,                               # Fixed core size of 50 nm
                        er.model	   = 4,                                  # Geiss track width parametrization
                        D.Gy		   = 0)                                  # For later use

ii                   <-  df$particle.no == 1001                              # Add particle names
df$particle.name     <-  "Carbon-12"
df$particle.name[ii] <-  "Protons"

for (i in 1:nrow(df)){                                                       # Loop through particles/energies
	df$D.Gy[i]	<-	AT.D.RDD.Gy( 	r.m		      = df$r.m[i],
                                    E.MeV.u	      = df$E.MeV.u[i],
                                    particle.no	  = df$particle.no[i],
                                    material.no	  = df$material.no[i],
                                    rdd.model	  = df$rdd.model[i],
                                    rdd.parameter = df$rdd.parameter[i],
                                    er.model	  = df$er.model[i])[[1]]
}

xyplot( log10(D.Gy) ~ log10(r.m)|particle.name,                              # Plot
		df,
		type      = 'l',
		groups    = E.MeV.u,
		auto.key  = TRUE)[c(2,1)]
		
#####################################################################
#################### 3. DETECTOR EFFICIENCY #########################
cat("Compute the relative efficiency of Alanine in 10 MeV protons\n")
cat("Waligorskis version of the Katz' model\n")

AT.run.IGK.method( particle.no                          = 1001,                        # namely protons with
                   E.MeV.u                              = 10,                          # 10 MeV/u
                   fluence.cm2.or.dose.Gy               = c(-1.0),                     # delivering 1 Gy
                   material.no                          = 5,                           # i.e. Alanine
                   rdd.model                            = 4,                           # Katz parametrization of radial dose distribution with simplified extended targets
                   rdd.parameter                        = c(5e-8,1e-10),               # with 50 nm target size and 1e-10 dose minimum
                   er.model                             = 2,                           # Butts&Katz parametrization of track radius
                   gamma.model                          = 2,                           # Use general target/hit model but here...
                   gamma.parameters                     = c(1,500,1,1,0),              # ...as exponential saturation with characteristic dose 500 Gy
                   saturation.cross.section.factor      = 1.4,                         # factor to take 'brush' around track into account
				   write.output                         = TRUE)                        # write a log file
}
