\name{libamtrack-package}
\alias{libamtrack-package}
\alias{libamtrack}
\docType{package}
\title{
libamtrack package
}
\description{
This package is the R interface to the open-source, ANSI C library libamtrack. 

libamtrack provides computational routines for the prediction of detector response and relative biological efficiency in heavy charged particle beams. It is designed for research in proton and ion dosimetry and radiotherapy. libamtrack also provides many auxiliary physics routines for proton and ion beams.

Please note that libamtrack is still under heavy development and so is the R interface. Function can be unstable especially when arguments are pushed out of their scope. If you experience any trouble your feed back is very appreciated: \email{s.greilich@dkfz.de}.
}

\details{
\tabular{ll}{
Package: \tab libamtrack\cr
Version: \tab 0.5.1 (Blue Wombat)\cr
Date: \tab 2011-08-05\cr
Depends: \tab R (>= 2.12.0)\cr
Suggests: \tab lattice\cr
SystemRequirements: \tab gsl\cr
License: \tab GPL (version 3 or later)
}

\bold{FUNCTION INDEX:}

\bold{Efficiency / RBE routines:}

These functions compute the relative efficiency / RBE of a mixed particle field according to a specific amorphous track model flavour and phyiscs

\tabular{ll}{
\code{\link{AT.run.GSM.method}} \tab Grid-summation ('checkerboard') method \cr
\code{\link{AT.run.IGK.method}} \tab Ion-gamma-kill ('Katz') method \cr
\code{\link{AT.run.CPPSC.method}} \tab Compund-Poison process with successive convolutions (CPP-SC, 'SPIFF') method \cr
\code{\link{AT.CPPSC.alpha.and.beta}} \tab Alpha, beta for ion beams as predicted by CPP-SC method
}

\bold{Track-structure routines:}

These functions handle underlying physic used in amorphous track modeling:

\tabular{ll}{
\code{\link{AT.D.RDD.Gy}} \tab Dose distribution around a particle track \cr
\code{\link{AT.r.RDD.m}} \tab Inverse dose distribution around a particle track \cr
\code{\link{AT.max.electron.ranges.m}} \tab Maximum electron range / track width \cr
\code{\link{AT.gamma.response}} \tab Gamma / X ray response of a system \cr
\code{\link{AT.inverse.gamma.response}} \tab Inverse Gamma / X ray response of a system
}

\bold{SPC routines:}

These functions provide an interface to spectral depth data:

\tabular{ll}{
\code{\link{AT.SPC.read}} \tab Read in spc data from file \cr
\code{\link{AT.SPC.spectrum.at.depth.step}} \tab Extract spectrum at given depth step\cr
\code{\link{AT.SPC.spectrum.at.depth.g.cm2}} \tab Extract spectrum at given depth (will use clostest step)\cr
\code{\link{AT.SPC.tapply}} \tab Applies a function to all depths of spc data
}

\bold{FLUKA routines:}

These function provide an interface to FLUKA output files:

\tabular{ll}{
\code{\link{AT.FLUKA.read.USRBIN.mesh}} \tab Reads USRBIN output (Cartesian mesh, also for multiple runs).\cr
\code{\link{AT.FLUKA.read.USRBIN.regs}} \tab Reads USRBIN output (regions, also for multiple runs).\cr
\code{\link{AT.FLUKA.read.USRTRACK}} \tab Reads USRTRACK output (energy spectra, also for multiple runs).\cr
\code{\link{AT.FLUKA.particle.name.to.libamtrack.particle.name}} \tab Converts FLUKA particle names to libamtrack conventions.
}

\bold{Physics routines:}

These functions handle the physics of proton and ion beams needed in libamtrack.

Stopping-power routines:

\tabular{ll}{
\code{\link{AT.Stopping.Power.keV.um}} \tab Unrestricted LET / Stopping power\cr
\code{\link{AT.Stopping.Power.MeV.cm2.g}} \tab Unrestricted LET per mass / mass stopping power\cr
\code{\link{AT.stopping.power.ratio}} \tab Computes stopping power ratios, also for mixed fields\cr

\code{\link{AT.CSDA.range.m}} \tab CSDA range in m\cr
\code{\link{AT.CSDA.range.g.cm2}} \tab CSDA range in mass thickness\cr

\code{\link{AT.Stopping.Power.Mass.Bethe.MeV.cm2.g}} \tab Restricted (!) LET / mass stopping power from the classical Bethe formula\cr
\code{\link{AT.CSDA.range.Bethe.g.cm2}} \tab CSDA range from Bethe formula\cr
\code{\link{AT.CSDA.energy.after.slab.E.MeV.u}} \tab Ion energy after slab of material according Bethe formula\cr
\code{\link{AT.WEPL.Bethe}} \tab Water equivalent path length from Bethe range
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
\code{\link{AT.Rutherford.SDCS}} \tab Rutherford cross section of energy transferred to sec. electrons\cr
\code{\link{AT.beta.from.E}} \tab Relativistic beta of an ion with \cr
\code{\link{AT.E.from.beta}} \tab Kinetic energy for given beta \cr
\code{\link{AT.gamma.from.E}} \tab Relativistic gamma of an ion
}


Other routines:

\tabular{ll}{
\code{\link{AT.particle.name.from.particle.no}} \tab Converts particle index numbers into particle names\cr
\code{\link{AT.particle.no.from.particle.name}} \tab Converts particle names into particle index numbers\cr
\code{\link{AT.particle.no.from.Z.and.A}} \tab Returns particle index number for given mass and atomic number\cr
\code{\link{AT.A.from.particle.no}} \tab Returns mass number for given particle number\cr
\code{\link{AT.Z.from.particle.no}} \tab Returns atomic number for given particle number\cr
\code{\link{AT.nuclear.spin.from.particle.no}} \tab Returns nuclear spin for given particle number\cr
\code{\link{AT.electron.density.m3}} \tab Returns electron density from average Z and A\cr
\code{\link{AT.electron.density.m3.from.material.no}} \tab Returns electron density for given material\cr
\code{\link{AT.material.name.from.material.no}} \tab Converts material index numbers into material names\cr
\code{\link{AT.material.no.from.material.name}} \tab Converts material names into material index numbers\cr
\code{\link{AT.electron.density.m3.from.composition}} \tab Computes electron density from material composition\cr
\code{\link{AT.average.A.from.composition}} \tab Computes average mass number from material composition\cr
\code{\link{AT.average.Z.from.composition}} \tab Computes average atomic number from material composition\cr
\code{\link{AT.effective.Z.from.composition}} \tab Computes effective atomic number from material composition\cr
\code{\link{AT.I.eV.from.composition}} \tab Computes average I value from material composition\cr
\code{\link{AT.set.user.material}} \tab Sets properties of user defined material (CAVE: only valid until library is freed)\cr
\code{\link{AT.set.user.material.from.composition}} \tab Sets properties of user defined material from elemental composition (CAVE: only valid until library is freed)\cr
\code{\link{AT.get.materials.data}} \tab Returns properties of pre-defined materials
}


Helper routines shipped with libamtrack:

\tabular{ll}{
\code{\link{AT.add.leading.zeros}} \tab Adds leading zeros to a character string representing a number\cr
\code{\link{AT.add.trailing.zeros}} \tab Adds trailing zeros to a character string representing a number
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
cat("Compute the LET (in MeV*cm2/g) of a 270 MeV/u carbon ion in Aluminum\n")
cat("using the PSTAR stopping power data:\n")

AT.Stopping.Power.MeV.cm2.g( E.MeV.u     = 270,
                  particle.no = AT.particle.no.from.particle.name("12C"),
                  material.no = AT.material.no.from.material.name("Aluminum"),
                  stopping.power.source.no             = 0)
                  
cat("... and in water:\n")

AT.Stopping.Power.MeV.cm2.g( E.MeV.u     = 270,
                  particle.no = AT.particle.no.from.particle.name("12C"),
                  material.no = AT.material.no.from.material.name("Water, Liquid"),
                  stopping.power.source.no             = 0)

########################################################################################
########################### 2. DOSE AROUND A TRACK #####################################
cat("Compare the Geiss parametrization for protons and Carbon at different energies:\n")

require(lattice)
df    <-    expand.grid(    E.MeV.u           = 10^seq(0, 3, length.out = 4),       # from 1 to 1000 MeV/u in 4 steps
                            particle.no       = c(1001,6012),                       # protons and carbons
                            r.m               = 10^seq(-9, -2, length.out = 100),   # from 1 nm to 1 cm in 100 steps
                            material.no       = 2,                                  # Aluminium Oxide
                            rdd.model         = 3,                                  # Geiss parametrization
                            rdd.parameter     = 5e-8,                               # Fixed core size of 50 nm
                            er.model          = 4,                                  # Geiss track width parametrization
                            D.Gy              = 0)                                  # For later use

ii                   <-  df$particle.no == 1001                                     # Add particle names
df$particle.name     <-  "Carbon-12"
df$particle.name[ii] <-  "Protons"

for (i in 1:nrow(df)){                                                              # Loop through particles/energies
    df$D.Gy[i]    <-    AT.D.RDD.Gy(     r.m              = df$r.m[i],
                                         E.MeV.u          = df$E.MeV.u[i],
                                         particle.no      = df$particle.no[i],
                                         material.no      = df$material.no[i],
                                         rdd.model        = df$rdd.model[i],
                                         rdd.parameter    = df$rdd.parameter[i],
                                         er.model         = df$er.model[i],
                                         stopping.power.source.no             = 0)[[1]]                           # use PSTAR data
}

xyplot( log10(D.Gy) ~ log10(r.m)|particle.name,                                     # Plot
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
                   write.output                         = TRUE,                        # write a log file
                   stopping.power.source.no             = 0)                           # use PSTAR data
}
