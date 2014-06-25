\name{AT.SPC.read}
\alias{AT.SPC.read}
\title{AT.SPC.read}
\description{Read a spc-formatted data file with energy-fluence in n depth steps. For original TRiP format definition by M. Kraemer, please see \url{http://bio.gsi.de/DOCS/TRiP98BEAM/DOCS/trip98fmtspc.html} and Kraemer and Scholz, Treatment planning for heavy-ion radiotherapy: calculation and optimization of biologically effective dose, Phys. Med. Biol. 45 (2000) 3319-3330. Please note that the user has to take care of picking the spc-file for the projectile and target material desired.

Presently, two versions of the reader exists: a slower but stable R version ('vanilla') and a faster but still buggy C version.

IMPORTANT: SPC files report DIFFERENTIAL fluences dN/dE per primary particle, i.e. normalized by the energy bin. This reader converts them to ABSOLUTE fluences by multiplying dE. This facilitates e.g. summation to get the total fluence etc. but produces funny results when ABSOLUTE fluences are plotted. 

Also, SPCs only contain the left limit of the energy bins. This reader also reports the energy bin midpoint which is more reasonable e.g. for the computation of mean dE/dx etc. 

Please note that the "compress" statement has been depreciated due to preparation for C translation of spc handling.
}
\usage{AT.SPC.read( file.name, flavour, endian, mean, header.only = FALSE)
}
\arguments{
  \item{file.name}{ name of spc-file to be read (with extension).}
  \item{flavour}{ 'vanilla' (default) to use stable R version.}
  \item{endian}{ byte-order of the spc file to read, "big" (AIX) or "little" (Linux, VMS, etc., default). Only necessary if flavour = 'vanilla'}
  \item{mean}{ method for computing bin midpoints, "geometric" or "arithmetic" (default).}
  \item{header.only}{ if true, only information on spc file but no data will be read.}
}
\value{
  A list with (1) the 'spc' data frame with the following columns
  \item{depth.step}{ Index number of depth step (one-based)}
  \item{depth.g.cm2}{ Depth in g/cm2}
  \item{particle.no}{ Particle index number (see also(\code{\link{particle.no}}).}
  \item{E.low.MeV.u}{ Energy bin lower limit in MeV/u (see also(\code{\link{E.MeV.u}}).}
  \item{E.mid.MeV.u}{ Energy bin mid in MeV/u (see also(\code{\link{E.MeV.u}}).}
  \item{E.high.MeV.u}{ Energy bin upper limit in MeV/u (see also(\code{\link{E.MeV.u}}).}
  \item{dE.MeV.u}{ Energy bin width in MeV/u}
  \item{dN.dE.per.MeV.u.per.primary}{ Fluence differential in energy (bin width) per primary.}
  \item{N.per.primary}{ Fluence per primary.}
  and (2) variables containing information as beam energy, peak position, projectile and target material.
}

\seealso{
\code{\link{AT.SPC.tapply}} for more example to derive quantaties from spc-data.
}

\examples{
\dontrun{
# Load required libraries
require(libamtrack)
require(lattice)

# Use example data set
file.name <- system.file("extdata", "libamtrack.12C.H2O.active3.MeV27000.zip", package = "libamtrack")
endian    <- c("big", "little")[2]


# Read in spc data, we use old style R reader, so endianess has to be given
spc                    <- AT.SPC.read ( file.name = file.name,
                                        endian    = endian,
                                        flavour   = "vanilla")$spc

# Get name and Z of particles in spectrum
spc$particle.name      <- AT.particle.name.from.particle.no(particle.no = spc$particle.no)
spc$particle.name      <- factor(spc$particle.name, levels = unique(spc$particle.name))
spc$Z                  <- AT.Z.from.particle.no(particle.no = spc$particle.no)$Z

                                   
#####################################################################
# Plot parts of spectra
ii      <- spc$depth.step%in%c(1,10,30,40,60,70)

#xyplot( log10(dN.dE.per.MeV.u.per.primary) ~ E.mid.MeV.u|sprintf("%4d", depth.step),
#        spc[ii,],
#        type        = 's',
#        subset      = log10(dN.dE.per.MeV.u.per.primary) >= -4,
#        groups      = particle.name,
#        as.table    = TRUE,
#        ylab        = 'log (fluence per primary per bin width / (1/MeV.u)',
#        xlab        = 'kinetic energy / MeV.u',
#        auto.key    = list(space = 'right', lines = TRUE, points = FALSE),
#        main        = 'Energy spectra at selected depth steps')


#####################################################################
# Process data, compute total dose and dose-weighted LET with depth
df.total                <- AT.SPC.tapply( spc                  = spc, 
                                          INDEX                = "depth.g.cm2", 
                                          FUN                  = AT.total.D.Gy, 
                                          additional.arguments = list(c("material.no", "AT.material.no.from.material.name('Water, Liquid')", FALSE),
                                                                      c("stopping.power.source.no", "0", FALSE)),
                                          names.results        = "D.Gy")
df.total$LET.MeV.cm2.g  <- AT.SPC.tapply( spc                  = spc, 
                                          INDEX                = "depth.g.cm2", 
                                          FUN                  = AT.dose.weighted.LET.MeV.cm2.g,
                                          additional.arguments = list(c("material.no", "AT.material.no.from.material.name('Water, Liquid')", FALSE),
                                                                      c("stopping.power.source.no", "0", FALSE)))$returnValue
df.total$fLET.MeV.cm2.g <- AT.SPC.tapply( spc                  = spc, 
                                          INDEX                = "depth.g.cm2", 
                                          FUN                  = AT.fluence.weighted.LET.MeV.cm2.g,
                                          additional.arguments = list(c("material.no", "AT.material.no.from.material.name('Water, Liquid')", FALSE),
                                                                      c("stopping.power.source.no", "0", FALSE)))$returnValue
df.total$Mean.E.C.MeV.u <- AT.SPC.tapply( spc                  = spc[spc$particle.no == 6012,], 
                                          INDEX                = "depth.g.cm2", 
                                          FUN                  = AT.dose.weighted.E.MeV.u,
                                          additional.arguments = list(c("material.no", "AT.material.no.from.material.name('Water, Liquid')", FALSE),
                                                                      c("stopping.power.source.no", "0", FALSE)))$returnValue

xyplot( D.Gy ~ depth.g.cm2,
        df.total,
        grid        = TRUE,
        type        = 'o',
        ylab        = 'total dose / Gy',
        xlab        = 'depth / (g/cm2)')

xyplot( LET.MeV.cm2.g + fLET.MeV.cm2.g ~ depth.g.cm2,
        df.total,
        type        = 'o',
        main        = 'Fluence (fLET) and dose (LET) weighted LET',
        grid        = TRUE,
        auto.key    = list(space = 'right', lines = TRUE, points = FALSE),
        ylab        = 'LET / (MeV*cm2/g)',
        xlab        = 'depth / (g/cm2)')

xyplot( Mean.E.C.MeV.u ~ depth.g.cm2,
        df.total,
        type        = 'o',
        grid        = TRUE,
        main        = file.name,
        ylab        = 'Average kinetic energy of C / (MeV/u)',
        xlab        = 'depth / (g/cm2)')

}
}

