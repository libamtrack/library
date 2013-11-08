\name{AT.SPC.read}
\alias{AT.SPC.read}
\title{AT.SPC.read}
\description{Read a spc-formatted data file with energy-fluence in n depth steps. For original TRiP format definition by M. Kraemer, please see \url{http://bio.gsi.de/DOCS/TRiP98BEAM/DOCS/trip98fmtspc.html} and Kraemer and Scholz, Treatment planning for heavy-ion radiotherapy: calculation and optimization of biologically effective dose, Phys. Med. Biol. 45 (2000) 3319–3330. Please note that the user has to take care of picking the spc-file for the projectile and target material desired.

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
  		names(df)    <- c("depth.step", "depth.g.cm2", "particle.no", "E.low.MeV.u", "E.mid.MeV.u", "E.high.MeV.u", "dE.MeV.u", "dN.dE.per.primary", "dN.E.MeV.u.per.primary")
  \item{E.low.MeV.u}{ Energy bin lower limit in MeV/u (see also(\code{\link{E.MeV.u}}).}
  \item{E.mid.MeV.u}{ Energy bin mid in MeV/u (see also(\code{\link{E.MeV.u}}).}
  \item{E.high.MeV.u}{ Energy bin upper limit in MeV/u (see also(\code{\link{E.MeV.u}}).}
  \item{dE.MeV.u}{ Energy bin width in MeV/u}
  \item{dN.dE.per.primary}{ Fluence in 1/cm2 for each bin (NOT normalized per bin width).}
  and (2) variables containing information as beam energy, peak position, projectile and target material.
}
\examples{
\dontrun{

rm(list = ls())

# Load required libraries, libamtrack > r900
require(libamtrack)
require(lattice)

filename <- system.file("extdata", "libamtrack.12C.H2O.active3.MeV27000.spc", package = "libamtrack")
endian                 <- c("big", "little")[2]

# Read in spc data, we use old style R reader, so endianess has to be given
spc                    <- AT.SPC.read ( file.name = file.name,
                                        endian    = endian,
                                        flavour   = "vanilla")$spc

# Get name and Z of particles in spectrum
spc$particle.name      <- AT.particle.name.from.particle.no(particle.no = spc$particle.no)
spc$Z                  <- AT.Z.from.particle.no(particle.no = spc$particle.no)$Z

# Plot parts of spectra

spc <- spc[,-4]

ii      <- spc$depth.step%%10 == 0 & spc$depth.step != 70

xyplot( DN.DE.MeV.u ~ E.MeV.u|sprintf("%4d", depth.step),
        spc[ii,],
        type        = 's',
        groups      = particle.no,
        as.table    = TRUE,
        auto.key    = list( space = 'right'),
        main        = "Energy-Fluence Correlation in various Depths" )

xyplot( log10(DN.DE.MeV.u) ~ E.MeV.u|sprintf("%4d", depth.step),
        spc[ii,],
        type        = 's',
		    subset      = log10(fluence.cm2) >= -6,
        groups      = particle.name,
        as.table    = TRUE,
        auto.key    = list(space = 'right', c("blue", "darkgreen", "red", "orange", "purple", "cyan")),
        main        = "Energy-Fluence Correlation in various Depths" )

xyplot( DN.DE.MeV.u ~ E.MeV.u,
        aggregate(DN.DE.MeV.u ~ E.MeV.u, subset(spc, particle.no==1002), sum),
        type        = 'o',
        col         = "darkgreen",
        xlim        = c(0,300),
        main        = "2H Energy Spectrum \n (cut-off at 294.68 MeV/u)" )

xyplot( DN.DE.MeV.u ~ E.MeV.u,
        aggregate(DN.DE.MeV.u ~ E.MeV.u, subset(spc, particle.no==2004), sum),
        type        = 'o',
        col         = "red",
        xlim        = c(0,300),
        main        = "4He Energy Spectrum \n (cut-off at 294.68 MeV/u)" )

xyplot( DN.DE.MeV.u ~ E.MeV.u,
        aggregate(DN.DE.MeV.u ~ E.MeV.u, subset(spc, particle.no==3006), sum),
        type        = 'o',
        col         = "orange",
        xlim        = c(0,300),
        main        = "6Li Energy Spectrum \n (cut-off at 294.68 MeV/u)" )


xyplot( DN.DE.MeV.u ~ E.MeV.u,
        aggregate(DN.DE.MeV.u ~ E.MeV.u, subset(spc, particle.no==4008), sum),
        type        = 'o',
        col         = "green",
        xlim        = c(0,300),
        main        = "8Be Energy Spectrum \n (cut-off at 294.68 MeV/u)" )

xyplot( DN.DE.MeV.u ~ E.MeV.u,
        aggregate(DN.DE.MeV.u ~ E.MeV.u, subset(spc, particle.no==5010), sum),
        type        = 'o',
        col         = "blue",
        xlim        = c(0,300),
        main        = "10B Energy Spectrum \n (cut-off at 294.68 MeV/u)" )

xyplot( DN.DE.MeV.u ~ E.MeV.u,
        aggregate(DN.DE.MeV.u ~ E.MeV.u, subset(spc, particle.no==6012), sum),
        type        = 'o',
        col         = "purple",
        xlim        = c(0,300),
        main        = "12C Energy Spectrum \n (cut-off at 296.71 MeV/u)" )


#########
# Process data, compute total dose and dose-weighted LET with depth
df.total                <- AT.SPC.tapply( spc                  = spc, 
                                          INDEX                = "depth.g.cm2", 
                                          FUN                  = AT.total.D.Gy, 
                                          additional.arguments = list(c("material.no", "AT.material.no.from.material.name('Water, Liquid')", F),
                                                                      c("stopping.power.source.no", "0", F)),
                                          names.results        = "D.Gy")
df.total$LET.MeV.cm2.g  <- AT.SPC.tapply( spc                  = spc, 
                                          INDEX                = "depth.g.cm2", 
                                          FUN                  = AT.dose.weighted.LET.MeV.cm2.g,
                                          additional.arguments = list(c("material.no", "AT.material.no.from.material.name('Water, Liquid')", F),
                                                                      c("stopping.power.source.no", "0", F)))$returnValue
df.total$fLET.MeV.cm2.g <- AT.SPC.tapply( spc                  = spc, 
                                          INDEX                = "depth.g.cm2", 
                                          FUN                  = AT.fluence.weighted.LET.MeV.cm2.g,
                                          additional.arguments = list(c("material.no", "AT.material.no.from.material.name('Water, Liquid')", F),
                                                                      c("stopping.power.source.no", "0", F)))$returnValue
df.total$Mean.E.C.MeV.u <- AT.SPC.tapply( spc                  = spc, #spc[spc$particle.no == 6012,], 
                                          INDEX                = "depth.g.cm2", 
                                          FUN                  = AT.dose.weighted.E.MeV.u,
                                          additional.arguments = list(c("material.no", "AT.material.no.from.material.name('Water, Liquid')", F),
                                                                      c("stopping.power.source.no", "0", F)))$returnValue

xyplot( D.Gy ~ depth.g.cm2,
        df.total,
        type   = 'o',
        ylab        = 'total dose / Gy',
        xlab        = 'depth / (g/cm2)')
        
xyplot( LET.MeV.cm2.g + fLET.MeV.cm2.g ~ depth.g.cm2,
        df.total,
        type   = 'o',
        main        = file.name,
        grid        = TRUE,
        auto.key    = list(space = 'right'),
        ylab        = 'LET / (MeV*cm2/g)',
        xlab        = 'depth / (g/cm2)')

xyplot( Mean.E.C.MeV.u ~ depth.g.cm2,
        df.total,
        type   = 'o',
        main        = file.name,
        ylab        = 'Mean energy / (MeV/u)',
        xlab        = 'depth / (g/cm2)')

#########
# Process data, differential in particle type

df.T                  <- AT.SPC.tapply( spc                  = spc, 
                                        INDEX                = c("depth.g.cm2", "particle.no"), 
                                        FUN                  = AT.total.D.Gy, 
                                        additional.arguments = list(c("material.no", "AT.material.no.from.material.name('Water, Liquid')", FALSE),
																	      c("stopping.power.source.no", "2", FALSE)),
                                        names.results        = "D.Gy")
df.T$fluence.cm2      <- AT.SPC.tapply( spc                  = spc, 
                                        INDEX                = c("depth.g.cm2", "particle.no"), 
                                        FUN                  = mean,
                                        additional.arguments = list(c("material.no", "AT.material.no.from.material.name('Water, Liquid')", FALSE),
                                                                    c("x", "fluence.cm2", TRUE)))$res
df.T$LET.MeV.cm2.g    <- AT.SPC.tapply( spc                  = spc, 
                                        INDEX                = c("depth.g.cm2", "particle.no"), 
                                        FUN                  = AT.dose.weighted.LET.MeV.cm2.g,
                                        additional.arguments = list(c("material.no", "AT.material.no.from.material.name('Water, Liquid')", FALSE),
																	                                  c("stopping.power.source.no", "2", FALSE)))$returnValue

spc$cv                <- paste(spc$depth.g.cm2, spc$particle.no)
df.T$depth.step       <- tapply( X                    = spc$depth.step,
                                 INDEX                = spc$cv,
                                 FUN                  = head,
                                 n                    = 1)
df.T$particle.name    <- tapply( X                    = spc$particle.name,
                                 INDEX                = spc$cv,
                                 FUN                  = head,
                                 n                    = 1)
df.T$Z                <- tapply( X                    = spc$Z,
                                 INDEX                = spc$cv,
                                 FUN                  = head,
                                 n                    = 1)
                                 
# normalize dose to maximum / normalize fluence to C-12 entrace channel
norm.dose.Gy          <- max(df.total$D.Gy)
#norm.dose.Gy          <- df.T$D.Gy[df.T$particle.no == 6012 & df.T$depth.step == 2]
norm.fluence.cm2      <- df.T$fluence.cm2[df.T$particle.no == 6012 & df.T$depth.step == 2]

df.T$D.Gy             <- df.T$D.Gy / norm.dose.Gy
df.T$fluence.cm2      <- df.T$fluence.cm2 / norm.fluence.cm2


# Plots
xyplot( log10(D.Gy) ~ depth.g.cm2,
        df.T,
        type       = 'o',
        groups     = particle.name,
		    subset     = log10(D.Gy) >= -6,
        auto.key   = list(space = 'right'),
        panel      = function(...){
                         panel.grid(-1, -1)
                         panel.xyplot(...)},
        xlab       = 'depth / (g/cm2)',
        ylab       = 'log dose relative to C-12 entrance channel')

xyplot( log10(fluence.cm2) ~ depth.g.cm2,
        df.T,
        type       = 'o',
		    subset     = log10(fluence.cm2) >= -6,
        groups     = particle.name,
        auto.key   = list(space = 'right'),
        panel      = function(...){
                         panel.grid(-1, -1)
                         panel.xyplot(...)},
        xlab       = 'depth / (g/cm2)',
        ylab       = 'log fluence relative to C-12 entrance channel')

xyplot( log10(LET.MeV.cm2.g) ~ depth.g.cm2,
        df.T,
        type       = 'o',
        groups     = particle.name,
        auto.key   = list(space = 'right'),
        panel      = function(...){
                         panel.grid(-1, -1)
                         panel.xyplot(...)},
        xlab       = 'depth / (g/cm2)',
        ylab       = 'log dose-weighted LET / (MeV*cm2/g)')



df.total$Z             <- tapply(        X                    = spc$Z, 
                                         INDEX                = spc$depth.g.cm2, 
                                         FUN                  = head,
                                         n                    = 1)

}
}

