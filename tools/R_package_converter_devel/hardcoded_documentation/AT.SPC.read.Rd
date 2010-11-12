\name{AT.SPC.read}
\alias{AT.SPC.read}
\title{AT.SPC.read}
\description{Read a spc-formatted data file with energy-fluence in n depth steps. For original TRiP format definition by M. Kraemer, please see \url{http://bio.gsi.de/DOCS/TRiP98BEAM/DOCS/trip98fmtspc.html}. Please note that the user has to take care of picking the spc-file for the projectile and target material desired.}
\usage{AT.SPC.read( file.name, endian, mean, compress)
}
\arguments{
  \item{file.name}{ name of spc-file to be read (with extension).}
  \item{endian}{ byte-order in spc file, "big" or "little" depending on the machine of origin.}
  \item{mean}{ method for computing bin midpoints, "geometric" or "arithmetic" (default).}
  \item{compress}{ if TRUE all entries with zero fluence will be removed to save memory.}
}
\value{
  A data frame with the following columns
  \item{depth.step}{ Index number of depth step (one-based)}
  \item{depth.g.cm2}{ Depth in g/cm2}
  \item{particle.no}{ Particle index number (see also(\code{\link{particle.no}}).}
  \item{E.MeV.u}{ Energy bin mid point in MeV/u (see also(\code{\link{E.MeV.u}}).}
  \item{DE.MeV.u}{ Energy bin width in MeV/u}
  \item{fluence.cm2}{ Fluence in 1/cm2 for each bin (NOT normalized per bin width).}
}
\examples{
# Read in TRiP data from 142 MeV protons and
# reformat/compress for use with libamtrack
# functions
# 
# spc <- AT.SPC.read( file.name = "TRIP_GSI_1H.H2O.MeV14260.spc",
#                     endian    = "little",
#                     mean      = 1,
#                     raw       = FALSE,
#                     compress  = TRUE)
#
#
# Read in same data but give complete overview on data content
# 
# spc.raw <- AT.SPC.read( file.name = "TRIP_GSI_1H.H2O.MeV14260.spc",
#                         endian    = "little",
#                         mean      = 1,
#                         raw       = TRUE,
#                         compress  = FALSE)
}
