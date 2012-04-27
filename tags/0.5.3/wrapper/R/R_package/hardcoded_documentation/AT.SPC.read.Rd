\name{AT.SPC.read}
\alias{AT.SPC.read}
\title{AT.SPC.read}
\description{Read a spc-formatted data file with energy-fluence in n depth steps. For original TRiP format definition by M. Kraemer, please see \url{http://bio.gsi.de/DOCS/TRiP98BEAM/DOCS/trip98fmtspc.html}. Please note that the user has to take care of picking the spc-file for the projectile and target material desired.

Presently, two versions of the reader exists: a slower but stable R version ('vanilla') and a faster but still buggy C version.

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
  \item{E.MeV.u}{ Energy bin mid point in MeV/u (see also(\code{\link{E.MeV.u}}).}
  \item{DE.MeV.u}{ Energy bin width in MeV/u}
  \item{fluence.cm2}{ Fluence in 1/cm2 for each bin (NOT normalized per bin width).}
  and (2) variables containing information as beam energy, peak position, projectile and target material.
}
\examples{
## Download data from libamtrack homepage 
##
## To read a simple SPC file, use AT.SPC.read
##
## To get SPC data for arbitrary energy, you first have to make libamtrack aware
## of the appropriate SPC files:
#
#  spc.list <- AT.SPC.get.list(path.to.files,...) # where 'path.to.files' should point to spc files
#                                                 # OF ONE KIND, i.e. same projectile, target, active/passive
#                                                 # but different energies
## N.B.: This routine can also be used to browser spc files and check their integrity
##
## Using the returned list, you can get spc data for any energy in between
#
#  spc <- AT.SPC.get(spc.list, energy = ...)      # This will can AT.SPC.read and AT.SPC.interpolate
#                                                
## returns a list in which the actual data are found in spc$spc
## which can then used in spc related routines, e.g.
## AT.SPC.tapply, AT.SPC.convert.to.DDD, etc.
##
## A spectrum at arbitraty can be obtain by using the data in spc$spc
## with AT.SPC.spectrum.at.depth.g.cm2
}

