\name{AT.SPC.interpolate}
\alias{AT.SPC.interpolate}
\title{AT.SPC.interpolate}
\description{Interpolates between two spc files.
}
\usage{AT.SPC.interpolate(spc.lower, spc.upper, energy.MeV.u)
}
\arguments{
  \item{spc.lower}{ spc data (read by AT.SPC.read) with lower bracketing energy.}
  \item{spc.upper}{ spc data (read by AT.SPC.read) with upper bracketing energy.}
  \item{energy.MeV.u}{ energy of interpolated spc data}
}
\value{
  list equivalent to those returned by AT.SPC.read.
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

