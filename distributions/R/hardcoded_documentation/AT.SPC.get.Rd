\name{AT.SPC.get}
\alias{AT.SPC.get}
\title{AT.SPC.get}
\description{Returns interpolated spc data. Needs list of spc files as input which is obtained by \link{AT.SPC.get.list}.
TODO: use path.to.files here alternative argument and call AT.SPC.get.list internally.}
\usage{AT.SPC.get(spc.list, energy.MeV.u)
}
\arguments{
  \item{spc.list}{ data frame with spc file information as returned by \link{AT.SPC.get.list}}
  \item{energy.MeV.u}{ requested beam energy, has to be within energy grid covered by spc.list}
}
\value{
  list with spc data equivalent to those returned by \link{AT.SPC.read}.
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

