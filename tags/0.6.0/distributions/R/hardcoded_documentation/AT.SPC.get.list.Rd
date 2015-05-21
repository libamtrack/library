\name{AT.SPC.get.list}
\alias{AT.SPC.get.list}
\title{AT.SPC.get.list}
\description{Returns filenames and essential information (energy, projectile etc.) of all spc files that are found at a user-proveded location. Needed by \link{AT.SPC.get}.
}
\usage{AT.SPC.get.list(path.to.files, endian = "little")
}
\arguments{
  \item{path.to.files}{ absolute path to spc files, including wildcard, e.g. *.active*.spc}
  \item{endian}{ byte-order of the spc file to read, "big" (AIX) or "little" (Linux, VMS, etc., default). Only necessary if flavour = 'vanilla'}
}
\value{
  A dataframe with path and file names, number of depth steps, projectile, target, beam energy, peak position.
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

