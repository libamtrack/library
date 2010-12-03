\name{AT.FLUKA.read.USRBIN.regs}
\alias{AT.FLUKA.read.USRBIN.regs}
\title{AT.FLUKA.read.USRBIN.regs}
\description{Reads USRBIN output for a series of regions, also for multiple output files from cluster runs. As USRBIN (regions, option 12.0) scores energy deposited per primary weight (GeV), both the region volume(s) and density(ies) have to be given for dose computation.}
\usage{AT.FLUKA.read.USRBIN.regs(exp.name, number.of.runs, unit, data.source = 'local', vol.cm3 = NULL, density.g.cm3 = NULL)
}
\arguments{
  \item{exp.name}{ Experiment name, i.e. name of input file (without '.inp' extension)}
  \item{number.of.runs}{Number of output files from parallel (cluster) runs.}
  \item{unit}{FLUKA output unit number}
  \item{data.source}{'local' if output files are from a local machine, 'condor' if from condor cluster, 'condor_cleaned' if from condor cluster with clean option (-c) in \code{rcfluka.py}}
  \item{vol.cm3}{Volume of regions, either single value (will be applied to all regions) or vector of length matching the number of regions (individual volume for each region).}
  \item{density.g.cm3}{Physical density of regions, either single value (will be applied to all regions) or vector of length matching the number of regions (individual density for each region).}
}
\value{
  Data frame with deposited energy, deposited dose (and in case of multiple runs estimate of their standard deviation) for each region.
}

\examples{
# None yet, requires FLUKA output file
}
