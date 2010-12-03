\name{AT.FLUKA.read.USRTRACK}
\alias{AT.FLUKA.read.USRTRACK}
\title{AT.FLUKA.read.USRTRACK}
\description{Reads USRTRACK output (for energy spectra scoring using fluscw.SPC.f) for a series of regions, also for multiple output files from cluster runs (using rcfluka.py).}
\usage{AT.FLUKA.read.USRTRACK(exp.name, number.of.runs, unit, data.source = 'local', compress = TRUE, vol.cm3 = NULL)
}
\arguments{
  \item{exp.name}{ Experiment name, i.e. name of input file (without '.inp' extension)}
  \item{number.of.runs}{Number of output files from parallel (cluster) runs.}
  \item{unit}{FLUKA output unit number}
  \item{data.source}{'local' if output files are from a local machine, 'condor' if from condor cluster, 'condor_cleaned' if from condor cluster with clean option (-c) in \code{rcfluka.py}}
  \item{compress}{If TRUE, all entries with zero fluence will be removed from resulting data frame in order to save memory.}
  \item{vol.cm3}{Volume of regions, either single value (will be applied to all regions) or vector of length matching the number of regions (individual volume for each region).}
}
\value{
  Data frame with midpoints and widths of energy bins, particle index no \code{\link{particle.no}}, and fluence (absolute, i.e. not normalized to bin width!) for each region.
}

\examples{
# None yet, requires FLUKA output file
}
