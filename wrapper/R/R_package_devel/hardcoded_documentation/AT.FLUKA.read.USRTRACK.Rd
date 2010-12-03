\name{AT.FLUKA.read.USRTRACK}
\alias{AT.FLUKA.read.USRTRACK}
\title{AT.FLUKA.read.USRTRACK}
\description{Reads USRTRACK output (for energy spectra scoring using fluscw.SPC.f) for a series of regions, also for multiple output files from cluster runs (using rcfluka.py). USRTRACK scores the fluence per primary weight per bin width(cm-2/GeV) but only if the correct volume of the corresponding regions is given in the FLUKA input card! AT.FLUKA.read.USRTRACK will however out the absolute fluence for each energy bin (in cm-2) as this can be directly used by libamtrack.}
\usage{AT.FLUKA.read.USRTRACK(exp.name, number.of.runs, unit, data.source = 'local', compress = TRUE)
}
\arguments{
  \item{exp.name}{ Experiment name, i.e. name of input file (without '.inp' extension)}
  \item{number.of.runs}{Number of output files from parallel (cluster) runs.}
  \item{unit}{FLUKA output unit number}
  \item{data.source}{'local' if output files are from a local machine, 'condor' if from condor cluster, 'condor_cleaned' if from condor cluster with clean option (-c) in \code{rcfluka.py}}
  \item{compress}{If TRUE, all entries with zero fluence will be removed from resulting data frame in order to save memory.}
}
\value{
  Data frame with midpoints and widths of energy bins, particle index no \code{\link{particle.no}}, and fluence (absolute [cm-2], i.e. not normalized to bin width!) for each region.
}

\examples{
# None yet, requires FLUKA output file
}
