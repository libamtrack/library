\name{AT.FLUKA.read.USRBIN.mesh}
\alias{AT.FLUKA.read.USRBIN.mesh}
\title{AT.FLUKA.read.USRBIN.mesh}
\description{Reads USRBIN output for Cartesian mesh, also for multiple output files from cluster runs. As USRBIN (mesh) scores energy deposited per unit volume per primary weight (GeV/cm3), the density as to be given for dose computation.}
\usage{AT.FLUKA.read.USRBIN.mesh(exp.name, number.of.runs, unit, data.source = 'local', density.g.cm3 = 1.0)
}
\arguments{
  \item{exp.name}{ Experiment name, i.e. name of input file (without '.inp' extension)}
  \item{number.of.runs}{Number of output files from parallel (cluster) runs.}
  \item{unit}{FLUKA output unit number}
  \item{data.source}{'local' if output files are from a local machine, 'condor' if from condor cluster, 'condor_cleaned' if from condor cluster with clean option (-c) in \code{rcfluka.py}}
  \item{density.g.cm3}{Physical density of material in mesh, needed for computation of dose.}
}
\value{
  Data frame with deposited energy, deposited dose (and in case of multiple runs estimate of their standard deviation) for each mesh cell.
}

\examples{
# None yet, requires FLUKA output file
}
