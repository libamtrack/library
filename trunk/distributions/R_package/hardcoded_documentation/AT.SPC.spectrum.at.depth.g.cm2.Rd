\name{AT.SPC.spectrum.at.depth.g.cm2}
\alias{AT.SPC.spectrum.at.depth.g.cm2}
\title{AT.SPC.spectrum.at.depth.g.cm2}
\description{Returns spectrum from spc data in given depth step.}
\usage{AT.SPC.spectrum.at.depth.g.cm2( spc, depth.g.cm2)
}
\arguments{
  \item{spc}{ spc data}
  \item{depth.g.cm2}{ depth the spectrum should be taken at (up to now, the data from the closest depth step will be taken)}
}
\value{
  A data frame with the following columns (ready to use in most libamtrack functions):
  \item{E.MeV.u}{ Energy [MeV/u] - mid point of bin}
  \item{particle.no}{ Particle index number}
  \item{fluence.cm2}{ Fluence for each bin}
}

\seealso{
Uses \code{\link{AT.SPC.spectrum.at.depth.step}}. spc data should be read in by \code{\link{AT.SPC.read}}.}
\examples{
# None yet.
}
