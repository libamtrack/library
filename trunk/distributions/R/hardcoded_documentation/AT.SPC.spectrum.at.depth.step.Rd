\name{AT.SPC.spectrum.at.depth.step}
\alias{AT.SPC.spectrum.at.depth.step}
\title{AT.SPC.spectrum.at.depth.step}
\description{Returns spectrum from spc data in given depth step.}
\usage{AT.SPC.spectrum.at.depth.step( spc, depth.step)
}
\arguments{
  \item{spc}{ spc data}
  \item{depth.step}{ depth step the spectrum should be taken from (if outside data scope an errormessage will be issued)}
}
\value{
  A data frame with the following columns (ready to use in most libamtrack functions):
  \item{E.MeV.u}{ Energy [MeV/u] - mid point of bin}
  \item{particle.no}{ Particle index number}
  \item{fluence.cm2}{ Fluence for each bin}
}

\seealso{
Use \code{\link{AT.SPC.spectrum.at.depth.g.cm2}} to get spectrum at given depth. spc data should be read in by \code{\link{AT.SPC.read}}.}
\examples{
# None yet.
}
